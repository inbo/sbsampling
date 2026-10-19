use extendr_api::prelude::*;
use geo::Geometry;
use geo_traits::to_geo::ToGeoGeometry;
use h3o::geom::{ContainmentMode, TilerBuilder};
use h3o::{CellIndex, Resolution};
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_pcg::Pcg64Mcg;
use std::collections::BinaryHeap;
use std::convert::TryInto;
use wkb::reader::Wkb;
// ---------------------------------------------------------------------------
// H3 bit-layout constants (from the H3 spec)
// ---------------------------------------------------------------------------

/// Bit position where the base cell index starts in an H3 index.
const H3_BASE_CELL_SHIFT: u32 = 45;

/// 7-bit mask to extract the base cell number (values 0–121).
const H3_BASE_CELL_MASK: u64 = 0b111_1111;

/// Maximum valid base cell index.
const H3_MAX_BASE_CELL: usize = 121;

/// Number of base cells in the H3 grid.
const H3_NUM_BASE_CELLS: usize = 122;

/// 3-bit mask to extract a single resolution digit (child index 0–6).
const H3_DIGIT_MASK: u64 = 0b111;

/// Number of bits used to store the base cell in the sort key.
const SORT_KEY_BASE_BITS: u32 = 7;

/// Number of bits used per resolution level in the sort key.
const SORT_KEY_BITS_PER_RES: u32 = 3;

// ---------------------------------------------------------------------------
// Type aliases
// ---------------------------------------------------------------------------

/// One entry in the GRTS reservoir heap.
///
/// Fields are kept as raw `u64` bits so the tuple is `Copy + Ord` and can
/// live directly in a `BinaryHeap` without a wrapper type.
///
/// * `effective_key_bits` - `f64::to_bits(sort_key / area)` (or bare
///   `sort_key`): drives eviction; positive finite f64 bit patterns
///   preserve numeric order.
/// * `raw_sort_key` - The un-weighted GRTS sort key, preserved for
///   display in the output data frame.
/// * `h3_index` - The raw `u64` H3 cell index.
/// * `area_bits`- `f64::to_bits(cell_area_m2)`, carried through
///   so areas are available after the loop without re-querying h3o.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct HeapEntry {
    effective_key_bits: u64,
    raw_sort_key: u64,
    h3_index: u64,
    area_bits: u64,
}

impl HeapEntry {
    #[inline]
    fn new(effective_key_bits: u64, raw_sort_key: u64, h3_index: u64, area_bits: u64) -> Self {
        Self {
            effective_key_bits,
            raw_sort_key,
            h3_index,
            area_bits,
        }
    }
}

/// Output of `grts_core`: the sorted reservoir, the total cell count, and
/// the sum of all cell areas in m² across the full coverage.
type GrtsCoreResult = (Vec<HeapEntry>, usize, f64);

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

/// geometry -> coverage
/// @noRd
#[extendr]
pub fn h3_coverage(wkb_bytes: &[u8], res: i32, containment: &str) -> extendr_api::Result<Vec<u8>> {
    let geometry = decode_wkb(wkb_bytes)?;
    let containment_mode = parse_containment(containment)?;
    let resolution = Resolution::try_from(res as u8)
        .map_err(|_| Error::Other(format!("Invalid H3 resolution: {}", res)))?;

    let mut tiler = TilerBuilder::new(resolution)
        .containment_mode(containment_mode)
        .build();
    feed_geometry(&mut tiler, geometry)?;

    // Write the raw 8-byte u64 chunks directly. No headers, no bincode!
    let bytes: Vec<u8> = tiler
        .into_coverage()
        .flat_map(|cell| u64::from(cell).to_le_bytes())
        .collect();

    Ok(bytes)
}
/// Generate a GRTS sample (Internal C-ABI function)
///
/// Geometry -> Coverage -> Sample
/// This function is wrapped by `grts_h3()` in R and should not be
/// called directly by the user.
/// Generate a GRTS (Generalized Random Tessellation Stratified) sample of H3
/// cell indices covering a polygon.
///
/// # Arguments
/// * `wkb_bytes`      – Study area geometry encoded as Well-Known Binary.
/// * `res`            – H3 resolution (0 = coarsest, 15 = finest).
/// * `containment`    – One of `"intersect"`, `"centroid"`, `"boundary"`, `"covers"`.
/// * `n`              – Desired sample size.
/// * `global_seed`    – Reproducibility seed (passed as f64 from R; values above
///                      2^53 lose precision — prefer seeds within that range).
/// * `area_correction`– If `false` (default), all H3 cells have equal inclusion
///                      probability pi_i = n / N. Cell areas vary by up to ~1.9x
///                      across the icosahedral projection; this is ignored in the
///                      sampling step but can be corrected at the estimation stage
///                      using the returned `area_m2` and `ip` columns.
///                      If `true`, sequential poisson sampling is applied, making the
///                      effective inclusion probability proportional to area:
///                      pi_i = n * area_i / sum(area).
///
/// # Returns
/// A `data.frame` with four columns ordered by ascending sort key (GRTS order):
/// * `cell`     – H3 cell identifier as a lowercase hexadecimal string.
/// * `sort_key` – GRTS sort key as a zero-padded 16-character lowercase hex
///                string. Zero-padding ensures correct lexicographic ordering
///                on the R side (e.g. for merging tile results).
/// * `area_m2` – True area of the sampled cell in m².
/// * `ip`      – Inclusion probability ip_i of the sampled cell:
///                  area_correction = false: ip_i = n / N  (constant)
///                  area_correction = true:  ip_i = n * area_i / sum(area)
///
/// Design-based estimators using the returned columns:
///   HT total:   T_hat   = sum(y_i * area_m2_i / ip_i)
///   Hajek mean: mu_Hajek = sum(y_i * area_m2_i / ip_i) / sum(area_m2_i / ip_i)
///
/// The returned `data.frame` also carries the following design attributes:
/// * `n_cells`        – Total H3 cells N covering the study area.
/// * `sum_area_m2`    – Sum of all cell areas in m².
/// * `n`              – Requested sample size.
/// * `resolution`     – H3 resolution used.
/// * `containment`    – Containment mode used.
/// * `area_correction`– Whether area-proportional sampling was applied.
/// @noRd
#[extendr]
pub fn grts_sample_from_wkb(
    wkb_bytes: &[u8],
    res: i32,
    containment: &str,
    n: i32,
    global_seed: f64,
    area_correction: bool,
) -> extendr_api::Result<Robj> {
    // NOTE: f64 -> u64 truncates; seeds above 2^53 lose precision.
    let seed_u64 = global_seed as u64;
    let target_n = n as usize;

    // 1. Decode WKB and resolve tiler settings ------------------------------
    let geometry = decode_wkb(wkb_bytes)?;
    let containment_mode = parse_containment(containment)?;
    let resolution = Resolution::try_from(res as u8)
        .map_err(|_| Error::Other(format!("Invalid H3 resolution: {}", res)))?;

    // 2. Build the main tiler ----------------------------------------------
    let mut tiler = TilerBuilder::new(resolution)
        .containment_mode(containment_mode)
        .build();
    feed_geometry(&mut tiler, geometry)?;

    // 3. Shuffle base cells (top level of the GRTS address) ----------------
    // This randomises which continental regions receive low sort keys,
    // ensuring global spatial balance across the entire Earth surface.
    let base_cells = shuffled_base_cells(seed_u64);

    // 4. Stream cells through the heap reservoir ----------------------------
    // The heap stores (sort_key, h3_index, area_bits) tuples. BinaryHeap is a
    // max-heap so peek() always returns the *largest* key — exactly what we
    // want to evict when a smaller key arrives. area_bits is the cell area
    // (m2) bit-cast to u64 so it can be stored without a separate Vec.
    let (sorted, total_cells, sum_all_areas) = grts_core(
        tiler.into_coverage(),
        target_n,
        res,
        seed_u64,
        area_correction,
        &base_cells,
    )?;

    build_result_df(
        sorted,
        total_cells,
        sum_all_areas,
        n,
        res,
        containment,
        area_correction,
    )
}

/// Pre-computed cells -> Sample
/// @noRd
#[extendr]
pub fn grts_sample_from_cells(
    cells_bytes: &[u8], // Pure raw bytes, no serialization headers
    n: i32,
    global_seed: f64,
    area_correction: bool,
) -> extendr_api::Result<Robj> {
    let seed_u64 = global_seed as u64;
    let target_n = n as usize;

    let (cells, n_dup) = parse_cells(cells_bytes)?;
    let res = cells[0].resolution() as i32;
    let base_cells = shuffled_base_cells(seed_u64);

    let (sorted, total_cells, sum_all_areas) = grts_core(
        cells.iter().copied(),
        target_n,
        res,
        seed_u64,
        area_correction,
        &base_cells,
    )?;

    let mut df = build_result_df(
        sorted,
        total_cells,
        sum_all_areas,
        n,
        res,
        "pre-computed",
        area_correction,
    )?;
    df.set_attrib("n_duplicates_removed", n_dup as i32)?;
    Ok(df)
}

// ---------------------------------------------------------------------------
// Shared Core Logic
// ---------------------------------------------------------------------------

fn grts_core(
    cell_iter: impl Iterator<Item = CellIndex>,
    target_n: usize,
    res: i32,
    seed_u64: u64,
    area_correction: bool,
    base_cells: &[u64; H3_NUM_BASE_CELLS],
) -> extendr_api::Result<GrtsCoreResult> {
    // HeapEntry stores effective_key_bits = f64::to_bits(effective_key), where:
    //   area_correction = false: effective_key = sort_key as f64
    //   area_correction = true:  effective_key = sort_key as f64 / cell_area
    //
    // For all positive finite f64, to_bits() preserves order, so the
    // max-heap correctly evicts the cell with the largest effective key.
    let mut heap: BinaryHeap<HeapEntry> = BinaryHeap::with_capacity(target_n + 1);
    // This 16-element array
    // stores (parent_id, [permutation]) for each of the 16 H3 resolutions.
    let mut last_seen_perm = [(u64::MAX, [0u64; 7]); 16];

    let mut total_cells: usize = 0;
    // Accumulate area over the FULL coverage regardless of acceptance, since
    // sum(area) — the denominator in pi_i — refers to the whole study area.
    let mut sum_all_areas: f64 = 0.0;

    for cell in cell_iter {
        total_cells += 1;
        let cell_area = cell.area_m2();
        sum_all_areas += cell_area;

        // Pass the tiny array
        let sort_key = compute_sort_key(cell, res, seed_u64, base_cells, &mut last_seen_perm)?;

        // Compute the effective key used for heap comparison.
        //
        // Equal probability: sort_key as f64 — identical ordering to the
        //   original u64 sort_key for all keys < 2^53, which holds up to
        //   res 15 (max key = 2^52, exactly representable).
        //
        // Proportional: sort_key as f64 / cell_area — cells with larger
        //   area get smaller effective keys and are thus more competitive.
        //   The constant factor (sum_area / n) cancels in all comparisons
        //   and need not be computed here.
        let effective_key: f64 = if area_correction {
            sort_key as f64 / cell_area
        } else {
            sort_key as f64
        };

        // to_bits() is safe here: effective_key is always positive and
        // finite (sort_key >= 0, cell_area > 0). For positive finite f64,
        // bit order == value order.
        let effective_key_bits = effective_key.to_bits();

        reservoir_push(
            &mut heap,
            target_n,
            effective_key_bits, // drives eviction
            sort_key,           // preserved for display
            u64::from(cell),
            cell_area,
        );
    }

    // Validate sample size -----------------------------------------------
    // If the heap has fewer items than target_n,
    // the study area simply doesn't contain enough cells.
    if heap.len() < target_n {
        return Err(Error::Other(format!(
            "Target n ({}) exceeds the total number of cells covering the geometry ({}). Try a higher resolution or larger area.",
            target_n,
            total_cells
        )));
    }

    Ok((heap.into_sorted_vec(), total_cells, sum_all_areas))
}

fn build_result_df(
    sorted: Vec<HeapEntry>,
    total_cells: usize,
    sum_all_areas: f64,
    n: i32,
    res: i32,
    containment: &str,
    area_correction: bool,
) -> extendr_api::Result<Robj> {
    let effective_keys: Vec<f64> = sorted
        .iter()
        .map(|e| f64::from_bits(e.effective_key_bits))
        .collect();

    let sort_keys: Vec<String> = sorted
        .iter()
        .map(|e| format!("{:016x}", e.raw_sort_key))
        .collect();

    let cells: Vec<String> = sorted.iter().map(|e| format!("{:x}", e.h3_index)).collect();

    let areas_m2: Vec<f64> = sorted.iter().map(|e| f64::from_bits(e.area_bits)).collect();

    let n_f64 = n as f64;

    let ip: Vec<f64> = if area_correction {
        sorted
            .iter()
            .map(|e| {
                let area = f64::from_bits(e.area_bits);
                n_f64 * area / sum_all_areas
            })
            .collect()
    } else {
        vec![n_f64 / total_cells as f64; n as usize]
    };

    let df = data_frame!(
        cell = cells,
        sort_key = sort_keys,
        effective_key = effective_keys,
        area_m2 = areas_m2,
        ip = ip
    );

    let mut df_robj: Robj = df;
    df_robj.set_attrib("n_cells", total_cells as i32)?;
    df_robj.set_attrib("sum_area_m2", sum_all_areas)?;
    df_robj.set_attrib("n", n)?;
    df_robj.set_attrib("resolution", res)?;
    df_robj.set_attrib("containment", containment)?;
    df_robj.set_attrib("area_correction", area_correction)?;

    Ok(df_robj)
}

// ---------------------------------------------------------------------------
// Sort key computation
// ---------------------------------------------------------------------------

/// Compute the GRTS sort key for a single H3 cell.
///
/// The sort key is a packed u64 built top-down through the H3 hierarchy:
///
/// ```text
/// bits 63..52  unused
/// bits 51..45  shuffled base cell  (7 bits, SORT_KEY_BASE_BITS)
/// bits 44..42  shuffled digit at res 1  (3 bits each)
/// bits 41..39  shuffled digit at res 2
/// ...
/// bits (44 - 3*(res-1))..(42 - 3*(res-1))  shuffled digit at target res
/// ```
///
/// Each digit is the child index (0–6) within its parent, remapped through a
/// per-parent pseudo-random permutation. Pentagon parents only have 6 valid
/// children (digits 0–5), but the same 7-element permutation is also used for
/// these (harmless).
fn compute_sort_key(
    cell: CellIndex,
    res: i32,
    seed_u64: u64,
    base_cells: &[u64; H3_NUM_BASE_CELLS],
    last_seen_perm: &mut [(u64, [u64; 7]); 16],
) -> extendr_api::Result<u64> {
    let h3_index = u64::from(cell);

    // -- Base cell (most significant part of the key) ----------------------
    let base_cell_raw = ((h3_index >> H3_BASE_CELL_SHIFT) & H3_BASE_CELL_MASK) as usize;
    let remapped_base = if base_cell_raw <= H3_MAX_BASE_CELL {
        base_cells[base_cell_raw]
    } else {
        base_cell_raw as u64
    };
    let mut sort_key: u64 = remapped_base;

    // -- Resolution digits (one per level from 1 to res) -------------------
    for r in 1..=res {
        // Use the h3o API to get the parent cell at resolution r-1, which is
        // more robust than manual bit-masking.
        let parent_res = Resolution::try_from((r - 1) as u8)
            .map_err(|_| Error::Other(format!("Invalid parent resolution: {}", r - 1)))?;
        let parent_cell = cell
            .parent(parent_res)
            .ok_or_else(|| Error::Other(format!("Could not get parent at res {}", r - 1)))?;
        let parent_id = u64::from(parent_cell);
        let r_idx = r as usize;

        // Cache hit logic: only shuffle if the parent is different from the last cell we saw
        if last_seen_perm[r_idx].0 != parent_id {
            let mut rng = Pcg64Mcg::seed_from_u64(parent_id.wrapping_add(seed_u64));
            let mut perm = [0u64, 1, 2, 3, 4, 5, 6];
            perm.shuffle(&mut rng);
            last_seen_perm[r_idx] = (parent_id, perm);
        }

        let perm = last_seen_perm[r_idx].1;

        let h3_offset = H3_BASE_CELL_SHIFT
            .checked_sub(3 * r as u32)
            .ok_or_else(|| Error::Other("Resolution digit offset underflow".into()))?;
        let digit = ((h3_index >> h3_offset) & H3_DIGIT_MASK) as usize;

        if digit < 7 {
            let remapped_digit = perm[digit];
            // Pack the remapped digit into the sort key.
            // Layout: 7 base-cell bits, then 3 bits per resolution level.
            let key_offset = SORT_KEY_BASE_BITS + SORT_KEY_BITS_PER_RES * (r as u32 - 1);
            sort_key |= remapped_digit << key_offset;
        }
    }

    // Sanity check: at res 15 we need 7 + 15*3 = 52 bits, well within u64.
    debug_assert!(
        SORT_KEY_BASE_BITS + SORT_KEY_BITS_PER_RES * res as u32 <= 64,
        "Sort key overflows u64 at resolution {}",
        res
    );

    Ok(sort_key)
}

// ---------------------------------------------------------------------------
// Heap reservoir helper
// ---------------------------------------------------------------------------

/// Push a [`HeapEntry`] into the reservoir heap, evicting the maximum element
/// if the heap is already full and the new key is smaller.
/// Cell area is carried through so weights can be computed after the loop
/// without re-querying h3o.
#[inline]
fn reservoir_push(
    heap: &mut BinaryHeap<HeapEntry>,
    target_n: usize,
    effective_key_bits: u64,
    raw_sort_key: u64,
    h3_index: u64,
    cell_area: f64,
) {
    let entry = HeapEntry::new(
        effective_key_bits,
        raw_sort_key,
        h3_index,
        cell_area.to_bits(),
    );
    if heap.len() < target_n {
        heap.push(entry);
    } else if let Some(&max) = heap.peek() {
        if effective_key_bits < max.effective_key_bits {
            heap.push(entry);
            heap.pop();
        }
    }
}
// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------

/// Decode a WKB byte slice into a `geo::Geometry`.
fn decode_wkb(wkb_bytes: &[u8]) -> extendr_api::Result<Geometry> {
    // Create a zero-copy WKB reader from the raw byte slice.
    let wkb_obj =
        Wkb::try_new(wkb_bytes).map_err(|e| Error::Other(format!("WKB parse failed: {}", e)))?;
    // Use the geo-traits extension to convert into an owned geo::Geometry
    Ok(wkb_obj.to_geometry())
}
/// Parse the containment mode string into an `h3o` `ContainmentMode`.
fn parse_containment(containment: &str) -> extendr_api::Result<ContainmentMode> {
    match containment {
        "intersect" => Ok(ContainmentMode::IntersectsBoundary),
        "centroid" => Ok(ContainmentMode::ContainsCentroid),
        "boundary" => Ok(ContainmentMode::ContainsBoundary),
        "covers" => Ok(ContainmentMode::Covers),
        other => Err(Error::Other(format!(
            "Unknown containment mode: '{}'. \
             Expected one of: intersect, centroid, boundary, covers.",
            other
        ))),
    }
}

/// Feed a `geo::Geometry` into an h3o tiler, recursing into collections.
/// Returns an error (rather than silently skipping) if unexpected geometry
/// types such as LineStrings are encountered, as these likely indicate an
/// upstream data pipeline issue.
fn feed_geometry(tiler: &mut h3o::geom::Tiler, geometry: Geometry) -> extendr_api::Result<()> {
    match geometry {
        Geometry::Polygon(p) => {
            tiler
                .add(p)
                .map_err(|e| Error::Other(format!("Invalid polygon: {}", e)))?;
        }
        Geometry::MultiPolygon(mp) => {
            for p in mp {
                tiler
                    .add(p)
                    .map_err(|e| Error::Other(format!("Invalid polygon in multipolygon: {}", e)))?;
            }
        }
        Geometry::GeometryCollection(gc) => {
            for geom in gc {
                match geom {
                    Geometry::Polygon(_) | Geometry::MultiPolygon(_) => {
                        feed_geometry(tiler, geom)?;
                    }
                    other => {
                        return Err(Error::Other(format!(
                            "Unsupported geometry type inside GeometryCollection: '{}'. \
                             Only Polygon and MultiPolygon are supported.",
                            geometry_type_name(&other)
                        )));
                    }
                }
            }
        }
        other => {
            return Err(Error::Other(format!(
                "Unsupported top-level geometry type: '{}'. \
                 Must be Polygon, MultiPolygon, or GeometryCollection.",
                geometry_type_name(&other)
            )));
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Base cell shuffle
// ---------------------------------------------------------------------------

/// Return a shuffled mapping of base cell indices 0–121.
/// This forms the most significant level of the GRTS address, ensuring that
/// sample coverage is balanced across the entire Earth surface.
fn shuffled_base_cells(seed_u64: u64) -> [u64; H3_NUM_BASE_CELLS] {
    let mut cells: [u64; H3_NUM_BASE_CELLS] = std::array::from_fn(|i| i as u64);
    let mut rng = Pcg64Mcg::seed_from_u64(seed_u64);
    cells.shuffle(&mut rng);
    cells
}

/// Return a human-readable name for a `geo::Geometry` variant.
fn geometry_type_name<T: geo::CoordNum>(geom: &Geometry<T>) -> &'static str {
    match geom {
        Geometry::Point(_) => "Point",
        Geometry::Line(_) => "Line",
        Geometry::LineString(_) => "LineString",
        Geometry::Polygon(_) => "Polygon",
        Geometry::MultiPoint(_) => "MultiPoint",
        Geometry::MultiLineString(_) => "MultiLineString",
        Geometry::MultiPolygon(_) => "MultiPolygon",
        Geometry::GeometryCollection(_) => "GeometryCollection",
        Geometry::Rect(_) => "Rect",
        Geometry::Triangle(_) => "Triangle",
    }
}

/// Parse, validate, and dedupe a raw little-endian u64 cell vector.
/// Returns the sorted unique cells and the number of duplicates removed.
fn parse_cells(cells_bytes: &[u8]) -> extendr_api::Result<(Vec<CellIndex>, usize)> {
    if cells_bytes.is_empty() || !cells_bytes.len().is_multiple_of(8) {
        return Err(Error::Other(
            "Provided cells vector is empty or its length is not a multiple of 8 bytes.".into(),
        ));
    }

    let total = cells_bytes.len() / 8;
    let mut cells: Vec<CellIndex> = Vec::with_capacity(total);
    let mut n_invalid = 0usize;
    let mut first_bad: Option<u64> = None;

    for chunk in cells_bytes.chunks_exact(8) {
        let raw = u64::from_le_bytes(chunk.try_into().unwrap());
        match CellIndex::try_from(raw) {
            Ok(c) => cells.push(c),
            Err(_) => {
                n_invalid += 1;
                first_bad.get_or_insert(raw);
            }
        }
    }
    if n_invalid > 0 {
        return Err(Error::Other(format!(
            "{n_invalid} of {total} values are not valid H3 cell indices (first: {:x}).",
            first_bad.unwrap_or(0)
        )));
    }

    // Non-empty here: the byte vector was non-empty and every value parsed.
    let res = cells[0].resolution();
    if let Some(other) = cells.iter().find(|c| c.resolution() != res) {
        return Err(Error::Other(format!(
            "All cells must share one resolution; found {} and {}.",
            u8::from(res),
            u8::from(other.resolution())
        )));
    }

    cells.sort_unstable();
    cells.dedup();
    let n_dup = total - cells.len();
    Ok((cells, n_dup))
}

#[cfg(test)]
mod tests {
    use super::*;
    use h3o::LatLng;

    fn cell(lat: f64, lng: f64, res: Resolution) -> u64 {
        u64::from(LatLng::new(lat, lng).unwrap().to_cell(res))
    }
    fn bytes(v: &[u64]) -> Vec<u8> {
        v.iter().flat_map(|c| c.to_le_bytes()).collect()
    }

    #[test]
    fn parse_cells_dedupes_and_sorts() {
        let a = cell(50.0, 4.0, Resolution::Five);
        let b = cell(51.0, 5.0, Resolution::Five);
        let (cells, dup) = parse_cells(&bytes(&[b, a, b])).unwrap();
        assert_eq!(dup, 1);
        assert_eq!(cells.len(), 2);
        assert!(cells.windows(2).all(|w| w[0] < w[1]));
    }

    #[test]
    fn parse_cells_rejects_invalid_index() {
        let a = cell(50.0, 4.0, Resolution::Five);
        assert!(parse_cells(&bytes(&[a, 0])).is_err());
    }

    #[test]
    fn parse_cells_rejects_mixed_resolution() {
        let a = cell(50.0, 4.0, Resolution::Five);
        let b = cell(50.0, 4.0, Resolution::Six);
        assert!(parse_cells(&bytes(&[a, b])).is_err());
    }

    #[test]
    fn parse_cells_rejects_empty_and_misaligned() {
        assert!(parse_cells(&[]).is_err());
        assert!(parse_cells(&[1, 2, 3]).is_err());
    }

    #[test]
    fn sample_is_independent_of_input_order() {
        let cs: Vec<u64> = (0..40)
            .map(|i| cell(48.0 + i as f64 * 0.1, 4.0, Resolution::Seven))
            .collect();
        let mut rev = cs.clone();
        rev.reverse();
        let run = |v: &[u64]| {
            let (cells, _) = parse_cells(&bytes(v)).unwrap();
            let base = shuffled_base_cells(7);
            grts_core(cells.iter().copied(), 10, 7, 7, false, &base)
                .unwrap()
                .0
                .iter()
                .map(|e| e.h3_index)
                .collect::<Vec<_>>()
        };
        assert_eq!(run(&cs), run(&rev));
    }
}

// ---------------------------------------------------------------------------
// extendr module registration for the `grts_h3` submodule
// ---------------------------------------------------------------------------

extendr_module! {
    mod grts_h3;
    fn h3_coverage;
    fn grts_sample_from_wkb;
    fn grts_sample_from_cells;
}
