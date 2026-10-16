// src/sobol.rs — Unified Spatially Balanced Sampling (Sobol Engine)
//
// Supports two complementary execution strategies:
// 1. "acceptance sampling" (BAS): Generates candidate points across the bounding
//    box and accepts those inside the polygon. Best for compact/dense shapes.
// 2. "routing" (Continuous GRTS): Maps the 1D Morton-interleaved sequence
//    down an exact-area quadtree. Zero rejection, ideal for sparse/narrow domains.

use extendr_api::prelude::*;
use geo::algorithm::bounding_rect::BoundingRect;
use geo::algorithm::contains::Contains;
use geo::{Coord, Geometry, LineString, Point, Rect};
use geo_traits::to_geo::ToGeoGeometry;
use rand::{RngCore, SeedableRng};
use rand_pcg::Pcg64Mcg;
use wkb::reader::Wkb;

const LON_MIN: f64 = -180.0;
const LON_MAX: f64 = 180.0;
const LEAF_DEPTH: u32 = 28;
const DENSIFY_DEG: f64 = 0.01;

// ---------------------------------------------------------------------------
// Coordinate Transforms (Equal-Area on Sphere)
// ---------------------------------------------------------------------------

#[inline]
fn lat_to_s(lat_deg: f64) -> f64 {
    lat_deg.to_radians().sin()
}

#[inline]
fn s_to_lat(s: f64) -> f64 {
    s.clamp(-1.0, 1.0).asin().to_degrees()
}

// ---------------------------------------------------------------------------
// Sobol Generation & Owen Scrambling
// ---------------------------------------------------------------------------

const DIRECTION_VECTORS: [[u32; 32]; 5] = [
    // dim 0: Van der Corput
    [
        0x80000000, 0x40000000, 0x20000000, 0x10000000, 0x08000000, 0x04000000, 0x02000000,
        0x01000000, 0x00800000, 0x00400000, 0x00200000, 0x00100000, 0x00080000, 0x00040000,
        0x00020000, 0x00010000, 0x00008000, 0x00004000, 0x00002000, 0x00001000, 0x00000800,
        0x00000400, 0x00000200, 0x00000100, 0x00000080, 0x00000040, 0x00000020, 0x00000010,
        0x00000008, 0x00000004, 0x00000002, 0x00000001,
    ],
    // dim 1: Joe-Kuo d=2 (s=1, a=0, m=[1])
    [
        0x80000000, 0xc0000000, 0xa0000000, 0xf0000000, 0x88000000, 0xcc000000, 0xaa000000,
        0xff000000, 0x80800000, 0xc0c00000, 0xa0a00000, 0xf0f00000, 0x88880000, 0xcccc0000,
        0xaaaa0000, 0xffff0000, 0x80008000, 0xc000c000, 0xa000a000, 0xf000f000, 0x88008800,
        0xcc00cc00, 0xaa00aa00, 0xff00ff00, 0x80808080, 0xc0c0c0c0, 0xa0a0a0a0, 0xf0f0f0f0,
        0x88888888, 0xcccccccc, 0xaaaaaaaa, 0xffffffff,
    ],
    // dim 2: Joe-Kuo d=3 (s=2, a=1, m=[1,3])
    [
        0x80000000, 0xc0000000, 0x60000000, 0x90000000, 0xe8000000, 0x5c000000, 0x8e000000,
        0xc5000000, 0x68800000, 0x9cc00000, 0xee600000, 0x55900000, 0x80680000, 0xc09c0000,
        0x60ee0000, 0x90550000, 0xe8808000, 0x5cc0c000, 0x8e606000, 0xc5909000, 0x6868e800,
        0x9c9c5c00, 0xeeee8e00, 0x5555c500, 0x8000e880, 0xc0005cc0, 0x60008e60, 0x9000c590,
        0xe8006868, 0x5c009c9c, 0x8e00eeee, 0xc5005555,
    ],
    // dim 3: Joe-Kuo d=4 (s=3, a=1, m=[1,3,1])
    [
        0x80000000, 0xc0000000, 0x20000000, 0x50000000, 0xf8000000, 0x74000000, 0xa2000000,
        0x93000000, 0xd8800000, 0x25400000, 0x59e00000, 0xe6d00000, 0x78080000, 0xb40c0000,
        0x82020000, 0xc3050000, 0x208f8000, 0x51474000, 0xfbea2000, 0x75d93000, 0xa0858800,
        0x914e5400, 0xdbe79e00, 0x25db6d00, 0x58800080, 0xe54000c0, 0x79e00020, 0xb6d00050,
        0x800800f8, 0xc00c0074, 0x200200a2, 0x50050093,
    ],
    // dim 4: Joe-Kuo d=5 (s=3, a=2, m=[1,1,1])
    [
        0x80000000, 0x40000000, 0x20000000, 0xb0000000, 0xf8000000, 0xdc000000, 0x7a000000,
        0x9d000000, 0x5a800000, 0x2fc00000, 0xa1600000, 0xf0b00000, 0xda880000, 0x6fc40000,
        0x81620000, 0x40bb0000, 0x22878000, 0xb3c9c000, 0xfb65a000, 0xddb2d000, 0x78022800,
        0x9c0b3c00, 0x5a0fb600, 0x2d0ddb00, 0xa2878080, 0xf3c9c040, 0xdb65a020, 0x6db2d0b0,
        0x800228f8, 0x400b3cdc, 0x200fb67a, 0xb00ddb9d,
    ],
];

#[inline]
fn sobol_raw(mut index: u32, dimension: usize) -> u32 {
    let v = &DIRECTION_VECTORS[dimension];
    let mut result = 0u32;
    let mut i = 0usize;
    while index > 0 {
        if index & 1 != 0 {
            result ^= v[i];
        }
        index >>= 1;
        i += 1;
    }
    result
}

#[inline]
fn lk_hash(mut x: u32, seed: u32) -> u32 {
    x ^= x.wrapping_mul(0x3d20_adea);
    x = x.wrapping_add(seed);
    x = x.wrapping_mul((seed >> 16) | 1);
    x ^= x.wrapping_mul(0x0552_6c56);
    x ^= x.wrapping_mul(0x53a2_2864);
    x
}

const N_DIMS: usize = DIRECTION_VECTORS.len();

#[inline]
fn dim_seeds(global_seed: u64) -> [u32; N_DIMS] {
    let mut rng = Pcg64Mcg::seed_from_u64(global_seed);
    std::array::from_fn(|_| (rng.next_u64() >> 32) as u32)
}

#[inline]
fn owen_scramble(x: u32, seed: u32) -> u32 {
    let x = x.reverse_bits();
    let x = lk_hash(x, seed);
    x.reverse_bits()
}

#[inline]
fn expand_bits_32(v: u32) -> u64 {
    let mut x = v as u64;
    x = (x | (x << 16)) & 0x0000FFFF0000FFFF;
    x = (x | (x << 8)) & 0x00FF00FF00FF00FF;
    x = (x | (x << 4)) & 0x0F0F0F0F0F0F0F0F;
    x = (x | (x << 2)) & 0x3333333333333333;
    x = (x | (x << 1)) & 0x5555555555555555;
    x
}

#[inline]
fn morton_interleave(x: u32, y: u32) -> u64 {
    expand_bits_32(x) | (expand_bits_32(y) << 1)
}

// ---------------------------------------------------------------------------
// Geometry Processing & Clipping
// ---------------------------------------------------------------------------

#[derive(Clone)]
struct SignedRing {
    coords: Vec<Coord>,
    sign: f64,
}

fn densify_and_transform_ring(ring: &LineString) -> Vec<Coord> {
    let coords = ring.0.as_slice();
    if coords.len() < 2 {
        return coords
            .iter()
            .map(|c| Coord {
                x: c.x,
                y: lat_to_s(c.y),
            })
            .collect();
    }

    let mut out: Vec<Coord> = Vec::with_capacity(coords.len() * 4);

    for window in coords.windows(2) {
        let (p0, p1) = (window[0], window[1]);
        let max_span = (p1.x - p0.x).abs().max((p1.y - p0.y).abs());
        let n_segs = if max_span > DENSIFY_DEG {
            (max_span / DENSIFY_DEG).ceil() as usize
        } else {
            1
        };

        for k in 0..n_segs {
            let t = k as f64 / n_segs as f64;
            let lon = p0.x + t * (p1.x - p0.x);
            let lat = p0.y + t * (p1.y - p0.y);
            out.push(Coord {
                x: lon,
                y: lat_to_s(lat),
            });
        }
    }

    let last = coords.last().unwrap();
    out.push(Coord {
        x: last.x,
        y: lat_to_s(last.y),
    });
    out
}

fn collect_rings(geometry: &Geometry, rings: &mut Vec<SignedRing>) -> extendr_api::Result<()> {
    match geometry {
        Geometry::Polygon(p) => {
            rings.push(SignedRing {
                coords: densify_and_transform_ring(p.exterior()),
                sign: 1.0,
            });
            for hole in p.interiors() {
                rings.push(SignedRing {
                    coords: densify_and_transform_ring(hole),
                    sign: -1.0,
                });
            }
        }
        Geometry::MultiPolygon(mp) => {
            for p in &mp.0 {
                collect_rings(&Geometry::Polygon(p.clone()), rings)?;
            }
        }
        Geometry::GeometryCollection(gc) => {
            for g in gc.iter() {
                collect_rings(g, rings)?;
            }
        }
        _ => {
            return Err(Error::Other(
                "Geometry must be POLYGON or MULTIPOLYGON".into(),
            ))
        }
    }
    Ok(())
}

fn extract_signed_rings(geometry: &Geometry) -> extendr_api::Result<Vec<SignedRing>> {
    let mut rings = Vec::new();
    collect_rings(geometry, &mut rings)?;
    Ok(rings)
}

#[inline]
fn inside(p: &Coord, edge: u8, val: f64) -> bool {
    match edge {
        0 => p.x >= val,
        1 => p.x <= val,
        2 => p.y >= val,
        3 => p.y <= val,
        _ => unreachable!(),
    }
}

#[inline]
fn intersect(p1: &Coord, p2: &Coord, edge: u8, val: f64) -> Coord {
    match edge {
        0 | 1 => Coord {
            x: val,
            y: p1.y + (val - p1.x) / (p2.x - p1.x) * (p2.y - p1.y),
        },
        2 | 3 => Coord {
            x: p1.x + (val - p1.y) / (p2.y - p1.y) * (p2.x - p1.x),
            y: val,
        },
        _ => unreachable!(),
    }
}

fn clip_edge(input: &[Coord], edge: u8, val: f64) -> Vec<Coord> {
    let mut out = Vec::with_capacity(input.len());
    if input.is_empty() {
        return out;
    }
    let mut prev = *input.last().unwrap();
    let mut prev_in = inside(&prev, edge, val);

    for &curr in input {
        let curr_in = inside(&curr, edge, val);
        if curr_in != prev_in {
            out.push(intersect(&prev, &curr, edge, val));
        }
        if curr_in {
            out.push(curr);
        }
        prev = curr;
        prev_in = curr_in;
    }
    out
}

fn clip_sh(input: &[Coord], xmin: f64, xmax: f64, ymin: f64, ymax: f64) -> Vec<Coord> {
    let p = clip_edge(input, 0, xmin);
    let p = clip_edge(&p, 1, xmax);
    let p = clip_edge(&p, 2, ymin);
    clip_edge(&p, 3, ymax)
}

fn shoelace_abs(coords: &[Coord]) -> f64 {
    if coords.len() < 3 {
        return 0.0;
    }
    let mut area = 0.0;
    for i in 0..coords.len() {
        let p1 = coords[i];
        let p2 = coords[(i + 1) % coords.len()];
        area += p1.x * p2.y - p2.x * p1.y;
    }
    (area / 2.0).abs()
}

fn calculate_clipped_area(rings: &[SignedRing], xmin: f64, xmax: f64, ymin: f64, ymax: f64) -> f64 {
    let mut total_area = 0.0;
    for ring in rings {
        let clipped = clip_sh(&ring.coords, xmin, xmax, ymin, ymax);
        if clipped.len() >= 3 {
            total_area += shoelace_abs(&clipped) * ring.sign;
        }
    }
    total_area.max(0.0)
}

fn bounding_rect_geom(geometry: &Geometry) -> Option<Rect> {
    match geometry {
        Geometry::Polygon(p) => p.bounding_rect(),
        Geometry::MultiPolygon(m) => m.bounding_rect(),
        Geometry::GeometryCollection(gc) => {
            gc.iter().filter_map(bounding_rect_geom).reduce(|a, b| {
                Rect::new(
                    Coord {
                        x: a.min().x.min(b.min().x),
                        y: a.min().y.min(b.min().y),
                    },
                    Coord {
                        x: a.max().x.max(b.max().x),
                        y: a.max().y.max(b.max().y),
                    },
                )
            })
        }
        _ => None,
    }
}

fn geometry_contains(geometry: &Geometry, point: &Point) -> bool {
    match geometry {
        Geometry::Polygon(p) => p.contains(point),
        Geometry::MultiPolygon(mp) => mp.contains(point),
        Geometry::GeometryCollection(gc) => gc.iter().any(|g| geometry_contains(g, point)),
        _ => false,
    }
}

// ---------------------------------------------------------------------------
// Routing Engine
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct CellBounds {
    xmin: f64,
    xmax: f64,
    ymin: f64,
    ymax: f64,
}

#[derive(Debug, Clone)]
struct TargetPoint {
    id: usize,
    area_target: f64,
    lon: f64,
    lat: f64,
}

fn route(
    pts: &mut [TargetPoint],
    rings: &[SignedRing],
    cell: CellBounds,
    depth: u32,
    base_area: f64,
) {
    if pts.is_empty() {
        return;
    }

    if depth == LEAF_DEPTH {
        let cx = (cell.xmin + cell.xmax) * 0.5;
        let cy = (cell.ymin + cell.ymax) * 0.5;
        for p in pts.iter_mut() {
            p.lon = cx;
            p.lat = cy;
        }
        return;
    }

    let mut local_rings = Vec::new();
    for ring in rings {
        let clipped = clip_sh(&ring.coords, cell.xmin, cell.xmax, cell.ymin, cell.ymax);
        if clipped.len() >= 3 {
            local_rings.push(SignedRing {
                coords: clipped,
                sign: ring.sign,
            });
        }
    }

    if local_rings.is_empty() {
        return;
    }

    let xmid = (cell.xmin + cell.xmax) * 0.5;
    let ymid = (cell.ymin + cell.ymax) * 0.5;

    let a0 = calculate_clipped_area(&local_rings, cell.xmin, xmid, cell.ymin, ymid);
    let a1 = calculate_clipped_area(&local_rings, xmid, cell.xmax, cell.ymin, ymid);
    let a2 = calculate_clipped_area(&local_rings, cell.xmin, xmid, ymid, cell.ymax);

    let t0 = base_area + a0;
    let t1 = t0 + a1;
    let t2 = t1 + a2;

    let mut i0 = 0;
    while i0 < pts.len() && pts[i0].area_target < t0 {
        i0 += 1;
    }
    let (pts0, rest) = pts.split_at_mut(i0);

    let mut i1 = 0;
    while i1 < rest.len() && rest[i1].area_target < t1 {
        i1 += 1;
    }
    let (pts1, rest) = rest.split_at_mut(i1);

    let mut i2 = 0;
    while i2 < rest.len() && rest[i2].area_target < t2 {
        i2 += 1;
    }
    let (pts2, pts3) = rest.split_at_mut(i2);

    route(
        pts0,
        &local_rings,
        CellBounds {
            xmax: xmid,
            ymax: ymid,
            ..cell
        },
        depth + 1,
        base_area,
    );
    route(
        pts1,
        &local_rings,
        CellBounds {
            xmin: xmid,
            ymax: ymid,
            ..cell
        },
        depth + 1,
        t0,
    );
    route(
        pts2,
        &local_rings,
        CellBounds {
            xmax: xmid,
            ymin: ymid,
            ..cell
        },
        depth + 1,
        t1,
    );
    route(
        pts3,
        &local_rings,
        CellBounds {
            xmin: xmid,
            ymin: ymid,
            ..cell
        },
        depth + 1,
        t2,
    );
}

// ---------------------------------------------------------------------------
// Unified Public Entry Point
// ---------------------------------------------------------------------------

#[extendr]
pub fn sobol_sample_from_wkb(
    wkb_bytes: &[u8],
    n: i32,
    seed: f64,
    master_bbox: Robj,
    algorithm: &str,
) -> extendr_api::Result<Robj> {
    let target_n = n as usize;
    if target_n == 0 {
        return Err(Error::Other("n must be >= 1".into()));
    }

    // Preferred XOR-folding seed representation
    if !(seed.is_finite() && seed >= 0.0 && seed.fract() == 0.0 && seed < 9_007_199_254_740_992.0) {
        return Err(Error::Other(
            "seed must be a whole number in [0, 2^53)".into(),
        ));
    }
    let seeds = dim_seeds(seed as u64);
    let (seed_x, seed_y) = (seeds[0], seeds[1]);

    // Decode WKB
    let wkb_obj =
        Wkb::try_new(wkb_bytes).map_err(|e| Error::Other(format!("WKB parse failed: {e}")))?;
    let geom = wkb_obj.to_geometry();

    let bbox_geom = bounding_rect_geom(&geom)
        .ok_or_else(|| Error::Other("Could not compute bounding box".into()))?;

    let (poly_xmin, poly_ymin, poly_xmax, poly_ymax) = (
        bbox_geom.min().x,
        bbox_geom.min().y,
        bbox_geom.max().x,
        bbox_geom.max().y,
    );

    // Sequence Domain resolution
    let (seq_xmin, seq_ymin, seq_xmax, seq_ymax) = if !master_bbox.is_null() {
        let v = master_bbox.as_real_vector().ok_or_else(|| {
            Error::Other("master_bbox must be a numeric vector of length 4".into())
        })?;
        if v.len() != 4 {
            return Err(Error::Other(
                "master_bbox must contain exactly 4 elements: c(xmin, ymin, xmax, ymax)".into(),
            ));
        }
        (v[0], v[1], v[2], v[3])
    } else {
        (poly_xmin, poly_ymin, poly_xmax, poly_ymax)
    };

    if !(seq_xmin < seq_xmax && seq_ymin < seq_ymax)
        || seq_xmin < LON_MIN
        || seq_xmax > LON_MAX
        || seq_ymin < -90.0
        || seq_ymax > 90.0
    {
        return Err(Error::Other(
            "master_bbox is invalid or outside WGS-84 bounds".into(),
        ));
    }

    // Check intersection between sequence box and polygon box
    if seq_xmax < poly_xmin || seq_xmin > poly_xmax || seq_ymax < poly_ymin || seq_ymin > poly_ymax
    {
        return Err(Error::Other(
            "master_bbox and geometry bounding box do not intersect".into(),
        ));
    }

    let root_cell = CellBounds {
        xmin: seq_xmin,
        xmax: seq_xmax,
        ymin: lat_to_s(seq_ymin),
        ymax: lat_to_s(seq_ymax),
    };

    let (lons, lats, indices, sobol_scanned, acceptance_rate) = match algorithm {
        "acceptance sampling" | "acceptance" => {
            let mut lons = Vec::with_capacity(target_n);
            let mut lats = Vec::with_capacity(target_n);
            let mut indices = Vec::with_capacity(target_n);

            let lon_range = root_cell.xmax - root_cell.xmin;
            let s_range = root_cell.ymax - root_cell.ymin;

            let mut scan_idx: u32 = 0;
            while lons.len() < target_n {
                let rx = sobol_raw(scan_idx, 0);
                let ry = sobol_raw(scan_idx, 1);
                let sx = owen_scramble(rx, seed_x);
                let sy = owen_scramble(ry, seed_y);

                // division by 2^32 for double precision, floating point representation
                let u = sx as f64 * (1.0 / 4_294_967_296.0);
                let v = sy as f64 * (1.0 / 4_294_967_296.0);

                let cand_lon = root_cell.xmin + u * lon_range;
                let cand_lat = s_to_lat(root_cell.ymin + v * s_range);

                scan_idx = scan_idx.checked_add(1).ok_or_else(|| {
                    Error::Other("Sobol index overflow (2^32): geometry is too sparse relative to bounding box".into())
                })?;

                // Bounding box pre-filter
                if cand_lon < poly_xmin
                    || cand_lon > poly_xmax
                    || cand_lat < poly_ymin
                    || cand_lat > poly_ymax
                {
                    continue;
                }

                if geometry_contains(
                    &geom,
                    &Point(Coord {
                        x: cand_lon,
                        y: cand_lat,
                    }),
                ) {
                    lons.push(cand_lon);
                    lats.push(cand_lat);
                    indices.push((scan_idx - 1) as i32);
                }
            }

            let scanned = scan_idx as f64;
            let acc_rate = target_n as f64 / scanned;
            (lons, lats, indices, scanned, acc_rate)
        }
        "routing" => {
            let signed_rings = extract_signed_rings(&geom)?;
            let total_area = calculate_clipped_area(
                &signed_rings,
                root_cell.xmin,
                root_cell.xmax,
                root_cell.ymin,
                root_cell.ymax,
            );
            if total_area <= 0.0 {
                return Err(Error::Other(
                    "Geometry has zero computable area inside master_bbox".into(),
                ));
            }

            let mut targets = Vec::with_capacity(target_n);
            for i in 0..target_n {
                let rx = sobol_raw(i as u32, 0);
                let ry = sobol_raw(i as u32, 1);
                let sx = owen_scramble(rx, seed_x);
                let sy = owen_scramble(ry, seed_y);
                let m64 = morton_interleave(sx, sy);

                let u = (m64 >> 11) as f64 / (1u64 << 53) as f64;
                targets.push(TargetPoint {
                    id: i,
                    area_target: u * total_area,
                    lon: 0.0,
                    lat: 0.0,
                });
            }

            targets.sort_by(|a, b| a.area_target.partial_cmp(&b.area_target).unwrap());
            route(&mut targets, &signed_rings, root_cell, 1, 0.0);
            targets.sort_by_key(|t| t.id);

            let lons: Vec<f64> = targets.iter().map(|t| t.lon).collect();
            let lats: Vec<f64> = targets.iter().map(|t| s_to_lat(t.lat)).collect();
            let indices: Vec<i32> = targets.iter().map(|t| t.id as i32).collect();

            (lons, lats, indices, target_n as f64, 1.0)
        }
        other => {
            return Err(Error::Other(format!(
                "Unknown algorithm: '{other}'. Expected 'acceptance sampling' or 'routing'."
            )));
        }
    };

    let df = data_frame!(lon = lons, lat = lats, rank = indices);
    let mut df_robj: Robj = df;
    df_robj.set_attrib("seed", seed)?;
    df_robj.set_attrib("n", n)?;
    df_robj.set_attrib("algorithm", algorithm)?;
    df_robj.set_attrib("bbox", r!([poly_xmin, poly_ymin, poly_xmax, poly_ymax]))?;
    df_robj.set_attrib(
        "sequence_bbox",
        r!([seq_xmin, seq_ymin, seq_xmax, seq_ymax]),
    )?;
    df_robj.set_attrib("sobol_scanned", sobol_scanned)?;
    df_robj.set_attrib("acceptance_rate", acceptance_rate)?;

    Ok(df_robj)
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[inline]
    fn sobol_sample(index: u32, dimension: usize, seed: u32) -> f64 {
        let raw = sobol_raw(index, dimension);
        let scrambled = owen_scramble(raw, seed);
        scrambled as f64 * (1.0 / 4_294_967_296.0_f64)
    }

    // ── Direction vector correctness ─────────────────────────────────────────

    #[test]
    fn dim0_is_van_der_corput() {
        // Van der Corput base-2: 0, 1/2, 1/4, 3/4, 1/8, 5/8, 3/8, 7/8
        let expected = [0.0, 0.5, 0.25, 0.75, 0.125, 0.625, 0.375, 0.875];
        for (i, &exp) in expected.iter().enumerate() {
            let got = sobol_raw(i as u32, 0) as f64 / 4_294_967_296.0;
            assert!(
                (got - exp).abs() < 1e-12,
                "dim 0 index {}: expected {}, got {}",
                i,
                exp,
                got
            );
        }
    }

    #[test]
    fn all_dims_pairwise_distinct() {
        // Every pair of dims must produce a different raw sequence over 64 points.
        for d1 in 0..5usize {
            for d2 in (d1 + 1)..5usize {
                let s1: Vec<u32> = (0..64).map(|i| sobol_raw(i, d1)).collect();
                let s2: Vec<u32> = (0..64).map(|i| sobol_raw(i, d2)).collect();
                assert_ne!(s1, s2, "dims {} and {} produce identical sequences", d1, d2);
            }
        }
    }

    // ── Scramble properties ──────────────────────────────────────────────────

    #[test]
    fn sobol_sample_in_unit_interval() {
        for dim in 0..5usize {
            for idx in 0..1024u32 {
                let v = sobol_sample(idx, dim, 0xdeadbeef);
                assert!(
                    (0.0..1.0).contains(&v),
                    "dim {} index {}: {} out of [0,1)",
                    dim,
                    idx,
                    v
                );
            }
        }
    }

    #[test]
    fn sobol_sample_deterministic() {
        for dim in 0..5usize {
            for idx in 0..64u32 {
                assert_eq!(
                    sobol_sample(idx, dim, 42),
                    sobol_sample(idx, dim, 42),
                    "dim {} index {}: not deterministic",
                    dim,
                    idx
                );
            }
        }
    }

    #[test]
    fn different_seeds_produce_different_streams() {
        for dim in 0..5usize {
            let s0: Vec<f64> = (0..32).map(|i| sobol_sample(i, dim, 0)).collect();
            let s1: Vec<f64> = (0..32).map(|i| sobol_sample(i, dim, 1)).collect();
            assert_ne!(s0, s1, "seed 0 and seed 1 gave same stream for dim {}", dim);
        }
    }

    // ── Seed independence ────────────────────────────────────────────────────
    #[test]
    fn dim_seeds_distinct_and_nonzero_for_seed_zero() {
        let s = dim_seeds(0);
        for i in 0..N_DIMS {
            assert_ne!(s[i], 0);
            for j in (i + 1)..N_DIMS {
                assert_ne!(s[i], s[j], "dims {i} and {j} share a scramble seed");
            }
        }
    }

    #[test]
    fn dim_seeds_vary_across_global_seeds() {
        for dim in 0..N_DIMS {
            assert_ne!(dim_seeds(0)[dim], dim_seeds(1)[dim]);
        }
    }

    /// Regression test for the old fold+lk_hash seeding, where small integer seeds
    /// shared the same coarse-level scramble and produced near-identical points.
    /// Old code: mean |du| ≈ 0.13; well-mixed seeds: ≈ 0.28.
    #[test]
    fn adjacent_seeds_are_decorrelated() {
        let mut total = 0.0;
        for s in 0..64u64 {
            let (a, b) = (dim_seeds(s)[0], dim_seeds(s + 1)[0]);
            total += (0..64u32)
                .map(|i| {
                    let u = |sd| owen_scramble(sobol_raw(i, 0), sd) as f64 / 4_294_967_296.0;
                    (u(a) - u(b)).abs()
                })
                .sum::<f64>()
                / 64.0;
        }
        let mean = total / 64.0;
        assert!(
            mean > 0.2,
            "adjacent seeds look correlated: mean |du| = {mean}"
        );
    }

    // ── Prefix stability ─────────────────────────────────────────────────────

    #[test]
    fn prefix_of_stream_is_stable() {
        let seed = 0xabcd_1234u32;
        let short: Vec<f64> = (0..100).map(|i| sobol_sample(i, 0, seed)).collect();
        let long: Vec<f64> = (0..200).map(|i| sobol_sample(i, 0, seed)).collect();
        assert_eq!(
            &short[..],
            &long[..100],
            "first 100 values differ between n=100 and n=200 runs"
        );
    }

    /// (s, a, initial m) for dims 1..=4, from new-joe-kuo-6.21201 (rows d=2..5).
    const JK_PARAMS: [(usize, u32, &[u32]); 4] = [
        (1, 0, &[1]),
        (2, 1, &[1, 3]),
        (3, 1, &[1, 3, 1]),
        (3, 2, &[1, 1, 1]),
    ];

    fn vectors_from_recurrence(s: usize, a: u32, m_init: &[u32]) -> [u32; 32] {
        let mut m = [0u32; 32];
        m[..s].copy_from_slice(m_init);
        for k in s..32 {
            let mut new = m[k - s] ^ (m[k - s] << s);
            for i in 1..s {
                if (a >> (s - 1 - i)) & 1 == 1 {
                    new ^= m[k - i] << i;
                }
            }
            m[k] = new;
        }
        std::array::from_fn(|k| m[k] << (31 - k))
    }

    #[test]
    fn direction_vectors_match_joe_kuo_recurrence() {
        for (d, &(s, a, m)) in JK_PARAMS.iter().enumerate() {
            assert_eq!(
                DIRECTION_VECTORS[d + 1],
                vectors_from_recurrence(s, a, m),
                "dim {} deviates from its Joe-Kuo recurrence",
                d + 1
            );
        }
    }

    #[test]
    fn every_dim_is_one_dimensionally_stratified() {
        // Any valid Sobol dim puts the first 2^10 points one per 1/1024 stratum.
        for dim in 0..N_DIMS {
            let mut seen = vec![false; 1024];
            for i in 0..1024u32 {
                let b = (sobol_raw(i, dim) >> 22) as usize;
                assert!(!seen[b], "dim {dim}: stratum {b} hit twice");
                seen[b] = true;
            }
        }
    }
}

extendr_module! {
    mod sobol;
    fn sobol_sample_from_wkb;
}
