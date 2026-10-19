#' Generate a Spatially Balanced Sample using a Sobol Sequence
#'
#' @description
#' Generates a spatially balanced sample from a `POLYGON` or `MULTIPOLYGON`.
#' Points are drawn from a 2-D Owen-scrambled Sobol sequence (Pharr et al 2023)
#' projected onto the study area via an area-equal coordinate transform,
#' so that equal fractions of the sequence correspond to equal surface areas on
#' the sphere.
#'
#' Unlike [grts_h3()], this function requires no H3 discretisation step and
#' therefore scales to very large study areas without memory pressure.
#' H3 cell indices can be obtained afterwards using [h3o::h3_from_xy()] if
#' needed.
#'
#' @param wkb Raw byte vector representing a Well-Known Binary (WKB) geometry
#'   in WGS84 (EPSG:4326). Must be a `POLYGON` or `MULTIPOLYGON`.
#' @param n Integer. Target number of sample locations to draw.
#' @param seed Numeric. Random seed (whole number `>= 0` and `< 2^53`).
#' @param master_bbox Optional numeric vector `c(xmin, ymin, xmax, ymax)` in
#'   decimal degrees (WGS84).
#'   Defines the domain of the sequence.
#'   If `NULL` (default), the bounding box of `wkb` is used.
#'   See the section "Master bounding box".
#' @param algorithm Character string specifying the sampling approach:
#'   \itemize{
#'     \item \code{"acceptance sampling"}: Balanced Acceptance Sampling (BAS) bounding-box rejection.
#'     \item \code{"routing"}: Continuous GRTS exact-area quadtree routing (zero rejection).
#'   }
#'
#' @section Master bounding box:
#' Points are only drawn from the part of the polygon that lies inside
#' `master_bbox`.
#' If the polygon extends beyond the box, that part is never
#' sampled and a warning is issued.
#' The `n` requested points are still returned, spread over the area inside the
#' box.
#' They are not an equal-probability sample of the whole polygon.
#'
#' If the polygon only touches the box without overlapping it, the call fails
#' (`"routing"`) or can scan a very long stretch of the sequence before
#' failing (`"acceptance sampling"`).
#'
#' With the same `master_bbox`, `seed` and `algorithm`, the sequence domain is
#' identical for different polygons.
#' Spatial consistency across polygons holds for `"acceptance sampling"`, but
#' not for `"routing"`, which depends on the exact area of each polygon
#' (see [grts_sobol()]).
#'
#' @return A `data.frame` with `n` rows and three columns:
#' * `lon`         – Longitude of the sample point (WGS84 decimal degrees).
#' * `lat`         – Latitude  of the sample point (WGS84 decimal degrees).
#' * `rank`        – Index into the raw Sobol stream at which this point was
#'                   accepted. Useful for diagnostics and for verifying prefix
#'                   stability across calls.
#'
#' The returned `data.frame` also carries the following attributes:
#' * `seed`          – The seed used (numeric).
#' * `n`             – Requested sample size (integer).
#' * `bbox` – Named numeric vector `c(xmin, ymin, xmax, ymax)` of the polygon
#'            bounding box in decimal degrees.
#' * `sequence_bbox` – Named numeric vector `c(xmin, ymin, xmax, ymax)` of the
#'                     geometry bounding box in decimal degrees.
#'                     Either from master bounding box or polygon bounding box.
#' * `sobol_scanned` – Total number of Sobol stream indices evaluated
#'                     (accepted + rejected). Equals the last `sobol_index` + 1.
#' * `fill_ratio`    – `n / sobol_scanned`. The fraction of bounding-box
#'                     candidates that fell inside the polygon. Low values
#'                     (< 0.05) indicate a sparse geometry relative to its
#'                     bounding box and imply higher computational cost.
#' @export
sobol_sample <- function(wkb,
                         n,
                         seed,
                         master_bbox = NULL,
                         algorithm = c("acceptance sampling", "routing")) {
  algorithm <- match.arg(algorithm)

  stopifnot(
    "Input 'wkb' must be a raw byte vector representing WKB geometry." = is.raw(wkb),
    "`n` must be a whole number, > 0, not `NA`." = is_valid_n(n),
    "`seed` must be a whole number, >= 0, and < 2^53." = is_valid_seed(seed),
    "`master_bbox` must be length 4 numeric in decimal degrees" = is_valid_bbox(master_bbox)
  )

  res_df <- sobol_sample_from_wkb(
    wkb_bytes = wkb,
    n = as.integer(n),
    seed = as.numeric(seed),
    master_bbox = master_bbox,
    algorithm = algorithm
  )

  # Attach standard names to bbox attributes
  for (attr_name in c("bbox", "sequence_bbox")) {
    bb <- attr(res_df, attr_name)
    if (!is.null(bb)) {
      names(bb) <- c("xmin", "ymin", "xmax", "ymax")
      attr(res_df, attr_name) <- bb
    }
  }

  if (!is.null(master_bbox)) {
    poly_bb <- attr(res_df, "bbox")
    seq_bb <- attr(res_df, "sequence_bbox")
    tol <- 1e-9
    crosses <- poly_bb[["xmin"]] < seq_bb[["xmin"]] - tol ||
      poly_bb[["ymin"]] < seq_bb[["ymin"]] - tol ||
      poly_bb[["xmax"]] > seq_bb[["xmax"]] + tol ||
      poly_bb[["ymax"]] > seq_bb[["ymax"]] + tol
    if (crosses) {
      warning(
        "The bounding box of `wkb` extends beyond `master_bbox`. ",
        "Only the part of the polygon inside `master_bbox` was sampled.",
        call. = FALSE
      )
    }
  }
  res_df
}

#' Backward-compatible wrapper for BAS acceptance sampling
#'
#' @description
#' Generates a spatially balanced sample from a `POLYGON` or `MULTIPOLYGON`
#' using a balanced acceptance sampling (BAS) approach.
#' Points are drawn from a 2-D Owen-scrambled Sobol sequence (Pharr et al 2023)
#' projected onto the study area via an area-equal coordinate transform,
#' so that equal fractions of the sequence correspond to equal surface areas on
#' the sphere.
#'
#' Unlike [grts_h3()], this function requires no H3 discretisation step and
#' therefore scales to very large study areas without memory pressure.
#' H3 cell indices can be obtained afterwards using [h3o::h3_from_xy()] if
#' needed.
#'
#' @param wkb Raw byte vector representing a Well-Known Binary (WKB) geometry
#'   in WGS84 (EPSG:4326). Must be a `POLYGON` or `MULTIPOLYGON`.
#' @param n Integer. The target number of sample locations to draw.
#' @param seed Numeric. A random seed passed to the Rust engine to ensure
#'   reproducible sampling. Must be a whole number `>= 0` and `< 2^53`.
#'   The seed controls the Owen scramble of the Sobol sequence; different seeds
#'   yield statistically independent samples.
#' @param master_bbox Numeric vector `c(xmin, ymin, xmax, ymax)`
#' (WGS-84 decimal degrees).
#'
#' @return A `data.frame` with `n` rows and three columns:
#' * `lon`         – Longitude of the sample point (WGS84 decimal degrees).
#' * `lat`         – Latitude  of the sample point (WGS84 decimal degrees).
#' * `sobol_index` – Index into the raw Sobol stream at which this point was
#'                   accepted. Useful for diagnostics and for verifying prefix
#'                   stability across calls.
#'
#' The returned `data.frame` also carries the following attributes:
#' * `seed`          – The seed used (numeric).
#' * `n`             – Requested sample size (integer).
#' * `sequence_bbox` – Named numeric vector `c(xmin, ymin, xmax, ymax)` of the
#'                     geometry bounding box in decimal degrees.
#'                     Either from master bounding box or polygon bounding box.
#' * `sobol_scanned` – Total number of Sobol stream indices evaluated
#'                     (accepted + rejected). Equals the last `sobol_index` + 1.
#' * `fill_ratio`    – `n / sobol_scanned`. The fraction of bounding-box
#'                     candidates that fell inside the polygon. Low values
#'                     (< 0.05) indicate a sparse geometry relative to its
#'                     bounding box and imply higher computational cost.
#'
#' @section Area-equal sampling:
#' Longitude is mapped linearly from the Sobol stream.
#' Latitude is mapped through the inverse CDF of `sin(lat)` restricted to the
#' bounding box latitude band, so that the expected number of candidate points
#' per unit surface area is uniform across the globe regardless of latitude.
#'
#' @section Coordinate system:
#' Input and output are in WGS84 (EPSG:4326) decimal degrees. No projection
#' is required.
#'
#' @section Prefix stability: calling with the same `seed` and a larger sample
#' size `n2 > n1` produces a sample whose first `n1` rows are identical to the
#' `n1` row result.
#'
#' @references
#' Burley, B. (2020). Practical Hash-based Owen Scrambling.
#' \emph{Journal of Computer Graphics Techniques (JCGT)}, 9(4), 1-20.
#' \url{https://jcgt.org/published/0009/04/01/}
#'
#' Joe, S., & Kuo, F. Y. (2008). Constructing Sobol sequences with better
#' two-dimensional projections.
#' \emph{SIAM Journal on Scientific Computing}, 30(5), 2635-2654.
#' \doi{10.1137/070709359}
#'
#' Pharr, M., Jakob, W., & Humphreys, G. (2023).
#' \emph{Physically Based Rendering: From Theory to Implementation} (4th ed.).
#' MIT Press.
#'
#' Vegdahl, N. (2021). Building a Better LK Hash.
#' \emph{Psychopath Renderer Blog}.
#' \url{https://psychopath.io/post/2021_01_30_building_a_better_lk_hash}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # Load the built-in North Carolina dataset from the sf package
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#'
#' # Use the first county (Ashe) as our study area
#' study_area <- nc[1, ] |> st_geometry() |> st_as_binary(EWKB = FALSE)
#'
#' # Draw a spatially balanced sample of 50 points
#' sample_df <- bas_sobol(
#'   wkb  = study_area[[1]],
#'   n    = 50,
#'   seed = 2021
#' )
#'
#' # Prefix stability: extend to 100 points — first 50 rows are identical
#' sample_df_100 <- bas_sobol(
#'   wkb  = study_area[[1]],
#'   n    = 100,
#'   seed = 2021
#' )
#' stopifnot(identical(sample_df[1:50, c("lon", "lat")],
#'                     sample_df_100[1:50, c("lon", "lat")]))
#'
#' # Convert to sf and visualise
#' sample_sf <- st_as_sf(sample_df, coords = c("lon", "lat"), crs = 4326)
#' plot(st_geometry(nc[1, ]), main = "BAS-Sobol Sample (Ashe County)")
#' plot(st_geometry(sample_sf), add = TRUE, pch = 20, col = "steelblue")
#' }
bas_sobol <- function(wkb, n, seed, master_bbox = NULL) {
  df <- sobol_sample(
    wkb = wkb,
    n = n,
    seed = seed,
    master_bbox = master_bbox,
    algorithm = "acceptance sampling"
  )
  # Legacy column name compatibility
  names(df)[names(df) == "rank"] <- "sobol_index"
  attr(df, "fill_ratio") <- attr(df, "acceptance_rate")
  df
}

#' Backward-compatible wrapper for GRTS routing
#'
#' @description
#' Generates a continuous Generalized Random Tessellation Stratified (GRTS)
#' sample from a `POLYGON` or `MULTIPOLYGON` without relying on discrete grid
#' approximations or rejection sampling.
#'
#' The algorithm generates a 1D spatial sequence (via interleaved 2D
#' Owen-scrambled Sobol sequences) that routes directly through the geometry's
#' exact area.
#' It traverses a lazy top-down quadtree in equal-area `(lon, sin(lat))` space,
#' dynamically clipping the polygon to compute exact intersection areas at every
#' branch.
#' This guarantees equal(-area) inclusion probabilities and (practically)
#' infinite continuous resolution.
#'
#' @param wkb Raw byte vector representing a Well-Known Binary (WKB) geometry
#'   in WGS84 (EPSG:4326). Must be a `POLYGON` or `MULTIPOLYGON`.
#' @param n Integer. The target number of sample locations to draw.
#' @param seed Numeric. A random seed passed to the Rust engine to ensure
#'   reproducible sampling. Must be a whole number `>= 0` and `< 2^53`.
#'   The seed controls the independent Owen scrambles applied to the spatial
#'   dimensions before Morton interleaving.
#' @param master_bbox Numeric vector `c(xmin, ymin, xmax, ymax)`
#'   (WGS-84 decimal degrees).
#'
#' @return A `data.frame` with `n` rows and three columns:
#' * `lon`         – Longitude of the sample point (WGS84 decimal degrees).
#' * `lat`         – Latitude  of the sample point (WGS84 decimal degrees).
#' * `grts_rank` – Spatial GRTS sequence order.
#'
#' The returned `data.frame` also carries the following attributes:
#' * `seed`            – The seed used (numeric).
#' * `n`               – Requested sample size (integer).
#' * `bbox`            – Named numeric vector `c(xmin, ymin, xmax, ymax)` of the
#'                       geometry bounding box in decimal degrees.
#' * `sequence_bbox`   - Named numeric vector `c(xmin, ymin, xmax, ymax)` of the
#'                       geometry bounding box in decimal degrees.
#' * `sobol_scanned`   – Because this is a zero-rejection algorithm, this will
#'                       always exactly equal `n`.
#' * `acceptance_rate` – Always `1.0`. The algorithm strictly routes points
#'                       into valid geometry, completely bypassing rejection.
#'
#' @section Exact Equal-Area Routing:
#' The algorithm maps the sequence to a quadtree rooted at either the
#' bounding box of `wkb` or the `master_bbox`.
#' By operating in `sin(lat)` space, all cells are exactly equal-area on the
#' planetary sphere, eliminating polar distortion.
#' Polygon edges are automatically densified prior to transformation to limit
#' any spherical curvature containment error to < 1.4 m.
#'
#' @section Coordinate system & Antimeridian:
#' Input and output are in WGS84 (EPSG:4326) decimal degrees.
#' No external projection is required.
#' However, polygons straddling the antimeridian (180°/-180°) must be
#' pre-split into a `MULTIPOLYGON` (e.g., using `sf::st_wrap_dateline()`) before
#' conversion to WKB.
#'
#' @section Local Prefix Stability & Master Samples:
#' Calling with the same `seed` and a larger sample size `n2 > n1` on the
#' **same polygon** produces a sample whose first `n1` rows are identical to the
#' `n1` row result.
#'
#' Note: Because the 1D routing relies on the exact local area of the provided
#' polygon to achieve zero rejection, spatial consistency is **not** preserved
#' across overlapping but differently shaped polygons.
#' For globally consistent master samples, users should generate points over a
#' shared bounding box and spatially filter them in R.
#' Alternatively, use `bas_sobol()` for this purpose.
#'
#' @references
#' Robertson, B. L., et al. (2013). Spatial-balanced sampling in continuous
#' space. \emph{Environmental and Ecological Statistics}, 20(1), 149-167.
#'
#' Burley, B. (2020). Practical Hash-based Owen Scrambling.
#' \emph{Journal of Computer Graphics Techniques (JCGT)}, 9(4), 1-20.
#' \url{https://jcgt.org/published/0009/04/01/}
#'
#' Joe, S., & Kuo, F. Y. (2008). Constructing Sobol sequences with better
#' two-dimensional projections.
#' \emph{SIAM Journal on Scientific Computing}, 30(5), 2635-2654.
#' \doi{10.1137/070709359}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # Load the built-in North Carolina dataset from the sf package
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#'
#' # Use the first county (Ashe) as our study area
#' study_area <- nc[1, ] |> st_geometry() |> st_as_binary(EWKB = FALSE)
#'
#' # Draw a spatially balanced sample of 50 points
#' sample_df <- sobol_sample(
#'   wkb  = study_area[[1]],
#'   n    = 50,
#'   seed = 2021,
#'   algorithm = "routing"
#' )
#'
#' # Prefix stability: extend to 100 points — first 50 rows are identical
#' sample_df_100 <- sobol_sample(
#'   wkb  = study_area[[1]],
#'   n    = 100,
#'   seed = 2021,
#'   algorithm = "routing"
#' )
#' stopifnot(identical(sample_df[1:50, c("lon", "lat")],
#'                     sample_df_100[1:50, c("lon", "lat")]))
#'
#' # Convert to sf and visualise
#' sample_sf <- st_as_sf(sample_df, coords = c("lon", "lat"), crs = 4326)
#' plot(st_geometry(nc[1, ]), main = "GRTS-Sobol Sample (Ashe County)")
#' plot(st_geometry(sample_sf), add = TRUE, pch = 20, col = "darkgreen")
#' }
#' @export
grts_sobol <- function(wkb, n, seed, master_bbox = NULL) {
  df <- sobol_sample(
    wkb = wkb,
    n = n,
    seed = seed,
    master_bbox = master_bbox,
    algorithm = "routing"
  )
  # Legacy column name compatibility
  names(df)[names(df) == "rank"] <- "grts_rank"
  df
}

is_valid_bbox <- function(bbox) {
  is.null(bbox) || (
    #Numeric vector c(xmin, ymin, xmax, ymax) (WGS-84 decimal degrees)
    is.double(bbox) &&
      length(bbox) == 4 &&
      bbox_is_decimal_degrees(bbox)
  )
}

bbox_is_decimal_degrees <- function(bbox) {
  all(bbox[c(1, 3)] <= 180) && all(bbox[c(1, 3)] >= -180) &&
    all(bbox[c(2, 4)] <= 90) && all(bbox[c(2, 4)] >= -90)
}
