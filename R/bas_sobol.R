#' Generate a Spatially Balanced BAS Sample using a Sobol Sequence
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
  # 0. Assertions
  stopifnot(
    "Input 'wkb' must be a raw byte vector representing WKB geometry." =
      is.raw(wkb),
    "`n` must be a whole number, > 0, not `NA`." =
      is_valid_n(n),
    "`seed` must be a whole number, >= 0, and < 2^53." =
      is_valid_seed(seed),
    "`master_bbox` must be length 4 numeric in decimal degrees" =
      is_valid_bbox(master_bbox)
  )

  # Call the Rust engine
  res_df <- bas_sample_from_wkb(
    wkb_bytes = wkb,
    n = as.integer(n),
    seed = as.numeric(seed),
    master_bbox = master_bbox
  )

  # Attach names to the bbox attribute
  bbox <- attr(res_df, "sequence_bbox")
  names(bbox) <- c("xmin", "ymin", "xmax", "ymax")
  attr(res_df, "sequence_bbox") <- bbox

  return(res_df)
}

is_valid_bbox <- function(bbox) {
  is.null(bbox) || (
    #Numeric vector c(xmin, ymin, xmax, ymax) (WGS-84 decimal degrees)
    is.numeric(bbox) &&
      length(bbox) == 4 &&
      bbox_is_decimal_degrees(bbox)
  )
}

bbox_is_decimal_degrees <- function(bbox) {
  all(bbox[c(1, 3)] <= 180) && all(bbox[c(1, 3)] >= -180) &&
    all(bbox[c(2, 4)] <= 90) && all(bbox[c(2, 4)] >= -90)
}
