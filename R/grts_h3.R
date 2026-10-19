#' Generate a Spatially Balanced GRTS Sample using H3
#'
#' @description
#' Generates a continuous Generalized Random-Tessellation Stratified (GRTS)
#' sample.
#' The function can accept either a geometry (`wkb`) to tessellate
#' on the fly, or a pre-computed raw byte vector of H3 cells (`cells`).
#'
#' @param wkb Raw byte vector representing Well-Known Binary (WKB) geometry in
#' WGS84 (EPSG:4326).
#' @param cells Raw byte vector of pre-computed H3 cell indices, all at the
#'  same resolution. Use `h3_get_coverage()` to generate them. Invalid indices
#'  or mixed resolutions are an error; duplicates are removed with a warning.
#'  If provided, `wkb`, `resolution`, and `containment` are ignored.
#' @param n Integer. The target number of sample locations to draw.
#' @param resolution The H3 resolution to use (0-15).
#' Ignored if `cells` is provided.
#' @param containment The containment mode to use. One of "centroid",
#'  "intersect", "boundary", or "covers". Ignored if `cells` is provided.
#' @param area_correction If `FALSE` (default), all H3 cells have equal
#'  inclusion probability pi_i = n / N.
#'  Cell areas vary by up to ~1.9x across the icosahedral projection;
#'  this is ignored in the sampling step but can be corrected at the estimation
#'  stage using the returned `area_m2` and `ip` columns.
#'  If `TRUE`, Sequential Poisson Sampling is applied, making the effective
#'  inclusion probability approximately proportional to area:
#'  pi_i = n * area_i / sum(area).
#' @param seed Numeric. A random seed passed to the Rust engine to ensure a
#'   reproducible hierarchical shuffle. The seed must be a whole number,
#'   `>= 0`, and `< 2^53`.
#'
#'
#' @return A `data.frame` with four columns and `n` rows, ordered by ascending
#'  effective key (GRTS visiting order):
#' * `cell`     – H3 cell identifier as a lowercase hexadecimal string.
#' * `sort_key` – GRTS sort key as a zero-padded 16-character lowercase hex
#'                string. Zero-padding ensures correct lexicographic ordering
#'                (e.g. for merging tile results).
#' * `effective_key` – Value used to rank cells for selection; the cells with
#'                     the smallest values are sampled. Equal to the numeric
#'                     value of `sort_key` if `area_correction = FALSE`, and
#'                     `sort_key / area_m2` if `TRUE`. The data frame is
#'                     ordered by this column.
#' * `area_m2` – True area of the cell in m².
#' * `ip`      – Inclusion probability pi_i of the sampled cell:
#'                  area_correction = FALSE: pi_i = n / N  (constant)
#'                  area_correction = TRUE:  pi_i = n * area_i / sum(area)
#'
#' Use `h3o::h3_to_points(h3o::h3_from_strings())` to convert the `cell` strings
#' back to `sf` geometries.
#'
#' The returned `data.frame` also carries the following attributes which
#' provide the quantities needed to reconstruct design-based estimators:
#' * `n_cells`     – Total number of H3 cells N covering the study area.
#'                   Used to compute pi_i = n / N.
#' * `sum_area_m2`– Sum of all cell areas in m² (i.e. the approximate
#'                   area of the study area as seen by the H3 grid).
#' * `n`           – Requested sample size.
#' * `resolution`  – H3 resolution used.
#' * `containment` – Containment mode used.
#' * `area_correction` – Whether area-proportional sampling was applied.
#' * `n_duplicates_removed` – (cells input only) number of duplicate cells
#'                          dropped before sampling.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(h3o)
#'
#' # Load the built-in North Carolina dataset from the sf package
#' nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#'
#' # Use the first county (Ashe) as our study area
#' study_area <- nc[1, ] |> st_geometry() |> st_as_binary(EWKB = FALSE)
#'
#' # Draw a sample of 25 points at resolution 7
#' sample_df <- grts_h3(
#'   wkb = study_area[[1]],
#'   n = 25,
#'   resolution = 7,
#'   containment = "centroid",
#'   seed = 123
#' )
#'
#' # Convert the H3 indices back into spatial points
#' sample_df$geom <- h3_to_points(h3_from_strings(sample_df$cell))
#' sample_df <- st_as_sf(sample_df)
#'
#' # Visualize the spatially balanced sample
#' plot(st_geometry(nc[1,]), main = "H3-GRTS Sample (Ashe County)")
#' plot(st_geometry(sample_df), add = TRUE, pch = 20, col = "red")
#' }
grts_h3 <- function(
  wkb,
  cells,
  n,
  resolution = 5,
  containment = c("centroid", "intersect", "boundary", "covers"),
  seed,
  area_correction = FALSE
) {
  # 0. Assertions
  stopifnot(
    "Either 'wkb' or 'cells' must be provided, not both." =
      xor(!missing(wkb), !missing(cells)),
    is.logical(area_correction) && !is.na(area_correction),
    "`n` must be a whole number, > 0, not `NA`" = is_valid_n(n),
    "`seed` must be a whole number, >= 0, and < 2^53." = is_valid_seed(seed)
  )

  containment <- match.arg(containment)

  # 1. Call the appropriate Rust Pipeline
  if (!missing(cells)) {
    stopifnot("Input 'cells' must be a raw byte vector." = is.raw(cells))

    res_df <- grts_sample_from_cells(
      cells_bytes = cells,
      n = as.integer(n),
      global_seed = as.numeric(seed),
      area_correction = area_correction
    )
    n_dup <- attr(res_df, "n_duplicates_removed")
    if (!is.null(n_dup) && n_dup > 0) {
      warning(
        n_dup,
        " duplicate H3 cell(s) in `cells` were removed.",
        call. = FALSE
      )
    }
    return(res_df)
  }

  # wkb path
  stopifnot(
    "Input 'wkb' must be a raw byte vector representing WKB geometry." =
      is.raw(wkb),
    "`resolution` must be a whole number, >= 0, and <= 15." =
      is_valid_resolution(resolution)
  )

  res_df <- grts_sample_from_wkb(
    wkb_bytes = wkb,
    res = as.integer(resolution),
    containment = containment,
    n = as.integer(n),
    global_seed = as.numeric(seed),
    area_correction = area_correction
  )
  return(res_df)
}

is_valid_resolution <- function(resolution) {
  is.numeric(resolution) &&
    resolution %% 1 == 0 &&
    resolution >= 0 &&
    resolution <= 15
}

is_valid_seed <- function(seed) {
  !missing(seed) &&
    is.numeric(seed) && seed %% 1 == 0 && seed >= 0 && seed < 2^53
}

is_valid_n <- function(n) {
  !missing(n) && is.numeric(n) && n %% 1 == 0 && n > 0
}
