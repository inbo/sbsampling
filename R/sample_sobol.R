#' Generate a Spatially Balanced Sample using a Sobol Sequence
#'
#' @description
#' Generates a spatially balanced sample from a `POLYGON` or `MULTIPOLYGON`.
#' Points are drawn from a 2-D Owen-scrambled Sobol sequence (Pharr et al 2023)
#' projected onto the study area via an area-equal coordinate transform,
#' so that equal fractions of the sequence correspond to equal surface areas on
#' the sphere.
#'
#' Unlike [sample_hex()], this function requires no H3 discretisation step and
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
#'     \item \code{"acceptance sampling"}: Balanced Acceptance Sampling (BAS)
#'      algorithm. Points are generated in the bounding box and are accepted
#'      when they fall inside the polygon.
#'     \item \code{"warp"}: (2-D hierarchical sample warping): Keeps the two
#'      Sobol coordinates (`u1`, `u2`) separate and pushes them down the same
#'      exact-area quadtree as a binary-tree Rosenblatt transform (x-split
#'      with `u1`, then conditional y-split with `u2`, rescaling both after
#'      every level).
#'      Zero rejection, prefix stable.
#'     \item \code{"routing"}: Maps a 1D Morton-interleaved sequence down an
#'      exact-area quadtree (zero rejection, prefix stable). Ideal for sparse
#'      and narrow domains.
#'      The algorithm generates a 1D spatial sequence (via interleaved 2D
#'      Owen-scrambled Sobol sequences) that routes directly through the
#'      geometry's exact area.
#'      It traverses a lazy top-down quadtree in equal-area `(lon, sin(lat))`
#'      space, dynamically clipping the polygon to compute exact intersection
#'      areas at every branch.
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
#' not for `"routing"` or `"warp"`, which depends on the exact area of each
#' polygon.
#'
#' @section Area-equal sampling:
#' Longitude is mapped linearly from the Sobol stream.
#' Authalic latitude is mapped through the inverse CDF of `sin(lat)` restricted
#' to the bounding box latitude band, so that the expected number of candidate
#' points per unit surface area is uniform across the globe regardless of
#' latitude.
#'
#' @section Coordinate system:
#' Input and output are in longitude, latitude decimal degrees without assuming
#' a specific datum.
#' Projection to planar coordinates is not required.
#' However, polygons straddling the antimeridian (180°/-180°) must be
#' pre-split into a `MULTIPOLYGON` (e.g., using `sf::st_wrap_dateline()`) before
#' conversion to WKB.
#'
#' @section Prefix stability:
#' Calling with the same `seed` and a larger sample size `n2 > n1` on the
#' **same polygon** produces a sample whose first `n1` rows are identical to the
#' `n1` row result.
#'
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
#' Robertson, B. L., et al. (2013). Spatial-balanced sampling in continuous
#' space. \emph{Environmental and Ecological Statistics}, 20(1), 149-167.
#'
#' Vegdahl, N. (2021). Building a Better LK Hash.
#' \emph{Psychopath Renderer Blog}.
#' \url{https://psychopath.io/post/2021_01_30_building_a_better_lk_hash}
#'
#' @return A `data.frame` with `n` rows and three columns:
#' * `lon`         – Longitude of the sample point (decimal degrees).
#' * `lat`         – Latitude  of the sample point (decimal degrees).
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
#' * `acceptance_rate` – `n / sobol_scanned`. The fraction of bounding-box
#'                     candidates that fell inside the polygon. Low values
#'                     (< 0.05) indicate a sparse geometry relative to its
#'                     bounding box.
#' @export
sample_sobol <- function(
  wkb,
  n,
  seed,
  master_bbox = NULL,
  algorithm = c("acceptance sampling", "warp", "routing")
) {
  algorithm <- match.arg(algorithm)

  stopifnot(
    "Input 'wkb' must be a raw byte vector representing WKB geometry." =
      is.raw(wkb),
    "`n` must be a whole number, > 0, not `NA`." = is_valid_n(n),
    "`seed` must be a whole number, >= 0, and < 2^53." = is_valid_seed(seed),
    "`master_bbox` must be length 4 numeric in decimal degrees" =
      is_valid_bbox(master_bbox)
  )

  res_df <- sample_sobol_from_wkb(
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
