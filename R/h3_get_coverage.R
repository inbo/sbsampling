#' Get H3 Coverage for a Polygon
#'
#' @description
#' Tessellates a study area polygon into an H3 hexagonal grid and returns
#' the vector of cell indices. Useful for caching grid covers to speed up
#' repeated simulations.
#'
#' @param wkb Raw byte vector representing Well-Known Binary (WKB) geometry.
#' @param resolution The H3 resolution to use (0-15).
#' @param containment The containment mode to use. One of "centroid",
#' "intersect", "boundary", or "covers".
#'
#' @return A character vector of H3 cell identifiers.
#' @export
h3_get_coverage <- function(
  wkb,
  resolution = 5,
  containment = c("centroid", "intersect", "boundary", "covers")
) {
  stopifnot(
    "Input 'wkb' must be a raw byte vector representing WKB geometry." =
      is.raw(wkb),
    "`resolution` must be a whole number, >= 0, and <= 15." =
      is_valid_resolution(resolution)
  )
  containment <- match.arg(containment)

  h3_coverage(
    wkb_bytes = wkb,
    res = as.integer(resolution),
    containment = containment
  )
}
