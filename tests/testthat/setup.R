# ---------------------------------------------------------------------------
# Setup: shared helpers for sbsampling tests
# ---------------------------------------------------------------------------

# Create a valid WKB polygon for testing.
# A 1x1 degree square near the equator — produces a tractable number of H3
# cells at resolutions 3-5. Requires sf (listed in Suggests).
create_test_wkb <- function() {
  poly <- sf::st_polygon(
    list(matrix(c(0, 0, 0, 1, 1, 1, 1, 0, 0, 0), ncol = 2, byrow = TRUE))
  )
  sf::st_as_binary(sf::st_sfc(poly, crs = 4326), EWKB = FALSE)[[1]]
}

# A tiny polygon that produces at most a handful of cells even at fine
# resolution. Used to test n > N error paths cheaply.
create_tiny_wkb <- function() {
  poly <- sf::st_polygon(
    list(matrix(
      c(0, 0, 0, 0.05, 0.05, 0.05, 0.05, 0, 0, 0),
      ncol = 2, byrow = TRUE
    ))
  )
  sf::st_as_binary(sf::st_sfc(poly, crs = 4326), EWKB = FALSE)[[1]]
}
