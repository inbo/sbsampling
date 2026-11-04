test_that(
  "acceptance sampling fails fast when geometry misses the master bbox", {
    skip_if_not_installed("sf")
    tri <- sf::st_as_binary(sf::st_polygon(list(
      rbind(c(0, 0), c(10, 0), c(0, 10), c(0, 0))
    )))
    t <- system.time(
      expect_error(
        sobol_sample(tri, n = 100, seed = 1,
                     master_bbox = c(8, 8, 10, 10),
                     algorithm = "acceptance sampling"),
        "zero computable area"
      )
    )[["elapsed"]]
    expect_lt(t, 5)
  }
)

test_that("non-polygon geometries give a clear error", {
  skip_if_not_installed("sf")
  line <- sf::st_as_binary(sf::st_linestring(rbind(c(0, 0), c(1, 1))))
  for (alg in c("acceptance sampling", "routing", "warp")) {
    expect_error(
      sobol_sample(line, n = 10, seed = 1, algorithm = alg),
      "Unsupported geometry type: LineString"
    )
  }
})
