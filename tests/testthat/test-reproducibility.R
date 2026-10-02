# test-reproducibility.R
# Tests that sampler output is deterministic for a given seed, differs
# across seeds and seeds are cross-platform stable

# ---------------------------------------------------------------------------
# 1. Snapshot — cross-OS portability anchor
# ---------------------------------------------------------------------------

test_that("output is strictly portable across OS architectures", {
  skip_if_not_installed("sf")
  test_wkb <- create_test_wkb()
  grts_h3_out <- grts_h3(
    wkb         = test_wkb,
    n           = 5,
    resolution  = 4,
    containment = "centroid",
    seed        = 42
  )
  # style = "serialize" captures exact binary equivalence including column
  # types, float bit patterns, and attribute values.
  expect_snapshot_value(grts_h3_out, style = "serialize")

  grts_sobol_out <- grts_sobol(
    wkb         = test_wkb,
    n           = 5,
    seed        = 42
  )
  expect_snapshot_value(grts_sobol_out, style = "serialize")

  bas_sobol_out <- bas_sobol(
    wkb         = test_wkb,
    n           = 5,
    seed        = 42
  )
  expect_snapshot_value(bas_sobol_out, style = "serialize")
})

# ---------------------------------------------------------------------------
# 2. Within session test
# ---------------------------------------------------------------------------

test_that("same seed produces identical results across calls", {
  skip_if_not_installed("sf")
  wkb <- create_test_wkb()
  r1 <- grts_h3(wkb, n = 5, resolution = 4, seed = 42)
  r2 <- grts_h3(wkb, n = 5, resolution = 4, seed = 42)
  expect_equal(r1, r2)

  r1 <- grts_h3(wkb, n = 100, resolution = 6, area_correction = TRUE, seed = 7)
  r2 <- grts_h3(wkb, n = 100, resolution = 6, area_correction = TRUE, seed = 7)
  expect_equal(r1, r2)

  r1 <- grts_sobol(wkb, n = 5, seed = 42)
  r2 <- grts_sobol(wkb, n = 5, seed = 42)
  expect_equal(r1, r2)

  r1 <- bas_sobol(wkb, n = 5, seed = 42)
  r2 <- bas_sobol(wkb, n = 5, seed = 42)
  expect_equal(r1, r2)
})

test_that("different seeds produce different selections", {
  skip_if_not_installed("sf")
  wkb <- create_test_wkb()
  r1 <- grts_h3(wkb, n = 5, resolution = 4, seed = 1)
  r2 <- grts_h3(wkb, n = 5, resolution = 4, seed = 2)
  expect_false(identical(r1$cell, r2$cell))

  r1 <- grts_h3(wkb, n = 5, resolution = 4, area_correction = TRUE, seed = 7)
  r2 <- grts_h3(wkb, n = 5, resolution = 4, area_correction = TRUE, seed = 8)
  expect_false(identical(r1$cell, r2$cell))

  r1 <- grts_sobol(wkb, n = 5, seed = 1)
  r2 <- grts_sobol(wkb, n = 5, seed = 2)
  expect_false(identical(r1$lon, r2$lon))

  r1 <- bas_sobol(wkb, n = 5, seed = 1)
  r2 <- bas_sobol(wkb, n = 5, seed = 2)
  expect_false(identical(r1$lon, r2$lon))
})


test_that(
  "master_bbox maintains spatial sequence across overlapping polygons", {
    # 1. Setup geometries
    # Polygon 1: [0, 1] x [0, 1] from setup.R
    wkb1 <- create_test_wkb()

    # Polygon 2: [0.5, 1.5] x [0.5, 1.5] - heavily overlapping Polygon 1
    poly2 <- sf::st_polygon(
      list(matrix(
        c(0.5, 0.5,
          0.5, 1.5,
          1.5, 1.5,
          1.5, 0.5,
          0.5, 0.5),
        ncol = 2, byrow = TRUE
      ))
    )
    wkb2 <- sf::st_as_binary(sf::st_sfc(poly2, crs = 4326), EWKB = FALSE)[[1]]

    # 2. Define a master_bbox that covers both polygons
    mbbox <- c(0, 0, 1.5, 1.5)

    # 3. Sample both geometries with the same seed and master_bbox
    seed <- 42
    target_n <- 50

    res1 <- bas_sobol(wkb1, n = target_n, seed = seed, master_bbox = mbbox)
    res2 <- bas_sobol(wkb2, n = target_n, seed = seed, master_bbox = mbbox)

    # 4. Verify the sequence bounds were correctly registered in attributes
    expect_equal(unname(attr(res1, "sequence_bbox")), mbbox)
    expect_equal(unname(attr(res2, "sequence_bbox")), mbbox)

    # 5. Extract points that were sampled in both runs
    # Because the sequence is master-bbox-relative, points falling in the
    # intersection of wkb1 and wkb2 MUST share the exact same Sobol index.
    shared_indices <- intersect(res1$sobol_index, res2$sobol_index)

    # Assert that overlapping points exist (with n=50 and 25% area overlap,
    # we confidently expect multiple shared points).
    expect_true(length(shared_indices) > 0)

    # 6. Verify identical coordinates for shared indices
    # Subset both result sets to just the shared indices
    shared_res1 <- res1[res1$sobol_index %in% shared_indices, ]
    shared_res2 <- res2[res2$sobol_index %in% shared_indices, ]
    attributes(shared_res1) <- NULL
    attributes(shared_res2) <- NULL

    # The coordinates mapped from the Sobol stream must be identical
    expect_equal(shared_res1, shared_res2)

    # 7. Verification of independent (non-master) sampling behavior
    # If we sample WITHOUT master_bbox, the domains shift, and sobol_index
    # collisions will yield completely different geographic coordinates.
    res_independent <- bas_sobol(
      wkb2, n = target_n, seed = seed, master_bbox = NULL
    )

    # Find an index shared between the master-bbox run and independent run
    control_indices <- intersect(res2$sobol_index, res_independent$sobol_index)
    expect_true(length(control_indices) > 0)

    idx <- control_indices[1]
    indep_point <- res_independent[res_independent$sobol_index == idx, ]
    master_point <- res2[res2$sobol_index == idx, ]

    # Assert they do NOT map to the same coordinates when domains differ
    expect_false(
      isTRUE(all.equal(indep_point$lon, master_point$lon)) &&
        isTRUE(all.equal(indep_point$lat, master_point$lat))
    )
  }
)
