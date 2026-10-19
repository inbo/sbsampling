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

test_that(
  "acceptance sampling and exact area routing are quasi identical for
  rectangular wkb inputs", {
    # Polygon 1: [0, 1] x [0, 1] from setup.R
    wkb1 <- create_test_wkb()
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

    ac1 <- sobol_sample(wkb1, n = 4, seed = 1, algorithm = "acceptance")
    ro1 <- sobol_sample(wkb1, n = 4, seed = 1, algorithm = "routing")
    expect_equal(
      ac1[, c("lon", "lat", "rank")],
      ro1[, c("lon", "lat", "rank")],
      tolerance = 1e-7
    )
    ac2 <- sobol_sample(wkb2, n = 4, seed = 1, algorithm = "acceptance")
    ro2 <- sobol_sample(wkb2, n = 4, seed = 1, algorithm = "routing")
    expect_equal(
      ac2[, c("lon", "lat", "rank")],
      ro2[, c("lon", "lat", "rank")],
      tolerance = 1e-7
    )

  }
)

test_that("master_bbox smaller than polygon bbox warns", {
  sq <- sf::st_polygon(list(rbind(c(0, 0), c(2, 0), c(2, 2), c(0, 2), c(0, 0))))
  wkb <- sf::st_as_binary(sq)
  expect_warning(
    sobol_sample(wkb, n = 10, seed = 1, master_bbox = c(0, 0, 1, 2)),
    "extends beyond"
  )
})

test_that("master_bbox covering polygon bbox does not warn", {
  sq <- sf::st_polygon(list(rbind(c(0, 0), c(2, 0), c(2, 2), c(0, 2), c(0, 0))))
  wkb <- sf::st_as_binary(sq)
  expect_no_warning(
    sobol_sample(wkb, n = 10, seed = 1, master_bbox = c(-1, -1, 3, 3))
  )
  expect_no_warning(
    sobol_sample(wkb, n = 10, seed = 1, master_bbox = c(0, 0, 2, 2))
  )
})

test_that("sort keys are pinned for the whole coverage", {
  skip_if_not_installed("sf")
  wkb <- create_test_wkb()
  n_all <- attr(grts_h3(wkb, n = 1, resolution = 4, seed = 1), "n_cells")
  r <- grts_h3(wkb, n = n_all, resolution = 4, seed = 215215)
  # Any change to the key pipeline changes this snapshot, i.e. every sample.
  expect_snapshot_value(setNames(r$sort_key, r$cell), style = "deparse")
})

# Hex string -> raw vector of little-endian u64 (the `cells` input format)
h3_hex_to_raw <- function(hex) {
  hex <- paste0(strrep("0", 16 - nchar(hex)), hex)
  bytes <- vapply(
    hex,
    function(h) strtoi(substring(h, seq(1, 15, by = 2), seq(2, 16, by = 2)), 16L),
    integer(8)
  )
  as.raw(as.vector(bytes[8:1, , drop = FALSE]))
}

test_that("res-15 sort keys are pinned across base cells and a pentagon", {
  # Expected keys come from an independent re-implementation of the key
  # (base-cell shuffle, per-parent permutations, bit layout), seed = 42.
  # Listed in ascending key order, which is the output order.
  expected <- c(
    "8f464b968480463" = "000206966a2dc89f", # Honolulu
    "8f7eb57221acba8" = "0003762b286a81b3", # near antimeridian
    "8ff2938565b506e" = "00052c403d32a73e", # near south pole
    "8f0326222a584e8" = "00085a58c61d4127", # near north pole
    "8f0800000000000" = "0008a4262a862d46", # pentagon (base cell 4)
    "8f66d33ab84b8db" = "0008f4d94631820c", # Quito
    "8f1fa4418c6905a" = "000a61af12d8a834", # Brussels
    "8fbe0e35cb15192" = "000d25a93602b308"  # Sydney
  )
  r <- grts_h3(
    cells = h3_hex_to_raw(names(expected)),
    n = length(expected),
    seed = 42
  )
  expect_equal(setNames(r$sort_key, r$cell), expected)
})

test_that("sort keys do not depend on the other cells or on input order", {
  skip_if_not_installed("sf")
  cov <- h3_get_coverage(create_test_wkb(), resolution = 5)
  m <- matrix(cov, nrow = 8)             # one H3 index (8 raw bytes) per column
  n_all <- ncol(m)
  full <- grts_h3(cells = cov, n = n_all, seed = 42)

  # same cells, reversed input order
  rev_cov <- as.vector(m[, rev(seq_len(n_all))])
  expect_equal(grts_h3(cells = rev_cov, n = n_all, seed = 42), full)

  # a subset of the cells: shared cells must carry identical keys
  idx <- seq(1, n_all, by = 3)
  sub <- grts_h3(cells = as.vector(m[, idx]), n = length(idx), seed = 42)
  expect_equal(sub$sort_key, full$sort_key[match(sub$cell, full$cell)])
})
