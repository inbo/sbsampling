# test-grts_h3-equivalence.R
# Tests for equivalence between the
# WKB and pre-computed cells input pathways.

# ---------------------------------------------------------------------------
# 1. WKB vs pre-computed cells pathway equivalence
# ---------------------------------------------------------------------------

test_that("WKB and cells pathways produce identical samples", {
  skip_if_not_installed("sf")
  wkb <- create_test_wkb()
  n <- 10
  res <- 5
  seed <- 123

  out_wkb <- grts_h3(
    wkb         = wkb,
    n           = n,
    resolution  = res,
    containment = "centroid",
    seed        = seed
  )

  coverage_bytes <- h3_get_coverage(
    wkb         = wkb,
    resolution  = res,
    containment = "centroid"
  )
  out_cells <- grts_h3(
    cells = coverage_bytes,
    n     = n,
    seed  = seed
  )

  # The cells path records containment as "pre-computed" — nullify before
  # comparing so the structural equality check is not confused by this
  # expected metadata difference.
  expect_equal(attr(out_wkb,   "containment"), "centroid")
  expect_equal(attr(out_cells, "containment"), "pre-computed")
  expect_equal(attr(out_cells, "n_duplicates_removed"), 0L)
  attr(out_wkb,   "containment") <- NULL
  attr(out_cells, "containment") <- NULL
  attr(out_cells, "n_duplicates_removed") <- NULL

  expect_equal(out_wkb, out_cells)
})

test_that("cells pathway: sort keys for shared cells match the WKB pathway", {
  skip_if_not_installed("sf")
  # The sort key for any H3 cell is a pure function of (cell, seed).
  # A cell appearing in both samples must carry the same sort_key.
  wkb  <- create_test_wkb()
  seed <- 42

  out_wkb <- grts_h3(
    wkb, n = 5, resolution = 4,
    containment = "centroid", seed = seed
  )
  coverage_bytes <- h3_get_coverage(
    wkb, resolution = 4, containment = "centroid"
  )
  out_cells <- grts_h3(cells = coverage_bytes, n = 5, seed = seed)

  shared <- intersect(out_wkb$cell, out_cells$cell)

  keys_wkb   <- out_wkb$sort_key[match(shared, out_wkb$cell)]
  keys_cells <- out_cells$sort_key[match(shared, out_cells$cell)]
  expect_equal(keys_wkb, keys_cells)
})

# ---------------------------------------------------------------------------
# 3. Input validation (dispatch and type checks on R side)
# ---------------------------------------------------------------------------

test_that("grts_h3 errors when neither wkb nor cells is provided", {
  expect_error(
    grts_h3(n = 5, seed = 42),
    "Either 'wkb' or 'cells' must be provided"
  )
})

test_that("grts_h3 errors when cells is a character vector instead of raw", {
  expect_error(
    grts_h3(cells = c("842a993ffffffff"), n = 5, seed = 42),
    "must be a raw byte vector"
  )
})

test_that("grts_h3 errors for invalid resolution", {
  skip_if_not_installed("sf")
  expect_error(grts_h3(create_test_wkb(), n = 3, resolution = -1, seed = 42))
  expect_error(grts_h3(create_test_wkb(), n = 3, resolution = 16, seed = 42))
})

test_that("grts_h3 errors for invalid containment mode", {
  skip_if_not_installed("sf")
  expect_error(
    grts_h3(create_test_wkb(), n = 3, containment = "overlap", seed = 42),
    "should be one of"
  )
})

test_that("grts_h3 containment ordering: intersect >= centroid cells", {
  skip_if_not_installed("sf")
  wkb <- create_test_wkb()
  n_intersect <- attr(
    grts_h3(wkb, n = 1, resolution = 4, containment = "intersect", seed = 42),
    "n_cells"
  )
  n_centroid <- attr(
    grts_h3(wkb, n = 1, resolution = 4, containment = "centroid", seed = 42),
    "n_cells"
  )
  expect_gte(n_intersect, n_centroid)
})
