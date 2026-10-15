# test-h3_get_coverage.R
# Tests for h3_get_coverage(): input validation, output format, and behaviour.

# ---------------------------------------------------------------------------
# 1. Input validation
# ---------------------------------------------------------------------------

test_that("h3_get_coverage rejects non-raw wkb input", {
  expect_error(
    h3_get_coverage(wkb = "not a raw vector"),
    "Input 'wkb' must be a raw byte vector"
  )
  expect_error(
    h3_get_coverage(wkb = NULL),
    "Input 'wkb' must be a raw byte vector"
  )
})

test_that("h3_get_coverage rejects invalid resolution", {
  skip_if_not_installed("sf")
  valid_wkb <- create_test_wkb()
  expect_error(h3_get_coverage(valid_wkb, resolution = -1))
  expect_error(h3_get_coverage(valid_wkb, resolution = 16))
  expect_error(h3_get_coverage(valid_wkb, resolution = "5"))
})

test_that("h3_get_coverage rejects invalid containment mode", {
  skip_if_not_installed("sf")
  valid_wkb <- create_test_wkb()
  expect_error(
    h3_get_coverage(valid_wkb, containment = "inside"),
    "should be one of"
  )
})

# ---------------------------------------------------------------------------
# 2. Output format
# ---------------------------------------------------------------------------

test_that("h3_get_coverage returns a raw bincode vector", {
  skip_if_not_installed("sf")
  cells_raw <- h3_get_coverage(create_test_wkb(), resolution = 4)
  expect_type(cells_raw, "raw")
})

# ---------------------------------------------------------------------------
# 3. Behaviour
# ---------------------------------------------------------------------------

test_that("h3_get_coverage accepts all valid containment modes", {
  skip_if_not_installed("sf")
  valid_wkb <- create_test_wkb()
  expect_no_error(h3_get_coverage(valid_wkb, containment = "centroid"))
  expect_no_error(h3_get_coverage(valid_wkb, containment = "intersect"))
  expect_no_error(h3_get_coverage(valid_wkb, containment = "boundary"))
  expect_no_error(h3_get_coverage(valid_wkb, containment = "covers"))
})

test_that("h3_get_coverage produces more bytes at higher resolution", {
  skip_if_not_installed("sf")
  valid_wkb <- create_test_wkb()
  res3 <- h3_get_coverage(valid_wkb, resolution = 3)
  res4 <- h3_get_coverage(valid_wkb, resolution = 4)
  # Each H3 cell encodes to 8 bytes; more cells → larger payload
  expect_true(length(res4) > length(res3))
})

# ---------------------------------------------------------------------------
# 4. Snapshot stability
# ---------------------------------------------------------------------------

test_that("h3_get_coverage produces stable byte payloads across builds", {
  skip_if_not_installed("sf")
  cells_raw <- h3_get_coverage(
    create_test_wkb(), resolution = 4, containment = "centroid"
  )
  expect_snapshot_value(cells_raw, style = "serialize")
})
