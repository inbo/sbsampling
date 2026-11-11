# test-get_coverage_hex.R
# Tests for get_coverage_hex(): input validation, output format, and behaviour.

# ---------------------------------------------------------------------------
# 1. Input validation
# ---------------------------------------------------------------------------

test_that("get_coverage_hex rejects non-raw wkb input", {
  expect_error(
    get_coverage_hex(wkb = "not a raw vector"),
    "Input 'wkb' must be a raw byte vector"
  )
  expect_error(
    get_coverage_hex(wkb = NULL),
    "Input 'wkb' must be a raw byte vector"
  )
})

test_that("get_coverage_hex rejects invalid resolution", {
  skip_if_not_installed("sf")
  valid_wkb <- create_test_wkb()
  expect_error(get_coverage_hex(valid_wkb, resolution = -1))
  expect_error(get_coverage_hex(valid_wkb, resolution = 16))
  expect_error(get_coverage_hex(valid_wkb, resolution = "5"))
})

test_that("get_coverage_hex rejects invalid containment mode", {
  skip_if_not_installed("sf")
  valid_wkb <- create_test_wkb()
  expect_error(
    get_coverage_hex(valid_wkb, containment = "inside"),
    "should be one of"
  )
})

# ---------------------------------------------------------------------------
# 2. Output format
# ---------------------------------------------------------------------------

test_that("get_coverage_hex returns a raw bincode vector", {
  skip_if_not_installed("sf")
  cells_raw <- get_coverage_hex(create_test_wkb(), resolution = 4)
  expect_type(cells_raw, "raw")
})

# ---------------------------------------------------------------------------
# 3. Behaviour
# ---------------------------------------------------------------------------

test_that("get_coverage_hex accepts all valid containment modes", {
  skip_if_not_installed("sf")
  valid_wkb <- create_test_wkb()
  expect_no_error(get_coverage_hex(valid_wkb, containment = "centroid"))
  expect_no_error(get_coverage_hex(valid_wkb, containment = "intersect"))
  expect_no_error(get_coverage_hex(valid_wkb, containment = "boundary"))
  expect_no_error(get_coverage_hex(valid_wkb, containment = "covers"))
})

test_that("get_coverage_hex produces more bytes at higher resolution", {
  skip_if_not_installed("sf")
  valid_wkb <- create_test_wkb()
  res3 <- get_coverage_hex(valid_wkb, resolution = 3)
  res4 <- get_coverage_hex(valid_wkb, resolution = 4)
  # Each H3 cell encodes to 8 bytes; more cells → larger payload
  expect_true(length(res4) > length(res3))
})

# ---------------------------------------------------------------------------
# 4. Snapshot stability
# ---------------------------------------------------------------------------

test_that("get_coverage_hex produces stable byte payloads across builds", {
  skip_if_not_installed("sf")
  cells_raw <- get_coverage_hex(
    create_test_wkb(), resolution = 4, containment = "centroid"
  )
  expect_snapshot_value(cells_raw, style = "serialize")
})
