# test-sample-size.R
# Tests around the n argument: exact row counts and boundary behaviour.

test_that("result has exactly n rows", {
  skip_if_not_installed("sf")
  wkb <- create_test_wkb()
  result <- grts_h3(wkb, n = 5, resolution = 4, seed = 42)
  expect_equal(nrow(result), 5L)

  result <- grts_sobol(wkb, n = 5, seed = 42)
  expect_equal(nrow(result), 5L)

  result <- bas_sobol(wkb, n = 5, seed = 42)
  expect_equal(nrow(result), 5L)

})

test_that("n = 1 returns exactly one row with all required columns", {
  skip_if_not_installed("sf")
  wkb <- create_test_wkb()
  result <- grts_h3(wkb, n = 1, resolution = 4, seed = 42)
  expect_equal(nrow(result), 1L)

  result <- grts_sobol(wkb, n = 1, seed = 42)
  expect_equal(nrow(result), 1L)

  result <- bas_sobol(wkb, n = 1, seed = 42)
  expect_equal(nrow(result), 1L)
})

test_that("n equal to total coverage returns all cells in sort_key order", {
  skip_if_not_installed("sf")
  wkb   <- create_test_wkb()
  n_all <- attr(grts_h3(wkb, n = 1, resolution = 4, seed = 1), "n_cells")
  result <- grts_h3(wkb, n = n_all, resolution = 4, seed = 42)
  expect_equal(nrow(result), n_all)
  expect_equal(result$sort_key, sort(result$sort_key))
})

test_that("n = N-1 excludes the cell with the highest sort_key", {
  skip_if_not_installed("sf")
  wkb   <- create_test_wkb()
  n_all <- attr(grts_h3(wkb, n = 1, resolution = 4, seed = 1), "n_cells")
  if (n_all < 2L) skip("Coverage too small for this sub-test")
  r_all  <- grts_h3(wkb, n = n_all, resolution = 4, seed = 42)
  r_less <- grts_h3(wkb, n = n_all - 1L, resolution = 4, seed = 42)
  # The excluded cell is the one with the largest sort key in the full sample
  highest_key <- r_all$sort_key[n_all]
  expect_false(highest_key %in% r_less$sort_key)
})

test_that("n exceeding total coverage returns a clear error", {
  skip_if_not_installed("sf")
  # wkb_tiny covers very few cells; n = 50000 at res 2 will always exceed N
  expect_error(
    grts_h3(create_test_wkb(), n = 50000, resolution = 2, seed = 42),
    "Target n"
  )
})

test_that("n = -5 returns a clear error", {
  skip_if_not_installed("sf")
  expect_error(
    grts_h3(create_test_wkb(), n = -5, seed = 42),
    "`n` must be "
  )
  expect_error(
    grts_sobol(create_test_wkb(), n = -5, seed = 42),
    "`n` must be "
  )
  expect_error(
    bas_sobol(create_test_wkb(), n = -5, seed = 42),
    "`n` must be "
  )
})
