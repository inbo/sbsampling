# test-grts_h3-ip.R
# Tests for inclusion probabilities (ip column), area_m2, and design
# attributes. Covers both equal-probability and area-proportional modes.

# ---------------------------------------------------------------------------
# 1. area_correction = FALSE (equal-probability)
# ---------------------------------------------------------------------------

test_that("ip unit tests when area_correction = FALSE", {
  skip_if_not_installed("sf")
  result <- grts_h3(
    create_test_wkb(), n = 5, resolution = 4, area_correction = FALSE, seed = 42
  )
  expect_true(all(result$ip == result$ip[1]))
  expected_ip <- attr(result, "n") / attr(result, "n_cells")
  expect_equal(result$ip, rep(expected_ip, nrow(result)))
  expect_true(all(result$ip > 0))
  expect_true(all(result$ip <= 1))
})

test_that("sum of ip over all N cells equals n when area_correction = FALSE", {
  skip_if_not_installed("sf")
  wkb   <- create_test_wkb()
  n_all <- attr(grts_h3(wkb, n = 1, resolution = 4, seed = 1), "n_cells")
  result <- grts_h3(wkb, n = n_all, resolution = 4,
                    area_correction = FALSE, seed = 42)
  expect_equal(sum(result$ip), as.numeric(n_all), tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 2. area_correction = TRUE (area-proportional)
# ---------------------------------------------------------------------------

test_that("ip varies across cells when area_correction = TRUE", {
  skip_if_not_installed("sf")
  result <- grts_h3(create_test_wkb(), n = 5, resolution = 4,
                    area_correction = TRUE, seed = 42)
  expect_all_false(result$ip == mean(result$ip))
})

test_that("ip formula holds when area_correction = TRUE", {
  skip_if_not_installed("sf")
  result   <- grts_h3(create_test_wkb(), n = 5, resolution = 4,
                      area_correction = TRUE, seed = 42)
  expected <- attr(result, "n") * result$area_m2 / attr(result, "sum_area_m2")
  expect_equal(result$ip, expected, tolerance = 1e-12)
})

test_that(
  "ip rank order matches area_m2 rank order when area_correction = TRUE", {
    skip_if_not_installed("sf")
    result <- grts_h3(create_test_wkb(), n = 5, resolution = 4,
                      area_correction = TRUE, seed = 42)
    expect_equal(rank(result$ip), rank(result$area_m2))
  }
)

test_that("ip values are in (0, 1] when area_correction = TRUE", {
  skip_if_not_installed("sf")
  result <- grts_h3(create_test_wkb(), n = 5, resolution = 4,
                    area_correction = TRUE, seed = 42)
  expect_true(all(result$ip > 0))
  expect_true(all(result$ip <= 1))
})

test_that("sum of ip over all N cells equals n when area_correction = TRUE", {
  skip_if_not_installed("sf")
  wkb   <- create_test_wkb()
  n_all <- attr(grts_h3(wkb, n = 1, resolution = 4, seed = 1), "n_cells")
  result <- grts_h3(wkb, n = n_all, resolution = 4,
                    area_correction = TRUE, seed = 42)
  expect_equal(sum(result$ip), as.numeric(n_all), tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 3. area_m2 and sum_area_m2 consistency
# ---------------------------------------------------------------------------

test_that("all area_m2 values are positive", {
  skip_if_not_installed("sf")
  result <- grts_h3(create_test_wkb(), n = 5, resolution = 4, seed = 42)
  expect_true(all(result$area_m2 > 0))
})

test_that("sum_area_m2 equals sum of area_m2 over the full coverage", {
  skip_if_not_installed("sf")
  wkb   <- create_test_wkb()
  n_all <- attr(grts_h3(wkb, n = 1, resolution = 4, seed = 1), "n_cells")
  result <- grts_h3(wkb, n = n_all, resolution = 4, seed = 42)
  expect_equal(
    sum(result$area_m2),
    attr(result, "sum_area_m2"),
    tolerance = 1e-9
  )
})

test_that("sum_area_m2 is identical for area_correction TRUE and FALSE", {
  skip_if_not_installed("sf")
  # sum_area_m2 covers the full coverage regardless of acceptance
  wkb <- create_test_wkb()
  r_no  <- grts_h3(
    wkb, n = 3, resolution = 4, area_correction = FALSE, seed = 1
  )
  r_yes <- grts_h3(
    wkb, n = 3, resolution = 4, area_correction = TRUE,  seed = 1
  )
  expect_equal(
    attr(r_no,  "sum_area_m2"),
    attr(r_yes, "sum_area_m2"),
    tolerance = 1e-9
  )
  expect_equal(attr(r_no, "n_cells"), attr(r_yes, "n_cells"))
})
