# test-grts_h3-area-correction.R
# Tests specific to the area_correction = TRUE path
# (sequential poisson sampling).
# Tests that are shared between modes (ip formula, reproducibility) live in
# test-grts_h3-ip.R and test-grts_h3-reproducibility.R respectively.

test_that("area_correction = TRUE executes and returns a data.frame", {
  skip_if_not_installed("sf")
  out <- grts_h3(
    wkb = create_test_wkb(), n = 10, resolution = 5, area_correction = TRUE,
    seed = 42
  )
  expect_s3_class(out, "data.frame")
  expect_true(attr(out, "area_correction"))
  expect_true(is.numeric(out$ip))
  expect_all_false(out$ip == mean(out$ip))
})

test_that(
  "area_correction changes the sample relative to area_correction = FALSE", {
    skip_if_not_installed("sf")
    wkb  <- create_test_wkb()
    seed <- 42
    r_no  <- grts_h3(
      wkb, n = 100, resolution = 6, area_correction = FALSE, seed = seed
    )
    r_yes <- grts_h3(
      wkb, n = 100, resolution = 6, area_correction = TRUE,  seed = seed
    )
    expect_false(identical(r_no$cell, r_yes$cell))
  }
)

test_that("area_correction = TRUE: sort keys for shared cells are identical to FALSE path", { # nolint
  skip_if_not_installed("sf")
  # The sort key is a pure function of (cell, seed) — independent of
  # the acceptance step — so shared cells must carry the same sort_key.
  wkb  <- create_test_wkb()
  seed <- 42
  r_no  <- grts_h3(
    wkb, n = 5, resolution = 4, area_correction = FALSE, seed = seed
  )
  r_yes <- grts_h3(
    wkb, n = 5, resolution = 4, area_correction = TRUE,  seed = seed
  )
  shared <- intersect(r_no$cell, r_yes$cell)
  if (length(shared) == 0L) skip("No shared cells between modes")
  expect_equal(
    r_no$sort_key[match(shared, r_no$cell)],
    r_yes$sort_key[match(shared, r_yes$cell)]
  )
})
