# test-output-structure.R
# Tests for the structure, column types, and attributes of sampler output.

# ---------------------------------------------------------------------------
# 1. Return type, columns
# ---------------------------------------------------------------------------

test_that("sampler returns a data.frame with correct columns", {
  skip_if_not_installed("sf")
  result <- grts_h3(create_test_wkb(), n = 5, resolution = 4, seed = 42)
  expect_s3_class(result, "data.frame")
  expect_equal(
    names(result),
    c("cell", "sort_key", "effective_key", "area_m2", "ip")
  )
  expect_type(result$cell, "character")
  expect_type(result$sort_key, "character")
  expect_type(result$area_m2, "double")
  expect_type(result$ip, "double")
  expect_type(attr(result, "n_cells"), "integer")
  expect_type(attr(result, "sum_area_m2"), "double")
  expect_type(attr(result, "n"), "integer")
  expect_type(attr(result, "resolution"), "integer")
  expect_type(attr(result, "containment"), "character")
  expect_type(attr(result, "area_correction"), "logical")
  expect_equal(attr(result, "n"), nrow(result))
  expect_true(all(grepl("^[0-9a-f]+$", result$cell)))
  # Canonical H3 format: no leading zeros (unless the value itself is zero)
  expect_false(any(grepl("^0[0-9a-f]", result$cell)))
  expect_true(all(nchar(result$sort_key) == 16L))
  expect_true(all(grepl("^[0-9a-f]{16}$", result$sort_key)))
  expect_equal(result$sort_key, sort(result$sort_key))
  expect_equal(length(unique(result$cell)), nrow(result))
  expect_equal(length(unique(result$sort_key)), nrow(result))

  result <- grts_sobol(create_test_wkb(), n = 5, seed = 42)
  expect_s3_class(result, "data.frame")
  expect_equal(
    names(result),
    c("lon", "lat", "grts_rank")
  )
  expect_type(result$lon, "double")
  expect_type(result$lat, "double")
  expect_type(result$grts_rank, "integer")
  expect_type(attr(result, "seed"), "double")
  expect_type(attr(result, "n"), "integer")
  expect_type(attr(result, "sobol_scanned"), "double")
  expect_type(attr(result, "acceptance_rate"), "double")
  expect_type(attr(result, "bbox"), "double")
  expect_equal(result$grts_rank, sort(result$grts_rank))
  expect_equal(length(unique(result$grts_rank)), nrow(result))
  expect_equal(attr(result, "n"), nrow(result))

  result <- bas_sobol(create_test_wkb(), n = 5, seed = 42)
  expect_s3_class(result, "data.frame")
  expect_equal(
    names(result),
    c("lon", "lat", "sobol_index")
  )
  expect_type(result$lon, "double")
  expect_type(result$lat, "double")
  expect_type(result$sobol_index, "integer")
  expect_type(attr(result, "n"), "integer")
  expect_type(attr(result, "seed"), "double")
  expect_type(attr(result, "n"), "integer")
  expect_type(attr(result, "sobol_scanned"), "double")
  expect_type(attr(result, "fill_ratio"), "double")
  expect_type(attr(result, "sequence_bbox"), "double")
  expect_equal(result$sobol_index, sort(result$sobol_index))
  expect_equal(length(unique(result$sobol_index)), nrow(result))
  expect_equal(attr(result, "n"), nrow(result))

})


test_that("attribute containment echoes the input argument", {
  skip_if_not_installed("sf")
  result <- grts_h3(create_test_wkb(), n = 5, resolution = 4,
                    containment = "centroid", seed = 42)
  expect_equal(attr(result, "containment"), "centroid")
})

test_that("area_correction attribute is FALSE by default", {
  skip_if_not_installed("sf")
  result <- grts_h3(create_test_wkb(), n = 5, resolution = 4, seed = 42)
  expect_false(attr(result, "area_correction"))
})

test_that("n_cells >= n", {
  skip_if_not_installed("sf")
  result <- grts_h3(create_test_wkb(), n = 5, resolution = 4, seed = 42)
  expect_gte(attr(result, "n_cells"), 5L)
})

test_that("sum_area_m2 is positive", {
  skip_if_not_installed("sf")
  result <- grts_h3(create_test_wkb(), n = 5, resolution = 4, seed = 42)
  expect_gt(attr(result, "sum_area_m2"), 0)
})
