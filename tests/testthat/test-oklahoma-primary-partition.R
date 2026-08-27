test_that("Oklahoma primary tessellation defaults to grid_1.0R and skips itself in sensitivity", {
  geom_candidates <- c(
    system.file("oklahoma", "oklahoma_geometry.R", package = "PPDisentangle"),
    file.path("inst", "oklahoma", "oklahoma_geometry.R"),
    file.path("..", "..", "inst", "oklahoma", "oklahoma_geometry.R")
  )
  geom_file <- geom_candidates[nzchar(geom_candidates) & file.exists(geom_candidates)][1]
  if (is.na(geom_file) || !nzchar(geom_file) || !file.exists(geom_file)) {
    skip("oklahoma_geometry.R is not available in this check")
  }
  sys.source(geom_file, envir = environment())

  expect_equal(oklahoma_normalize_primary_partition(""), "grid_1.0R")
  expect_equal(oklahoma_normalize_primary_partition("grid_1.0R"), "grid_1.0R")
  expect_equal(oklahoma_normalize_primary_partition("1r"), "grid_1.0R")
  expect_equal(oklahoma_normalize_primary_partition("county"), "county")
  expect_equal(
    oklahoma_normalize_primary_partition("grid_1.0R", quick_check = TRUE),
    "grid_coarse"
  )
  expect_error(oklahoma_normalize_primary_partition("grid_2.0R"))

  expect_equal(
    oklahoma_sensitivity_partition_labels("grid_1.0R"),
    c("county", "grid_2.0R", "grid_5.0R", "aoi_region")
  )
  expect_equal(
    oklahoma_sensitivity_partition_labels("county"),
    c("grid_1.0R", "grid_2.0R", "grid_5.0R", "aoi_region")
  )
  expect_equal(
    oklahoma_sensitivity_partition_labels("grid_coarse"),
    c("county", "aoi_region")
  )
  expect_false("grid_1.0R" %in% oklahoma_sensitivity_partition_labels("grid_1.0R"))
})
