source(testthat::test_path("..", "..", "inst", "sim_study",
                          "structured_robustness_utils.R"))

test_that("effect-modification design is fixed and exactly balanced", {
  omega <- spatstat.geom::owin(c(0, 100), c(0, 100))
  partition <- spatstat.geom::quadrats(omega, nx = 10, ny = 10)
  design_a <- make_effect_modification_design(partition, seed = 101L)
  design_b <- make_effect_modification_design(partition, seed = 101L)

  expect_identical(design_a$z_obs, design_b$z_obs)
  expect_equal(sum(design_a$z_obs), 50L)
  expect_equal(sum(design_a$z_obs & design_a$grid$X == 1), 25L)
  expect_equal(sum(design_a$z_obs & design_a$grid$X == -1), 25L)
  expect_equal(mean(effect_multiplier(design_a$grid$X[design_a$z_obs], 0.6)), 1)
  expect_equal(
    effect_multiplier(1, 0.6) / effect_multiplier(-1, 0.6),
    exp(1.2)
  )
  expect_equal(design_a$grid$X[design_a$flip_plus_cell], 1)
  expect_equal(design_a$grid$X[design_a$flip_minus_cell], -1)
  expect_equal(sum(design_a$z_flip_plus), 1L)
  expect_equal(sum(design_a$z_flip_minus), 1L)
  expect_true(design_a$z_flip_plus[design_a$flip_plus_cell])
  expect_true(design_a$z_flip_minus[design_a$flip_minus_cell])
})

test_that("geometry path has the requested endpoints and invariants", {
  omega <- spatstat.geom::owin(c(0, 100), c(0, 100))
  partition <- spatstat.geom::quadrats(omega, nx = 10, ny = 10)
  design <- make_geometry_transport_design(
    partition, path_seed = 202L, observed_seed = 203L
  )

  expect_equal(nrow(design$edges), 180L)
  expect_equal(nrow(design$discrepant_pairs), 25L)
  expect_identical(unname(design$focal_cells), c(5L, 25L, 45L, 65L, 85L))
  expect_true(all(vapply(design$allocations, sum, integer(1L)) == 50L))
  expect_identical(design$allocations[[1L]], design$z_checkerboard)
  expect_identical(design$allocations[[6L]], design$z_block)
  expect_equal(unname(design$coarseness[c(1L, 6L)]), c(0, 1))
  expect_true(all(vapply(
    design$allocations,
    function(z) all(z[design$focal_cells]),
    logical(1L)
  )))
})
