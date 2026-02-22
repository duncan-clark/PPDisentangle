test_that("ATE_estim returns non-parametric estimate", {
  obs <- data.frame(
    x = runif(100), y = runif(100), t = runif(100),
    inferred_process = sample(c("control", "treated"), 100, replace = TRUE),
    location_process = sample(c("control", "treated"), 100, replace = TRUE)
  )
  Omega <- spatstat.geom::owin(xrange = c(0, 1), yrange = c(0, 1))
  partition <- spatstat.geom::quadrats(Omega, nx = 2, ny = 2)
  treated_parts <- tilenames(partition)[c(3, 4)]

  result <- ATE_estim(Omega, partition, obs, treated_parts, windowT = c(0, 1))
  expect_true(is.numeric(result$ATE_non_param))
  expect_true(is.finite(result$ATE_non_param))
})

test_that("generate_poisson_data returns valid structure", {
  Omega <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  partition <- spatstat.geom::quadrats(Omega, nx = 3, ny = 3)
  set.seed(42)
  result <- generate_poisson_data(Omega, partition, treat_prop = 0.5, lambda_1 = 50, lambda_2 = 80)

  expect_true(nrow(result$obs_data) > 0)
  expect_true("process" %in% names(result$obs_data))
  expect_true("location_process" %in% names(result$obs_data))
  expect_true(all(result$obs_data$process %in% c("control", "treated")))
})
