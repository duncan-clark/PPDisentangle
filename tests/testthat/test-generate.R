test_that("generate_inhomogeneous_hawkes produces valid output", {
  Omega <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  partition <- spatstat.geom::quadrats(Omega, nx = 2, ny = 2)
  partition_processes <- c("control", "treated", "control", "treated")
  time_window <- c(0, 10)
  hawkes_params <- list(
    control = list(mu = 20, alpha = 0.5, beta = 1, K = 0.2),
    treated = list(mu = 30, alpha = 0.5, beta = 1, K = 0.3)
  )
  set.seed(42)
  result <- generate_inhomogeneous_hawkes(
    Omega, partition, time_window, partition_processes, hawkes_params
  )

  expect_true(nrow(result) > 0)
  expect_true("x" %in% names(result))
  expect_true("y" %in% names(result))
  expect_true("t" %in% names(result))
  expect_true("process" %in% names(result))
  expect_true(all(result$process %in% c("control", "treated")))
})

test_that("fit_hawkes returns sensible Poisson MLE", {
  params_init <- list(mu = 50, alpha = 1, beta = 1, K = 0)
  windowT <- c(0, 10)
  windowS <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  set.seed(42)
  sim <- sim_hawkes(list(mu = 50, alpha = 1, beta = 1, K = 0), windowT, windowS)
  realiz <- data.frame(x = sim$x, y = sim$y, t = sim$t)

  result <- fit_hawkes(params_init, realiz, windowT, windowS, maxit = 100, poisson_flag = TRUE)
  expect_true(!is.null(result$par))
  expect_equal(result$par$K, 0)
  expect_true(result$par$mu > 0)
})
