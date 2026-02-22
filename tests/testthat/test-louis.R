test_that("complete_data_loglik returns finite for valid labelling", {
  Omega <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  partition <- spatstat.geom::quadrats(Omega, nx = 2, ny = 2)
  partition_processes <- c("control", "treated", "control", "treated")

  theta <- c(20, 0.5, 1, 0.2, 30, 0.5, 1, 0.3)

  set.seed(42)
  sim <- generate_inhomogeneous_hawkes(
    Omega, partition, c(0, 10), partition_processes,
    list(control = list(mu = 20, alpha = 0.5, beta = 1, K = 0.2),
         treated = list(mu = 30, alpha = 0.5, beta = 1, K = 0.3))
  )
  sim$inferred_process <- sim$process
  sim$W <- rep(1, nrow(sim))

  ll <- complete_data_loglik(theta, sim, treatment_time = 0,
                              statespace = Omega, partition = partition,
                              partition_processes = partition_processes)
  expect_true(is.finite(ll))
  expect_true(ll < 0)
})

test_that("complete_data_loglik penalizes bad parameters", {
  Omega <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  partition <- spatstat.geom::quadrats(Omega, nx = 2, ny = 2)
  partition_processes <- c("control", "treated", "control", "treated")

  set.seed(42)
  sim <- generate_inhomogeneous_hawkes(
    Omega, partition, c(0, 10), partition_processes,
    list(control = list(mu = 20, alpha = 0.5, beta = 1, K = 0.2),
         treated = list(mu = 30, alpha = 0.5, beta = 1, K = 0.3))
  )
  sim$inferred_process <- sim$process
  sim$W <- rep(1, nrow(sim))

  theta_good <- c(20, 0.5, 1, 0.2, 30, 0.5, 1, 0.3)
  theta_bad <- c(20, 0.5, 1, 0.2, 30, 0.5, 1, 1.5)

  ll_good <- complete_data_loglik(theta_good, sim, treatment_time = 0,
                                   statespace = Omega, partition = partition,
                                   partition_processes = partition_processes)
  ll_bad <- complete_data_loglik(theta_bad, sim, treatment_time = 0,
                                  statespace = Omega, partition = partition,
                                  partition_processes = partition_processes)
  expect_true(ll_good > ll_bad)
})
