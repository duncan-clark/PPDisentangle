test_that("batched bivariate ETAS loglik matches per-labelling calls", {
  skip_if_not_installed("spatstat.geom")
  set.seed(20260810)
  n <- 40
  t <- sort(runif(n, 0, 20))
  x <- runif(n, 0, 10)
  y <- runif(n, 0, 10)
  mag <- 2.5 + rexp(n, 1.5)
  windowT <- c(0, 20)
  windowS <- spatstat.geom::owin(c(0, 10), c(0, 10))
  ctrl <- spatstat.geom::owin(c(0, 5), c(0, 10))
  treat <- spatstat.geom::owin(c(5, 10), c(0, 10))
  W0 <- as.numeric(!spatstat.geom::inside.owin(x, y, treat))
  W1 <- as.numeric(!spatstat.geom::inside.owin(x, y, ctrl))
  area0 <- spatstat.geom::area(ctrl)
  area1 <- spatstat.geom::area(treat)

  params <- c(
    mu_0 = 0.8, mu_1 = 0.6,
    A_00 = 0.2, alpha_m_00 = 0.5,
    A_11 = 0.25, alpha_m_11 = 0.6,
    A_01 = 0.02, alpha_m_01 = 0.4,
    A_10 = 0.02, alpha_m_10 = 0.4,
    c = 0.5, p = 2.1, D = 1.0, gamma = 0, q = 1.6
  )
  K <- 5L
  pid_mat <- matrix(rbinom(n * K, 1, 0.5), n, K)
  W0s <- matrix(W0, n, K)
  W1s <- matrix(W1, n, K)

  batch <- loglik_etas_bivariate_batch(
    params = params, t = t, x = x, y = y, mag = mag,
    process_ids = pid_mat, W0s = W0s, W1s = W1s,
    areaS_0 = area0, areaS_1 = area1, windowT = windowT,
    t_already_shifted = FALSE, m0 = 2.5, beta_gr = 1.5,
    max_branching_radius = 0.98, t_trunc = 5
  )
  single <- vapply(seq_len(K), function(k) {
    realiz <- data.frame(
      x = x, y = y, t = t, mag = mag,
      inferred_process = ifelse(pid_mat[, k] == 1L, "treated", "control")
    )
    loglik_etas_bivariate(
      params = params, realiz = realiz, windowT = windowT, windowS = windowS,
      control_state_space = ctrl, treated_state_space = treat,
      precomp = list(W_0 = W0, W_1 = W1, areaS_0 = area0, areaS_1 = area1),
      m0 = 2.5, beta_gr = 1.5, max_branching_radius = 0.98, t_trunc = 5
    )
  }, 0)

  expect_equal(batch, single, tolerance = 1e-9)
})

test_that("threaded K=1 batch likelihood matches scalar likelihood", {
  skip_if_not_installed("spatstat.geom")
  set.seed(20260811)
  n <- 160L
  t <- sort(runif(n, -10, 20))
  x <- runif(n, 0, 10)
  y <- runif(n, 0, 10)
  mag <- 2.5 + rexp(n, 1.5)
  windowT <- c(-10, 20)
  windowS <- spatstat.geom::owin(c(0, 10), c(0, 10))
  ctrl <- spatstat.geom::owin(c(0, 5), c(0, 10))
  treat <- spatstat.geom::owin(c(5, 10), c(0, 10))
  pid <- as.integer(x >= 5 & t >= 0)
  W0 <- as.numeric(!(x >= 5 & t >= 0))
  W1 <- as.numeric(x >= 5 & t >= 0)
  params <- c(
    mu_0 = 0.8, mu_1 = 0.6,
    A_00 = 0.2, alpha_m_00 = 0.5,
    A_11 = 0.25, alpha_m_11 = 0.6,
    A_01 = 0.02, alpha_m_01 = 0.4,
    A_10 = 0.02, alpha_m_10 = 0.4,
    c = 0.5, p = 2.1, D = 1.0, gamma = 0, q = 1.6
  )
  realiz <- data.frame(
    x = x, y = y, t = t, mag = mag,
    inferred_process = ifelse(pid == 1L, "treated", "control")
  )
  precomp <- list(W_0 = W0, W_1 = W1, areaS_0 = 50, areaS_1 = 50)
  scalar <- loglik_etas_bivariate(
    params, realiz, windowT, windowS,
    control_state_space = ctrl, treated_state_space = treat,
    precomp = precomp, m0 = 2.5, beta_gr = 1.5,
    max_branching_radius = 0.98, t_trunc = 5
  )
  batch_args <- list(
    params = params, t = t - windowT[1], x = x, y = y, mag = mag,
    process_ids = matrix(pid, ncol = 1),
    W0s = matrix(W0, ncol = 1), W1s = matrix(W1, ncol = 1),
    areaS_0 = 50, areaS_1 = 50, t_max = diff(windowT),
    windowT = windowT, t_already_shifted = TRUE,
    m0 = 2.5, beta_gr = 1.5,
    max_branching_radius = 0.98, t_trunc = 5
  )
  serial_batch <- do.call(
    loglik_etas_bivariate_batch,
    c(batch_args, list(n_threads = 1L))
  )
  threaded_batch <- do.call(
    loglik_etas_bivariate_batch,
    c(batch_args, list(n_threads = 4L))
  )

  expect_equal(serial_batch[[1]], scalar, tolerance = 1e-9)
  expect_equal(threaded_batch, serial_batch, tolerance = 1e-9)
})

test_that("closed-form bivariate spectral radius matches eigen", {
  pv <- c(
    mu_0 = 1, mu_1 = 1,
    A_00 = 0.3, alpha_m_00 = 0.5,
    A_11 = 0.25, alpha_m_11 = 0.6,
    A_01 = 0.05, alpha_m_01 = 0.4,
    A_10 = 0.04, alpha_m_10 = 0.3,
    c = 1, p = 2.1, D = 1, gamma = 0, q = 1.6
  )
  beta <- 2.0
  eta <- function(A, a) A * beta / (beta - a)
  M <- matrix(c(
    eta(pv[["A_00"]], pv[["alpha_m_00"]]),
    eta(pv[["A_01"]], pv[["alpha_m_01"]]),
    eta(pv[["A_10"]], pv[["alpha_m_10"]]),
    eta(pv[["A_11"]], pv[["alpha_m_11"]])
  ), 2, byrow = TRUE)
  expect_equal(
    PPDisentangle:::.etas_biv_spectral_radius(pv, beta),
    max(Re(eigen(M, only.values = TRUE)$values)),
    tolerance = 1e-12
  )
})
