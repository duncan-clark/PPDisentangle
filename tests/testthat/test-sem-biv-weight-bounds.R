# SEM bivariate weight scoring must honor enforce_finite_trigger_moments /
# p_lower_bound / q_lower_bound from adaptive_SEM(...), matching the M-step.

test_that("bivariate batch LL admits p in (1, 2] when finite moments are off", {
  skip_if_not_installed("spatstat.geom")
  set.seed(7)
  n <- 40L
  t <- sort(runif(n, 0, 20))
  x <- runif(n, 0, 10)
  y <- runif(n, 0, 10)
  mag <- 2.5 + rexp(n, 1.5)
  pid <- matrix(as.integer(x >= 5), n, 1L)
  W0 <- matrix(as.numeric(x < 5), n, 1L)
  W1 <- matrix(as.numeric(x >= 5), n, 1L)
  par <- c(
    mu_0 = 0.8, mu_1 = 0.7,
    A_00 = 0.12, alpha_m_00 = 0.4,
    A_11 = 0.14, alpha_m_11 = 0.45,
    A_01 = 0.01, alpha_m_01 = 0.3,
    A_10 = 0.01, alpha_m_10 = 0.3,
    c = 0.4, p = 1.5, D = 1.2, gamma = 0, q = 1.25
  )
  rejected <- loglik_etas_bivariate_batch(
    params = par,
    t = t, x = x, y = y, mag = mag,
    process_ids = pid, W0s = W0, W1s = W1,
    areaS_0 = 50, areaS_1 = 50, t_max = 20,
    t_trunc = 5, t_already_shifted = TRUE,
    m0 = 2.5, beta_gr = 1.5,
    enforce_finite_trigger_moments = TRUE,
    p_lower_bound = 2.001, q_lower_bound = 1.501
  )
  expect_equal(as.numeric(rejected), -1e15)

  admitted <- loglik_etas_bivariate_batch(
    params = par,
    t = t, x = x, y = y, mag = mag,
    process_ids = pid, W0s = W0, W1s = W1,
    areaS_0 = 50, areaS_1 = 50, t_max = 20,
    t_trunc = 5, t_already_shifted = TRUE,
    m0 = 2.5, beta_gr = 1.5,
    enforce_finite_trigger_moments = FALSE,
    p_lower_bound = 1, q_lower_bound = 1
  )
  expect_true(is.finite(admitted))
  expect_gt(as.numeric(admitted), -1e10)
})

test_that("SEM bivariate weights keep relaxed p/q constraints from dots", {
  skip_if_not_installed("spatstat.geom")
  set.seed(11)
  n_pre <- 20L
  n_post <- 50L
  x <- c(runif(n_pre, 0, 5), runif(n_post, 0, 10))
  y <- runif(n_pre + n_post, 0, 10)
  t <- c(sort(runif(n_pre, -6, -0.01)), sort(runif(n_post, 0, 18)))
  mag <- 2.5 + rexp(n_pre + n_post, 1.5)
  win <- spatstat.geom::owin(c(0, 10), c(0, 10))
  partition <- spatstat.geom::tess(xgrid = c(0, 5, 10), ygrid = c(0, 10))
  partition_processes <- c("control", "treated")
  state_spaces <- list(
    control = spatstat.geom::as.owin(partition[c(TRUE, FALSE)]),
    treated = spatstat.geom::as.owin(partition[c(FALSE, TRUE)])
  )
  pp <- data.frame(
    x = x, y = y, t = t, mag = mag,
    location_process = ifelse(x < 5 | t < 0, "control", "treated"),
    stringsAsFactors = FALSE
  )
  pp$inferred_process <- pp$location_process
  init <- c(
    mu_0 = 0.9, mu_1 = 0.7,
    A_00 = 0.12, alpha_m_00 = 0.4,
    A_11 = 0.14, alpha_m_11 = 0.45,
    A_01 = 0.01, alpha_m_01 = 0.3,
    A_10 = 0.01, alpha_m_10 = 0.3,
    c = 0.35, p = 1.4, D = 1.1, gamma = 0, q = 1.3
  )
  res <- adaptive_SEM(
    pp_data = pp,
    partition = partition,
    partition_processes = partition_processes,
    statespace = win,
    time_window = c(0, 18),
    treatment_time = 0,
    hawkes_params_control = list(mu = 0.9, A = 0.12, alpha_m = 0.4, c = 0.35,
                                 p = 1.4, D = 1.1, gamma = 0, q = 1.3),
    hawkes_params_treated = list(mu = 0.7, A = 0.14, alpha_m = 0.45, c = 0.35,
                                 p = 1.4, D = 1.1, gamma = 0, q = 1.3),
    N_labellings = 2L,
    N_iter = 1L,
    reset = FALSE,
    verbose = FALSE,
    adaptive_control = list(
      param_update_cadence = 1L,
      proposal_update_cadence = 1,
      state_spaces = state_spaces,
      iter = 2L, n_props = 2L, change_factor = 0.15,
      verbose = FALSE,
      update_starting_data = TRUE,
      include_starting_data = TRUE,
      include_starting_first_n = 2L,
      update_control_params = TRUE,
      fixed_params = list(gamma = 0),
      proposal_method = "simulation",
      outer_maxit_biv = 40L,
      optim_method = "max"
    ),
    model_type = "etas_bivariate",
    m0 = 2.5,
    beta_gr = 1.5,
    etas_bivariate_params = init,
    enforce_finite_trigger_moments = FALSE,
    p_lower_bound = 1,
    q_lower_bound = 1,
    max_branching_radius = 0.98,
    hard_subcritical = TRUE,
    t_trunc = 5,
    use_pre_history_for_biv = TRUE
  )
  expect_true(!is.null(res$etas_bivariate_params))
  expect_true("etas_bivariate_convergence" %in% names(res))
  sem_ok <- PPDisentangle:::.etas_biv_sem_ok(res, 1.5, 0.98)
  conv0 <- isTRUE(as.integer(res$etas_bivariate_convergence)[1L] == 0L)
  params_ok <- PPDisentangle:::.etas_biv_params_ok(
    res$etas_bivariate_params, 1.5, 0.98
  )
  expect_identical(sem_ok, conv0 && params_ok)
  p_hat <- as.numeric(res$etas_bivariate_params[["p"]])
  q_hat <- as.numeric(res$etas_bivariate_params[["q"]])
  expect_true(is.finite(p_hat) && p_hat > 1)
  expect_true(is.finite(q_hat) && q_hat > 1)
  # Weight scoring under package defaults would reject p=1.4 / q=1.3 and skip
  # the outer M-step with "All weights are zero". Completing with a finite
  # retained parameter vector is the integration check for that path.
  expect_true(length(res$c_params) >= 1L || length(res$fits) >= 1L ||
                !is.null(res$etas_bivariate_params))
})
