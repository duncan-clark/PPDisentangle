# Bivariate SEM must retain param updates in the local em_style_labelling scope.

test_that("em_style_labelling updates etas_bivariate_params after cadence refits", {
  skip_if_not_installed("spatstat.geom")
  set.seed(42)
  n_pre <- 30L
  n_post <- 80L
  x <- c(runif(n_pre, 0, 5), runif(n_post, 0, 10))
  y <- runif(n_pre + n_post, 0, 10)
  t <- c(sort(runif(n_pre, -8, -0.01)), sort(runif(n_post, 0, 25)))
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
    A_00 = 0.15, alpha_m_00 = 0.4,
    A_11 = 0.18, alpha_m_11 = 0.45,
    A_01 = 0.01, alpha_m_01 = 0.3,
    A_10 = 0.01, alpha_m_10 = 0.3,
    c = 0.4, p = 2.05, D = 1.2, gamma = 0, q = 1.55
  )
  res <- em_style_labelling(
    pp_data = pp,
    partition = partition,
    partition_processes = partition_processes,
    statespace = win,
    state_spaces = state_spaces,
    time_window = c(0, 25),
    treatment_time = 0,
    hawkes_params_control = list(mu = 0.9, A = 0.15, alpha_m = 0.4, c = 0.4,
                                 p = 2.05, D = 1.2, gamma = 0, q = 1.55),
    hawkes_params_treated = list(mu = 0.7, A = 0.18, alpha_m = 0.45, c = 0.4,
                                 p = 2.05, D = 1.2, gamma = 0, q = 1.55),
    update_control_params = TRUE,
    param_update_cadence = 1L,
    proposal_update_cadence = 1,
    update_starting_data = TRUE,
    include_starting_data = TRUE,
    include_starting_first_n = 3L,
    optim_method = "max",
    metric_name = "post_likelihood",
    iter = 3L, n_props = 3L, change_factor = 0.2,
    stagnation_trigger_every = 50,
    proposal_method = "simulation",
    fixed_params = list(gamma = 0),
    verbose = FALSE,
    model_type = "etas_bivariate",
    outer_maxit_biv = 80L,
    m0 = 2.5, beta_gr = 1.5,
    etas_bivariate_params = init,
    hard_subcritical = TRUE,
    max_branching_radius = 0.98,
    alpha_beta_gap_min = 1e-4,
    t_trunc = 5
  )
  expect_true(!is.null(res$etas_bivariate_params))
  expect_true(length(res$fits) >= 1L)
  expect_equal(
    as.integer(res$etas_bivariate_convergence),
    as.integer(res$fits[[length(res$fits)]]$convergence)
  )
  last_fit_par <- res$fits[[length(res$fits)]]$par
  expect_equal(
    unname(as.numeric(res$etas_bivariate_params[names(last_fit_par)])),
    unname(as.numeric(last_fit_par)),
    tolerance = 1e-10
  )
  # The retained params must not remain exactly at initialization after refits.
  expect_false(isTRUE(all.equal(
    unname(as.numeric(res$etas_bivariate_params[names(init)])),
    unname(as.numeric(init)),
    tolerance = 0
  )))
})
