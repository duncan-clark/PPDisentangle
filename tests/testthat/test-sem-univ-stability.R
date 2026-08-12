# Univariate SEM M-steps should use the soft-barrier / interior-init stack.

test_that("em_style_labelling univ M-step stays inside the eta margin from explosive starts", {
  skip_if_not_installed("spatstat.geom")
  set.seed(7)
  n_pre <- 24L
  n_post <- 60L
  x <- c(runif(n_pre, 0, 5), runif(n_post, 0, 10))
  y <- runif(n_pre + n_post, 0, 10)
  t <- c(sort(runif(n_pre, -6, -0.01)), sort(runif(n_post, 0, 20)))
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
  # Explosive productivity start (eta >> 0.98 under beta_gr=2.3).
  init_ctrl <- list(
    mu = 0.8, A = 1.1, alpha_m = 0.6,
    c = 0.4, p = 2.05, D = 1.2, gamma = 0, q = 1.55
  )
  init_treat <- list(
    mu = 0.9, A = 1.2, alpha_m = 0.7,
    c = 0.4, p = 2.05, D = 1.2, gamma = 0, q = 1.55
  )
  expect_gt(PPDisentangle:::.etas_univ_branching_ratio(unlist(init_treat), 2.3), 0.98)

  res <- em_style_labelling(
    pp_data = pp,
    partition = partition,
    partition_processes = partition_processes,
    statespace = win,
    state_spaces = state_spaces,
    time_window = c(0, 20),
    treatment_time = 0,
    hawkes_params_control = init_ctrl,
    hawkes_params_treated = init_treat,
    update_control_params = TRUE,
    param_update_cadence = 1L,
    proposal_update_cadence = 1,
    update_starting_data = TRUE,
    include_starting_data = TRUE,
    include_starting_first_n = 2L,
    optim_method = "max",
    metric_name = "post_likelihood",
    iter = 2L, n_props = 2L, change_factor = 0.2,
    stagnation_trigger_every = 50,
    proposal_method = "simulation",
    fixed_params = list(gamma = 0),
    verbose = FALSE,
    model_type = "etas",
    outer_maxit = 60L,
    m0 = 2.5, beta_gr = 2.3,
    hard_subcritical = TRUE,
    soft_branching_barrier = TRUE,
    log_transform = TRUE,
    init_branching_margin = 0.9,
    max_branching_ratio = 0.98,
    alpha_beta_gap_min = 1e-4,
    t_trunc = 5
  )

  expect_true(length(res$treated_par) >= 2L)
  expect_true(length(res$control_par) >= 2L)
  eta_t <- PPDisentangle:::.etas_univ_branching_ratio(
    unlist(res$treated_par[[length(res$treated_par)]]), 2.3
  )
  eta_c <- PPDisentangle:::.etas_univ_branching_ratio(
    unlist(res$control_par[[length(res$control_par)]]), 2.3
  )
  expect_lt(eta_t, 0.98)
  expect_lt(eta_c, 0.98)
  expect_true(is.finite(eta_t) && is.finite(eta_c))
  last_fit <- res$fits[[length(res$fits)]]
  expect_true(!is.null(last_fit$branching_ratio))
  expect_lt(last_fit$branching_ratio, 0.98)
})
