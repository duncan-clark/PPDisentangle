test_that("SEM packs named ETAS parameters in canonical order", {
  scrambled <- list(
    q = 1.6, gamma = 0, D = 1, p = 2.1,
    c = 0.4, alpha_m = 0.4, A = 0.2, mu = 0.8
  )
  packed <- PPDisentangle:::.etas_param_vector(scrambled)
  expect_identical(
    names(packed),
    c("mu", "A", "alpha_m", "c", "p", "D", "gamma", "q")
  )
  expect_equal(
    unname(packed),
    c(0.8, 0.2, 0.4, 0.4, 2.1, 1, 0, 1.6)
  )
  expect_error(
    PPDisentangle:::.etas_param_vector(scrambled[names(scrambled) != "A"]),
    "missing: A"
  )
})

test_that("univariate ETAS SEM scores match the conditional direct likelihood", {
  skip_if_not_installed("spatstat.geom")
  set.seed(20260811)
  n_pre <- 12L
  n_each_post <- 10L
  x <- c(
    runif(n_pre, 0.5, 4.5),
    runif(n_each_post, 0.5, 4.5),
    runif(n_each_post, 5.5, 9.5)
  )
  y <- runif(length(x), 0.5, 9.5)
  t <- c(
    sort(runif(n_pre, -5, -0.05)),
    sort(runif(2L * n_each_post, 0.05, 10))
  )
  win <- spatstat.geom::owin(c(0, 10), c(0, 10))
  partition <- spatstat.geom::tess(
    xgrid = c(0, 5, 10), ygrid = c(0, 10)
  )
  state_spaces <- list(
    control = spatstat.geom::as.owin(partition[c(TRUE, FALSE)]),
    treated = spatstat.geom::as.owin(partition[c(FALSE, TRUE)])
  )
  labels <- ifelse(x < 5 | t < 0, "control", "treated")
  pp <- data.frame(
    x = x, y = y, t = t,
    mag = 2.5 + rexp(length(x), 1.8),
    W = 1,
    process = labels,
    location_process = labels,
    inferred_process = labels,
    stringsAsFactors = FALSE
  )
  params <- list(
    mu = 0.8, A = 0.2, alpha_m = 0.4, c = 0.4,
    p = 2.1, D = 1, gamma = 0, q = 1.6
  )

  res <- em_style_labelling(
    pp_data = pp,
    partition = partition,
    partition_processes = c("control", "treated"),
    statespace = win,
    state_spaces = state_spaces,
    time_window = c(0, 10),
    treatment_time = 0,
    hawkes_params_control = params,
    hawkes_params_treated = params,
    update_control_params = TRUE,
    param_update_cadence = NULL,
    proposal_update_cadence = 1,
    update_starting_data = TRUE,
    include_starting_data = TRUE,
    include_starting_first_n = 1L,
    optim_method = "max",
    metric_name = "post_likelihood",
    iter = 1L,
    n_props = 2L,
    change_factor = 0.2,
    proposal_method = "single_flip",
    verbose = FALSE,
    model_type = "etas",
    etas_use_filtration_history = TRUE,
    m0 = 2.5,
    beta_gr = 1.8,
    max_branching_ratio = 0.98,
    t_trunc = 5
  )

  proposal_scores <- unlist(res$all_metrics, use.names = FALSE)
  expect_true(any(is.finite(proposal_scores)))
  expect_true(all(proposal_scores[is.finite(proposal_scores)] > -1e14))
  expect_gt(length(unique(round(proposal_scores[is.finite(proposal_scores)], 8))), 1L)

  selected <- res$labelling
  post <- selected[selected$t >= 0, , drop = FALSE]
  history <- selected[
    selected$t < 0 & selected$inferred_process == "control",
    c("x", "y", "t", "mag"),
    drop = FALSE
  ]
  control_post <- post[post$inferred_process == "control", , drop = FALSE]
  treated_post <- post[post$inferred_process == "treated", , drop = FALSE]
  expected <- loglik_etas(
    params, control_post, c(0, 10), win,
    m0 = 2.5, beta_gr = 1.8, max_branching_ratio = 0.98,
    zero_background_region = state_spaces$treated,
    t_trunc = 5, history = history
  ) + loglik_etas(
    params, treated_post, c(0, 10), win,
    m0 = 2.5, beta_gr = 1.8, max_branching_ratio = 0.98,
    zero_background_region = state_spaces$control,
    t_trunc = 5,
    history = history[0, , drop = FALSE]
  )
  expect_equal(res$metrics[[1]], expected, tolerance = 1e-8)
})

test_that("adaptive univariate ETAS keeps conditional SEM weights valid", {
  skip_if_not_installed("spatstat.geom")
  set.seed(20260812)
  n_pre <- 8L
  n_each_post <- 6L
  x <- c(
    runif(n_pre, 0.5, 4.5),
    runif(n_each_post, 0.5, 4.5),
    runif(n_each_post, 5.5, 9.5)
  )
  t <- c(
    sort(runif(n_pre, -4, -0.05)),
    sort(runif(2L * n_each_post, 0.05, 8))
  )
  labels <- ifelse(x < 5 | t < 0, "control", "treated")
  pp <- data.frame(
    x = x,
    y = runif(length(x), 0.5, 9.5),
    t = t,
    mag = 2.5 + rexp(length(x), 1.8),
    W = 1,
    process = labels,
    location_process = labels,
    inferred_process = labels,
    stringsAsFactors = FALSE
  )
  win <- spatstat.geom::owin(c(0, 10), c(0, 10))
  partition <- spatstat.geom::tess(
    xgrid = c(0, 5, 10), ygrid = c(0, 10)
  )
  state_spaces <- list(
    control = spatstat.geom::as.owin(partition[c(TRUE, FALSE)]),
    treated = spatstat.geom::as.owin(partition[c(FALSE, TRUE)])
  )
  params <- list(
    mu = 0.8, A = 0.2, alpha_m = 0.4, c = 0.4,
    p = 2.1, D = 1, gamma = 0, q = 1.6
  )

  fit <- adaptive_SEM(
    pp_data = pp,
    partition = partition,
    partition_processes = c("control", "treated"),
    statespace = win,
    time_window = c(-4, 8),
    treatment_time = 0,
    hawkes_params_control = params,
    hawkes_params_treated = params,
    N_labellings = 0L,
    N_iter = 1L,
    reset = FALSE,
    verbose = FALSE,
    adaptive_control = list(
      update_control_params = TRUE,
      param_update_cadence = 1000L,
      proposal_update_cadence = 1L,
      state_spaces = state_spaces,
      iter = 1L,
      n_props = 1L,
      change_factor = 0.2,
      verbose = FALSE,
      update_starting_data = TRUE,
      include_starting_data = TRUE,
      include_starting_first_n = 1L,
      proposal_method = "single_flip",
      fixed_params = params,
      outer_maxit = 5L
    ),
    model_type = "etas",
    etas_use_filtration_history = TRUE,
    m0 = 2.5,
    beta_gr = 1.8,
    max_branching_ratio = 0.98,
    t_trunc = 5
  )

  adaptive_scores <- unlist(fit$adaptive$all_metrics, use.names = FALSE)
  expect_true(all(adaptive_scores[is.finite(adaptive_scores)] > -1e14))
  expect_equal(fit$weights, 1)
  expect_equal(
    unname(unlist(fit$hawkes_params_control)[names(params)]),
    unname(unlist(params)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(unlist(fit$hawkes_params_treated)[names(params)]),
    unname(unlist(params)),
    tolerance = 1e-12
  )
})
