# Inner-loop optimizations: hoisted bivariate metric precompute (L1),
# threaded metric kernel (L2), and cached nncross proposal weights (S1).

.make_inner_workload <- function() {
  set.seed(31415)
  n <- 120
  n_pre <- 40
  x <- c(runif(n_pre, 0, 5), runif(n, 0, 10))       # pre events on control side
  y <- runif(n_pre + n, 0, 10)
  t <- c(sort(runif(n_pre, -10, -0.01)), sort(runif(n, 0, 30)))
  mag <- 2.5 + rexp(n_pre + n, 1.5)
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
  pp$process <- pp$location_process
  pp$inferred_process <- pp$location_process
  params <- c(
    mu_0 = 0.8, mu_1 = 0.6,
    A_00 = 0.2, alpha_m_00 = 0.5,
    A_11 = 0.25, alpha_m_11 = 0.6,
    A_01 = 0.02, alpha_m_01 = 0.4,
    A_10 = 0.02, alpha_m_10 = 0.4,
    c = 0.5, p = 2.1, D = 1.0, gamma = 0, q = 1.6
  )
  list(pp = pp, win = win, partition = partition,
       partition_processes = partition_processes,
       state_spaces = state_spaces, params = params,
       windowT = c(0, 30))
}

.run_inner_small <- function(w, n_threads = 1L, iter = 3L, seed = 2718,
                             verbose = FALSE) {
  set.seed(seed)
  em_style_labelling(
    pp_data = w$pp,
    partition = w$partition,
    partition_processes = w$partition_processes,
    statespace = w$win,
    state_spaces = w$state_spaces,
    time_window = w$windowT,
    treatment_time = 0,
    hawkes_params_control = list(mu = 0.8, A = 0.2, alpha_m = 0.5, c = 0.5,
                                 p = 2.1, D = 1.0, gamma = 0, q = 1.6),
    hawkes_params_treated = list(mu = 0.6, A = 0.25, alpha_m = 0.6, c = 0.5,
                                 p = 2.1, D = 1.0, gamma = 0, q = 1.6),
    update_control_params = TRUE,
    param_update_cadence = 100000L,
    proposal_update_cadence = 1,
    update_starting_data = TRUE,
    include_starting_data = TRUE,
    include_starting_first_n = iter,
    optim_method = "max",
    metric_name = "post_likelihood",
    iter = iter, n_props = 4L, change_factor = 0.5,
    stagnation_trigger_every = 50,
    proposal_method = "simulation",
    fixed_params = as.list(w$params),   # freeze all params: no refits
    verbose = verbose,
    model_type = "etas_bivariate",
    temporal_weight = 0,
    biv_n_threads = n_threads,
    m0 = 2.5, beta_gr = 1.5,
    etas_bivariate_params = w$params,
    hard_subcritical = TRUE,
    max_branching_radius = 0.98,
    alpha_beta_gap_min = 1e-4,
    t_trunc = 5
  )
}

test_that("hoisted bivariate metric matches a direct batch likelihood", {
  skip_if_not_installed("spatstat.geom")
  w <- .make_inner_workload()
  res <- .run_inner_small(w)

  final <- res$labelling
  final <- final[order(final$t), , drop = FALSE]
  # The inner-loop metric conditions on the pre-treatment history:
  # times are shifted to the history window start.
  history_window <- c(min(final$t), w$windowT[2])
  pid <- matrix(as.integer(final$inferred_process == "treated"), ncol = 1)
  W0 <- as.numeric(!spatstat.geom::inside.owin(final$x, final$y,
                                               w$state_spaces$treated))
  W1 <- as.numeric(!spatstat.geom::inside.owin(final$x, final$y,
                                               w$state_spaces$control))
  expected <- loglik_etas_bivariate_batch(
    params = res$etas_bivariate_params,
    t = final$t - history_window[1],
    x = final$x, y = final$y, mag = final$mag,
    process_ids = pid,
    W0s = matrix(W0, ncol = 1), W1s = matrix(W1, ncol = 1),
    areaS_0 = spatstat.geom::area(w$state_spaces$control),
    areaS_1 = spatstat.geom::area(w$state_spaces$treated),
    t_max = diff(history_window), windowT = history_window,
    t_already_shifted = TRUE, t_trunc = 5,
    m0 = 2.5, beta_gr = 1.5,
    alpha_beta_gap_min = 1e-4, max_branching_radius = 0.98,
    alpha_m_lower_bound = 0
  )
  # metrics[iter] is the score of the selected (returned) labelling
  expect_equal(res$metrics[length(res$metrics)], expected[1], tolerance = 1e-8)
  expect_length(res$class_results, sum(lengths(res$all_metrics)))
})

test_that("biv_n_threads does not change the inner-loop results", {
  skip_if_not_installed("spatstat.geom")
  w <- .make_inner_workload()
  r1 <- .run_inner_small(w, n_threads = 1L)
  r4 <- .run_inner_small(w, n_threads = 4L)
  expect_equal(r1$metrics, r4$metrics, tolerance = 1e-8)
  expect_identical(r1$labelling$inferred_process, r4$labelling$inferred_process)
})

test_that("SEM verbosity does not change stochastic results", {
  skip_if_not_installed("spatstat.geom")
  w <- .make_inner_workload()
  quiet <- .run_inner_small(w, verbose = FALSE)
  traced <- capture.output(
    loud <- .run_inner_small(w, verbose = TRUE)
  )

  expect_gt(length(traced), 0L)
  expect_equal(quiet$metrics, loud$metrics, tolerance = 1e-10)
  expect_identical(
    quiet$labelling$inferred_process,
    loud$labelling$inferred_process
  )
  expect_equal(
    quiet$etas_bivariate_params,
    loud$etas_bivariate_params,
    tolerance = 1e-10
  )
  expect_identical(quiet$class_results, loud$class_results)
})

test_that("cached nncross weights leave proposals bitwise identical", {
  skip_if_not_installed("spatstat.geom")
  w <- .make_inner_workload()
  sim_prop <- PPDisentangle:::simulation_labeling_hawkes_hawkes_fast
  args <- list(
    w$pp, partition = w$partition, partition_process = w$partition_processes,
    statespace = w$win, state_spaces = w$state_spaces, windowT = w$windowT,
    hawkes_params_control = list(mu = 0.8), hawkes_params_treated = list(mu = 0.6),
    change_factor = 1.0, filtration = w$pp[0, , drop = FALSE],
    verbose = FALSE, temporal_weight = 0,
    model_type = "etas_bivariate",
    return_proposal_sim_cache = TRUE,
    m0 = 2.5, beta_gr = 1.5, etas_bivariate_params = w$params, t_trunc = 5
  )
  set.seed(777)
  first <- do.call(sim_prop, c(args, list(proposal_sim_cache = NULL)))
  cache_full <- first$proposal_sim_cache

  # A later proposal in the same batch must compute the identical weights,
  # whether it recomputes nncross (cache stripped) or reuses the cache.
  cache_nogeom <- cache_full
  cache_nogeom$geom_weights <- NULL
  set.seed(999)
  recompute <- do.call(sim_prop, c(args, list(proposal_sim_cache = cache_nogeom)))
  set.seed(999)
  reuse <- do.call(sim_prop, c(args, list(proposal_sim_cache = cache_full)))

  expect_identical(recompute$data, reuse$data)
  # the recomputing call must also have produced/stored the weights
  if (!is.null(cache_full$geom_weights)) {
    expect_identical(recompute$proposal_sim_cache$geom_weights,
                     cache_full$geom_weights)
  }
})
