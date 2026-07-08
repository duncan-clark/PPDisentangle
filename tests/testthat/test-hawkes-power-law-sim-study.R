make_power_law_sim_study_fixture <- function(seed = 20260706, sim_kernel = "power_law") {
  set.seed(seed)
  Omega <- spatstat.geom::owin(xrange = c(0, 20), yrange = c(0, 20))
  partition <- spatstat.geom::quadrats(Omega, nx = 2, ny = 2)
  partition_processes <- c("control", "treated", "treated", "control")
  treated_idx <- partition_processes == "treated"
  control_ss <- spatstat.geom::as.owin(partition[!treated_idx])
  treated_ss <- spatstat.geom::as.owin(partition[treated_idx])
  state_spaces <- list(control_ss, treated_ss)

  power_control <- list(mu = 12, alpha = 0.08, c = 0.35, p = 2.2,
                        beta = 0.35, K = 0.25, kernel = "power_law")
  power_treated <- list(mu = 12, alpha = 0.08, c = 0.35, p = 2.2,
                        beta = 0.35, K = 0.55, kernel = "power_law")
  exp_control <- list(mu = 12, alpha = 0.08, beta = 2.0, K = 0.25,
                      kernel = "exponential")
  exp_treated <- list(mu = 12, alpha = 0.08, beta = 2.0, K = 0.55,
                      kernel = "exponential")

  sim_control <- if (identical(sim_kernel, "power_law")) power_control else exp_control
  sim_treated <- if (identical(sim_kernel, "power_law")) power_treated else exp_treated

  treatment_time <- 2
  end_time <- 6
  pre <- sim_hawkes_fast(sim_control, c(0, treatment_time), Omega)
  pre_df <- data.frame(
    x = pre$x, y = pre$y, t = pre$t,
    background = pre$background,
    process = rep("control", length(pre$t)),
    location_process = rep("control", length(pre$t)),
    inferred_process = rep("control", length(pre$t)),
    W = rep(1, length(pre$t))
  )

  post <- generate_inhomogeneous_hawkes(
    Omega = Omega,
    partition = partition,
    time_window = c(treatment_time, end_time),
    partition_processes = partition_processes,
    hawkes_params = list(control = sim_control, treated = sim_treated),
    state_spaces = state_spaces,
    filtration = pre_df
  )
  post_df <- as.data.frame(post)
  post_df$inferred_process <- post_df$location_process
  post_df$W <- rep(1, nrow(post_df))

  list(
    Omega = Omega,
    partition = partition,
    partition_processes = partition_processes,
    state_spaces = state_spaces,
    treatment_time = treatment_time,
    end_time = end_time,
    pre = pre_df,
    post = post_df,
    observed = rbind(
      pre_df[, c("x", "y", "t", "background", "process", "location_process", "inferred_process", "W")],
      post_df[, c("x", "y", "t", "background", "process", "location_process", "inferred_process", "W")]
    ),
    power_control = power_control,
    power_treated = power_treated
  )
}

test_that("power-law Hawkes works through simulation-study labelling paths", {
  old_timing <- Sys.getenv("OK_SEM_TIMING_VERBOSE", unset = NA_character_)
  Sys.setenv(OK_SEM_TIMING_VERBOSE = "false", OK_SEM_PROPOSAL_TIMING = "false")
  on.exit({
    if (is.na(old_timing)) Sys.unsetenv("OK_SEM_TIMING_VERBOSE") else Sys.setenv(OK_SEM_TIMING_VERBOSE = old_timing)
  }, add = TRUE)

  d <- make_power_law_sim_study_fixture()
  expect_gt(nrow(d$post), 20)
  expect_true(all(d$post$inferred_process %in% c("control", "treated")))

  set.seed(20260707)
  proposal <- simulation_labeling_hawkes_hawkes_fast(
    d$post,
    partition = d$partition,
    partition_process = d$partition_processes,
    statespace = d$Omega,
    state_spaces = d$state_spaces,
    windowT = c(d$treatment_time, d$end_time),
    hawkes_params_control = d$power_control,
    hawkes_params_treated = d$power_treated,
    change_factor = 0.02,
    filtration = d$pre,
    proximity_weight = 0,
    kernel = "power_law",
    verbose = FALSE
  )
  expect_equal(nrow(proposal), nrow(d$post))
  expect_true(all(proposal$inferred_process %in% c("control", "treated")))

  set.seed(20260708)
  sem <- em_style_labelling(
    pp_data = d$observed,
    partition = d$partition,
    partition_processes = d$partition_processes,
    statespace = d$Omega,
    time_window = c(d$treatment_time, d$end_time),
    treatment_time = d$treatment_time,
    hawkes_params_control = d$power_control,
    hawkes_params_treated = d$power_control,
    update_control_params = FALSE,
    param_update_cadence = 1,
    proposal_update_cadence = 1,
    iter = 1,
    n_props = 1,
    change_factor = 0.02,
    state_spaces = d$state_spaces,
    proximity_weight = 0,
    outer_maxit = 5,
    kernel = "power_law",
    verbose = FALSE
  )
  expect_length(sem$accuracies, 1)
  expect_true(all(c("c", "p", "K") %in% names(sem$treated_par[[length(sem$treated_par)]])))
})

test_that("simulation-study ATE fitting supports exponential simulation with power-law fit", {
  d <- make_power_law_sim_study_fixture(sim_kernel = "exponential")

  set.seed(20260709)
  ate <- ATE_estim_hawkes(
    statespace = d$Omega,
    partition = d$partition,
    observed_data = d$post,
    treated_partitions = spatstat.geom::tilenames(d$partition)[d$partition_processes == "treated"],
    hawkes_params = NULL,
    n_tau_sims = 1,
    n_tau_i = 1,
    n_sims = 1,
    windowT = c(d$treatment_time, d$end_time),
    windowS = d$Omega,
    maxit = 5,
    poisson_flags = list(control = FALSE, treated = FALSE),
    filtration_data = d$pre,
    control_filtration_aware = TRUE,
    treated_params_init = d$power_control,
    keep_all_nothing_sim = FALSE,
    compute_tau = FALSE,
    kernel = "power_law"
  )

  expect_true(is.list(ate$control_pp))
  expect_true(is.list(ate$treated_pp))
  expect_identical(ate$control_pp$kernel, "power_law")
  expect_identical(ate$treated_pp$kernel, "power_law")
  expect_true(all(c("c", "p", "K") %in% names(ate$control_pp)))
  expect_true(is.finite(ate$all_nothing_theory$ATE))
})

test_that("power-law spatial Hawkes kernel is mean matched to exponential scale", {
  alpha <- 0.08
  q <- 2.2
  d <- PPDisentangle:::hawkes_power_law_spatial_d(alpha, q)
  exp_mean <- sqrt(pi) / (2 * sqrt(alpha))
  pl_mean <- sqrt(d) * (q - 1) * sqrt(pi) / 2 * gamma(q - 1.5) / gamma(q)

  expect_true(is.finite(d))
  expect_equal(pl_mean, exp_mean, tolerance = 1e-8)

  params <- PPDisentangle:::as_hawkes_params(
    list(mu = 5, alpha = alpha, beta = 1, K = 0.2, kernel = "exponential",
         spatial_kernel = "power_law", spatial_q = q)
  )
  expect_identical(params$spatial_kernel, "power_law")
  expect_equal(params$spatial_d, d, tolerance = 1e-8)
})

test_that("power-law spatial Hawkes works through C++ simulation and likelihood", {
  set.seed(20260710)
  params <- PPDisentangle:::as_hawkes_params(
    list(mu = 15, alpha = 0.08, beta = 1.5, K = 0.35,
         kernel = "exponential", spatial_kernel = "power_law", spatial_q = 2.2)
  )
  win <- spatstat.geom::owin(xrange = c(0, 20), yrange = c(0, 20))
  sim <- sim_hawkes_fast(params, c(0, 5), win)
  expect_gt(length(sim$t), 0)

  ll <- loglik_hawk_fast(
    params = params,
    realiz = data.frame(x = sim$x, y = sim$y, t = sim$t, W = rep(1, length(sim$t))),
    windowT = c(0, 5),
    windowS = win
  )
  expect_true(is.finite(ll))
})
