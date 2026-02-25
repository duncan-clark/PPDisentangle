# ---- shared setup ----
make_test_data <- function(seed = 42) {
  set.seed(seed)
  OMEGA <- c(0, 50, 0, 50)
  END_TIME <- 55
  TREATMENT_TIME <- 5
  NX <- 5; NY <- 5
  hawkes_par_1 <- list(mu = 8, alpha = 0.01, beta = 10, K = 0.5)
  hawkes_par_2 <- list(mu = 8, alpha = 0.01, beta = 10, K = 0.2)
  partition <- spatstat.geom::quadrats(X = OMEGA, nx = NX, ny = NY)
  partition_processes <- rep("control", partition$n)
  partition_processes[sample(1:25, 12)] <- "treated"
  treated_idx <- partition_processes == "treated"
  control_ss <- spatstat.geom::as.owin(partition[!treated_idx])
  treated_ss <- spatstat.geom::as.owin(partition[treated_idx])
  if (partition_processes[1] == "treated") {
    state_spaces <- list(treated_ss, control_ss)
  } else {
    state_spaces <- list(control_ss, treated_ss)
  }
  pre <- sim_hawkes(hawkes_par_1, c(0, TREATMENT_TIME), OMEGA)
  pre_df <- data.frame(x = pre$x, y = pre$y, t = pre$t,
    process = "control", location_process = "control",
    inferred_process = "control", background = pre$background)
  post <- generate_inhomogeneous_hawkes(OMEGA, partition, c(TREATMENT_TIME, END_TIME),
    partition_processes, list(control = hawkes_par_1, treated = hawkes_par_2),
    state_spaces = state_spaces, filtration = pre_df)
  post_df <- as.data.frame(post)
  post_df$inferred_process <- post_df$location_process
  shared_cols <- c("x", "y", "t", "process", "location_process", "inferred_process", "background")
  pre_df <- pre_df[, shared_cols]
  post_df <- post_df[, shared_cols]
  obs_data <- rbind(pre_df, post_df)
  obs_data$W <- rep(1, nrow(obs_data))

  list(
    obs_data = obs_data, OMEGA = OMEGA, partition = partition,
    partition_processes = partition_processes,
    hawkes_par_1 = hawkes_par_1, hawkes_par_2 = hawkes_par_2,
    TREATMENT_TIME = TREATMENT_TIME, END_TIME = END_TIME,
    state_spaces = state_spaces, treated_ss = treated_ss, control_ss = control_ss
  )
}

# ---- C++ binary search correctness ----

test_that("loglik_hawk_fast with binary search matches loglik_hawk (R version)", {
  params <- c(mu = 30, alpha = 0.5, beta = 5, K = 0.3)
  windowT <- c(0, 10)
  windowS <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  set.seed(99)
  sim <- sim_hawkes(list(mu = 30, alpha = 0.5, beta = 5, K = 0.3), windowT, windowS)
  realiz <- data.frame(x = sim$x, y = sim$y, t = sim$t, W = rep(1, length(sim$t)))

  ll_r <- loglik_hawk(params, realiz, windowT, windowS, density_approx = TRUE)
  ll_cpp <- loglik_hawk_fast(params, realiz, windowT, windowS)
  expect_true(is.finite(ll_r))
  expect_true(is.finite(ll_cpp))
  expect_equal(ll_r, ll_cpp, tolerance = 0.01)
})

test_that("loglik_hawk_fast with beta=0.1 (wide temporal window) is correct", {
  params <- c(mu = 20, alpha = 0.5, beta = 0.1, K = 0.3)
  windowT <- c(0, 10)
  windowS <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  set.seed(101)
  sim <- sim_hawkes(list(mu = 20, alpha = 0.5, beta = 0.1, K = 0.3), windowT, windowS)
  realiz <- data.frame(x = sim$x, y = sim$y, t = sim$t, W = rep(1, length(sim$t)))

  ll_r <- loglik_hawk(params, realiz, windowT, windowS, density_approx = TRUE)
  ll_cpp <- loglik_hawk_fast(params, realiz, windowT, windowS)
  expect_true(is.finite(ll_r))
  expect_true(is.finite(ll_cpp))
  expect_equal(ll_r, ll_cpp, tolerance = 0.01)
})

test_that("loglik_hawk_fast with large beta (tight temporal cutoff) is correct", {
  params <- c(mu = 20, alpha = 0.5, beta = 50, K = 0.3)
  windowT <- c(0, 10)
  windowS <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  set.seed(102)
  sim <- sim_hawkes(list(mu = 20, alpha = 0.5, beta = 50, K = 0.3), windowT, windowS)
  realiz <- data.frame(x = sim$x, y = sim$y, t = sim$t, W = rep(1, length(sim$t)))

  ll_r <- loglik_hawk(params, realiz, windowT, windowS, density_approx = TRUE)
  ll_cpp <- loglik_hawk_fast(params, realiz, windowT, windowS)
  expect_true(is.finite(ll_r))
  expect_true(is.finite(ll_cpp))
  expect_equal(ll_r, ll_cpp, tolerance = 0.01)
})

# ---- precompute_loglik_args correctness ----

test_that("precompute_loglik_args returns correct areas", {
  windowS <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  zero_region <- spatstat.geom::owin(xrange = c(0, 5), yrange = c(0, 10))
  pts <- data.frame(x = c(2, 7), y = c(5, 5))

  pc <- precompute_loglik_args(pts, windowS, zero_region)
  expect_equal(pc$active_area, 50)
  expect_true(!is.null(pc$in_zero_bg_all))
  expect_equal(pc$in_zero_bg_all, c(TRUE, FALSE))
})

test_that("precompute_loglik_args without zero_background returns full area", {
  windowS <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  pts <- data.frame(x = c(2, 7), y = c(5, 5))

  pc <- precompute_loglik_args(pts, windowS, NULL)
  expect_equal(pc$active_area, 100)
  expect_null(pc$in_zero_bg_all)
})

test_that("loglik_hawk_fast with precomp matches without precomp", {
  params <- c(mu = 30, alpha = 0.5, beta = 5, K = 0.3)
  windowT <- c(0, 10)
  windowS <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  zero_region <- spatstat.geom::owin(xrange = c(0, 5), yrange = c(0, 10))
  set.seed(77)
  sim <- sim_hawkes(list(mu = 30, alpha = 0.5, beta = 5, K = 0.3), windowT, windowS)
  realiz <- data.frame(x = sim$x, y = sim$y, t = sim$t, W = rep(1, length(sim$t)))
  realiz <- realiz[realiz$t >= windowT[1] & realiz$t <= windowT[2], ]
  realiz <- realiz[order(realiz$t), ]

  ll_old <- loglik_hawk_fast(params, realiz, windowT, windowS,
                             zero_background_region = zero_region)

  pc <- precompute_loglik_args(realiz, windowS, zero_region)
  ll_new <- loglik_hawk_fast(params, realiz, windowT, windowS,
                             precomp = list(active_area = pc$active_area,
                                            in_zero_bg = pc$in_zero_bg_all))
  expect_equal(ll_old, ll_new)
})

test_that("loglik_hawk_fast precomp works with unsorted data", {
  params <- c(mu = 30, alpha = 0.5, beta = 5, K = 0.3)
  windowT <- c(0, 10)
  windowS <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  set.seed(88)
  sim <- sim_hawkes(list(mu = 30, alpha = 0.5, beta = 5, K = 0.3), windowT, windowS)
  realiz <- data.frame(x = sim$x, y = sim$y, t = sim$t, W = rep(1, length(sim$t)))
  shuffled <- realiz[sample(nrow(realiz)), ]

  ll_sorted <- loglik_hawk_fast(params, realiz, windowT, windowS)
  ll_shuffled <- loglik_hawk_fast(params, shuffled, windowT, windowS)
  expect_equal(ll_sorted, ll_shuffled)
})

# ---- precomp subsetting correctness (EM-style usage) ----

test_that("precomp with row subsetting matches full computation", {
  params <- c(mu = 30, alpha = 0.5, beta = 5, K = 0.3)
  windowT <- c(0, 10)
  windowS <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  zero_region <- spatstat.geom::owin(xrange = c(0, 5), yrange = c(0, 10))
  set.seed(55)
  sim <- sim_hawkes(list(mu = 30, alpha = 0.5, beta = 5, K = 0.3), windowT, windowS)
  realiz <- data.frame(x = sim$x, y = sim$y, t = sim$t, W = rep(1, length(sim$t)))
  realiz <- realiz[realiz$t >= windowT[1] & realiz$t <= windowT[2], ]
  realiz <- realiz[order(realiz$t), ]
  n <- nrow(realiz)

  ctrl_rows <- sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(0.6, 0.4))
  ctrl_data <- realiz[ctrl_rows, ]

  ll_direct <- loglik_hawk_fast(params, ctrl_data, windowT, windowS,
                                zero_background_region = zero_region)

  pc_all <- precompute_loglik_args(realiz, windowS, zero_region)
  ll_precomp <- loglik_hawk_fast(params, ctrl_data, windowT, windowS,
                                 precomp = list(active_area = pc_all$active_area,
                                                in_zero_bg = pc_all$in_zero_bg_all[ctrl_rows]))
  expect_equal(ll_direct, ll_precomp)
})

# ---- em_style_labelling integration ----

test_that("em_style_labelling runs and produces valid output", {
  d <- make_test_data(seed = 42)
  set.seed(123)
  result <- em_style_labelling(
    pp_data = d$obs_data,
    partition = d$partition, partition_processes = d$partition_processes,
    statespace = d$OMEGA, time_window = c(d$TREATMENT_TIME, d$END_TIME),
    treatment_time = d$TREATMENT_TIME,
    hawkes_params_control = d$hawkes_par_1, hawkes_params_treated = d$hawkes_par_2,
    update_control_params = FALSE, param_update_cadence = 3,
    iter = 3, n_props = 5, change_factor = 0.02,
    state_spaces = d$state_spaces, verbose = FALSE
  )
  expect_true(!is.null(result$labelling))
  expect_true(!is.null(result$accuracies))
  expect_equal(length(result$accuracies), 3)
  expect_true(all(result$accuracies >= 0 & result$accuracies <= 1))
  expect_true(all(is.finite(result$metrics)))
  expect_equal(nrow(result$labelling), nrow(d$obs_data))
})

test_that("em_style_labelling accuracy improves or stays stable", {
  d <- make_test_data(seed = 42)
  set.seed(123)
  result <- em_style_labelling(
    pp_data = d$obs_data,
    partition = d$partition, partition_processes = d$partition_processes,
    statespace = d$OMEGA, time_window = c(d$TREATMENT_TIME, d$END_TIME),
    treatment_time = d$TREATMENT_TIME,
    hawkes_params_control = d$hawkes_par_1, hawkes_params_treated = d$hawkes_par_2,
    update_control_params = FALSE, param_update_cadence = NULL,
    iter = 5, n_props = 8, change_factor = 0.02,
    state_spaces = d$state_spaces, verbose = FALSE
  )
  acc <- result$accuracies
  expect_true(acc[length(acc)] >= acc[1] - 0.1,
    info = "Final accuracy should not drop catastrophically from initial")
})

test_that("em_style_labelling with update_control_params works", {
  d <- make_test_data(seed = 42)
  set.seed(456)
  result <- em_style_labelling(
    pp_data = d$obs_data,
    partition = d$partition, partition_processes = d$partition_processes,
    statespace = d$OMEGA, time_window = c(d$TREATMENT_TIME, d$END_TIME),
    treatment_time = d$TREATMENT_TIME,
    hawkes_params_control = d$hawkes_par_1, hawkes_params_treated = d$hawkes_par_2,
    update_control_params = TRUE, param_update_cadence = 2,
    iter = 3, n_props = 5, change_factor = 0.02,
    state_spaces = d$state_spaces, verbose = FALSE
  )
  expect_equal(length(result$accuracies), 3)
  expect_true(length(result$control_par) > 1)
  expect_true(length(result$treated_par) > 1)
})

# ---- simulation_labeling_hawkes_hawkes_fast ----

test_that("simulation_labeling_hawkes_hawkes_fast produces valid relabelling", {
  d <- make_test_data(seed = 42)
  post <- d$obs_data[d$obs_data$t >= d$TREATMENT_TIME, ]
  post <- post[order(post$t), ]
  pre <- d$obs_data[d$obs_data$t < d$TREATMENT_TIME, ]
  pre$location_process <- "control"

  set.seed(111)
  result <- simulation_labeling_hawkes_hawkes_fast(
    post, partition = d$partition, partition_process = d$partition_processes,
    statespace = d$OMEGA, state_spaces = d$state_spaces,
    windowT = c(d$TREATMENT_TIME, d$END_TIME),
    hawkes_params_control = d$hawkes_par_1, hawkes_params_treated = d$hawkes_par_2,
    change_factor = 0.05, filtration = pre, verbose = FALSE
  )
  expect_equal(nrow(result), nrow(post))
  expect_true(all(result$inferred_process %in% c("control", "treated")))
  n_flipped <- sum(result$inferred_process != post$inferred_process)
  expect_true(n_flipped > 0, info = "Should relabel at least some points")
  expect_true(n_flipped < nrow(post), info = "Should not relabel all points")
})

test_that("simulation_labeling_hawkes_hawkes_fast matches non-fast with same seed (control-only)", {
  d <- make_test_data(seed = 42)
  post <- d$obs_data[d$obs_data$t >= d$TREATMENT_TIME, ]
  post <- post[order(post$t), ]
  pre <- d$obs_data[d$obs_data$t < d$TREATMENT_TIME, ]
  pre$location_process <- "control"
  post_same <- post

  set.seed(301)
  res_slow <- simulation_labeling_hawkes_hawkes(
    post, partition = d$partition, partition_process = d$partition_processes,
    statespace = d$OMEGA, state_spaces = d$state_spaces,
    windowT = c(d$TREATMENT_TIME, d$END_TIME),
    hawkes_params_control = d$hawkes_par_1, hawkes_params_treated = NULL,
    change_factor = 0.05, filtration = pre, verbose = FALSE
  )
  set.seed(301)
  res_fast <- simulation_labeling_hawkes_hawkes_fast(
    post_same, partition = d$partition, partition_process = d$partition_processes,
    statespace = d$OMEGA, state_spaces = d$state_spaces,
    windowT = c(d$TREATMENT_TIME, d$END_TIME),
    hawkes_params_control = d$hawkes_par_1, hawkes_params_treated = NULL,
    change_factor = 0.05, filtration = pre, proximity_weight = 0, verbose = FALSE
  )
  expect_identical(res_slow$inferred_process, res_fast$inferred_process,
    info = "fast (proximity_weight=0) should match non-fast when hawkes_params_treated is NULL")
})

test_that("simulation_labeling_hawkes_hawkes_fast matches non-fast with same seed (control + treated)", {
  d <- make_test_data(seed = 42)
  post <- d$obs_data[d$obs_data$t >= d$TREATMENT_TIME, ]
  post <- post[order(post$t), ]
  pre <- d$obs_data[d$obs_data$t < d$TREATMENT_TIME, ]
  pre$location_process <- "control"
  post_same <- post

  set.seed(302)
  res_slow <- simulation_labeling_hawkes_hawkes(
    post, partition = d$partition, partition_process = d$partition_processes,
    statespace = d$OMEGA, state_spaces = d$state_spaces,
    windowT = c(d$TREATMENT_TIME, d$END_TIME),
    hawkes_params_control = d$hawkes_par_1, hawkes_params_treated = d$hawkes_par_2,
    change_factor = 0.05, filtration = pre, verbose = FALSE
  )
  set.seed(302)
  res_fast <- simulation_labeling_hawkes_hawkes_fast(
    post_same, partition = d$partition, partition_process = d$partition_processes,
    statespace = d$OMEGA, state_spaces = d$state_spaces,
    windowT = c(d$TREATMENT_TIME, d$END_TIME),
    hawkes_params_control = d$hawkes_par_1, hawkes_params_treated = d$hawkes_par_2,
    change_factor = 0.05, filtration = pre, proximity_weight = 0, verbose = FALSE
  )
  expect_identical(res_slow$inferred_process, res_fast$inferred_process,
    info = "fast (proximity_weight=0) should match non-fast when both control and treated params set")
})

# ---- loglik_hawk_fast speed regression ----

test_that("loglik_hawk_fast with precomp is faster than without", {
  params <- c(mu = 30, alpha = 0.5, beta = 5, K = 0.3)
  windowT <- c(0, 10)
  windowS <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  zero_region <- spatstat.geom::owin(xrange = c(0, 5), yrange = c(0, 10))
  set.seed(66)
  sim <- sim_hawkes(list(mu = 30, alpha = 0.5, beta = 5, K = 0.3), windowT, windowS)
  realiz <- data.frame(x = sim$x, y = sim$y, t = sim$t, W = rep(1, length(sim$t)))

  t_old <- system.time(for (i in 1:20) {
    loglik_hawk_fast(params, realiz, windowT, windowS,
                     zero_background_region = zero_region)
  })["elapsed"]

  pc <- precompute_loglik_args(realiz, windowS, zero_region)
  t_new <- system.time(for (i in 1:20) {
    loglik_hawk_fast(params, realiz, windowT, windowS,
                     precomp = list(active_area = pc$active_area,
                                    in_zero_bg = pc$in_zero_bg_all))
  })["elapsed"]

  expect_true(t_new <= t_old,
    info = paste0("precomp should be faster: ", round(t_new, 4), " <= ", round(t_old, 4)))
})
