test_that("ATE_estim returns non-parametric estimate", {
  obs <- data.frame(
    x = runif(100), y = runif(100), t = runif(100),
    inferred_process = sample(c("control", "treated"), 100, replace = TRUE),
    location_process = sample(c("control", "treated"), 100, replace = TRUE)
  )
  Omega <- spatstat.geom::owin(xrange = c(0, 1), yrange = c(0, 1))
  partition <- spatstat.geom::quadrats(Omega, nx = 2, ny = 2)
  treated_parts <- tilenames(partition)[c(3, 4)]

  result <- ATE_estim(Omega, partition, obs, treated_parts, windowT = c(0, 1))
  expect_true(is.numeric(result$ATE_non_param))
  expect_true(is.finite(result$ATE_non_param))
})

test_that("generate_poisson_data returns valid structure", {
  Omega <- spatstat.geom::owin(xrange = c(0, 10), yrange = c(0, 10))
  partition <- spatstat.geom::quadrats(Omega, nx = 3, ny = 3)
  set.seed(42)
  result <- generate_poisson_data(Omega, partition, treat_prop = 0.5, lambda_1 = 50, lambda_2 = 80)

  expect_true(nrow(result$obs_data) > 0)
  expect_true("process" %in% names(result$obs_data))
  expect_true("location_process" %in% names(result$obs_data))
  expect_true(all(result$obs_data$process %in% c("control", "treated")))
})

test_that("oracle recovers superposition Hawkes means in sim-study setup", {
  set.seed(20260324)
  Omega <- c(0, 50, 0, 50)
  partition <- spatstat.geom::quadrats(X = Omega, nx = 6, ny = 6)

  partition_processes <- rep("control", partition$n)
  partition_processes[sample(seq_len(partition$n), partition$n / 2)] <- "treated"
  treated_partitions <- tilenames(partition[partition_processes == "treated"])

  control_pp <- list(mu = 3, alpha = 0.02, beta = 8, K = 0.35)
  treated_pp <- list(mu = 4, alpha = 0.02, beta = 8, K = 0.20)
  dt <- 80
  n_ctrl_tiles <- sum(partition_processes == "control")
  n_treat_tiles <- sum(partition_processes == "treated")

  one_rep <- function(seed) {
    set.seed(seed)
    sim <- generate_hawkes_superposition(
      Omega = Omega,
      partition = partition,
      time_window = c(0, dt),
      partition_processes = partition_processes,
      hawkes_params = list(control = control_pp, treated = treated_pp)
    )
    obs <- as.data.frame(sim$superposed)
    tile_idx <- as.integer(obs$tile)
    obs$location_process <- partition_processes[pmin(pmax(tile_idx, 1L), partition$n)]
    obs$process <- obs$location_process
    obs <- oracle_labeling(obs)

    fit <- ATE_estim_hawkes(
      statespace = Omega,
      partition = partition,
      observed_data = obs,
      treated_partitions = treated_partitions,
      hawkes_params = NULL,
      n_tau_sims = 2,
      n_tau_i = 2,
      n_sims = 6,
      windowT = c(0, dt),
      windowS = Omega,
      maxit = 700,
      poisson_flags = list(control = FALSE, treated = FALSE),
      filtration_data = obs[0, c("x", "y", "t", "inferred_process", "location_process"), drop = FALSE],
      control_filtration_aware = FALSE
    )

    est_ctrl_pts <- fit$control_pp$mu * dt * (1 / (1 - fit$control_pp$K)) / n_ctrl_tiles
    est_treat_pts <- fit$treated_pp$mu * dt * (1 / (1 - fit$treated_pp$K)) / n_treat_tiles
    emp_ctrl_pts <- sum(obs$process == "control") / n_ctrl_tiles
    emp_treat_pts <- sum(obs$process == "treated") / n_treat_tiles
    c(est_ctrl_pts = est_ctrl_pts, emp_ctrl_pts = emp_ctrl_pts,
      est_treat_pts = est_treat_pts, emp_treat_pts = emp_treat_pts)
  }

  out <- do.call(rbind, lapply(seq_len(8), function(i) one_rep(93000 + i)))
  mae_ctrl <- mean(abs(out[, "est_ctrl_pts"] - out[, "emp_ctrl_pts"]))
  mae_treat <- mean(abs(out[, "est_treat_pts"] - out[, "emp_treat_pts"]))

  expect_true(mae_ctrl < 0.5)
  expect_true(mae_treat < 0.5)
})
