#' Estimate the one-flip treatment effect for cell i
#'
#' Simulates under z_plus (cell i treated) and z_minus (cell i control),
#' returns the difference in expected counts in cell i.
#'
#' @param i Cell index
#' @param partition A spatstat tess object
#' @param treated_partitions Tilenames of treated partitions
#' @param statespace Full observation window
#' @param windowT Numeric vector c(start, end)
#' @param control_pp Control Hawkes parameters
#' @param treated_pp Treated Hawkes parameters
#' @param n_sim Number of simulations
#' @return Scalar estimate of tau_i
#' @export
tau_i <- function(i, partition, treated_partitions, statespace,
                  windowT, control_pp, treated_pp, n_sim) {
  treated_idx <- tilenames(partition) %in% treated_partitions

  z_plus <- treated_idx
  z_plus[i] <- TRUE
  z_minus <- treated_idx
  z_minus[i] <- FALSE

  partition_plus <- rep("control", partition$n)
  partition_plus[z_plus] <- "treated"
  partition_minus <- rep("control", partition$n)
  partition_minus[z_minus] <- "treated"

  plus_control_state_space <- as.owin(partition[!z_plus])
  plus_treated_state_space <- as.owin(partition[z_plus])
  minus_control_state_space <- as.owin(partition[!z_minus])
  minus_treated_state_space <- as.owin(partition[z_minus])

  if (partition_plus[1] == "treated") {
    state_spaces_plus <- list(plus_treated_state_space, plus_control_state_space)
  } else {
    state_spaces_plus <- list(plus_control_state_space, plus_treated_state_space)
  }

  if (partition_minus[1] == "treated") {
    state_spaces_minus <- list(minus_treated_state_space, minus_control_state_space)
  } else {
    state_spaces_minus <- list(minus_control_state_space, minus_treated_state_space)
  }

  plus_total <- 0
  minus_total <- 0
  for (j in seq_len(n_sim)) {
    pp <- generate_inhomogeneous_hawkes(
      Omega = statespace, partition = partition, time_window = windowT,
      partition_processes = partition_plus,
      hawkes_params = list(control = control_pp, treated = treated_pp),
      state_spaces = state_spaces_plus, space_triggering = FALSE
    )
    plus_total <- plus_total + sum(as.numeric(tileindex(pp$x, pp$y, partition)) == i)

    pp <- generate_inhomogeneous_hawkes(
      Omega = statespace, partition = partition, time_window = windowT,
      partition_processes = partition_minus,
      hawkes_params = list(control = control_pp, treated = treated_pp),
      state_spaces = state_spaces_minus, space_triggering = FALSE
    )
    minus_total <- minus_total + sum(as.numeric(tileindex(pp$x, pp$y, partition)) == i)
  }
  return(plus_total / n_sim - minus_total / n_sim)
}

#' Non-parametric ATE estimation from labeled data
#'
#' @param statespace Full observation window
#' @param partition A spatstat tess object
#' @param observed_data Data frame with inferred_process column
#' @param treated_partitions Tilenames of treated partitions
#' @param pp_rhs_form Formula for intensity (unused, kept for API compat)
#' @param windowT Numeric vector c(start, end)
#' @return List with ATE_non_param
#' @export
ATE_estim <- function(statespace, partition, observed_data, treated_partitions,
                      pp_rhs_form = "~1", windowT = c(0, 1)) {
  treated_idx <- tilenames(partition) %in% treated_partitions
  control_points_per_tile <- sum(observed_data$inferred_process == "control") / sum(!treated_idx)
  treated_points_per_tile <- sum(observed_data$inferred_process == "treated") / sum(treated_idx)
  ATE_naive <- treated_points_per_tile - control_points_per_tile
  return(list(ATE_non_param = ATE_naive, ATE_non_param_z = ATE_naive))
}

#' Parametric ATE estimation using fitted Hawkes processes
#'
#' Fits Hawkes models to control/treated subsets, then estimates ATE via
#' simulation (all-or-nothing and one-flip).
#'
#' @param statespace Full observation window
#' @param partition A spatstat tess object
#' @param observed_data Data frame with labeled point process data
#' @param treated_partitions Tilenames of treated partitions
#' @param hawkes_params Optional pre-fitted params (list with control, treated)
#' @param n_tau_sims Number of sims per tau_i estimate
#' @param n_tau_i Number of tau_i estimates
#' @param n_sims Number of all-or-nothing simulations
#' @param windowT Numeric vector c(start, end)
#' @param windowS Observation window for fitting
#' @param maxit Max optim iterations
#' @param poisson_flags List with control and treated Poisson flags
#' @return List with all-or-nothing sims, tau estimates, and decomposed ATEs
#' @export
ATE_estim_hawkes <- function(statespace, partition, observed_data, treated_partitions,
                             hawkes_params = NULL, n_tau_sims = 10, n_tau_i = 10,
                             n_sims = 100, windowT = c(0, 1), windowS = c(0, 1, 0, 1),
                             maxit = 1000, poisson_flags = list(control = FALSE, treated = FALSE)) {
  treated_idx <- tilenames(partition) %in% treated_partitions
  control_state_space <- as.owin(partition[!treated_idx])
  treated_state_space <- as.owin(partition[treated_idx])

  control_mean <- sum(observed_data$location_process == "control") / sum(!treated_idx)
  treated_mean <- sum(observed_data$location_process == "treated") / sum(treated_idx)
  control_points_in_control_cells <- sum(observed_data$inferred_process == "control" & observed_data$location_process == "control") / sum(!treated_idx)
  treated_points_in_treated_cells <- sum(observed_data$inferred_process == "treated" & observed_data$location_process == "treated") / sum(treated_idx)

  ATE_total <- treated_mean - control_points_in_control_cells
  ATE_treatment <- treated_mean - control_mean
  ATE_spillover <- control_mean - control_points_in_control_cells
  ATE_naive <- ATE_treatment

  if (is.null(hawkes_params)) {
    control_pp <- fit_hawkes(
      unlist(list(mu = (sum(observed_data$location_process == "control")) / (windowT[2] - windowT[1]),
                  alpha = 0, beta = (windowT[2] - windowT[1]) / 100, K = 0.01)),
      realiz = observed_data[observed_data$inferred_process == "control", ],
      zero_background_region = treated_state_space,
      windowT = windowT, windowS = windowS, trace = 0, maxit = maxit,
      density_approx = FALSE, numeric_integral = FALSE,
      poisson_flag = poisson_flags$control
    )$par
    control_pp <- as.list(control_pp)
    names(control_pp) <- c("mu", "alpha", "beta", "K")

    treated_pp <- fit_hawkes(
      unlist(list(mu = (sum(observed_data$location_process == "treated")) / (windowT[2] - windowT[1]),
                  alpha = 0, beta = (windowT[2] - windowT[1]) / 100, K = 0.01)),
      realiz = observed_data[observed_data$inferred_process == "treated", ],
      zero_background_region = control_state_space,
      windowT = windowT, windowS = windowS, trace = 0, maxit = maxit,
      density_approx = FALSE, numeric_integral = FALSE,
      poisson_flag = poisson_flags$treated
    )$par
    treated_pp <- as.list(treated_pp)
    names(treated_pp) <- c("mu", "alpha", "beta", "K")
  } else {
    control_pp <- hawkes_params$control
    treated_pp <- hawkes_params$treated
  }

  partition_process <- rep("control", partition$n)
  partition_process[treated_idx] <- "treated"

  tau_i_estim <- vapply(seq_len(n_tau_i), function(j) {
    tau_i(
      sample(length(partition_process), 1),
      partition = partition, treated_partitions = treated_partitions,
      statespace = statespace, windowT = windowT,
      control_pp = control_pp, treated_pp = treated_pp, n_sim = n_tau_sims
    )
  }, numeric(1))
  tau_1_estim <- mean(tau_i_estim)

  c_counts <- numeric(n_sims)
  t_counts <- numeric(n_sims)
  for (s in seq_len(n_sims)) {
    c_counts[s] <- sim_hawkes(control_pp, windowT, windowS = statespace)$n[1]
    t_counts[s] <- sim_hawkes(treated_pp, windowT, windowS = statespace)$n[1]
  }
  all_nothing_sim <- data.frame(
    c_mean = c_counts / partition$n,
    t_mean = t_counts / partition$n,
    ATE = (t_counts - c_counts) / partition$n
  )

  all_nothing_theory <- data.frame(
    c_mean = control_pp$mu * (windowT[2] - windowT[1]) * (1 / (1 - control_pp$K)) / partition$n,
    t_mean = treated_pp$mu * (windowT[2] - windowT[1]) * (1 / (1 - treated_pp$K)) / partition$n,
    ATE = (treated_pp$mu * (windowT[2] - windowT[1]) * (1 / (1 - treated_pp$K)) -
             control_pp$mu * (windowT[2] - windowT[1]) * (1 / (1 - control_pp$K))) / partition$n
  )

  return(list(
    all_nothing_sim = all_nothing_sim,
    all_nothing_theory = all_nothing_theory,
    tau_1_estim = tau_1_estim,
    ATE_total = ATE_total, ATE_treatment = ATE_treatment,
    ATE_spillover = ATE_spillover, ATE_naive = ATE_naive,
    treated_pp = treated_pp, control_pp = control_pp
  ))
}


#' Estimate one-flip treatment effect for cell i (ETAS model)
#'
#' Analogue of \code{tau_i} for the ETAS model.
#'
#' @inheritParams tau_i
#' @param control_pp  Control ETAS parameter list.
#' @param treated_pp  Treated ETAS parameter list.
#' @param m0  Reference magnitude.
#' @param beta_gr  Gutenberg-Richter beta (or NULL for resampling).
#' @param mag_pool  Magnitude pool for resampling.
#' @return Scalar estimate of tau_i.
#' @export
tau_i_etas <- function(i, partition, treated_partitions, statespace,
                       windowT, control_pp, treated_pp, n_sim,
                       m0, beta_gr = NULL, mag_pool = NULL) {
  treated_idx <- tilenames(partition) %in% treated_partitions

  z_plus <- treated_idx;  z_plus[i]  <- TRUE
  z_minus <- treated_idx; z_minus[i] <- FALSE

  partition_plus  <- ifelse(z_plus, "treated", "control")
  partition_minus <- ifelse(z_minus, "treated", "control")

  make_ss <- function(z) {
    pr <- unique(ifelse(z, "treated", "control"))
    ss <- list()
    for (p in pr) ss[[p]] <- as.owin(partition[ifelse(z, "treated", "control") == p])
    ss
  }

  plus_total <- 0; minus_total <- 0
  for (j in seq_len(n_sim)) {
    pp <- generate_inhomogeneous_etas(
      Omega = statespace, partition = partition, time_window = windowT,
      partition_processes = partition_plus,
      etas_params = list(control = control_pp, treated = treated_pp),
      m0 = m0, beta_gr = beta_gr, mag_pool = mag_pool
    )
    plus_total <- plus_total + sum(as.numeric(tileindex(pp$x, pp$y, partition)) == i)

    pp <- generate_inhomogeneous_etas(
      Omega = statespace, partition = partition, time_window = windowT,
      partition_processes = partition_minus,
      etas_params = list(control = control_pp, treated = treated_pp),
      m0 = m0, beta_gr = beta_gr, mag_pool = mag_pool
    )
    minus_total <- minus_total + sum(as.numeric(tileindex(pp$x, pp$y, partition)) == i)
  }
  return(plus_total / n_sim - minus_total / n_sim)
}


#' Parametric ATE estimation using fitted ETAS processes
#'
#' ETAS analogue of \code{ATE_estim_hawkes}.  Fits ETAS models to
#' control/treated subsets, then estimates ATE via simulation.
#'
#' @inheritParams ATE_estim_hawkes
#' @param etas_params  Optional pre-fitted ETAS params
#'   (list with control, treated).
#' @param m0  Reference magnitude.
#' @param beta_gr  Gutenberg-Richter beta (or NULL for resampling).
#' @param mag_pool  Numeric vector of magnitudes for resampling.
#' @param fixed_params  Named list of ETAS params to hold fixed.
#' @return List with all-or-nothing sims, tau estimates, decomposed ATEs,
#'   and fitted parameter lists.
#' @export
ATE_estim_etas <- function(statespace, partition, observed_data,
                           treated_partitions,
                           etas_params = NULL,
                           n_tau_sims = 10, n_tau_i = 10,
                           n_sims = 100,
                           windowT = c(0, 1),
                           windowS = NULL,
                           maxit = 3000,
                           m0 = NULL,
                           beta_gr = NULL,
                           mag_pool = NULL,
                           fixed_params = NULL) {
  if (is.null(windowS)) windowS <- statespace
  treated_idx <- tilenames(partition) %in% treated_partitions
  control_state_space <- as.owin(partition[!treated_idx])
  treated_state_space <- as.owin(partition[treated_idx])

  if (is.null(m0)) m0 <- min(observed_data$mag, na.rm = TRUE)

  control_mean <- sum(observed_data$location_process == "control") / sum(!treated_idx)
  treated_mean <- sum(observed_data$location_process == "treated") / sum(treated_idx)
  control_in_ctrl <- sum(observed_data$inferred_process == "control" &
                         observed_data$location_process == "control") / sum(!treated_idx)
  treated_in_treat <- sum(observed_data$inferred_process == "treated" &
                          observed_data$location_process == "treated") / sum(treated_idx)

  ATE_total     <- treated_mean - control_in_ctrl
  ATE_treatment <- treated_mean - control_mean
  ATE_spillover <- control_mean - control_in_ctrl
  ATE_naive     <- ATE_treatment

  if (is.null(etas_params)) {
    dt <- windowT[2] - windowT[1]
    init_ctrl <- list(mu = sum(observed_data$inferred_process == "control") / dt,
                      A = 0.1, alpha_m = 0.5, c = 0.05, p = 1.2,
                      D = 1.0, gamma = 0.3, q = 1.5)
    init_treat <- list(mu = sum(observed_data$inferred_process == "treated") / dt,
                       A = 0.1, alpha_m = 0.5, c = 0.05, p = 1.2,
                       D = 1.0, gamma = 0.3, q = 1.5)

    ctrl_fit <- fit_etas(
      params_init = init_ctrl,
      realiz = observed_data[observed_data$inferred_process == "control", ],
      windowT = windowT, windowS = windowS, m0 = m0,
      maxit = maxit, fixed_params = fixed_params,
      zero_background_region = treated_state_space
    )
    control_pp <- as.list(ctrl_fit$par)
    names(control_pp) <- .etas_par_names

    treat_fit <- fit_etas(
      params_init = init_treat,
      realiz = observed_data[observed_data$inferred_process == "treated", ],
      windowT = windowT, windowS = windowS, m0 = m0,
      maxit = maxit, fixed_params = fixed_params,
      zero_background_region = control_state_space
    )
    treated_pp <- as.list(treat_fit$par)
    names(treated_pp) <- .etas_par_names
  } else {
    control_pp <- etas_params$control
    treated_pp <- etas_params$treated
  }

  partition_process <- rep("control", partition$n)
  partition_process[treated_idx] <- "treated"

  tau_i_estim <- vapply(seq_len(n_tau_i), function(j) {
    tau_i_etas(
      sample(length(partition_process), 1),
      partition = partition, treated_partitions = treated_partitions,
      statespace = statespace, windowT = windowT,
      control_pp = control_pp, treated_pp = treated_pp, n_sim = n_tau_sims,
      m0 = m0, beta_gr = beta_gr, mag_pool = mag_pool
    )
  }, numeric(1))
  tau_1_estim <- mean(tau_i_estim)

  c_counts <- numeric(n_sims)
  t_counts <- numeric(n_sims)
  for (s in seq_len(n_sims)) {
    c_sim <- sim_etas(control_pp, windowT, windowS = statespace,
                      m0 = m0, beta_gr = beta_gr, mag_pool = mag_pool)
    t_sim <- sim_etas(treated_pp, windowT, windowS = statespace,
                      m0 = m0, beta_gr = beta_gr, mag_pool = mag_pool)
    c_counts[s] <- length(c_sim$t)
    t_counts[s] <- length(t_sim$t)
  }
  all_nothing_sim <- data.frame(
    c_mean = c_counts / partition$n,
    t_mean = t_counts / partition$n,
    ATE    = (t_counts - c_counts) / partition$n
  )

  return(list(
    all_nothing_sim = all_nothing_sim,
    tau_1_estim = tau_1_estim,
    ATE_total = ATE_total, ATE_treatment = ATE_treatment,
    ATE_spillover = ATE_spillover, ATE_naive = ATE_naive,
    treated_pp = treated_pp, control_pp = control_pp
  ))
}
