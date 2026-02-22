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
    state_spaces <- list(plus_treated_state_space, plus_control_state_space)
  } else {
    state_spaces <- list(plus_control_state_space, plus_treated_state_space)
  }

  plus_sims <- lapply(1:n_sim, function(j) {
    generate_inhomogeneous_hawkes(
      Omega = statespace, partition = partition, time_window = windowT,
      partition_processes = partition_plus,
      hawkes_params = list(control = control_pp, treated = treated_pp),
      state_spaces = state_spaces, space_triggering = FALSE
    )
  })

  if (partition_minus[1] == "treated") {
    state_spaces <- list(plus_treated_state_space, plus_control_state_space)
  } else {
    state_spaces <- list(plus_control_state_space, plus_treated_state_space)
  }

  minus_sims <- lapply(1:n_sim, function(j) {
    generate_inhomogeneous_hawkes(
      Omega = statespace, partition = partition, time_window = windowT,
      partition_processes = partition_minus,
      hawkes_params = list(control = control_pp, treated = treated_pp),
      state_spaces = state_spaces, space_triggering = FALSE
    )
  })

  plus_counts <- sapply(plus_sims, function(pp) {
    inds <- as.numeric(tileindex(pp$x, pp$y, partition))
    sum(inds == i)
  })
  minus_counts <- sapply(minus_sims, function(pp) {
    inds <- as.numeric(tileindex(pp$x, pp$y, partition))
    sum(inds == i)
  })
  return(mean(plus_counts) - mean(minus_counts))
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

  tau_i_estim <- sapply(1:n_tau_i, function(j) {
    tau_i(
      sample(length(partition_process), 1),
      partition = partition, treated_partitions = treated_partitions,
      statespace = statespace, windowT = windowT,
      control_pp = control_pp, treated_pp = treated_pp, n_sim = n_tau_sims
    )
  })
  tau_1_estim <- mean(tau_i_estim)

  all_nothing_sim <- lapply(1:n_sims, function(n) {
    control_sim <- sim_hawkes(control_pp, windowT, windowS = statespace)
    treated_sim <- sim_hawkes(treated_pp, windowT, windowS = statespace)
    c_mean <- control_sim$n[1] / partition$n
    t_mean <- treated_sim$n[1] / partition$n
    data.frame(c_mean = c_mean, t_mean = t_mean, ATE = t_mean - c_mean)
  })
  all_nothing_sim <- do.call(rbind, all_nothing_sim)

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
