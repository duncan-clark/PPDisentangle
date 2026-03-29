#' @importFrom dplyr filter arrange pull group_by summarize bind_rows
NULL

#' Oracle labeling: use the true process labels
#' @param pp_data Data frame with a `process` column
#' @return pp_data with `inferred_process` set to `process`
#' @export
oracle_labeling <- function(pp_data) {
  pp_data$inferred_process <- pp_data$process
  return(pp_data)
}

#' Naive labeling: use the location-based process labels
#' @param pp_data Data frame with a `location_process` column
#' @return pp_data with `inferred_process` set to `location_process`
#' @export
naive_labeling <- function(pp_data) {
  pp_data$inferred_process <- pp_data$location_process
  return(pp_data)
}

#' Dumb labeling: probabilistic relabeling based on expected counts
#' @param pp_data Data frame with columns x, y, process, location_process, background
#' @param partition A spatstat tess object
#' @param tiles_to_thin Integer vector of tile indices to relabel
#' @param expected_points Expected number of points per tile
#' @return pp_data with updated inferred_process
#' @export
dumb_labeling <- function(pp_data, partition, tiles_to_thin, expected_points) {
  inds <- as.numeric(tileindex(pp_data$x, pp_data$y, partition))
  labels <- unique(pp_data$process)
  pp_data$inferred_process <- pp_data$location_process
  if (length(expected_points) == 1) {
    expected_points_vec <- rep(0, max(tiles_to_thin))
    expected_points_vec[tiles_to_thin] <- expected_points
    expected_points <- expected_points_vec
  }
  for (i in tiles_to_thin) {
    n_expected <- expected_points[i]
    n_actual <- sum(inds == i)
    if (n_actual == 0) next
    correct_process_prob <- min(1, n_expected / n_actual)
    inferred <- rbinom(sum(inds == i), 1, correct_process_prob)
    current_label <- unique(pp_data$location_process[inds == i])
    pp_data$inferred_process[inds == i] <- c(labels[labels != current_label], current_label)[inferred + 1]
    pp_data$background[inds == i][inferred == 0] <- TRUE
  }
  pp_data$n <- dim(pp_data)[1]
  return(pp_data)
}

#' Semi-dumb labeling: boundary-weighted probabilistic relabeling
#' @param pp_data Data frame with point process data
#' @param partition A spatstat tess object
#' @param tiles_to_thin Integer vector of tile indices
#' @param expected_points Expected number of points per tile
#' @param beta Distance weighting parameter
#' @return pp_data with updated inferred_process
#' @export
semi_dumb_labeling <- function(pp_data, partition, tiles_to_thin, expected_points, beta = 1) {
  inds <- as.numeric(tileindex(pp_data$x, pp_data$y, partition))
  labels <- unique(pp_data$process)
  pp_data$inferred_process <- pp_data$location_process
  x_diff <- outer(pp_data$x, pp_data$x, "-")
  y_diff <- outer(pp_data$y, pp_data$y, "-")
  space_dist <- x_diff^2 + y_diff^2
  treated_points <- which(pp_data$location_process == "treated")

  if (length(expected_points) == 1) {
    expected_points_vec <- rep(0, max(tiles_to_thin))
    expected_points_vec[tiles_to_thin] <- expected_points
    expected_points <- expected_points_vec
  }
  for (i in tiles_to_thin) {
    n_expected <- expected_points[i]
    n_actual <- sum(inds == i)
    if (n_actual == 0 | n_actual <= n_expected) next
    dists <- space_dist[inds == i, treated_points]
    min_dists <- apply(dists, 1, min)
    probs <- exp(-beta * min_dists)
    probs <- probs / sum(probs)
    inferred <- sample(1:sum(inds == i), size = n_actual - n_expected, replace = FALSE, prob = probs)
    current_label <- unique(pp_data$location_process[inds == i])
    label <- rep(current_label, sum(inds == i))
    label[inferred] <- labels[labels != current_label]
    pp_data$inferred_process[inds == i] <- label
    pp_data$background[inds == i][inferred] <- TRUE
  }
  pp_data$n <- dim(pp_data)[1]
  return(pp_data)
}

#' Simulation-based labeling for Hawkes-Hawkes processes
#'
#' Relabels points by comparing observed vs simulated counts per tile.
#'
#' @param pp_data Data frame with columns x, y, t, inferred_process, location_process
#' @param partition A spatstat tess object
#' @param partition_process Character vector of process names per tile
#' @param statespace The full observation window (owin)
#' @param windowT Numeric vector c(start, end)
#' @param state_spaces Optional precomputed state spaces
#' @param hawkes_params_control Control Hawkes parameters
#' @param hawkes_params_treated Treated Hawkes parameters
#' @param change_factor Scaling factor for number of relabelings
#' @param filtration Optional data frame of pre-treatment history
#' @param verbose Print progress
#' @param ... Additional arguments passed to generate_inhomogeneous_hawkes
#' @return pp_data with updated inferred_process
#' @export
simulation_labeling_hawkes_hawkes <- function(pp_data,
                                              partition,
                                              partition_process,
                                              statespace,
                                              windowT,
                                              state_spaces = NULL,
                                              hawkes_params_control = NULL,
                                              hawkes_params_treated = NULL,
                                              change_factor = 1,
                                              filtration = NULL,
                                              verbose = FALSE,
                                              ...) {
  no_points_hawkes <- list(mu = 0, alpha = 0, beta = 0, K = 0)
  treated_idx <- partition_process == "treated"
  control_idx <- partition_process == "control"
  control_tiles <- which(control_idx)
  treated_tiles <- which(treated_idx)
  inds <- as.numeric(tileindex(pp_data$x, pp_data$y, partition))
  control_inds <- inds[pp_data$inferred_process == "control"]
  treated_inds <- inds[pp_data$inferred_process == "treated"]

  dat <- pp_data
  if (is.null(dat$inferred_process)) {
    dat$inferred_process <- dat$location_process
  }
  if (is.null(hawkes_params_control)) {
    thin_control <- 0
  } else {
    sim_data <- generate_inhomogeneous_hawkes(
      Omega = statespace, partition = partition, time_window = windowT,
      partition_processes = partition_process,
      hawkes_params = list(control = hawkes_params_control, treated = no_points_hawkes),
      filtration = filtration, state_spaces = state_spaces, space_triggering = FALSE, ...
    )
    sim_data <- sim_data[sim_data$t > windowT[1], ]
    sim_inds <- as.numeric(tileindex(sim_data$x, sim_data$y, partition))
    where_to_thin <- tabulate(control_inds, nbins = partition$n) - tabulate(sim_inds, nbins = partition$n)
    thin_control <- length(control_inds) - length(sim_inds)
    if (verbose) {
      print(paste0("Total in data ", sum(dat$inferred_process == "control")))
      print(paste0("Total in sim ", length(sim_inds)))
    }
  }
  thin_control <- thin_control * change_factor
  if (!is.finite(thin_control) || is.na(thin_control)) thin_control <- 0
  if (thin_control < 0) {
    thin_control <- rpois(1, max(-thin_control, 1))
    where_to_thin <- -where_to_thin
    relabel_label <- "treated"
  } else {
    thin_control <- rpois(1, max(thin_control, 1))
    relabel_label <- "control"
  }
  if (verbose) {
    print(paste0("Relabelling ", thin_control, " points from ", relabel_label, " inferred_points"))
  }
  if (thin_control > 0) {
    probs <- where_to_thin
    probs[!is.finite(probs)] <- 0
    probs <- probs - min(0, min(probs))
    if (sum(probs, na.rm = TRUE) <= 0) probs <- rep(1 / length(probs), length(probs))
    probs <- probs / sum(probs, na.rm = TRUE)
    partition_thins <- rmultinom(1, thin_control, prob = probs)
    total_changes <- 0
    for (i in 1:length(partition_thins)) {
      n_thin <- partition_thins[i]
      if (length(n_thin) == 0 || is.null(n_thin) || is.na(n_thin) || n_thin <= 0) next
      points <- which(inds == i)
      if (length(points) == 1) {
        changes <- points
      } else {
        possible_labels <- which(inds == i & dat$inferred_process == relabel_label)
        changes <- sample(possible_labels, size = min(n_thin, length(possible_labels)), replace = FALSE)
      }
      total_changes <- total_changes + length(changes)
      current_labels <- dat$inferred_process[changes]
      new_labels <- c("control", "treated")[1 * (current_labels == "control") + 1]
      dat$inferred_process[changes] <- new_labels
    }
  }

  treated_inds <- inds[dat$inferred_process == "treated"]

  if (is.null(hawkes_params_treated)) {
    thin_treated <- 0
    where_to_thin <- rep(1, length(inds))
  }
  if (!is.null(hawkes_params_treated)) {
    sim_data <- generate_inhomogeneous_hawkes(
      Omega = statespace, partition = partition, time_window = windowT,
      partition_processes = partition_process,
      hawkes_params = list(control = no_points_hawkes, treated = hawkes_params_treated),
      filtration = filtration, state_spaces = state_spaces, space_triggering = FALSE, ...
    )
    sim_data <- sim_data[sim_data$t > windowT[1], ]
    sim_inds <- as.numeric(tileindex(sim_data$x, sim_data$y, partition))
    where_to_thin <- tabulate(treated_inds, nbins = partition$n) - tabulate(sim_inds, nbins = partition$n)
    thin_treated <- length(treated_inds) - length(sim_inds)
    if (verbose) {
      print(paste0("Total in data ", sum(dat$inferred_process == "treated")))
      print(paste0("Total in sim ", length(sim_inds)))
    }
  }
  thin_treated <- thin_treated * change_factor
  if (!is.finite(thin_treated) || is.na(thin_treated)) thin_treated <- 0
  if (thin_treated < 0) {
    thin_treated <- rpois(1, max(-thin_treated, 1))
    where_to_thin <- -where_to_thin
    relabel_label <- "control"
  } else {
    thin_treated <- rpois(1, max(thin_treated, 1))
    relabel_label <- "treated"
  }
  if (verbose) {
    print(paste0("Relabelling ", thin_treated, " points from ", relabel_label, " inferred_points"))
  }
  if (thin_treated > 0) {
    probs <- where_to_thin
    probs[!is.finite(probs)] <- 0
    probs <- probs - min(0, min(probs))
    if (sum(probs, na.rm = TRUE) <= 0) probs <- rep(1 / length(probs), length(probs))
    probs <- probs / sum(probs, na.rm = TRUE)
    probs[control_tiles] <- 0
    if (sum(probs, na.rm = TRUE) <= 0) {
      thin_treated <- 0
      partition_thins <- c()
    } else {
      probs <- probs / sum(probs, na.rm = TRUE)
      partition_thins <- rmultinom(1, thin_treated, prob = probs)
    }
    for (i in 1:length(partition_thins)) {
      n_thin <- partition_thins[i]
      if (length(n_thin) == 0 || is.null(n_thin) || is.na(n_thin) || n_thin <= 0) next
      points <- which(inds == i)
      if (length(points) == 1) {
        changes <- points
      } else {
        possible_labels <- which(inds == i & dat$inferred_process == relabel_label)
        changes <- sample(possible_labels, size = min(n_thin, length(possible_labels)), replace = FALSE)
      }
      current_labels <- dat$location_process[changes]
      new_labels <- c("control", "treated")[1 * (current_labels == "control") + 1]
      dat$inferred_process[changes] <- new_labels
    }
  }
  return(dat)
}

#' Fast labeling proposal by random flipping
#' @param pp_data Data frame with inferred_process column
#' @param partition A spatstat tess object
#' @param partition_process Character vector of process names per tile
#' @param n_thin Number of points to flip
#' @return pp_data with updated inferred_process
#' @export
fast_labelling_proposal <- function(pp_data, partition, partition_process, n_thin) {
  dat <- pp_data
  n_change_treated <- floor(n_thin * sum(pp_data$inferred_process == "treated") / dim(pp_data)[1])
  n_change_control <- n_thin - n_change_treated

  treated_idx <- partition_process == "treated"
  control_idx <- partition_process == "control"
  control_tiles <- which(control_idx)
  treated_tiles <- which(treated_idx)
  inds <- as.numeric(tileindex(pp_data$x, pp_data$y, partition))

  changes <- c(
    sample(which(inds %in% control_tiles), size = n_change_control, replace = FALSE),
    sample(which(inds %in% treated_tiles), size = n_change_treated, replace = FALSE)
  )
  current_labels <- dat$inferred_process[changes]
  new_labels <- c("control", "treated")[1 * (current_labels == "control") + 1]
  dat$inferred_process[changes] <- new_labels
  return(dat)
}

#' Global likelihood-ratio labeling
#' @param pp_data Data frame with point process data
#' @param partition A spatstat tess object
#' @param tiles_to_thin Integer vector of control tile indices
#' @param statespace Full observation window
#' @param windowT Numeric vector c(start, end)
#' @param hawkes_params Hawkes parameter list
#' @param n_traversal Number of labeling proposals to evaluate
#' @return List with labellings and likelihoods
#' @export
global_hawkes_likelihood_ratio_labeling <- function(pp_data, partition, tiles_to_thin,
                                                    statespace, windowT, hawkes_params,
                                                    n_traversal = 1000) {
  no_points_hawkes <- list(mu = 0, alpha = 0, beta = 0, K = 0)
  partition_process <- c("control", "treated")[(!(1:partition$n %in% tiles_to_thin)) * 1 + 1]
  treated_idx <- partition_process == "treated"
  control_state_space <- as.owin(partition[!treated_idx])
  treated_state_space <- as.owin(partition[treated_idx])

  inds <- as.numeric(tileindex(pp_data$x, pp_data$y, partition))
  labellings <- lapply(1:n_traversal, function(i) {
    dat <- pp_data
    dat$inferred_process <- dat$location_process
    thin_total <- sum(dat$location_process == "control") - (hawkes_params$mu * windowT[2] * mean(partition_process == "control"))
    thin_total <- rpois(1, max(thin_total, round(sum(partition_process == "control") / 2)))
    sim_data <- generate_inhomogeneous_hawkes(
      Omega = statespace, partition = partition, time_window = windowT,
      partition_processes = partition_process,
      hawkes_params = list(control = hawkes_params, treated = no_points_hawkes),
      space_triggering = FALSE
    )
    sim_inds <- as.numeric(tileindex(sim_data$x, sim_data$y, partition))
    where_to_thin <- summary(factor(inds, levels = 1:partition$n)) - summary(factor(sim_inds, levels = 1:partition$n))
    probs <- where_to_thin[which(!treated_idx)]
    probs <- probs - min(0, min(probs))
    probs <- probs / sum(probs)
    partition_thins <- rmultinom(1, thin_total, prob = probs)

    for (ii in 1:length(tiles_to_thin)) {
      j <- tiles_to_thin[ii]
      n_thin <- partition_thins[ii]
      if (length(n_thin) == 0 || is.null(n_thin) || is.na(n_thin) || n_thin <= 0) next
      points <- which(inds == j)
      if (length(points) == 1) {
        changes <- points
      } else {
        changes <- sample(which(inds == j), size = min(n_thin, sum(inds == j)), replace = FALSE)
      }
      current_labels <- dat$location_process[changes]
      new_labels <- c("control", "treated")[1 * (current_labels == "control") + 1]
      dat$inferred_process[changes] <- new_labels
    }
    return(dat)
  })
  likelihoods <- sapply(1:n_traversal, function(i) {
    data <- labellings[[i]][labellings[[i]]$inferred_process == "control", ]
    lik <- loglik_hawk_fast(
      params = unlist(hawkes_params), realiz = data,
      windowT = windowT, windowS = statespace,
      zero_background_region = treated_state_space,
      density_approx = FALSE, numeric_integral = FALSE
    )
    return(lik)
  })
  return(list(labellings = labellings, likelihoods = likelihoods))
}

#' Local likelihood-ratio labeling
#' @param pp_data Data frame with point process data
#' @param partition A spatstat tess object
#' @param tiles_to_thin Integer vector of control tile indices
#' @param statespace Full observation window
#' @param windowT Numeric vector c(start, end)
#' @param hawkes_params Hawkes parameter list
#' @param n_traversal Number of labeling proposals
#' @param sample_size Number of partitions for local likelihood
#' @return List with labellings and likelihoods
#' @export
local_hawkes_likelihood_ratio_labeling <- function(pp_data, partition, tiles_to_thin,
                                                   statespace, windowT, hawkes_params,
                                                   n_traversal = 1000, sample_size = 10) {
  no_points_hawkes <- list(mu = 0, alpha = 0, beta = 0, K = 0)
  partition_process <- c("control", "treated")[(!(1:partition$n %in% tiles_to_thin)) * 1 + 1]
  treated_idx <- partition_process == "treated"
  control_state_space <- as.owin(partition[!treated_idx])
  treated_state_space <- as.owin(partition[treated_idx])

  inds <- as.numeric(tileindex(pp_data$x, pp_data$y, partition))
  labellings <- lapply(1:n_traversal, function(i) {
    dat <- pp_data
    dat$inferred_process <- dat$location_process
    thin_total <- sum(dat$location_process == "control") - (hawkes_params$mu * windowT[2] * mean(partition_process == "control"))
    thin_total <- rpois(1, max(thin_total, round(sum(partition_process == "control") / 2)))
    sim_data <- generate_inhomogeneous_hawkes(
      Omega = statespace, partition = partition, time_window = windowT,
      partition_processes = partition_process,
      hawkes_params = list(control = hawkes_params, treated = no_points_hawkes),
      space_triggering = FALSE
    )
    sim_inds <- as.numeric(tileindex(sim_data$x, sim_data$y, partition))
    where_to_thin <- summary(factor(inds, levels = 1:partition$n)) - summary(factor(sim_inds, levels = 1:partition$n))
    probs <- where_to_thin[which(!treated_idx)]
    probs <- probs - min(0, min(probs))
    probs <- probs / sum(probs)
    partition_thins <- rmultinom(1, thin_total, prob = probs)

    for (ii in 1:length(tiles_to_thin)) {
      j <- tiles_to_thin[ii]
      n_thin <- partition_thins[ii]
      if (length(n_thin) == 0 || is.null(n_thin) || is.na(n_thin) || n_thin <= 0) next
      points <- which(inds == j)
      if (length(points) == 1) {
        changes <- points
      } else {
        changes <- sample(which(inds == j), size = min(n_thin, sum(inds == j)), replace = FALSE)
      }
      current_labels <- dat$location_process[changes]
      new_labels <- c("control", "treated")[1 * (current_labels == "control") + 1]
      dat$inferred_process[changes] <- new_labels
    }
    return(dat)
  })
  likelihoods <- sapply(1:n_traversal, function(i) {
    data <- labellings[[i]][labellings[[i]]$inferred_process == "control", ]
    local_inds <- as.numeric(tileindex(data$x, data$y, partition))
    lik <- sapply(sample(unique(local_inds), min(length(unique(local_inds)), sample_size)), function(ii) {
      dat_local <- data[local_inds == ii, ]
      loglik_hawk_fast(
        params = unlist(hawkes_params), realiz = dat_local,
        windowT = windowT, windowS = as.owin(partition[ii]),
        zero_background_region = NULL,
        density_approx = FALSE, numeric_integral = FALSE
      )
    })
    return(sum(lik))
  })
  return(list(labellings = labellings, likelihoods = likelihoods))
}

#' Simulation-based fast labeling with proximity weighting
#'
#' Uses discrepancy-guided proposals with distance-based weighting for
#' selecting which points to relabel near process boundaries.
#'
#' @param pp_data Data frame with columns x, y, t, inferred_process
#' @param partition A spatstat tess object
#' @param partition_process Character vector of process names per tile
#' @param statespace Full observation window (owin)
#' @param windowT Numeric vector c(start, end)
#' @param state_spaces Optional precomputed state spaces
#' @param hawkes_params_control Control Hawkes parameters
#' @param hawkes_params_treated Treated Hawkes parameters
#' @param change_factor Scaling factor for relabeling count
#' @param filtration Optional pre-treatment history
#' @param verbose Print progress
#' @param partition_mask Optional mask for fast tile lookup
#' @param proximity_weight Weight for distance-based vs uniform sampling (0-1)
#' @param points_tile_index Optional precomputed tile index for \code{pp_data} (same length as rows); skips \code{tileindex()} when provided
#' @param ... Additional arguments passed to generate_inhomogeneous_hawkes
#' @return pp_data with updated inferred_process
#' @export
#'
#' Generate all single-flip labelling proposals
#'
#' For each post-treatment point, create one proposal where that point's
#' label is toggled (control <-> treated). Returns N proposals for N post
#' points -- exhaustive over all single-point moves.
#'
#' @param post Data frame of post-treatment points (with inferred_process)
#' @param pre Data frame of pre-treatment points (all labelled "control")
#' @return List of N data frames, each with one point's label flipped
#' @export
single_flip_proposals <- function(post, pre) {
  n <- nrow(post)
  lapply(seq_len(n), function(j) {
    proposed <- post
    proposed$inferred_process[j] <- if (post$inferred_process[j] == "control") {
      "treated"
    } else {
      "control"
    }
    rbind(pre, proposed)
  })
}

simulation_labeling_hawkes_hawkes_fast <- function(pp_data,
                                                   partition,
                                                   partition_process,
                                                   statespace,
                                                   windowT,
                                                   state_spaces = NULL,
                                                   hawkes_params_control = NULL,
                                                   hawkes_params_treated = NULL,
                                                   change_factor = 1,
                                                   filtration = NULL,
                                                   verbose = FALSE,
                                                   partition_mask = NULL,
                                                   proximity_weight = 0.9,
                                                   temporal_weight = 0,
                                                   temporal_scale_days = NULL,
                                                   points_tile_index = NULL,
                                                   model_type = "hawkes",
                                                   ...) {
  proposal_trace <- isTRUE(verbose) &&
    (tolower(Sys.getenv("OK_SEM_PROPOSAL_TIMING", "true")) %in% c("1", "true", "yes", "y"))
  dat <- pp_data
  if (is.null(dat$inferred_process)) {
    dat$inferred_process <- dat$location_process
  }

  n_pts <- nrow(dat)
  if (!is.null(points_tile_index) && length(points_tile_index) == n_pts) {
    inds <- as.numeric(points_tile_index)
  } else if (!is.null(partition_mask)) {
    inds <- as.integer(tileindex(dat$x, dat$y, partition))
  } else {
    inds <- as.numeric(tileindex(dat$x, dat$y, partition))
  }

  control_inds <- inds[dat$inferred_process == "control"]
  treated_inds <- inds[dat$inferred_process == "treated"]

  dots <- list(...)
  is_etas <- identical(model_type, "etas")
  is_biv_etas <- identical(model_type, "etas_bivariate")
  sem_prop_branching_vals <- suppressWarnings(as.numeric(
    Sys.getenv("OK_SEM_PROPOSAL_BRANCHING_MAX", "0.98")
  ))
  sem_prop_branching_max <- sem_prop_branching_vals[[1]]
  if (!is.finite(sem_prop_branching_max) || is.na(sem_prop_branching_max) || sem_prop_branching_max <= 0) {
    sem_prop_branching_max <- NA_real_
  }
  beta_guard_vals <- suppressWarnings(as.numeric(dots$beta_gr))
  beta_guard <- beta_guard_vals[which(is.finite(beta_guard_vals) & !is.na(beta_guard_vals) & beta_guard_vals > 0)][1]
  if (!is.finite(beta_guard) || is.na(beta_guard) || beta_guard <= 0) beta_guard <- 1.0
  spectral_radius_proxy <- function(par_obj, beta_gr) {
    get_num <- function(nm, default = NA_real_) {
      v <- suppressWarnings(as.numeric(par_obj[[nm]]))
      if (length(v) < 1L || !is.finite(v[[1]])) return(default)
      v[[1]]
    }
    b00 <- get_num("A_00", 0) * exp(get_num("alpha_m_00", 0) / beta_gr)
    b11 <- get_num("A_11", 0) * exp(get_num("alpha_m_11", 0) / beta_gr)
    b01 <- get_num("A_01", 0) * exp(get_num("alpha_m_01", 0) / beta_gr)
    b10 <- get_num("A_10", 0) * exp(get_num("alpha_m_10", 0) / beta_gr)
    if (!all(is.finite(c(b00, b11, b01, b10)))) return(Inf)
    disc <- (b00 - b11)^2 + 4 * b01 * b10
    if (!is.finite(disc) || disc < 0) return(Inf)
    0.5 * ((b00 + b11) + sqrt(disc))
  }
  if (!is.finite(temporal_weight)) temporal_weight <- 0
  temporal_weight <- min(max(temporal_weight, 0), 1)
  if (is.null(temporal_scale_days) || !is.finite(temporal_scale_days) || temporal_scale_days <= 0) {
    temporal_scale_days <- max(1, as.numeric(windowT[2] - windowT[1]) / 5)
  }
  use_precompute <- (proximity_weight > 0 || temporal_weight > 0)
  inds_int <- as.integer(inds)
  tile_members <- NULL
  if (use_precompute) {
    valid_idx <- which(is.finite(inds_int) & inds_int >= 1L & inds_int <= partition$n)
    tile_members <- split(valid_idx, factor(inds_int[valid_idx], levels = seq_len(partition$n)))
  }
  time_weight_base <- NULL
  if (use_precompute && temporal_weight > 0) {
    dt_all <- pmax(dat$t - windowT[1], 0)
    tw_all <- exp(-dt_all / temporal_scale_days)
    if (all(is.finite(tw_all)) && sum(tw_all) > 0) {
      time_weight_base <- tw_all
    }
  }

  # For bivariate ETAS, simulate jointly and compare per-process counts
  if (is_biv_etas) {
    if (proposal_trace) {
      cat(sprintf("    [proposal-sim] start model=etas_bivariate n_pts=%d t_window=[%.3f, %.3f]\n",
                  nrow(dat), windowT[1], windowT[2]))
    }
    biv_params <- dots$etas_bivariate_params
    if (is.null(biv_params)) {
      biv_params <- init_bivariate_from_independent(
        hawkes_params_control, hawkes_params_treated)
    }
    if (is.finite(sem_prop_branching_max)) {
      rho_hat <- spectral_radius_proxy(biv_params, beta_guard)
      if (is.finite(rho_hat) && rho_hat > sem_prop_branching_max) {
        scale_fac <- sem_prop_branching_max / rho_hat
        for (nm in c("A_00", "A_11", "A_01", "A_10")) {
          biv_params[[nm]] <- as.numeric(biv_params[[nm]]) * scale_fac
        }
        if (proposal_trace) {
          cat(sprintf(
            "    [proposal-sim] branching guard: rho=%.3f > %.3f, scaling A by %.4f\n",
            rho_hat, sem_prop_branching_max, scale_fac
          ))
        }
      }
    }
    t_sim <- proc.time()[3]
    sim_data <- generate_inhomogeneous_etas_bivariate(
      Omega = statespace, partition = partition, time_window = windowT,
      partition_processes = partition_process,
      etas_bivariate_params = biv_params,
      m0 = dots$m0, beta_gr = dots$beta_gr, mag_pool = dots$mag_pool,
      state_spaces = state_spaces, filtration = filtration,
      t_trunc = dots$t_trunc
    )
    if (proposal_trace) {
      cat(sprintf("    [proposal-sim] generated n=%d in %.2fs (t_trunc=%s)\n",
                  nrow(sim_data), proc.time()[3] - t_sim,
                  ifelse(is.null(dots$t_trunc), "none", as.character(signif(as.numeric(dots$t_trunc), 6)))))
    }
    if (nrow(sim_data) > 0L && any(sim_data$t <= windowT[1])) {
      sim_data <- sim_data[sim_data$t > windowT[1], , drop = FALSE]
    }
    sim_inds <- if (!is.null(sim_data$tile_index)) {
      as.numeric(sim_data$tile_index)
    } else {
      as.numeric(tileindex(sim_data$x, sim_data$y, partition))
    }
    if ("process_id" %in% names(sim_data)) {
      sim_pid <- as.integer(sim_data$process_id)
      sim_ctrl_inds <- sim_inds[sim_pid == 0L]
      sim_treat_inds <- sim_inds[sim_pid == 1L]
    } else {
      sim_proc <- if ("process" %in% names(sim_data)) sim_data$process else
                    ifelse(sim_data$process_id == 0, "control", "treated")
      sim_ctrl_inds <- sim_inds[sim_proc == "control"]
      sim_treat_inds <- sim_inds[sim_proc == "treated"]
    }

    where_to_thin <- tabulate(control_inds, nbins = partition$n) -
                     tabulate(sim_ctrl_inds, nbins = partition$n)
    thin_control <- length(control_inds) - length(sim_ctrl_inds)
  } else {

  if (is_etas) {
    no_points <- list(mu = 0, A = 0, alpha_m = 0, c = 0.1, p = 1.2,
                      D = 1, gamma = 0, q = 1.5)
    gen_fn <- function(...) generate_inhomogeneous_etas(...)
  } else {
    no_points <- list(mu = 0, alpha = 0, beta = 0, K = 0)
    gen_fn <- function(...) generate_inhomogeneous_hawkes(...)
  }
  params_key <- if (is_etas) "etas_params" else "hawkes_params"

  if (is.null(hawkes_params_control)) {
    thin_control <- 0
    sim_inds <- numeric(0)
  } else {
    if (proposal_trace) {
      cat(sprintf("    [proposal-sim] start model=%s(control-only)\n",
                  ifelse(is_etas, "etas", "hawkes")))
    }
    gen_args <- list(
      Omega = statespace, partition = partition, time_window = windowT,
      partition_processes = partition_process,
      filtration = filtration, state_spaces = state_spaces
    )
    gen_args[[params_key]] <- list(control = hawkes_params_control,
                                   treated = no_points)
    if (!is_etas) gen_args$space_triggering <- FALSE
    if (is_etas) {
      gen_args$m0 <- dots$m0
      gen_args$beta_gr <- dots$beta_gr
      gen_args$mag_pool <- dots$mag_pool
    }
    extra <- dots[!names(dots) %in% c("m0", "beta_gr", "mag_pool")]
    gen_args <- c(gen_args, extra)
    t_sim_ctrl <- proc.time()[3]
    sim_data <- do.call(gen_fn, gen_args)
    if (proposal_trace) {
      cat(sprintf("    [proposal-sim] control-side generated n=%d in %.2fs\n",
                  nrow(sim_data), proc.time()[3] - t_sim_ctrl))
    }
    sim_data <- sim_data[sim_data$t > windowT[1], ]
    sim_inds <- if (!is.null(sim_data$tile_index)) as.numeric(sim_data$tile_index) else as.numeric(tileindex(sim_data$x, sim_data$y, partition))
    where_to_thin <- tabulate(control_inds, nbins = partition$n) - tabulate(sim_inds, nbins = partition$n)
    thin_control <- length(control_inds) - length(sim_inds)
  }

  } # end non-bivariate branch

  thin_control <- thin_control * change_factor
  if (!is.finite(thin_control) || is.na(thin_control)) thin_control <- 0

  if (thin_control < 0) {
    n_total_change <- rpois(1, max(-thin_control, 1))
    where_to_thin <- -where_to_thin
    relabel_label <- "treated"
    target_label <- "control"
  } else {
    n_total_change <- rpois(1, max(thin_control, 1))
    relabel_label <- "control"
    target_label <- "treated"
  }

  if (verbose) {
    print(paste0("Relabelling approx ", n_total_change, " points from ", relabel_label, " to ", target_label))
  }

  if (n_total_change > 0) {
    probs <- where_to_thin
    probs[!is.finite(probs)] <- 0
    probs <- probs - min(0, min(probs))
    if (sum(probs, na.rm = TRUE) <= 0) probs <- rep(1 / length(probs), length(probs))
    probs <- probs / sum(probs, na.rm = TRUE)
    partition_thins <- rmultinom(1, n_total_change, prob = probs)

    attractor_idx <- which(dat$inferred_process == target_label)
    has_attractors <- length(attractor_idx) > 0
    geom_weight_raw <- NULL

    if (has_attractors) {
      win_box <- owin(range(dat$x), range(dat$y))
      X_attract <- ppp(dat$x[attractor_idx], dat$y[attractor_idx], window = win_box, check = FALSE)
      if (use_precompute && proximity_weight > 0) {
        relabel_idx <- which(dat$inferred_process == relabel_label)
        if (length(relabel_idx) > 0) {
          X_relabel <- ppp(dat$x[relabel_idx], dat$y[relabel_idx], window = win_box, check = FALSE)
          dists_all <- nncross(X_relabel, X_attract)$dist
          geom_weight_raw <- numeric(n_pts)
          geom_weight_raw[relabel_idx] <- 1 / (dists_all + 1e-4)
        }
      }
    }

    for (i in seq_len(length(partition_thins))) {
      n_thin <- partition_thins[i]
      if (length(n_thin) == 0 || is.na(n_thin) || n_thin <= 0) next

      in_tile_idx <- if (use_precompute) tile_members[[i]] else which(inds == i)
      if (is.null(in_tile_idx) || length(in_tile_idx) == 0) next
      candidates <- in_tile_idx[dat$inferred_process[in_tile_idx] == relabel_label]
      if (length(in_tile_idx) == 1L && n_thin > 0) {
        changes <- in_tile_idx
        dat$inferred_process[changes] <- if (dat$inferred_process[changes] == "control") "treated" else "control"
        next
      }
      if (length(candidates) == 0) next

      sampling_weights <- NULL
      base_weights <- NULL
      if (has_attractors && length(candidates) > 0 && proximity_weight > 0) {
        if (use_precompute && !is.null(geom_weight_raw)) {
          geom_weights <- geom_weight_raw[candidates]
          if (all(is.finite(geom_weights)) && sum(geom_weights) > 0) {
            geom_weights <- geom_weights / sum(geom_weights)
          } else {
            geom_weights <- rep(1 / length(candidates), length(candidates))
          }
        } else {
          X_cand <- ppp(dat$x[candidates], dat$y[candidates], window = win_box, check = FALSE)
          dists <- nncross(X_cand, X_attract)$dist
          geom_weights <- 1 / (dists + 1e-4)
          if (sum(geom_weights) > 0) {
            geom_weights <- geom_weights / sum(geom_weights)
          } else {
            geom_weights <- rep(1 / length(candidates), length(candidates))
          }
        }
        uni_weights <- rep(1 / length(candidates), length(candidates))
        base_weights <- (proximity_weight * geom_weights) + ((1 - proximity_weight) * uni_weights)
      } else if (length(candidates) > 0) {
        base_weights <- rep(1 / length(candidates), length(candidates))
      }

      if (length(candidates) > 0 && temporal_weight > 0) {
        time_weights <- if (use_precompute && !is.null(time_weight_base)) {
          time_weight_base[candidates]
        } else {
          dt <- pmax(dat$t[candidates] - windowT[1], 0)
          exp(-dt / temporal_scale_days)
        }
        if (!any(is.finite(time_weights)) || sum(time_weights, na.rm = TRUE) <= 0) {
          time_weights <- rep(1 / length(candidates), length(candidates))
        } else {
          time_weights[!is.finite(time_weights)] <- 0
          time_weights <- time_weights / sum(time_weights)
        }
        if (is.null(base_weights)) {
          base_weights <- rep(1 / length(candidates), length(candidates))
        }
        sampling_weights <- ((1 - temporal_weight) * base_weights) + (temporal_weight * time_weights)
        sampling_weights <- sampling_weights / sum(sampling_weights)
      } else {
        sampling_weights <- base_weights
      }

      changes_local <- sample(
        length(candidates),
        size = min(n_thin, length(candidates)),
        replace = FALSE, prob = sampling_weights
      )
      changes <- candidates[changes_local]
      dat$inferred_process[changes] <- target_label
    }
  }

  if (is.null(hawkes_params_treated) && !is_biv_etas) return(dat)

  treated_inds <- inds[dat$inferred_process == "treated"]
  control_tiles <- which(partition_process == "control")

  if (is_biv_etas) {
    sim_inds_t <- sim_treat_inds
    where_to_thin_t <- tabulate(treated_inds, nbins = partition$n) -
                       tabulate(sim_inds_t, nbins = partition$n)
  } else {
  gen_args_t <- list(
    Omega = statespace, partition = partition, time_window = windowT,
    partition_processes = partition_process,
    filtration = filtration, state_spaces = state_spaces
  )
  gen_args_t[[params_key]] <- list(control = no_points,
                                    treated = hawkes_params_treated)
  if (!is_etas) gen_args_t$space_triggering <- FALSE
  if (is_etas) {
    gen_args_t$m0 <- dots$m0
    gen_args_t$beta_gr <- dots$beta_gr
    gen_args_t$mag_pool <- dots$mag_pool
  }
  extra_t <- dots[!names(dots) %in% c("m0", "beta_gr", "mag_pool")]
  gen_args_t <- c(gen_args_t, extra_t)
  if (proposal_trace) {
    cat(sprintf("    [proposal-sim] start model=%s(treated-only)\n",
                ifelse(is_etas, "etas", "hawkes")))
  }
  t_sim_treat <- proc.time()[3]
  sim_data_t <- do.call(gen_fn, gen_args_t)
  if (proposal_trace) {
    cat(sprintf("    [proposal-sim] treated-side generated n=%d in %.2fs\n",
                nrow(sim_data_t), proc.time()[3] - t_sim_treat))
  }
  sim_data_t <- sim_data_t[sim_data_t$t > windowT[1], ]
  sim_inds_t <- if (!is.null(sim_data_t$tile_index)) as.numeric(sim_data_t$tile_index) else as.numeric(tileindex(sim_data_t$x, sim_data_t$y, partition))
  where_to_thin_t <- tabulate(treated_inds, nbins = partition$n) - tabulate(sim_inds_t, nbins = partition$n)
  } # end non-bivariate branch
  thin_treated <- (length(treated_inds) - length(sim_inds_t)) * change_factor
  if (!is.finite(thin_treated) || is.na(thin_treated)) thin_treated <- 0
  if (thin_treated < 0) {
    n_total_t <- rpois(1, max(-thin_treated, 1))
    where_to_thin_t <- -where_to_thin_t
    relabel_label_t <- "control"
    target_label_t <- "treated"
  } else {
    n_total_t <- rpois(1, max(thin_treated, 1))
    relabel_label_t <- "treated"
    target_label_t <- "control"
  }
  if (n_total_t > 0) {
    probs_t <- where_to_thin_t
    probs_t[!is.finite(probs_t)] <- 0
    probs_t <- probs_t - min(0, min(probs_t))
    if (sum(probs_t, na.rm = TRUE) <= 0) probs_t <- rep(1 / length(probs_t), length(probs_t))
    probs_t <- probs_t / sum(probs_t, na.rm = TRUE)
    probs_t[control_tiles] <- 0
    if (sum(probs_t, na.rm = TRUE) > 0) {
      probs_t <- probs_t / sum(probs_t, na.rm = TRUE)
      partition_thins_t <- rmultinom(1, n_total_t, prob = probs_t)
      for (i in seq_len(length(partition_thins_t))) {
        n_thin <- partition_thins_t[i]
        if (length(n_thin) == 0 || is.na(n_thin) || n_thin <= 0) next
        in_tile_idx <- if (use_precompute) tile_members[[i]] else which(inds == i)
        if (is.null(in_tile_idx) || length(in_tile_idx) == 0) next
        candidates_t <- in_tile_idx[dat$inferred_process[in_tile_idx] == relabel_label_t]
        if (length(in_tile_idx) == 1L && n_thin > 0) {
          changes <- in_tile_idx
          dat$inferred_process[changes] <- if (dat$location_process[changes] == "control") "treated" else "control"
          next
        }
        if (length(candidates_t) == 0) next
        sampling_weights_t <- NULL
        if (temporal_weight > 0) {
          tw_t <- if (use_precompute && !is.null(time_weight_base)) {
            time_weight_base[candidates_t]
          } else {
            dt_t <- pmax(dat$t[candidates_t] - windowT[1], 0)
            exp(-dt_t / temporal_scale_days)
          }
          if (all(is.finite(tw_t)) && sum(tw_t) > 0) {
            sampling_weights_t <- tw_t / sum(tw_t)
          }
        }
        changes_local <- sample(length(candidates_t),
                                size = min(n_thin, length(candidates_t)),
                                replace = FALSE, prob = sampling_weights_t)
        changes <- candidates_t[changes_local]
        dat$inferred_process[changes] <- ifelse(dat$location_process[changes] == "control", "treated", "control")
      }
    }
  }
  return(dat)
}

#' EM-style iterative labeling with parameter updates
#'
#' Core inner loop of the SEM algorithm: generates labeling proposals,
#' selects the best by likelihood, and optionally re-estimates Hawkes parameters.
#'
#' @param pp_data Data frame with columns x, y, t, process, location_process
#' @param partition A spatstat tess object
#' @param partition_processes Character vector of process names per tile
#' @param statespace Full observation window (owin)
#' @param time_window Numeric vector c(start, end)
#' @param treatment_time Scalar treatment time
#' @param hawkes_params_control Control Hawkes parameters
#' @param hawkes_params_treated Treated Hawkes parameters
#' @param update_control_params Logical; update control params too
#' @param param_update_cadence How often to refit params (every N iterations)
#' @param proposal_update_cadence How often to regenerate proposals
#' @param update_starting_data Logical; update starting data each iteration
#' @param include_starting_data Logical; include current data in proposals
#' @param include_starting_first_n Integer; force-include starting data as a
#'   candidate proposal for the first N inner iterations.
#' @param max_relabel_step_frac Maximum fraction of post-treatment points that
#'   can change labels in a single inner iteration.
#' @param force_param_update_flip_frac Cumulative accepted-flip fraction
#'   threshold that forces a parameter update once reached.
#' @param optim_method One of "max", "sample_weighted", "mean", "truncated_mean"
#' @param selection_temperature Positive temperature for likelihood-weighted
#'   sampling when \code{optim_method = "sample_weighted"}. Lower values are
#'   more concentrated near the maximum-likelihood proposal.
#' @param change_factor_min_mult Lower bound multiplier for adaptive
#'   \code{change_factor} relative to the initial value.
#' @param change_factor_max_mult Upper bound multiplier for adaptive
#'   \code{change_factor} relative to the initial value.
#' @param state_spaces Optional precomputed state spaces
#' @param metric_name Metric for selecting best labeling
#' @param iter Number of iterations
#' @param n_props Number of proposals per iteration
#' @param change_factor Scaling factor for relabeling
#' @param stagnation_trigger_every Trigger exploration every N consecutive
#'   no-flip iterations.
#' @param MCMC_style Use MCMC acceptance/rejection
#' @param verbose Print progress
#' @param ... Additional arguments
#' @return List with labelling, treated_par, control_par, accuracies, etc.
#' @export
em_style_labelling <- function(pp_data,
                               partition,
                               partition_processes,
                               statespace,
                               time_window,
                               treatment_time,
                               hawkes_params_control,
                               hawkes_params_treated,
                               update_control_params = FALSE,
                               param_update_cadence = 20,
                               proposal_update_cadence = 1,
                               update_starting_data = TRUE,
                               include_starting_data = TRUE,
                               include_starting_first_n = 50,
                               max_relabel_step_frac = 1.0,
                               force_param_update_flip_frac = 1.0,
                               optim_method = "max",
                               selection_temperature = 0.15,
                               change_factor_min_mult = 0.2,
                               change_factor_max_mult = 2.0,
                               state_spaces = NULL,
                               metric_name = "post_likelihood",
                               iter = 100,
                               n_props = 100,
                               change_factor = 0.1,
                               stagnation_trigger_every = 10,
                               MCMC_style = FALSE,
                               proposal_method = "simulation",
                               fixed_params = NULL,
                               verbose = FALSE,
                               model_type = "hawkes",
                               temporal_weight = 0,
                               temporal_scale_days = NULL,
                               ...) {
  sem_timing_verbose <- tolower(Sys.getenv("OK_SEM_TIMING_VERBOSE", "true")) %in% c("1", "true", "yes", "y")
  sem_proposal_verbose <- tolower(Sys.getenv("OK_SEM_PROPOSAL_VERBOSE", "true")) %in% c("1", "true", "yes", "y")
  dots <- list(...)
  base_change_factor <- as.numeric(change_factor)
  if (!is.finite(base_change_factor) || base_change_factor <= 0) {
    stop("change_factor must be a positive finite number.")
  }
  stagnation_trigger_every <- as.integer(stagnation_trigger_every)
  if (is.na(stagnation_trigger_every) || stagnation_trigger_every < 1L) {
    stop("stagnation_trigger_every must be an integer >= 1.")
  }
  include_starting_first_n <- suppressWarnings(as.integer(include_starting_first_n))
  if (!is.finite(include_starting_first_n) || is.na(include_starting_first_n) || include_starting_first_n < 0L) {
    include_starting_first_n <- 0L
  }
  max_relabel_step_frac <- suppressWarnings(as.numeric(max_relabel_step_frac))
  if (!is.finite(max_relabel_step_frac) || is.na(max_relabel_step_frac) || max_relabel_step_frac <= 0) {
    max_relabel_step_frac <- 1.0
  }
  max_relabel_step_frac <- min(max_relabel_step_frac, 1.0)
  force_param_update_flip_frac <- suppressWarnings(as.numeric(force_param_update_flip_frac))
  if (!is.finite(force_param_update_flip_frac) || is.na(force_param_update_flip_frac) || force_param_update_flip_frac <= 0) {
    force_param_update_flip_frac <- 1.0
  }
  force_param_update_flip_frac <- min(force_param_update_flip_frac, 1.0)
  selection_temperature <- suppressWarnings(as.numeric(selection_temperature))
  if (!is.finite(selection_temperature) || is.na(selection_temperature) || selection_temperature <= 0) {
    selection_temperature <- 0.15
  }
  change_factor_min_mult <- suppressWarnings(as.numeric(change_factor_min_mult))
  change_factor_max_mult <- suppressWarnings(as.numeric(change_factor_max_mult))
  if (!is.finite(change_factor_min_mult) || is.na(change_factor_min_mult) || change_factor_min_mult <= 0) {
    change_factor_min_mult <- 0.2
  }
  if (!is.finite(change_factor_max_mult) || is.na(change_factor_max_mult) || change_factor_max_mult < change_factor_min_mult) {
    change_factor_max_mult <- max(2.0, change_factor_min_mult)
  }
  # Keep adaptive proposal size changes moderate around the initial setting.
  change_factor_min <- change_factor_min_mult * base_change_factor
  change_factor_max <- change_factor_max_mult * base_change_factor
  background_rate_var <- if ("background_rate_var" %in% names(dots)) dots$background_rate_var else NULL
  hawkes_use_filtration_history <- if ("hawkes_use_filtration_history" %in% names(dots)) {
    isTRUE(dots$hawkes_use_filtration_history)
  } else {
    TRUE
  }
  treated_background_zero_before <- if ("treated_background_zero_before" %in% names(dots)) {
    as.numeric(dots$treated_background_zero_before)
  } else {
    NULL
  }
  t_trunc <- if ("t_trunc" %in% names(dots)) dots$t_trunc else NULL
  sem_outer_maxit <- if ("outer_maxit" %in% names(dots)) {
    suppressWarnings(as.integer(dots$outer_maxit))
  } else {
    NA_integer_
  }
  if (!is.finite(sem_outer_maxit) || is.na(sem_outer_maxit) || sem_outer_maxit < 1L) {
    sem_outer_maxit <- 1000L
  }
  sem_outer_maxit_biv <- if ("outer_maxit_biv" %in% names(dots)) {
    suppressWarnings(as.integer(dots$outer_maxit_biv))
  } else {
    NA_integer_
  }
  if (!is.finite(sem_outer_maxit_biv) || is.na(sem_outer_maxit_biv) || sem_outer_maxit_biv < 1L) {
    sem_outer_maxit_biv <- sem_outer_maxit
  }
  is_etas <- identical(model_type, "etas")
  is_biv_etas <- identical(model_type, "etas_bivariate")
  biv_etas_params <- dots$etas_bivariate_params

  all_names <- if (is_biv_etas) .etas_bivariate_par_names
               else if (is_etas) .etas_par_names
               else c("mu", "alpha", "beta", "K")
  fixed_idx <- if (!is.null(fixed_params)) {
    idx <- match(names(fixed_params), all_names)
    idx[!is.na(idx)]
  } else integer(0)
  free_idx  <- setdiff(seq_along(all_names), fixed_idx)

  loglik_fn <- if (is_biv_etas) loglik_etas_bivariate
               else if (is_etas) loglik_etas
               else loglik_hawk_fast

  class_func <- function(labelling) {
    if (is.null(labelling$process)) {
      return(list())
    }
    truth_factor <- factor(labelling$process, levels = c("control", "treated"))
    inferred_factor <- factor(labelling$inferred_process, levels = c("control", "treated"))
    if (length(truth_factor) != length(inferred_factor)) return(list())
    tab <- table(truth_factor, inferred_factor)
    l <- as.list(tab)
    names(l) <- paste0(rep(rownames(tab), times = 2), "_as_", rep(colnames(tab), each = 2))
    l
  }
  class_results <- list()

  if (!inherits(statespace, "owin")) statespace <- as.owin(statespace)
  treated_idx <- (partition_processes == "treated")
  treated_state_space <- as.owin(partition[treated_idx])
  control_state_space <- as.owin(partition[!treated_idx])

  if (is.null(pp_data$inferred_process)) {
    pp_data$inferred_process <- pp_data$location_process
  }

  starting_data <- as.data.frame(pp_data)
  treated_par <- list(hawkes_params_treated)
  control_par <- list(hawkes_params_control)
  accuracies <- numeric(iter)
  average_flips <- numeric(iter)
  max_metric_flips <- numeric(iter)
  metric_vec <- numeric(iter)
  change_factor_trace <- numeric(iter)
  retained_starting_trace <- logical(iter)
  all_accuracies <- vector("list", iter)
  all_metrics <- vector("list", iter)

  is_pre <- starting_data$t < treatment_time
  pre_data <- starting_data[is_pre, , drop = FALSE]
  post_data <- starting_data[!is_pre, , drop = FALSE]
  post_data <- post_data[order(post_data$t), , drop = FALSE]
  n_post_total <- nrow(post_data)
  max_flips_per_step <- max(1L, as.integer(ceiling(max_relabel_step_frac * max(1L, n_post_total))))
  force_param_update_flip_n <- max(1L, as.integer(ceiling(force_param_update_flip_frac * max(1L, n_post_total))))
  accepted_flips_cum <- 0L
  next_forced_param_update_at <- force_param_update_flip_n
  pre_data$location_process <- "control"
  pre_data$inferred_process <- "control"
  history_window_start <- if (nrow(pre_data) > 0) min(pre_data$t) else time_window[1]
  history_window <- c(history_window_start, time_window[2])
  current_metric_cache <- NA_real_
  select_pre_history_by_label <- function(label) {
    if (nrow(pre_data) < 1L) return(pre_data[0, c("x", "y", "t"), drop = FALSE])
    src <- if ("inferred_process" %in% names(pre_data)) pre_data$inferred_process else pre_data$location_process
    out <- pre_data[src == label, c("x", "y", "t"), drop = FALSE]
    out <- out[out$t < time_window[1], , drop = FALSE]
    out[order(out$t), , drop = FALSE]
  }
  hawkes_conditional_loglik <- function(params_vec, post_realiz, zero_background_region, pre_hist) {
    if (nrow(post_realiz) < 1L) return(-Inf)
    post_realiz <- post_realiz[order(post_realiz$t), , drop = FALSE]
    p <- suppressWarnings(as.numeric(params_vec[c("mu", "alpha", "beta", "K")]))
    names(p) <- c("mu", "alpha", "beta", "K")
    if (any(!is.finite(p))) return(-Inf)
    if (p[["mu"]] < 0 || p[["alpha"]] < 0 || p[["beta"]] <= 0 || p[["K"]] < 0 || p[["K"]] >= 1) return(-Inf)
    if (is.null(pre_hist)) pre_hist <- post_realiz[0, c("x", "y", "t"), drop = FALSE]
    pre_hist <- pre_hist[pre_hist$t < time_window[1], c("x", "y", "t"), drop = FALSE]
    pre_hist <- pre_hist[order(pre_hist$t), , drop = FALSE]
    W_post <- if (!is.null(background_rate_var) && background_rate_var %in% names(post_realiz)) {
      as.numeric(post_realiz[[background_rate_var]])
    } else {
      rep(1, nrow(post_realiz))
    }
    W_post[!is.finite(W_post)] <- 0
    total_area <- spatstat.geom::area(statespace)
    active_area <- total_area
    if (!is.null(zero_background_region)) {
      zero_area <- spatstat.geom::area(zero_background_region)
      active_area <- max(1e-12, total_area - zero_area)
      in_zero <- inside.owin(post_realiz$x, post_realiz$y, w = zero_background_region)
      W_post[in_zero] <- 0
    }
    parent_x <- c(pre_hist$x, post_realiz$x)
    parent_y <- c(pre_hist$y, post_realiz$y)
    parent_t <- c(pre_hist$t, post_realiz$t)
    loglik <- hawkes_loglik_inhom_filtration_cpp(
      post_t = as.numeric(post_realiz$t),
      post_x = as.numeric(post_realiz$x),
      post_y = as.numeric(post_realiz$y),
      W_val = as.numeric(W_post),
      parent_t = as.numeric(parent_t),
      parent_x = as.numeric(parent_x),
      parent_y = as.numeric(parent_y),
      mu = p[["mu"]],
      alpha = p[["alpha"]],
      beta = p[["beta"]],
      K = p[["K"]],
      areaS = active_area,
      t_start = time_window[1],
      t_end = time_window[2],
      adjust_factor = 1.0,
      t_trunc = if (!is.null(t_trunc)) t_trunc else -1.0
    )
    if (!is.finite(loglik)) return(-Inf)
    loglik
  }

  fits <- list()
  labelling_proposals <- list()
  total_sampling <- 0
  total_likelihood <- 0
  total_param_update <- 0
  current_change_factor <- change_factor
  no_flip_streak <- 0L
  progress_every <- suppressWarnings(as.integer(Sys.getenv("OK_SEM_PROGRESS_EVERY", "25")))
  if (!is.finite(progress_every) || is.na(progress_every) || progress_every < 1L) progress_every <- 0L

  for (i in 1:iter) {
    t_iter <- proc.time()[3]
    t_samp <- proc.time()[3]
    if (verbose && sem_timing_verbose) {
      cat(sprintf("  [Iter %d/%d] start: n_post=%d change_factor=%.4f\n",
                  i, iter, nrow(post_data), current_change_factor))
    }
    trigger_explore <- (no_flip_streak > 0L) && ((no_flip_streak %% stagnation_trigger_every) == 0L)
    proposal_change_factor <- if (trigger_explore) {
      min(base_change_factor, current_change_factor * 2)
    } else {
      current_change_factor
    }
    if (proposal_method == "single_flip") {
      n_post_pts <- nrow(post_data)
      labelling_proposals <- lapply(seq_len(n_post_pts), function(j) {
        proposed <- post_data
        proposed$inferred_process[j] <- if (proposed$inferred_process[j] == "control") {
          "treated"
        } else {
          "control"
        }
        proposed
      })
      if (verbose && i == 1) {
        cat(sprintf("  single_flip: %d proposals (one per post-treatment point)\n", length(labelling_proposals)))
      }
    } else if (!is.null(proposal_update_cadence)) {
      if ((i %% proposal_update_cadence) == 0 | i == iter | i == 1) {
        if (verbose) print("Updating labelling proposals")
        post_inds <- as.numeric(tileindex(post_data$x, post_data$y, partition))
        pre_for_proposals <- if (is_biv_etas || hawkes_use_filtration_history) {
          pre_data
        } else {
          pre_data[0, , drop = FALSE]
        }
        filt_by_proc <- if (!is.null(pre_for_proposals$location_process)) {
          split(pre_for_proposals, pre_for_proposals$location_process)
        } else {
          NULL
        }
        post_proposals <- lapply(1:n_props, function(j) {
          t_prop <- proc.time()[3]
          prop_out <- simulation_labeling_hawkes_hawkes_fast(
            post_data, partition = partition, partition_process = partition_processes,
            statespace = statespace, state_spaces = state_spaces,
            windowT = time_window,
            hawkes_params_control = control_par[[length(control_par)]],
            hawkes_params_treated = treated_par[[length(treated_par)]],
            change_factor = proposal_change_factor, filtration = pre_for_proposals,
            verbose = isTRUE(verbose) && isTRUE(sem_proposal_verbose),
            temporal_weight = temporal_weight,
            temporal_scale_days = temporal_scale_days,
            points_tile_index = post_inds, filt_by_proc = filt_by_proc,
            model_type = model_type, ...
          )
          if (verbose && sem_timing_verbose) {
            flips_j <- sum(post_data$inferred_process != prop_out$inferred_process, na.rm = TRUE)
            cat(sprintf("    [proposal %d/%d] done in %.2fs flips=%d\n",
                        j, n_props, proc.time()[3] - t_prop, flips_j))
          }
          prop_out
        })
        # Proposals preserve post_data ordering; avoid redundant per-proposal sort.
        labelling_proposals <- lapply(post_proposals, as.data.frame)
      }
    }
    force_include_starting <- isTRUE(include_starting_data) && (i <= include_starting_first_n)
    include_starting_this_iter <- (isTRUE(include_starting_data) && !isTRUE(trigger_explore)) || force_include_starting
    if (include_starting_this_iter) {
      if (length(labelling_proposals) == 0) {
        labelling_proposals <- list()
      }
      labelling_proposals[[length(labelling_proposals) + 1]] <- post_data
    }

    if (length(labelling_proposals) == 0) {
      stop("No labelling proposals generated in iteration ", i)
    }

    current_post_ip <- post_data$inferred_process
    flips_per_proposal <- vapply(labelling_proposals, function(y) {
      sum(current_post_ip != y$inferred_process)
    }, numeric(1))
    proposal_ctrl_idx <- lapply(labelling_proposals, function(y) which(y$inferred_process == "control"))
    proposal_treat_idx <- lapply(labelling_proposals, function(y) which(y$inferred_process == "treated"))

    class_results <- c(class_results, lapply(labelling_proposals, function(d) {
      class_func(d)
    }))
    t_samp <- proc.time()[3] - t_samp
    total_sampling <- total_sampling + t_samp
    if (verbose && sem_timing_verbose) {
      cat(sprintf("  [Iter %d/%d] proposals ready: n=%d sampling=%.2fs\n",
                  i, iter, length(labelling_proposals), t_samp))
      cat(sprintf("  [Iter %d/%d] likelihood evaluation starting\n", i, iter))
    }

    t_lik <- proc.time()[3]
    if (metric_name == "post_likelihood") {
      ref_post <- post_data
      ctrl_params_vec <- unlist(hawkes_params_control)
      treat_params_vec <- unlist(treated_par[[length(treated_par)]])
      metric <- rep(NA_real_, length(labelling_proposals))
      unchanged_idx <- which(flips_per_proposal == 0)
      changed_idx <- which(flips_per_proposal != 0)

      if (is_biv_etas) {
        biv_par <- if (!is.null(biv_etas_params)) {
          biv_etas_params
        } else {
          init_bivariate_from_independent(ctrl_params_vec, treat_params_vec)
        }
        eval_biv <- function(realiz) {
          # Include pre-treatment control history so post-treatment scoring
          # reflects carryover/triggering from pre events.
          realiz_full <- rbind(pre_data, realiz)
          realiz_full <- realiz_full[order(realiz_full$t), , drop = FALSE]
          loglik_etas_bivariate(
            params = biv_par, realiz = realiz_full,
            windowT = history_window, windowS = statespace,
            control_state_space = control_state_space,
            treated_state_space = treated_state_space,
            background_rate_var = background_rate_var,
            treated_background_zero_before = treated_background_zero_before,
            t_trunc = t_trunc
          )
        }
        if (length(unchanged_idx) > 0) {
          if (!is.finite(current_metric_cache)) current_metric_cache <- eval_biv(ref_post)
          metric[unchanged_idx] <- current_metric_cache
        }
        if (length(changed_idx) > 0) {
          metric[changed_idx] <- vapply(changed_idx, function(j) {
            y <- labelling_proposals[[j]]
            eval_biv(y)
          }, numeric(1))
        }
      } else {
      pc_ctrl_all <- precompute_loglik_args(ref_post, statespace, treated_state_space)
      pc_treat_all <- precompute_loglik_args(ref_post, statespace, control_state_space)
      printed_metric_diag <- FALSE
      pre_ctrl_hist <- if (hawkes_use_filtration_history && !is_etas) select_pre_history_by_label("control") else NULL
      pre_treat_hist <- if (hawkes_use_filtration_history && !is_etas) select_pre_history_by_label("treated") else NULL
      eval_nonbiv <- function(realiz, ctrl_idx, treat_idx) {
        if (length(ctrl_idx) < 1L) return(-Inf)
        if (length(treat_idx) < 1L) return(-Inf)
        if (hawkes_use_filtration_history && !is_etas) {
          control_lik <- hawkes_conditional_loglik(
            params_vec = ctrl_params_vec,
            post_realiz = realiz[ctrl_idx, , drop = FALSE],
            zero_background_region = treated_state_space,
            pre_hist = pre_ctrl_hist
          )
          treat_lik <- hawkes_conditional_loglik(
            params_vec = treat_params_vec,
            post_realiz = realiz[treat_idx, , drop = FALSE],
            zero_background_region = control_state_space,
            pre_hist = pre_treat_hist
          )
        } else {
          control_lik <- loglik_fn(
            params = ctrl_params_vec, realiz = realiz[ctrl_idx, , drop = FALSE],
            windowT = time_window, windowS = statespace,
            precomp = list(active_area = pc_ctrl_all$active_area,
                           in_zero_bg = pc_ctrl_all$in_zero_bg_all[ctrl_idx]), ...
          )
          treat_lik <- loglik_fn(
            params = treat_params_vec, realiz = realiz[treat_idx, , drop = FALSE],
            windowT = time_window, windowS = statespace,
            precomp = list(active_area = pc_treat_all$active_area,
                           in_zero_bg = pc_treat_all$in_zero_bg_all[treat_idx]), ...
          )
        }
        if (verbose && !printed_metric_diag) {
          cat(sprintf("  [metric diag] proposal 1: n_ctrl=%d n_treat=%d ctrl_lik=%s treat_lik=%s\n",
                      length(ctrl_idx), length(treat_idx), signif(control_lik, 6), signif(treat_lik, 6)))
          printed_metric_diag <<- TRUE
        }
        control_lik + treat_lik
      }
      if (length(unchanged_idx) > 0) {
        if (!is.finite(current_metric_cache)) {
          ref_ctrl_idx <- which(ref_post$inferred_process == "control")
          ref_treat_idx <- which(ref_post$inferred_process == "treated")
          current_metric_cache <- eval_nonbiv(ref_post, ref_ctrl_idx, ref_treat_idx)
        }
        metric[unchanged_idx] <- current_metric_cache
      }
      if (length(changed_idx) > 0) {
        metric[changed_idx] <- vapply(changed_idx, function(j) {
          y <- labelling_proposals[[j]]
          eval_nonbiv(y, proposal_ctrl_idx[[j]], proposal_treat_idx[[j]])
        }, numeric(1))
      }
      } # end non-bivariate metric
    }
    if (metric_name == "post_likelihood_control") {
      metric <- vapply(labelling_proposals, function(y) {
        post_ctrl <- y[y$inferred_process == "control", ]
        if (hawkes_use_filtration_history && !is_etas && !is_biv_etas) {
          hawkes_conditional_loglik(
            params_vec = unlist(hawkes_params_control),
            post_realiz = post_ctrl,
            zero_background_region = treated_state_space,
            pre_hist = select_pre_history_by_label("control")
          )
        } else {
          loglik_hawk_fast(
            params = unlist(hawkes_params_control), realiz = post_ctrl,
            windowT = time_window, windowS = statespace,
            zero_background_region = treated_state_space, ...
          )
        }
      }, numeric(1))
    }
    if (metric_name == "cheating") {
      metric <- vapply(labelling_proposals, function(y) {
        y_post <- y[is_post, ]
        mean(y_post$inferred_process == y_post$process)
      }, numeric(1))
    }

    sample_metric_idx <- function(metric_values, temp) {
      vals <- as.numeric(metric_values)
      finite_idx <- which(is.finite(vals))
      if (length(finite_idx) < 1L) return(which.max(vals))
      shifted <- (vals[finite_idx] - max(vals[finite_idx])) / temp
      shifted <- pmax(shifted, -700) # avoid underflow warnings in exp()
      w <- exp(shifted)
      if (!all(is.finite(w)) || sum(w) <= 0) return(finite_idx[[which.max(vals[finite_idx])]])
      finite_idx[[sample.int(length(finite_idx), size = 1L, prob = w)]]
    }
    if (identical(optim_method, "sample_weighted") && metric_name == "post_likelihood") {
      best_metric_idx <- sample_metric_idx(metric, selection_temperature)
      max_metric_labelling <- labelling_proposals[[best_metric_idx]]
    } else if (MCMC_style & metric_name == "post_likelihood" & i != 1) {
      tmp <- exp(metric[1] - metric[2])
      tmp <- runif(1) < min(tmp, 1)
      if (tmp) {
        best_metric_idx <- which.max(metric)
        max_metric_labelling <- labelling_proposals[[best_metric_idx]]
      } else {
        best_metric_idx <- length(labelling_proposals)
        max_metric_labelling <- labelling_proposals[[best_metric_idx]]
      }
    } else {
      best_metric_idx <- which.max(metric)
      max_metric_labelling <- labelling_proposals[[best_metric_idx]]
    }
    t_lik <- proc.time()[3] - t_lik
    total_likelihood <- total_likelihood + t_lik
    if (verbose && sem_timing_verbose) {
      cat(sprintf("  [Iter %d/%d] likelihood evaluation done: %.2fs\n", i, iter, t_lik))
    }

    proposed_best <- as.data.frame(max_metric_labelling)
    accepted_labelling <- proposed_best
    proposed_best_flips <- sum(post_data$inferred_process != proposed_best$inferred_process, na.rm = TRUE)
    if (proposed_best_flips > max_flips_per_step) {
      flip_idx <- which(post_data$inferred_process != proposed_best$inferred_process)
      keep_idx <- sample(flip_idx, size = max_flips_per_step, replace = FALSE)
      accepted_labelling <- post_data
      accepted_labelling$inferred_process[keep_idx] <- proposed_best$inferred_process[keep_idx]
    }
    max_diff <- sum(post_data$inferred_process != accepted_labelling$inferred_process, na.rm = TRUE)
    average <- mean(flips_per_proposal)
    accepted_flips_cum_next <- accepted_flips_cum + max_diff
    force_param_update_due <- accepted_flips_cum_next >= next_forced_param_update_at
    do_param_update <- FALSE
    if (!is.null(param_update_cadence)) {
      do_param_update <- ((i %% param_update_cadence) == 0L) || i == iter || force_param_update_due
    }
    if (do_param_update) {
        t_param_start <- proc.time()[3]
        if (verbose) print("Updating Parameters")
        if (verbose && sem_timing_verbose) {
          cat(sprintf("  [Iter %d/%d] parameter update starting\n", i, iter))
        }

        mml_post <- accepted_labelling
        mml_post_treated <- mml_post[mml_post$inferred_process == "treated", ]
        mml_post_control <- mml_post[mml_post$inferred_process == "control", ]

        if (is_biv_etas) {
          # Joint bivariate ETAS parameter update
          biv_par <- if (!is.null(biv_etas_params)) {
            unlist(biv_etas_params)
          } else {
            init_bivariate_from_independent(
              control_par[[length(control_par)]],
              treated_par[[length(treated_par)]])
          }
          if (is.null(names(biv_par))) names(biv_par) <- .etas_bivariate_par_names
          # Parameter updates should also account for pre-treatment control history.
          mml_full <- rbind(pre_data, mml_post)
          mml_full <- mml_full[order(mml_full$t), , drop = FALSE]
          process_id_full <- if ("inferred_process" %in% names(mml_full)) {
            as.integer(mml_full$inferred_process == "treated")
          } else if ("process" %in% names(mml_full)) {
            as.integer(mml_full$process == "treated")
          } else if ("location_process" %in% names(mml_full)) {
            as.integer(mml_full$location_process == "treated")
          } else {
            rep(0L, nrow(mml_full))
          }
          areaS_0 <- spatstat.geom::area(control_state_space)
          areaS_1 <- spatstat.geom::area(treated_state_space)
          if (!is.finite(areaS_0) || areaS_0 <= 0) areaS_0 <- 1
          if (!is.finite(areaS_1) || areaS_1 <= 0) areaS_1 <- 1
          W_0 <- rep(1.0, nrow(mml_full))
          W_1 <- rep(1.0, nrow(mml_full))
          W_0[inside.owin(mml_full$x, mml_full$y, treated_state_space)] <- 0
          W_1[inside.owin(mml_full$x, mml_full$y, control_state_space)] <- 0
          if (!is.null(background_rate_var) && background_rate_var %in% names(mml_full)) {
            W_cov <- as.numeric(mml_full[[background_rate_var]])
            W_cov[!is.finite(W_cov)] <- 0
            min_pos <- suppressWarnings(min(W_cov[W_cov > 0], na.rm = TRUE))
            if (!is.finite(min_pos)) min_pos <- 1e-12
            W_cov[W_cov <= 0] <- min_pos
            W_0 <- W_0 * W_cov
            W_1 <- W_1 * W_cov
          }
          if (!is.null(treated_background_zero_before)) {
            W_1[mml_full$t < as.numeric(treated_background_zero_before)] <- 0
          }
          biv_precomp <- list(
            W_0 = W_0, W_1 = W_1,
            areaS_0 = areaS_0, areaS_1 = areaS_1,
            process_id = process_id_full
          )
          m0_full <- if ("mag" %in% names(mml_full)) min(mml_full$mag, na.rm = TRUE) else NULL
          biv_obj <- function(par15) {
            loglik_etas_bivariate(
              params = par15, realiz = mml_full,
              windowT = history_window, windowS = statespace,
              control_state_space = control_state_space,
              treated_state_space = treated_state_space,
              background_rate_var = background_rate_var,
              treated_background_zero_before = treated_background_zero_before,
              m0 = m0_full,
              t_trunc = t_trunc,
              precomp = biv_precomp
            )
          }
          if (length(fixed_idx) > 0) {
            biv_free <- biv_par[free_idx]
            biv_wrap <- function(fp) {
              p15 <- biv_par; p15[free_idx] <- fp; biv_obj(p15)
            }
            biv_res <- tryCatch(
              optim(par = biv_free, fn = biv_wrap, method = "Nelder-Mead",
                    control = list(fnscale = -1, trace = 0, maxit = sem_outer_maxit_biv)),
              error = function(e) { cat("  bivariate optim error:", e$message, "\n"); list(par = biv_free) }
            )
            biv_par[free_idx] <- biv_res$par
          } else {
            biv_res <- tryCatch(
              optim(par = biv_par, fn = biv_obj, method = "Nelder-Mead",
                    control = list(fnscale = -1, trace = 0, maxit = sem_outer_maxit_biv)),
              error = function(e) { cat("  bivariate optim error:", e$message, "\n"); list(par = biv_par) }
            )
            biv_par <- biv_res$par
          }
          names(biv_par) <- .etas_bivariate_par_names
          biv_etas_params <<- biv_par
          fits[[i]] <- biv_res
          # Extract marginal params
          treated_par[[length(treated_par) + 1]] <- as.list(c(
            mu = biv_par[["mu_1"]], A = biv_par[["A_11"]],
            alpha_m = biv_par[["alpha_m_11"]],
            c = biv_par[["c"]], p = biv_par[["p"]],
            D = biv_par[["D"]], gamma = biv_par[["gamma"]], q = biv_par[["q"]]))
          control_par[[length(control_par) + 1]] <- as.list(c(
            mu = biv_par[["mu_0"]], A = biv_par[["A_00"]],
            alpha_m = biv_par[["alpha_m_00"]],
            c = biv_par[["c"]], p = biv_par[["p"]],
            D = biv_par[["D"]], gamma = biv_par[["gamma"]], q = biv_par[["q"]]))
        } else {
        profile_optim <- function(full_par, obj_fn, label) {
          full_vec <- unlist(full_par)
          if (length(fixed_idx) > 0) {
            free_par <- full_vec[free_idx]
            wrap_fn <- function(fp, ...) {
              p4 <- full_vec; p4[free_idx] <- fp
              obj_fn(p4, ...)
            }
            res <- tryCatch(
              optim(par = free_par, fn = wrap_fn, method = "Nelder-Mead",
                    control = list(fnscale = -1, trace = 0, maxit = sem_outer_maxit), ...),
              error = function(e) { cat("  error fitting", label, ":", e$message, "\n"); list(par = free_par) }
            )
            out <- full_vec; out[free_idx] <- res$par
            res$par <- out
          } else {
            res <- tryCatch(
              optim(par = full_vec, fn = obj_fn, method = "Nelder-Mead",
                    control = list(fnscale = -1, trace = 0, maxit = sem_outer_maxit), ...),
              error = function(e) { cat("  error fitting", label, ":", e$message, "\n"); list(par = full_vec) }
            )
          }
          pv <- as.numeric(res$par)
          par_list <- as.list(pv)
          names(par_list) <- all_names
          list(fit = res, par_list = par_list)
        }

        if (update_control_params) {
          optim_func_treat <- function(params, ...) {
            if (hawkes_use_filtration_history && !is_etas) {
              hawkes_conditional_loglik(
                params_vec = params,
                post_realiz = mml_post_treated,
                zero_background_region = control_state_space,
                pre_hist = select_pre_history_by_label("treated")
              )
            } else {
              loglik_fn(
                params = params, realiz = mml_post_treated,
                windowT = time_window, windowS = statespace,
                zero_background_region = control_state_space, ...
              )
            }
          }
          optim_func_control <- function(params, ...) {
            if (hawkes_use_filtration_history && !is_etas) {
              hawkes_conditional_loglik(
                params_vec = params,
                post_realiz = mml_post_control,
                zero_background_region = treated_state_space,
                pre_hist = select_pre_history_by_label("control")
              )
            } else {
              loglik_fn(
                params = params, realiz = mml_post_control,
                windowT = time_window, windowS = statespace,
                zero_background_region = treated_state_space, ...
              )
            }
          }
          res_t <- profile_optim(treated_par[[length(treated_par)]], optim_func_treat, "treated")
          fits[[i]] <- res_t$fit
          treated_par[[length(treated_par) + 1]] <- res_t$par_list

          res_c <- profile_optim(control_par[[length(control_par)]], optim_func_control, "control")
          control_par[[length(control_par) + 1]] <- res_c$par_list
        } else {
          optim_func <- function(params, ...) {
            if (hawkes_use_filtration_history && !is_etas) {
              hawkes_conditional_loglik(
                params_vec = params,
                post_realiz = mml_post_treated,
                zero_background_region = control_state_space,
                pre_hist = select_pre_history_by_label("treated")
              )
            } else {
              loglik_fn(
                params = params, realiz = mml_post_treated,
                windowT = time_window, windowS = statespace,
                zero_background_region = control_state_space, ...
              )
            }
          }
          res_t <- profile_optim(treated_par[[length(treated_par)]], optim_func, "treated")
          fits[[i]] <- res_t$fit
          treated_par[[length(treated_par) + 1]] <- res_t$par_list
        }
        } # end non-bivariate param update
        total_param_update <- total_param_update + (proc.time()[3] - t_param_start)
        if (verbose) {
          print("Estimated model params")
          print("treated:"); print(unlist(treated_par[[length(treated_par)]]))
          print("control:"); print(unlist(control_par[[length(control_par)]]))
          print(paste0("Updating model params on iteration ", i, " took ", signif(proc.time()[3] - t_param_start, 2)))
          if (force_param_update_due) {
            print(paste0("Forced parameter update: cumulative accepted flips reached ",
                         accepted_flips_cum_next, " (next threshold=", next_forced_param_update_at, ")."))
          }
        }
        current_metric_cache <- NA_real_
        accepted_flips_cum <- accepted_flips_cum_next
        while (accepted_flips_cum >= next_forced_param_update_at) {
          next_forced_param_update_at <- next_forced_param_update_at + force_param_update_flip_n
        }
    } else {
      accepted_flips_cum <- accepted_flips_cum_next
    }

    retained_starting <- isTRUE(max_diff == 0)
    if (retained_starting) {
      no_flip_streak <- no_flip_streak + 1L
    } else {
      no_flip_streak <- 0L
    }
    retained_starting_trace[i] <- retained_starting
    if (retained_starting) {
      current_change_factor <- max(change_factor_min, current_change_factor / 2)
    } else {
      current_change_factor <- min(change_factor_max, current_change_factor * 2)
    }
    change_factor_trace[i] <- current_change_factor

    if (update_starting_data) {
      post_data <- as.data.frame(accepted_labelling)
      post_data <- post_data[order(post_data$t), , drop = FALSE]
    }

    mml_post_acc <- accepted_labelling
    accuracy <- mean(mml_post_acc$inferred_process == mml_post_acc$process)
    accuracies[i] <- accuracy
    metric_vec[i] <- metric[best_metric_idx]
    current_metric_cache <- metric_vec[i]
    average_flips[i] <- average
    max_metric_flips[i] <- max_diff

    all_accuracies[[i]] <- vapply(labelling_proposals, function(y) {
      mean(y$inferred_process == y$process)
    }, numeric(1))
    all_metrics[[i]] <- metric
    if (verbose) {
      n_post <- nrow(post_data)
      cat(sprintf("  [Iter %d] Proposals: %d | Proposed flips: min=%d avg=%.1f max=%d (of %d post pts) | Accepted (best): %d\n",
                  i, length(labelling_proposals),
                  min(flips_per_proposal), average, max(flips_per_proposal),
                  n_post, max_diff))
      cat("  Metric summary:\n")
      print(summary(metric))
      print(mml_post_acc %>% dplyr::group_by(.data$location_process, .data$process, .data$inferred_process) %>% dplyr::summarize(n = dplyr::n()))
      print(paste0("Accuracy: ", accuracy))
      print(paste0("Iteration ", i, " took ", signif(proc.time()[3] - t_iter, 2)))
      print(paste0("sampling took: ", signif(t_samp, 2)))
      print(paste0("Likelihood eval took: ", signif(t_lik, 2)))
      print(paste0("change_factor(next) = ", signif(current_change_factor, 4),
                   " [retained_starting=", retained_starting, "]"))
      if (trigger_explore) {
        print(paste0("stagnation trigger active (streak=", no_flip_streak,
                     "): excluded starting-data proposal this iter, proposal change_factor=",
                     signif(proposal_change_factor, 4)))
      }
    } else if (progress_every > 0L && (i == 1L || i == iter || (i %% progress_every) == 0L)) {
      # Compact heartbeat for long runs when full SEM tracing is disabled.
      cat(sprintf("  [SEM progress] iter=%d/%d accepted_flips=%d avg_flips=%.1f metric=%.4f change_factor=%.4f\n",
                  i, iter, max_diff, average, metric_vec[i], current_change_factor))
    }
  }
  final_labelling <- rbind(pre_data, max_metric_labelling)
  final_labelling <- final_labelling[order(final_labelling$t), , drop = FALSE]
  return(list(
    labelling = final_labelling,
    treated_par = treated_par, control_par = control_par,
    accuracies = accuracies, average_flips = average_flips,
    max_metric_flips = max_metric_flips, metrics = metric_vec,
    all_accuracies = all_accuracies, all_metrics = all_metrics,
    change_factor_trace = change_factor_trace,
    retained_starting_trace = retained_starting_trace,
    class_results = class_results, fits = fits,
    timing = list(
      n_iter = iter,
      sampling_s = total_sampling,
      likelihood_s = total_likelihood,
      param_update_s = total_param_update
    )
  ))
}

#' Plot flips over iterations for tuning
#' @param res Result from em_style_labelling or adaptive_SEM (adaptive list)
#' @return ggplot object
#' @export
plot_flips <- function(res) {
  if ("adaptive" %in% names(res)) res <- res$adaptive
  df <- data.frame(
    iteration = 1:length(res$average_flips),
    average_flips = res$average_flips,
    max_metric_flips = res$max_metric_flips
  )
  ggplot(df, aes(x = iteration)) +
    geom_line(aes(y = average_flips, color = "Average flips")) +
    geom_line(aes(y = max_metric_flips, color = "Max metric flips")) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(title = "Flips per iteration", x = "Iteration", y = "Number of flips") +
    theme_minimal() +
    scale_color_manual(name = "Metric", values = c("Average flips" = "blue", "Max metric flips" = "red"))
}

