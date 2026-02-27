#' @importFrom dplyr filter arrange pull group_by summarize bind_rows
#' @importFrom caret confusionMatrix
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
    probs <- probs - min(0, min(probs))
    if (sum(probs) == 0) probs <- rep(1 / length(probs), length(probs))
    probs <- probs / sum(probs)
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
    probs <- probs - min(0, min(probs))
    if (sum(probs) == 0) probs <- rep(1 / length(probs), length(probs))
    probs <- probs / sum(probs)
    probs[control_tiles] <- 0
    if (sum(probs) <= 0) {
      thin_treated <- 0
      partition_thins <- c()
    } else {
      probs <- probs / sum(probs)
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
                                                   points_tile_index = NULL,
                                                   ...) {
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

  no_points_hawkes <- list(mu = 0, alpha = 0, beta = 0, K = 0)

  if (is.null(hawkes_params_control)) {
    thin_control <- 0
    sim_inds <- numeric(0)
  } else {
    sim_data <- generate_inhomogeneous_hawkes(
      Omega = statespace, partition = partition, time_window = windowT,
      partition_processes = partition_process,
      hawkes_params = list(control = hawkes_params_control, treated = no_points_hawkes),
      filtration = filtration, state_spaces = state_spaces,
      space_triggering = FALSE, ...
    )
    sim_data <- sim_data[sim_data$t > windowT[1], ]
    sim_inds <- if (!is.null(sim_data$tile_index)) as.numeric(sim_data$tile_index) else as.numeric(tileindex(sim_data$x, sim_data$y, partition))
    where_to_thin <- tabulate(control_inds, nbins = partition$n) - tabulate(sim_inds, nbins = partition$n)
    thin_control <- length(control_inds) - length(sim_inds)
  }

  thin_control <- thin_control * change_factor

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
    probs <- probs - min(0, min(probs))
    if (sum(probs) == 0) probs <- rep(1 / length(probs), length(probs))
    probs <- probs / sum(probs)
    partition_thins <- rmultinom(1, n_total_change, prob = probs)

    attractor_idx <- which(dat$inferred_process == target_label)
    has_attractors <- length(attractor_idx) > 0

    if (has_attractors) {
      win_box <- owin(range(dat$x), range(dat$y))
      X_attract <- ppp(dat$x[attractor_idx], dat$y[attractor_idx], window = win_box, check = FALSE)
    }

    for (i in 1:length(partition_thins)) {
      n_thin <- partition_thins[i]
      if (length(n_thin) == 0 || is.na(n_thin) || n_thin <= 0) next

      in_tile_idx <- which(inds == i)
      candidates <- in_tile_idx[dat$inferred_process[in_tile_idx] == relabel_label]
      if (length(in_tile_idx) == 1L && n_thin > 0) {
        changes <- in_tile_idx
        dat$inferred_process[changes] <- if (dat$inferred_process[changes] == "control") "treated" else "control"
        next
      }
      if (length(candidates) == 0) next

      sampling_weights <- NULL
      if (has_attractors && length(candidates) > 0 && proximity_weight > 0) {
        X_cand <- ppp(dat$x[candidates], dat$y[candidates], window = win_box, check = FALSE)
        dists <- nncross(X_cand, X_attract)$dist
        geom_weights <- 1 / (dists + 1e-4)
        if (sum(geom_weights) > 0) {
          geom_weights <- geom_weights / sum(geom_weights)
        } else {
          geom_weights <- rep(1 / length(candidates), length(candidates))
        }
        uni_weights <- rep(1 / length(candidates), length(candidates))
        sampling_weights <- (proximity_weight * geom_weights) + ((1 - proximity_weight) * uni_weights)
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

  if (is.null(hawkes_params_treated)) return(dat)

  treated_inds <- inds[dat$inferred_process == "treated"]
  control_tiles <- which(partition_process == "control")
  sim_data_t <- generate_inhomogeneous_hawkes(
    Omega = statespace, partition = partition, time_window = windowT,
    partition_processes = partition_process,
    hawkes_params = list(control = no_points_hawkes, treated = hawkes_params_treated),
    filtration = filtration, state_spaces = state_spaces,
    space_triggering = FALSE, ...
  )
  sim_data_t <- sim_data_t[sim_data_t$t > windowT[1], ]
  sim_inds_t <- if (!is.null(sim_data_t$tile_index)) as.numeric(sim_data_t$tile_index) else as.numeric(tileindex(sim_data_t$x, sim_data_t$y, partition))
  where_to_thin_t <- tabulate(treated_inds, nbins = partition$n) - tabulate(sim_inds_t, nbins = partition$n)
  thin_treated <- (length(treated_inds) - length(sim_inds_t)) * change_factor
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
    probs_t <- probs_t - min(0, min(probs_t))
    if (sum(probs_t) == 0) probs_t <- rep(1 / length(probs_t), length(probs_t))
    probs_t <- probs_t / sum(probs_t)
    probs_t[control_tiles] <- 0
    if (sum(probs_t) > 0) {
      probs_t <- probs_t / sum(probs_t)
      partition_thins_t <- rmultinom(1, n_total_t, prob = probs_t)
      for (i in 1:length(partition_thins_t)) {
        n_thin <- partition_thins_t[i]
        if (length(n_thin) == 0 || is.na(n_thin) || n_thin <= 0) next
        in_tile_idx <- which(inds == i)
        candidates_t <- in_tile_idx[dat$inferred_process[in_tile_idx] == relabel_label_t]
        if (length(in_tile_idx) == 1L && n_thin > 0) {
          changes <- in_tile_idx
          dat$inferred_process[changes] <- if (dat$location_process[changes] == "control") "treated" else "control"
          next
        }
        if (length(candidates_t) == 0) next
        changes_local <- sample(length(candidates_t), size = min(n_thin, length(candidates_t)), replace = FALSE)
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
#' @param optim_method One of "max", "mean", "truncated_mean"
#' @param state_spaces Optional precomputed state spaces
#' @param metric_name Metric for selecting best labeling
#' @param iter Number of iterations
#' @param n_props Number of proposals per iteration
#' @param change_factor Scaling factor for relabeling
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
                               optim_method = "max",
                               state_spaces = NULL,
                               metric_name = "post_likelihood",
                               iter = 100,
                               n_props = 100,
                               change_factor = 0.1,
                               MCMC_style = FALSE,
                               proposal_method = "simulation",
                               fixed_params = NULL,
                               verbose = FALSE,
                               ...) {
  dots <- list(...)
  background_rate_var <- if ("background_rate_var" %in% names(dots)) dots$background_rate_var else NULL

  all_names <- c("mu", "alpha", "beta", "K")
  fixed_idx <- if (!is.null(fixed_params)) match(names(fixed_params), all_names) else integer(0)
  free_idx  <- setdiff(seq_along(all_names), fixed_idx)

  class_func <- function(labelling) {
    truth_factor <- factor(labelling$process, levels = c("control", "treated"))
    inferred_factor <- factor(labelling$inferred_process, levels = c("control", "treated"))
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
  orig_starting_data <- starting_data
  treated_par <- list(hawkes_params_treated)
  control_par <- list(hawkes_params_control)
  accuracies <- numeric(iter)
  average_flips <- numeric(iter)
  max_metric_flips <- numeric(iter)
  metric_vec <- numeric(iter)
  all_accuracies <- vector("list", iter)
  all_metrics <- vector("list", iter)

  is_pre <- starting_data$t < treatment_time
  is_post <- !is_pre

  fits <- list()
  labelling_proposals <- list()

  for (i in 1:iter) {
    t_iter <- proc.time()[3]
    pre <- starting_data[is_pre, ]
    post <- starting_data[is_post, ]
    post <- post[order(post$t), ]
    pre$location_process <- "control"
    pre$inferred_process <- NULL
    t_samp <- proc.time()[3]
    if (proposal_method == "single_flip") {
      pre$inferred_process <- "control"
      labelling_proposals <- single_flip_proposals(post, pre)
      if (verbose && i == 1) {
        cat(sprintf("  single_flip: %d proposals (one per post-treatment point)\n", length(labelling_proposals)))
      }
    } else if (!is.null(proposal_update_cadence)) {
      if ((i %% proposal_update_cadence) == 0 | i == iter | i == 1) {
        if (verbose) print("Updating labelling proposals")
        post_inds <- as.numeric(tileindex(post$x, post$y, partition))
        filt_by_proc <- if (!is.null(pre$location_process)) split(pre, pre$location_process) else NULL
        post_proposals <- lapply(1:n_props, function(j) {
          simulation_labeling_hawkes_hawkes_fast(
            post, partition = partition, partition_process = partition_processes,
            statespace = statespace, state_spaces = state_spaces,
            windowT = time_window,
            hawkes_params_control = control_par[[length(control_par)]],
            hawkes_params_treated = treated_par[[length(treated_par)]],
            change_factor = change_factor, filtration = pre, verbose = FALSE,
            points_tile_index = post_inds, filt_by_proc = filt_by_proc, ...
          )
        })
        pre$inferred_process <- "control"
        labelling_proposals <- lapply(post_proposals, function(tmp) rbind(pre, tmp))
      }
    }
    if (i != 1 && isTRUE(include_starting_data)) {
      if (length(labelling_proposals) == 0) {
        labelling_proposals <- list()
      }
      pre$inferred_process <- "control"
      labelling_proposals[[length(labelling_proposals) + 1]] <- rbind(pre, post)
    }

    if (length(labelling_proposals) == 0) {
      stop("No labelling proposals generated in iteration ", i)
    }

    class_results <- c(class_results, lapply(labelling_proposals, function(d) {
      d_post <- d[d$t > treatment_time, ]
      class_func(d_post)
    }))
    t_samp <- proc.time()[3] - t_samp

    t_lik <- proc.time()[3]
    if (metric_name == "post_likelihood") {
      ref_post <- post
      pc_ctrl_all <- precompute_loglik_args(ref_post, statespace, treated_state_space)
      pc_treat_all <- precompute_loglik_args(ref_post, statespace, control_state_space)
      ctrl_params_vec <- unlist(hawkes_params_control)
      treat_params_vec <- unlist(treated_par[[length(treated_par)]])

      metric <- vapply(seq_along(labelling_proposals), function(j) {
        y <- labelling_proposals[[j]]
        realiz <- y[is_post, ]

        ctrl_rows <- realiz$inferred_process == "control"
        if (!any(ctrl_rows)) return(-Inf)
        control_lik <- loglik_hawk_fast(
          params = ctrl_params_vec, realiz = realiz[ctrl_rows, ],
          windowT = time_window, windowS = statespace,
          precomp = list(active_area = pc_ctrl_all$active_area,
                         in_zero_bg = pc_ctrl_all$in_zero_bg_all[ctrl_rows]), ...
        )

        treat_rows <- !ctrl_rows
        if (!any(treat_rows)) return(-Inf)
        treat_lik <- loglik_hawk_fast(
          params = treat_params_vec, realiz = realiz[treat_rows, ],
          windowT = time_window, windowS = statespace,
          precomp = list(active_area = pc_treat_all$active_area,
                         in_zero_bg = pc_treat_all$in_zero_bg_all[treat_rows]), ...
        )
        if (verbose && j == 1 && i == 1) {
          cat(sprintf("  [metric diag] proposal 1: n_ctrl=%d n_treat=%d ctrl_lik=%s treat_lik=%s\n",
                      sum(ctrl_rows), sum(treat_rows), signif(control_lik, 6), signif(treat_lik, 6)))
          cat(sprintf("    ctrl params: %s\n", paste(names(ctrl_params_vec), signif(ctrl_params_vec, 5), sep="=", collapse="  ")))
          cat(sprintf("    treat params: %s\n", paste(names(treat_params_vec), signif(treat_params_vec, 5), sep="=", collapse="  ")))
          cat(sprintf("    pc_ctrl active_area=%.1f  pc_treat active_area=%.1f\n",
                      pc_ctrl_all$active_area, pc_treat_all$active_area))
          cat(sprintf("    in_zero_bg ctrl: %d of %d TRUE  treat: %d of %d TRUE\n",
                      sum(pc_ctrl_all$in_zero_bg_all[ctrl_rows]), sum(ctrl_rows),
                      sum(pc_treat_all$in_zero_bg_all[treat_rows]), sum(treat_rows)))
          ctrl_sub <- realiz[ctrl_rows, ]
          cat(sprintf("    ctrl W range: [%s, %s]  treat W range: [%s, %s]\n",
                      signif(min(ctrl_sub$W, na.rm=TRUE), 4), signif(max(ctrl_sub$W, na.rm=TRUE), 4),
                      signif(min(realiz[treat_rows, "W"], na.rm=TRUE), 4),
                      signif(max(realiz[treat_rows, "W"], na.rm=TRUE), 4)))
        }
        control_lik + treat_lik
      }, numeric(1))
    }
    if (metric_name == "post_likelihood_control") {
      metric <- vapply(labelling_proposals, function(y) {
        post_ctrl <- y[is_post & y$inferred_process == "control", ]
        loglik_hawk_fast(
          params = unlist(hawkes_params_control), realiz = post_ctrl,
          windowT = time_window, windowS = statespace,
          zero_background_region = treated_state_space, ...
        )
      }, numeric(1))
    }
    if (metric_name == "cheating") {
      metric <- vapply(labelling_proposals, function(y) {
        y_post <- y[is_post, ]
        mean(y_post$inferred_process == y_post$process)
      }, numeric(1))
    }

    if (MCMC_style & metric_name == "post_likelihood" & i != 1) {
      tmp <- exp(metric[1] - metric[2])
      tmp <- runif(1) < min(tmp, 1)
      if (tmp) {
        max_metric_labelling <- labelling_proposals[[which.max(metric)]]
      } else {
        max_metric_labelling <- labelling_proposals[[length(labelling_proposals)]]
      }
    } else {
      max_metric_labelling <- labelling_proposals[[which.max(metric)]]
    }
    t_lik <- proc.time()[3] - t_lik

    if (!is.null(param_update_cadence)) {
      if ((i %% param_update_cadence) == 0 | i == iter) {
        if (verbose) {
          print("Updating Hawkes Parameters")
          t_hawkes <- proc.time()[3]
        }

        mml_post <- max_metric_labelling[is_post, ]
        mml_post_treated <- mml_post[mml_post$inferred_process == "treated", ]
        mml_post_control <- mml_post[mml_post$inferred_process == "control", ]

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
                    control = list(fnscale = -1, trace = 0, maxit = 1000), ...),
              error = function(e) { cat("  error fitting", label, ":", e$message, "\n"); list(par = free_par) }
            )
            out <- full_vec; out[free_idx] <- res$par
            res$par <- out
          } else {
            res <- tryCatch(
              optim(par = full_vec, fn = obj_fn, method = "Nelder-Mead",
                    control = list(fnscale = -1, trace = 0, maxit = 1000), ...),
              error = function(e) { cat("  error fitting", label, ":", e$message, "\n"); list(par = full_vec) }
            )
          }
          pv <- as.numeric(res$par)
          list(fit = res, par_list = list(mu = pv[1], alpha = pv[2], beta = pv[3], K = pv[4]))
        }

        if (update_control_params) {
          optim_func_treat <- function(params, ...) {
            loglik_hawk_fast(
              params = params, realiz = mml_post_treated,
              windowT = time_window, windowS = statespace,
              zero_background_region = control_state_space, ...
            )
          }
          optim_func_control <- function(params, ...) {
            loglik_hawk_fast(
              params = params, realiz = mml_post_control,
              windowT = time_window, windowS = statespace,
              zero_background_region = treated_state_space, ...
            )
          }
          res_t <- profile_optim(treated_par[[length(treated_par)]], optim_func_treat, "treated")
          fits[[i]] <- res_t$fit
          treated_par[[length(treated_par) + 1]] <- res_t$par_list

          res_c <- profile_optim(control_par[[length(control_par)]], optim_func_control, "control")
          control_par[[length(control_par) + 1]] <- res_c$par_list
        } else {
          optim_func <- function(params, ...) {
            loglik_hawk_fast(
              params = params, realiz = mml_post_treated,
              windowT = time_window, windowS = statespace,
              zero_background_region = control_state_space, ...
            )
          }
          res_t <- profile_optim(treated_par[[length(treated_par)]], optim_func, "treated")
          fits[[i]] <- res_t$fit
          treated_par[[length(treated_par) + 1]] <- res_t$par_list
        }
        if (verbose) {
          print("Estimated Hawkes Params")
          print("treated:"); print(unlist(treated_par[[length(treated_par)]]))
          print("control:"); print(unlist(control_par[[length(control_par)]]))
          print(paste0("Updating hawkes params on iteration ", i, " took ", signif(proc.time()[3] - t_hawkes, 2)))
        }
      }
    }

    sd_ip <- starting_data$inferred_process[order(starting_data$t)]
    flips_per_proposal <- vapply(labelling_proposals, function(y) {
      sum(sd_ip != y$inferred_process[order(y$t)])
    }, numeric(1))
    mml_ip <- max_metric_labelling$inferred_process[order(max_metric_labelling$t)]
    max_diff <- sum(sd_ip != mml_ip)
    average <- mean(flips_per_proposal)

    if (update_starting_data) {
      starting_data <- as.data.frame(max_metric_labelling)
      is_pre <- starting_data$t < treatment_time
      is_post <- !is_pre
    }

    mml_post_acc <- max_metric_labelling[is_post, ]
    accuracy <- mean(mml_post_acc$inferred_process == mml_post_acc$process)
    accuracies[i] <- accuracy
    metric_vec[i] <- metric[which.max(metric)]
    average_flips[i] <- average
    max_metric_flips[i] <- max_diff

    all_accuracies[[i]] <- vapply(labelling_proposals, function(y) {
      y_post <- y[is_post, ]
      mean(y_post$inferred_process == y_post$process)
    }, numeric(1))
    all_metrics[[i]] <- metric
    if (verbose) {
      n_post <- sum(is_post)
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
    }
  }
  return(list(
    labelling = max_metric_labelling,
    treated_par = treated_par, control_par = control_par,
    accuracies = accuracies, average_flips = average_flips,
    max_metric_flips = max_metric_flips, metrics = metric_vec,
    all_accuracies = all_accuracies, all_metrics = all_metrics,
    class_results = class_results, fits = fits
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

