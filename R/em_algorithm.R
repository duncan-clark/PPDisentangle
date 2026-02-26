#' Adaptive Stochastic EM algorithm for point process disentangling
#'
#' Outer loop of the SEM: runs adaptive labeling proposals, computes
#' importance weights, and iteratively re-estimates Hawkes parameters.
#'
#' @param pp_data Data frame with columns x, y, t, process, location_process
#' @param partition A spatstat tess object
#' @param partition_processes Character vector of process names per tile
#' @param statespace Full observation window (owin)
#' @param time_window Numeric vector c(start, end) (used for internal naming)
#' @param treatment_time Scalar treatment time
#' @param hawkes_params_control Initial control Hawkes parameters
#' @param hawkes_params_treated Initial treated Hawkes parameters
#' @param N_labellings Number of labeling proposals per outer iteration
#' @param N_iter Number of outer EM iterations
#' @param verbose Print progress
#' @param adaptive_control List of adaptive control parameters
#' @param ... Additional arguments passed through (e.g. background_rate_var, covariate_lookup)
#' @return List with hawkes_params_control, hawkes_params_treated, t_params,
#'   labellings, adaptive results, and timing
#' @export
adaptive_SEM <- function(pp_data,
                         partition,
                         partition_processes,
                         statespace,
                         time_window,
                         treatment_time,
                         hawkes_params_control,
                         hawkes_params_treated,
                         N_labellings,
                         N_iter,
                         verbose = TRUE,
                         adaptive_control = list(
                           param_update_cadence = 20,
                           proposal_update_cadence = 1,
                           state_spaces = NULL,
                           iter = 100,
                           n_props = 10,
                           change_factor = 0.1,
                           verbose = FALSE,
                           update_starting_data = TRUE,
                           include_starting_data = FALSE
                         ),
                         ...) {
  t_global <- proc.time()[3]
  dots <- list(...)
  background_rate_var <- if ("background_rate_var" %in% names(dots)) dots$background_rate_var else NULL

  if (partition$type != "mask") {
    if (verbose) message("Converting partition to raster mask for speed...")
    statespace <- spatstat.geom::as.mask(as.owin(statespace), dimyx = c(512, 512))
    tiles_mask <- lapply(tiles(partition), function(w) spatstat.geom::as.mask(w, xy = statespace))
    partition <- tess(tiles = tiles_mask)
    if (!is.null(adaptive_control$state_spaces)) {
      adaptive_control$state_spaces <- lapply(adaptive_control$state_spaces,
        function(w) spatstat.geom::as.mask(w, xy = statespace))
      if (is.null(names(adaptive_control$state_spaces))) {
        procs <- unique(partition_processes)
        if (length(procs) == length(adaptive_control$state_spaces)) {
          names(adaptive_control$state_spaces) <- procs
        }
      }
    }
  }

  dots <- list(...)
  background_rate_var <- if ("background_rate_var" %in% names(dots)) dots$background_rate_var else NULL

  treated_idx <- (partition_processes == "treated")
  treated_state_space <- as.owin(partition[treated_idx])
  control_state_space <- as.owin(partition[!treated_idx])

  if (is.null(pp_data$inferred_process)) {
    pp_data$inferred_process <- pp_data$location_process
  }
  starting_data <- pp_data

  if (is.null(hawkes_params_treated)) {
    hawkes_params_treated <- hawkes_params_control
  }

  adaptive_step <- function(starting_data) {
    t_adapt_start <- proc.time()[3]
    result <- em_style_labelling(
      pp_data = starting_data,
      partition = partition,
      partition_processes = partition_processes,
      statespace = statespace,
      state_spaces = adaptive_control$state_spaces,
      time_window = c(treatment_time, max(starting_data$t)),
      treatment_time = treatment_time,
      hawkes_params_control = hawkes_params_control,
      hawkes_params_treated = hawkes_params_treated,
      update_control_params = adaptive_control$update_control_params,
      param_update_cadence = adaptive_control$param_update_cadence,
      proposal_update_cadence = adaptive_control$proposal_update_cadence,
      update_starting_data = adaptive_control$update_starting_data,
      include_starting_data = adaptive_control$include_starting_data,
      metric_name = "post_likelihood",
      optim_method = "max",
      iter = adaptive_control$iter,
      n_props = adaptive_control$n_props,
      change_factor = adaptive_control$change_factor,
      MCMC_style = FALSE,
      proposal_method = if (!is.null(adaptive_control$proposal_method)) adaptive_control$proposal_method else "simulation",
      fixed_params = if (!is.null(adaptive_control$fixed_params)) adaptive_control$fixed_params else NULL,
      verbose = TRUE,
      ...
    )
    t_adapt_end <- proc.time()[3]
    return(list(
      adaptive_labelling = result$labelling,
      treated_par = result$treated_par,
      control_par = result$control_par,
      accuracies = result$accuracies,
      average_flips = result$average_flips,
      max_metric_flips = result$max_metric_flips,
      metrics = result$metric_vec,
      all_accuracies = result$all_accuracies,
      all_metrics = result$all_metrics,
      class_results = result$class_results,
      fits = result$fits,
      time_taken = t_adapt_end - t_adapt_start
    ))
  }

  check_weights <- function(w) {
    w <- w / sum(w)
    w <- sort(w)
    top_5_perc <- sum(w[floor(length(w) * 0.95):length(w)])
    return(top_5_perc > 0.95)
  }

  weights <- rep(1, N_labellings)
  em_iter <- 0
  counter <- 0
  adaptive_counter <- 0

  pre <- as.data.frame(starting_data) %>% dplyr::filter(.data$t < treatment_time)
  post <- as.data.frame(starting_data) %>% dplyr::filter(.data$t >= treatment_time)
  post <- post[order(post$t), ]

  calculate_weights <- function(labellings, treat_par, control_par, ...) {
    sapply(labellings, function(y) {
      realiz <- y %>% dplyr::filter(.data$t >= treatment_time)
      include <- which(realiz$inferred_process == "control")
      if (length(include) == 0) return(-Inf)
      control_lik <- loglik_hawk_fast(
        params = control_par, realiz = realiz[include, ],
        windowT = c(treatment_time, max(starting_data$t)),
        windowS = statespace, zero_background_region = treated_state_space,
        density_approx = FALSE, numeric_integral = FALSE, ...
      )
      include <- which(realiz$inferred_process == "treated")
      if (length(include) == 0) return(-Inf)
      treat_lik <- loglik_hawk_fast(
        params = treat_par, realiz = realiz[include, ],
        windowT = c(treatment_time, max(starting_data$t)),
        windowS = statespace, zero_background_region = control_state_space,
        density_approx = FALSE, numeric_integral = FALSE, ...
      )
      return(control_lik + treat_lik)
    })
  }

  t_main_sem_start <- proc.time()[3]
  while (counter < N_iter) {
    if (verbose) print(paste0("doing iteration ", counter))
    if (check_weights(weights) == TRUE | adaptive_counter == 0) {
      print("running adaptive step")
      if (adaptive_counter != 0) {
        adapt <- adaptive_step(starting_data = labellings[[which.max(weights)]])
      } else {
        adapt <- adaptive_step(starting_data = starting_data)
      }
      print("adaptive step complete!")
      c_params <- adapt$control_par
      t_params <- adapt$treated_par
      print("generating labellings")
      labellings <- lapply(1:N_labellings, function(i) {
        simulation_labeling_hawkes_hawkes_fast(
          adapt$adaptive_labelling,
          partition = partition, partition_process = partition_processes,
          statespace = statespace,
          state_spaces = adaptive_control$state_spaces,
          windowT = c(treatment_time, max(starting_data$t)),
          hawkes_params_control = hawkes_params_control,
          hawkes_params_treated = t_params[[length(t_params)]],
          change_factor = adaptive_control$change_factor,
          filtration = pre, proximity_weight = 0, verbose = FALSE, ...
        )
      })
      counter <- 0
      adaptive_counter <- adaptive_counter + 1
    }
    raw_weights <- calculate_weights(
      labellings,
      treat_par = unlist(t_params[[length(t_params)]]),
      control_par = unlist(c_params[[length(c_params)]]),
      ...
    )
    keepers <- which(raw_weights != 0 & !is.na(raw_weights) & is.finite(raw_weights))
    if (verbose) {
      cat(sprintf("\n[Outer iter %d] raw weights: %s\n", counter,
                  paste(signif(raw_weights, 4), collapse = ", ")))
      cat(sprintf("  keepers: %d / %d labellings\n", length(keepers), length(labellings)))
    }
    if (length(keepers) == 0) {
      warning("All weights are zero/NA/Inf — skipping outer optim")
      counter <- counter + 1
      next
    }
    weights <- normalize_weights(raw_weights[keepers])
    if (verbose) {
      cat("  normalized weights: ", paste(signif(weights, 4), collapse = ", "), "\n")
    }

    fp <- if (!is.null(adaptive_control$fixed_params)) adaptive_control$fixed_params else NULL
    all_names <- c("mu", "alpha", "beta", "K")
    fp_idx <- if (!is.null(fp)) match(names(fp), all_names) else integer(0)
    fr_idx <- setdiff(seq_along(all_names), fp_idx)

    run_outer_optim <- function(process_label, zero_bg_region, par_list) {
      obj_fn <- function(params) {
        liks <- sapply(labellings[keepers], function(y) {
          loglik_hawk_fast(
            params = params,
            realiz = y %>% dplyr::filter(.data$t >= treatment_time,
                                         .data$inferred_process == process_label),
            windowT = c(treatment_time, max(starting_data$t)),
            windowS = statespace, zero_background_region = zero_bg_region,
            density_approx = FALSE, numeric_integral = FALSE,
            background_rate_var = "W"
          )
        })
        return(sum(liks * weights))
      }

      full_vec <- unlist(par_list)
      if (verbose) {
        cat(sprintf("  [%s] starting par: %s\n", process_label,
                    paste(all_names, signif(full_vec, 5), sep = "=", collapse = "  ")))
        cat(sprintf("  [%s] starting lik: %s\n", process_label, signif(obj_fn(full_vec), 6)))
      }

      t0 <- proc.time()[3]
      if (length(fp_idx) > 0) {
        wrap_fn <- function(free_par) {
          p4 <- full_vec; p4[fr_idx] <- free_par
          obj_fn(p4)
        }
        res <- tryCatch(
          optim(par = full_vec[fr_idx], fn = wrap_fn, method = "Nelder-Mead",
                control = list(fnscale = -1, trace = 1 * verbose, maxit = 1000)),
          error = function(e) { cat(sprintf("  [%s] OPTIM ERROR: %s\n", process_label, e$message)); list(par = full_vec[fr_idx], convergence = -99) }
        )
        out_par <- full_vec; out_par[fr_idx] <- res$par
        res$par <- out_par
      } else {
        res <- tryCatch(
          optim(par = full_vec, fn = obj_fn, method = "Nelder-Mead",
                control = list(fnscale = -1, trace = 1 * verbose, maxit = 1000)),
          error = function(e) { cat(sprintf("  [%s] OPTIM ERROR: %s\n", process_label, e$message)); list(par = full_vec, convergence = -99) }
        )
      }

      if (verbose) {
        cat(sprintf("  [%s] convergence: %s\n", process_label,
                    if (!is.null(res$convergence)) res$convergence else "NULL"))
        cat(sprintf("  [%s] final par:   %s\n", process_label,
                    paste(all_names, signif(as.numeric(res$par), 5), sep = "=", collapse = "  ")))
        cat(sprintf("  [%s] final lik: %s  (took %.1fs)\n", process_label,
                    signif(obj_fn(res$par), 6), proc.time()[3] - t0))
      }
      return(res$par)
    }

    if (verbose) cat("\n--- Outer optim: treated ---\n")
    fit_t <- run_outer_optim("treated", control_state_space, t_params[[length(t_params)]])
    t_params[[length(t_params) + 1]] <- fit_t

    if (verbose) cat("\n--- Outer optim: control ---\n")
    fit_c <- run_outer_optim("control", treated_state_space, c_params[[length(c_params)]])
    c_params[[length(c_params) + 1]] <- fit_c

    counter <- counter + 1
  }
  t_main_sem_end <- proc.time()[3]

  print(paste0("Total time taken for SEM: ", signif(proc.time()[3] - t_global, 2)))
  return(list(
    hawkes_params_control = c_params[[length(c_params)]],
    hawkes_params_treated = t_params[[length(t_params)]],
    t_params = t_params,
    labellings = labellings,
    adaptive = adapt,
    time = proc.time()[3] - t_global,
    time_main_sem = t_main_sem_end - t_main_sem_start
  ))
}

#' Check convergence of the SEM algorithm
#'
#' Plots the proportion of points in each labeling category over iterations.
#'
#' @param result Output from adaptive_SEM
#' @return Data frame of proportions per iteration
#' @export
convergence_check <- function(result) {
  check <- do.call(rbind, lapply(result$adaptive$class_results, as.data.frame))
  check$total <- rowSums(check)
  for (i in names(check)) {
    check[[i]] <- check[[i]] / check$total
  }
  check$iter <- 1:nrow(check)
  check <- check[, which(names(check) != "total")]
  check_long <- check %>% reshape2::melt(id.vars = "iter")
  p <- ggplot2::ggplot(check_long, ggplot2::aes(x = .data$iter, y = .data$value, color = .data$variable)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      title = "Proportion of Points in Each Category Over SEM Iterations",
      x = "Iteration", y = "Proportion"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.title = ggplot2::element_blank())
  print(p)
  return(check)
}
