# Shared helpers for the two structured robustness studies.
#
# The studies deliberately live in inst/sim_study rather than in the package
# API: they reuse the package's Hawkes simulator and SEM implementation while
# keeping study-specific designs, estimands, and reporting out of the public
# interface.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) < 1L) y else x
}

structured_stage_seed <- function(base_seed, stage, replication = 0L, substage = 0L) {
  modulus <- .Machine$integer.max - 1L
  as.integer((as.double(base_seed) + 1000003 * stage +
    1009 * replication + 37 * substage) %% modulus + 1)
}

grid_design_table <- function(partition) {
  centers <- t(vapply(seq_len(partition$n), function(j) {
    w <- spatstat.geom::as.owin(partition[j])
    c(x = mean(w$xrange), y = mean(w$yrange))
  }, numeric(2L)))
  x_levels <- sort(unique(round(centers[, "x"], 10)))
  # Match spatstat's rectangular tile index convention: row 1 is the
  # uppermost row and columns run left-to-right.
  y_levels <- sort(unique(round(centers[, "y"], 10)), decreasing = TRUE)
  out <- data.frame(
    cell = seq_len(partition$n),
    x = centers[, "x"],
    y = centers[, "y"],
    col = match(round(centers[, "x"], 10), x_levels),
    row = match(round(centers[, "y"], 10), y_levels),
    stringsAsFactors = FALSE
  )
  if (length(x_levels) != 10L || length(y_levels) != 10L || nrow(out) != 100L) {
    stop("Structured robustness studies require the existing 10x10 grid.")
  }
  out$X <- ifelse((out$row + out$col) %% 2L == 0L, 1, -1)
  out
}

allocation_to_processes <- function(z) {
  ifelse(as.logical(z), "treated", "control")
}

make_effect_modification_design <- function(partition, seed = 271828L) {
  grid <- grid_design_table(partition)
  set.seed(as.integer(seed))
  treated_plus <- sample(grid$cell[grid$X == 1], 25L)
  treated_minus <- sample(grid$cell[grid$X == -1], 25L)
  z_obs <- grid$cell %in% c(treated_plus, treated_minus)
  z_plus <- grid$X == 1
  z_minus <- grid$X == -1
  # Fixed single-cell interventions from the all-control baseline:
  # one X=+1 cell and one X=-1 cell (chosen from the observed treated set).
  flip_plus_cell <- as.integer(treated_plus[[1L]])
  flip_minus_cell <- as.integer(treated_minus[[1L]])
  z_flip_plus <- rep(FALSE, partition$n)
  z_flip_minus <- rep(FALSE, partition$n)
  z_flip_plus[flip_plus_cell] <- TRUE
  z_flip_minus[flip_minus_cell] <- TRUE
  if (sum(z_obs & grid$X == 1) != 25L ||
      sum(z_obs & grid$X == -1) != 25L) {
    stop("Observed effect-modification allocation is not balanced by X.")
  }
  if (grid$X[flip_plus_cell] != 1 || grid$X[flip_minus_cell] != -1) {
    stop("Single-cell flip targets are not in the requested X strata.")
  }
  list(
    grid = grid,
    seed = as.integer(seed),
    z_obs = z_obs,
    z_plus = z_plus,
    z_minus = z_minus,
    flip_plus_cell = flip_plus_cell,
    flip_minus_cell = flip_minus_cell,
    z_flip_plus = z_flip_plus,
    z_flip_minus = z_flip_minus,
    z_all = rep(TRUE, partition$n),
    z_none = rep(FALSE, partition$n)
  )
}

rook_edge_table <- function(grid) {
  by_rc <- setNames(grid$cell, paste(grid$row, grid$col, sep = ":"))
  edges <- list()
  k <- 0L
  for (r in seq_len(10L)) {
    for (c in seq_len(10L)) {
      from <- unname(by_rc[[paste(r, c, sep = ":")]])
      if (c < 10L) {
        k <- k + 1L
        edges[[k]] <- c(from, unname(by_rc[[paste(r, c + 1L, sep = ":")]]))
      }
      if (r < 10L) {
        k <- k + 1L
        edges[[k]] <- c(from, unname(by_rc[[paste(r + 1L, c, sep = ":")]]))
      }
    }
  }
  out <- do.call(rbind, edges)
  colnames(out) <- c("j", "ell")
  if (nrow(out) != 180L) stop("Expected 180 rook-adjacency edges.")
  out
}

allocation_cut_count <- function(z, edges) {
  sum(as.logical(z[edges[, "j"]]) != as.logical(z[edges[, "ell"]]))
}

allocation_coarseness <- function(z, edges) {
  (180 - allocation_cut_count(z, edges)) / 170
}

make_geometry_transport_design <- function(partition, path_seed = 314159L,
                                           observed_seed = 271828L,
                                           m_values = c(0L, 5L, 10L, 15L, 20L, 25L)) {
  grid <- grid_design_table(partition)
  edges <- rook_edge_table(grid)
  z_cb <- (grid$row + grid$col) %% 2L == 0L
  z_block <- grid$col <= 5L
  lookup <- setNames(grid$cell, paste(grid$row, grid$col, sep = ":"))
  left <- grid[grid$col <= 5L, , drop = FALSE]
  pairs <- data.frame(
    row = left$row,
    left_cell = left$cell,
    right_cell = unname(lookup[paste(left$row, 11L - left$col, sep = ":")]),
    stringsAsFactors = FALSE
  )
  discrepant <- pairs[
    z_cb[pairs$left_cell] != z_block[pairs$left_cell] |
      z_cb[pairs$right_cell] != z_block[pairs$right_cell],
    , drop = FALSE
  ]
  if (nrow(discrepant) != 25L) stop("Expected 25 discrepant reflected pairs.")
  set.seed(as.integer(path_seed))
  discrepant <- discrepant[sample.int(nrow(discrepant)), , drop = FALSE]

  allocations <- lapply(as.integer(m_values), function(m) {
    z <- z_cb
    if (m > 0L) {
      idx <- seq_len(m)
      cells <- c(discrepant$left_cell[idx], discrepant$right_cell[idx])
      z[cells] <- !z[cells]
    }
    z
  })
  names(allocations) <- paste0("m", as.integer(m_values))
  if (!identical(as.logical(allocations[[1L]]), as.logical(z_cb)) ||
      !identical(as.logical(allocations[[length(allocations)]]), as.logical(z_block))) {
    stop("Geometry path does not recover its checkerboard and block endpoints.")
  }
  if (any(vapply(allocations, sum, integer(1L)) != 50L)) {
    stop("Geometry path failed to preserve 50 percent treatment.")
  }

  focal_cells <- unname(lookup[paste(c(1L, 3L, 5L, 7L, 9L), 5L, sep = ":")])
  if (any(!vapply(allocations, function(z) all(z[focal_cells]), logical(1L)))) {
    stop("Focal-band treatment status is not fixed along the geometry path.")
  }

  set.seed(as.integer(observed_seed))
  z_obs <- rep(FALSE, partition$n)
  z_obs[sample.int(partition$n, 50L)] <- TRUE
  coarseness <- vapply(allocations, allocation_coarseness, numeric(1L), edges = edges)
  list(
    grid = grid,
    edges = edges,
    path_seed = as.integer(path_seed),
    observed_seed = as.integer(observed_seed),
    m = as.integer(m_values),
    allocations = allocations,
    z_checkerboard = z_cb,
    z_block = z_block,
    z_obs = z_obs,
    coarseness = coarseness,
    observed_coarseness = allocation_coarseness(z_obs, edges),
    focal_cells = focal_cells,
    discrepant_pairs = discrepant
  )
}

tile_index_for_points <- function(x, y, partition) {
  if (length(x) < 1L) return(integer(0))
  if (identical(partition$type, "rect") &&
      !is.null(partition$xgrid) && !is.null(partition$ygrid)) {
    return(as.integer(PPDisentangle:::tile_index_rect_cpp(
      x, y, partition$xgrid, partition$ygrid
    )))
  }
  as.integer(spatstat.geom::tileindex(x, y, partition))
}

effect_multiplier <- function(X, h) {
  exp(as.numeric(h) * as.numeric(X)) / cosh(as.numeric(h))
}

make_source_lookup <- function(partition, grid, z, h = 0) {
  force(partition); force(grid); force(z); force(h)
  function(x, y) {
    tile <- tile_index_for_points(x, y, partition)
    out <- rep(1, length(tile))
    valid <- is.finite(tile) & tile >= 1L & tile <= length(z)
    treated <- valid & z[tile]
    out[treated] <- effect_multiplier(grid$X[tile[treated]], h)
    out[!valid] <- 0
    out
  }
}

state_spaces_for_allocation <- function(partition, z) {
  z <- as.logical(z)
  processes <- allocation_to_processes(z)
  unique_processes <- unique(processes)
  spaces <- lapply(unique_processes, function(p) {
    spatstat.geom::as.owin(partition[processes == p])
  })
  names(spaces) <- unique_processes
  spaces
}

prepare_pre_history <- function(pre, partition) {
  pre <- as.data.frame(pre)
  if (nrow(pre) < 1L) {
    return(data.frame(
      x = numeric(), y = numeric(), t = numeric(), n = integer(),
      background = logical(), process = character(),
      location_process = character(), inferred_process = character()
    ))
  }
  pre$process <- "control"
  pre$location_process <- "control"
  pre$inferred_process <- "control"
  pre
}

simulate_structured_catalogue <- function(partition, omega, windowT, z,
                                          control_params, treated_params,
                                          filtration = NULL, grid = NULL, h = 0,
                                          t_trunc = NULL) {
  if (is.null(grid)) grid <- grid_design_table(partition)
  processes <- allocation_to_processes(z)
  lookup <- make_source_lookup(partition, grid, z, h)
  out <- PPDisentangle::generate_inhomogeneous_hawkes(
    Omega = omega,
    partition = partition,
    time_window = windowT,
    partition_processes = processes,
    hawkes_params = list(control = control_params, treated = treated_params),
    state_spaces = state_spaces_for_allocation(partition, z),
    filtration = filtration,
    covariate_lookup = lookup,
    t_trunc = t_trunc
  )
  out <- as.data.frame(out)
  if (nrow(out) > 0L) {
    out$tile_index <- tile_index_for_points(out$x, out$y, partition)
    out$X <- grid$X[out$tile_index]
  }
  out
}

attach_design_columns <- function(dat, partition, grid, z, h = 0) {
  dat <- as.data.frame(dat)
  tile <- tile_index_for_points(dat$x, dat$y, partition)
  dat$tile_index <- tile
  dat$X <- grid$X[tile]
  dat$location_process <- allocation_to_processes(z)[tile]
  dat$W_effect <- ifelse(
    z[tile],
    effect_multiplier(dat$X, h),
    1
  )
  dat
}

make_oracle_and_naive_labels <- function(catalogue, treatment_time) {
  oracle <- catalogue
  oracle$inferred_process <- ifelse(
    oracle$t < treatment_time, "control", as.character(oracle$process)
  )
  naive <- catalogue
  naive$inferred_process <- ifelse(
    naive$t < treatment_time, "control", as.character(naive$location_process)
  )
  list(oracle = oracle, naive = naive)
}

fit_effect_modified_hawkes <- function(labelled, partition, grid, z_obs,
                                       windowT, init_params, h_init = 0,
                                       homogeneous = FALSE, maxit = 1000L,
                                       spatial_kernel = "exponential",
                                       spatial_q = 2, t_trunc = NULL,
                                       compute_hessian = TRUE,
                                       fixed_h = NULL) {
  if (isTRUE(homogeneous)) fixed_h <- 0
  post <- as.data.frame(labelled)
  post <- post[
    post$t >= windowT[1L] & post$t <= windowT[2L] &
      post$inferred_process == "treated",
    , drop = FALSE
  ]
  if (nrow(post) < 2L) {
    return(list(
      par = c(mu = NA, alpha = NA, beta = NA, K = NA,
              h = if (!is.null(fixed_h)) fixed_h else NA),
      convergence = 99L, converged = FALSE, h_se = NA_real_,
      hessian_ok = FALSE, value = NA_real_, message = "fewer than two treated events"
    ))
  }
  if (!"tile_index" %in% names(post)) {
    post$tile_index <- tile_index_for_points(post$x, post$y, partition)
  }
  post$X <- grid$X[post$tile_index]
  treated_space <- spatstat.geom::as.owin(partition[as.logical(z_obs)])
  active_area <- spatstat.geom::area(treated_space)
  treated_X_area_mean <- mean(grid$X[as.logical(z_obs)])
  estimates_nonzero_h <- is.null(fixed_h) || abs(as.numeric(fixed_h)) > 1e-12
  if (estimates_nonzero_h && abs(treated_X_area_mean) > 1e-12) {
    stop("The h likelihood requires exact X balance in observed treated cells.")
  }

  p0 <- PPDisentangle:::as_hawkes_params(
    init_params, "exponential", spatial_kernel, spatial_q
  )
  raw0 <- c(
    log_mu = log(max(as.numeric(p0$mu), 1e-8)),
    log_alpha = log(max(as.numeric(p0$alpha), 1e-8)),
    log_beta = log(max(as.numeric(p0$beta), 1e-8)),
    logit_K = stats::qlogis(min(max(as.numeric(p0$K), 1e-6), 1 - 1e-6))
  )
  if (is.null(fixed_h)) raw0 <- c(raw0, h = as.numeric(h_init))
  cpp_ll <- getFromNamespace("hawkes_loglik_inhom_filtration_cpp", "PPDisentangle")
  kernel_type <- getFromNamespace("hawkes_kernel_type", "PPDisentangle")
  spatial_type <- getFromNamespace("hawkes_spatial_kernel_type", "PPDisentangle")
  dt <- diff(windowT)
  parent_t <- post$t
  parent_x <- post$x
  parent_y <- post$y
  spatial_d <- p0$spatial_d %||% NA_real_

  unpack <- function(raw) {
    c(
      mu = exp(raw[["log_mu"]]),
      alpha = exp(raw[["log_alpha"]]),
      beta = exp(raw[["log_beta"]]),
      K = stats::plogis(raw[["logit_K"]]),
      h = if (!is.null(fixed_h)) as.numeric(fixed_h) else raw[["h"]]
    )
  }
  objective <- function(raw) {
    pars <- unpack(raw)
    if (any(!is.finite(pars)) || abs(pars[["h"]]) > 25) return(1e15)
    W <- effect_multiplier(post$X, pars[["h"]])
    ll <- cpp_ll(
      post_t = post$t,
      post_x = post$x,
      post_y = post$y,
      W_val = W,
      parent_t = parent_t,
      parent_x = parent_x,
      parent_y = parent_y,
      mu = pars[["mu"]],
      alpha = pars[["alpha"]],
      beta = pars[["beta"]],
      K = pars[["K"]],
      areaS = active_area,
      t_start = windowT[1L],
      t_end = windowT[2L],
      adjust_factor = 1,
      t_trunc = t_trunc %||% -1,
      kernel_type = kernel_type("exponential"),
      cc = 1,
      p = 2,
      spatial_kernel_type = spatial_type(spatial_kernel),
      spatial_q = spatial_q,
      spatial_d = spatial_d
    )
    if (!is.finite(ll)) 1e15 else -ll
  }
  fit <- tryCatch(
    stats::optim(
      raw0, objective, method = "BFGS",
      hessian = isTRUE(compute_hessian),
      control = list(maxit = as.integer(maxit), reltol = 1e-9)
    ),
    error = function(e) structure(list(message = conditionMessage(e)), class = "fit_error")
  )
  if (inherits(fit, "fit_error")) {
    return(list(
      par = c(mu = NA, alpha = NA, beta = NA, K = NA,
              h = if (!is.null(fixed_h)) fixed_h else NA),
      convergence = 98L, converged = FALSE, h_se = NA_real_,
      hessian_ok = FALSE, value = NA_real_, message = fit$message
    ))
  }
  pars <- unpack(fit$par)
  covariance <- NULL
  hessian_ok <- FALSE
  h_se <- NA_real_
  if (isTRUE(compute_hessian) && !is.null(fit$hessian) &&
      all(is.finite(fit$hessian))) {
    covariance <- tryCatch(solve(fit$hessian), error = function(e) NULL)
    hessian_ok <- !is.null(covariance) && all(is.finite(covariance)) &&
      all(diag(covariance) > 0)
    if (hessian_ok && is.null(fixed_h)) {
      h_se <- sqrt(covariance["h", "h"])
    }
  }
  list(
    par = pars,
    raw_par = fit$par,
    convergence = fit$convergence,
    converged = identical(fit$convergence, 0L) && all(is.finite(pars)),
    value = fit$value,
    hessian = fit$hessian,
    covariance_raw = covariance,
    hessian_ok = hessian_ok,
    h_se = h_se,
    message = fit$message %||% ""
  )
}

profile_effect_h <- function(fit, labelled, partition, grid, z_obs, windowT,
                             h_grid, maxit = 400L) {
  rows <- lapply(h_grid, function(h0) {
    prof <- fit_effect_modified_hawkes(
      labelled = labelled,
      partition = partition,
      grid = grid,
      z_obs = z_obs,
      windowT = windowT,
      init_params = as.list(fit$par[c("mu", "alpha", "beta", "K")]),
      h_init = h0,
      fixed_h = h0,
      maxit = maxit,
      compute_hessian = FALSE
    )
    data.frame(h_grid = h0, h_attained = prof$par[["h"]],
               loglik = -prof$value, converged = prof$converged)
  })
  do.call(rbind, rows)
}

run_effect_aware_sem <- function(catalogue, partition, grid, z_obs,
                                 omega, windowT, treatment_time,
                                 control_params, treated_init,
                                 h_init = 0, interaction_rounds = 2L,
                                 sem_control = list(), seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  h_current <- as.numeric(h_init)
  labelled <- catalogue
  sem_fit <- NULL
  h_fit <- NULL
  defaults <- list(
    N_labellings = 10L,
    N_iter = 10L,
    iter = 100L,
    n_props = 10L,
    change_factor = 0.01,
    param_update_cadence = 10L,
    proposal_update_cadence = 1L,
    update_control_params = FALSE,
    include_starting_data = FALSE,
    update_starting_data = TRUE,
    verbose = FALSE,
    maxit = 1000L
  )
  for (nm in names(defaults)) {
    if (is.null(sem_control[[nm]])) sem_control[[nm]] <- defaults[[nm]]
  }
  for (round_id in seq_len(max(1L, as.integer(interaction_rounds)))) {
    labelled <- attach_design_columns(labelled, partition, grid, z_obs, h_current)
    sem_fit <- PPDisentangle::adaptive_SEM(
      pp_data = labelled,
      partition = partition,
      partition_processes = allocation_to_processes(z_obs),
      statespace = omega,
      time_window = windowT,
      treatment_time = treatment_time,
      hawkes_params_control = control_params,
      hawkes_params_treated = treated_init,
      N_labellings = sem_control$N_labellings,
      N_iter = sem_control$N_iter,
      verbose = isTRUE(sem_control$verbose),
      hawkes_use_filtration_history = TRUE,
      kernel = "exponential",
      spatial_kernel = "exponential",
      background_rate_var = "W_effect",
      adaptive_control = list(
        iter = sem_control$iter,
        n_props = sem_control$n_props,
        change_factor = sem_control$change_factor,
        param_update_cadence = sem_control$param_update_cadence,
        proposal_update_cadence = sem_control$proposal_update_cadence,
        update_control_params = sem_control$update_control_params,
        include_starting_data = sem_control$include_starting_data,
        update_starting_data = sem_control$update_starting_data,
        verbose = isTRUE(sem_control$verbose)
      )
    )
    labelled <- sem_fit$adaptive$adaptive_labelling
    h_fit <- fit_effect_modified_hawkes(
      labelled, partition, grid, z_obs, windowT,
      init_params = sem_fit$hawkes_params_treated %||% treated_init,
      h_init = h_current, homogeneous = FALSE,
      maxit = sem_control$maxit
    )
    if (!isTRUE(h_fit$converged)) break
    h_current <- h_fit$par[["h"]]
    treated_init <- as.list(h_fit$par[c("mu", "alpha", "beta", "K")])
  }
  list(sem = sem_fit, h_fit = h_fit, labelling = labelled, h = h_current)
}

run_homogeneous_sem <- function(catalogue, partition, z_obs, omega,
                                windowT, treatment_time, control_params,
                                treated_init, sem_control = list(), seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  defaults <- list(
    N_labellings = 10L, N_iter = 10L, iter = 100L, n_props = 10L,
    change_factor = 0.01, param_update_cadence = 10L,
    proposal_update_cadence = 1L, update_control_params = FALSE,
    include_starting_data = FALSE, update_starting_data = TRUE,
    verbose = FALSE
  )
  for (nm in names(defaults)) {
    if (is.null(sem_control[[nm]])) sem_control[[nm]] <- defaults[[nm]]
  }
  out <- PPDisentangle::adaptive_SEM(
    pp_data = catalogue,
    partition = partition,
    partition_processes = allocation_to_processes(z_obs),
    statespace = omega,
    time_window = windowT,
    treatment_time = treatment_time,
    hawkes_params_control = control_params,
    hawkes_params_treated = treated_init,
    N_labellings = sem_control$N_labellings,
    N_iter = sem_control$N_iter,
    verbose = isTRUE(sem_control$verbose),
    hawkes_use_filtration_history = TRUE,
    kernel = "exponential",
    spatial_kernel = "exponential",
    adaptive_control = list(
      iter = sem_control$iter,
      n_props = sem_control$n_props,
      change_factor = sem_control$change_factor,
      param_update_cadence = sem_control$param_update_cadence,
      proposal_update_cadence = sem_control$proposal_update_cadence,
      update_control_params = sem_control$update_control_params,
      include_starting_data = sem_control$include_starting_data,
      update_starting_data = sem_control$update_starting_data,
      verbose = isTRUE(sem_control$verbose)
    )
  )
  list(
    sem = out,
    labelling = out$adaptive$adaptive_labelling,
    treated_params = out$hawkes_params_treated
  )
}

label_recovery_by_stratum <- function(labelled, treatment_time) {
  dat <- as.data.frame(labelled)
  dat <- dat[dat$t >= treatment_time, , drop = FALSE]
  if (nrow(dat) < 1L) return(NULL)
  metric_row <- function(d, stratum) {
    truth <- factor(d$process, levels = c("control", "treated"))
    pred <- factor(d$inferred_process, levels = c("control", "treated"))
    tab <- table(truth, pred)
    tn <- unname(tab["control", "control"])
    fp <- unname(tab["control", "treated"])
    fn <- unname(tab["treated", "control"])
    tp <- unname(tab["treated", "treated"])
    recall_treated <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
    recall_control <- if (tn + fp > 0) tn / (tn + fp) else NA_real_
    accuracy <- if (sum(tab) > 0) (tp + tn) / sum(tab) else NA_real_
    data.frame(
      stratum = stratum,
      n = nrow(d),
      accuracy = accuracy,
      balanced_accuracy = mean(c(recall_treated, recall_control), na.rm = TRUE),
      recall_control = recall_control,
      recall_treated = recall_treated,
      stringsAsFactors = FALSE
    )
  }
  rows <- list(metric_row(dat, "all"))
  if ("X" %in% names(dat)) {
    rows <- c(rows, lapply(c(-1, 1), function(x) {
      metric_row(dat[dat$X == x, , drop = FALSE], paste0("X", ifelse(x > 0, "+1", "-1")))
    }))
  }
  do.call(rbind, rows)
}

simulate_regime_counts <- function(partition, omega, windowT, allocations,
                                   control_params, treated_params, filtration,
                                   grid, h = 0, region_cells = seq_len(partition$n),
                                   n_sim = 100L, seed = 1L, t_trunc = NULL) {
  allocation_names <- names(allocations)
  if (is.null(allocation_names)) allocation_names <- paste0("z", seq_along(allocations))
  counts <- matrix(NA_real_, nrow = as.integer(n_sim), ncol = length(allocations),
                   dimnames = list(NULL, allocation_names))
  for (s in seq_len(as.integer(n_sim))) {
    common_seed <- structured_stage_seed(seed, 1L, s)
    for (a in seq_along(allocations)) {
      set.seed(common_seed)
      sim <- simulate_structured_catalogue(
        partition, omega, windowT, allocations[[a]],
        control_params, treated_params, filtration, grid, h, t_trunc
      )
      counts[s, a] <- if (nrow(sim) < 1L) 0 else
        sum(sim$tile_index %in% region_cells)
    }
  }
  means <- colMeans(counts)
  ses <- apply(counts, 2L, stats::sd) / sqrt(nrow(counts))
  data.frame(
    allocation = allocation_names,
    mean = as.numeric(means),
    mc_se = as.numeric(ses),
    stringsAsFactors = FALSE
  )
}

causal_error_summary <- function(df, estimate_col = "estimate", truth_col = "truth",
                                 group_cols) {
  if (nrow(df) < 1L) return(df)
  estimate <- df[[estimate_col]]
  truth <- df[[truth_col]]
  df$.error <- estimate - truth
  df$.ok <- is.finite(estimate) & is.finite(truth)
  dplyr::group_by(df, dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      n_rep = dplyr::n(),
      n_ok = sum(.data$.ok),
      bias = mean(.data$.error[.data$.ok]),
      rmse = sqrt(mean(.data$.error[.data$.ok]^2)),
      mae = mean(abs(.data$.error[.data$.ok])),
      empirical_sd = stats::sd(.data[[estimate_col]][.data$.ok]),
      mc_se_summary = mean(.data$mc_se[.data$.ok], na.rm = TRUE),
      mcse_bias = stats::sd(.data$.error[.data$.ok]) / sqrt(max(1L, sum(.data$.ok))),
      mcse_rmse = if (sum(.data$.ok) > 1L) {
        stats::sd(.data$.error[.data$.ok]^2) /
          (2 * sqrt(mean(.data$.error[.data$.ok]^2)) * sqrt(sum(.data$.ok)))
      } else NA_real_,
      failure_rate = mean(!.data$.ok),
      .groups = "drop"
    )
}
