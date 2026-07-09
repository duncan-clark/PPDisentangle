# Forward-simulation decay-validation helpers shared by sim_study.R and
# sim_study_robustness.R.

normalize_hawkes_kernel <- getFromNamespace("normalize_hawkes_kernel", "PPDisentangle")

choose_decay_flip_cell <- function(partition, partition_processes, omega, requested_cell = NA_integer_) {
  if (is.finite(requested_cell) && !is.na(requested_cell) &&
      requested_cell >= 1L && requested_cell <= partition$n) {
    return(as.integer(requested_cell))
  }
  centers <- do.call(rbind, lapply(seq_len(partition$n), function(i) {
    wi <- as.owin(partition[i])
    c(x = mean(wi$xrange), y = mean(wi$yrange))
  }))
  omega_center <- c(mean(omega[1:2]), mean(omega[3:4]))
  d2 <- (centers[, "x"] - omega_center[1])^2 + (centers[, "y"] - omega_center[2])^2
  treated_cells <- which(partition_processes == "treated")
  candidates <- if (length(treated_cells) > 0L) treated_cells else seq_len(partition$n)
  candidates[which.min(d2[candidates])]
}

point_distance_to_rect <- function(x, y, win) {
  dx <- pmax(win$xrange[1] - x, 0, x - win$xrange[2])
  dy <- pmax(win$yrange[1] - y, 0, y - win$yrange[2])
  sqrt(dx^2 + dy^2)
}

point_distance_to_point <- function(x, y, x0, y0) {
  sqrt((x - x0)^2 + (y - y0)^2)
}

sample_flip_point_in_cell <- function(cell_win, flip_time) {
  pp <- spatstat.random::runifpoint(n = 1L, win = cell_win)
  data.frame(x = pp$x, y = pp$y, t = flip_time)
}

resolve_single_point_flip_context <- function(
    partition,
    partition_processes,
    omega,
    hawkes_par_1,
    hawkes_par_2,
    flip_cell = NA_integer_) {
  flip_cell <- choose_decay_flip_cell(partition, partition_processes, omega, flip_cell)
  cell_win <- as.owin(partition[flip_cell])
  original_process <- partition_processes[flip_cell]
  flipped_process <- if (identical(original_process, "treated")) "control" else "treated"
  params_original <- if (identical(original_process, "treated")) hawkes_par_2 else hawkes_par_1
  params_flipped <- if (identical(flipped_process, "treated")) hawkes_par_2 else hawkes_par_1
  list(
    flip_cell = flip_cell,
    cell_win = cell_win,
    original_process = original_process,
    flipped_process = flipped_process,
    params_original = params_original,
    params_flipped = params_flipped
  )
}

simulate_single_point_catalogue <- function(params, flip_point, end_time, omega_win) {
  t0 <- flip_point$t[1]
  if (!is.finite(end_time) || !is.finite(t0) || end_time <= t0) {
    return(data.frame(x = numeric(0), y = numeric(0), t = numeric(0)))
  }
  ev <- getFromNamespace("sim_hawkes_fast", "PPDisentangle")(
    params = params,
    windowT = c(t0, end_time),
    windowS = omega_win,
    background_realization = flip_point
  )
  data.frame(x = ev$x, y = ev$y, t = ev$t)
}

annulus_counts_from_point <- function(catalogue, x0, y0, annulus_width, max_bin) {
  if (nrow(catalogue) < 1L) {
    bins <- integer(0)
  } else {
    dist <- point_distance_to_point(catalogue$x, catalogue$y, x0, y0)
    bins <- pmin(max_bin, floor(dist / annulus_width))
  }
  as.numeric(tabulate(bins + 1L, nbins = max_bin + 1L))
}

time_bin_counts_since_event <- function(catalogue, event_time, time_bin_width, max_bin) {
  counts <- rep(0, max_bin + 1L)
  if (nrow(catalogue) < 1L) return(counts)
  elapsed <- catalogue$t - event_time
  keep <- is.finite(elapsed) & elapsed >= 0
  if (!any(keep)) return(counts)
  bins <- pmin(max_bin, floor(elapsed[keep] / time_bin_width))
  as.numeric(tabulate(bins + 1L, nbins = max_bin + 1L))
}

run_single_point_label_flip_reps <- function(
    flip_ctx,
    omega,
    omega_win,
    treatment_time,
    end_time,
    decay_reps,
    count_fn,
    bins_df,
    stage_seed_fn,
    stage_offset,
    log_fn = function(...) invisible(NULL)) {
  cell_win <- flip_ctx$cell_win
  cell_center_x <- mean(cell_win$xrange)
  cell_center_y <- mean(cell_win$yrange)
  log_fn(
    "Running single-point label-flip diagnostic: reps=", decay_reps,
    " | flip_cell=", flip_ctx$flip_cell,
    " | ", flip_ctx$original_process, "->", flip_ctx$flipped_process
  )
  t0 <- proc.time()[3]
  rep_rows <- lapply(seq_len(decay_reps), function(rep_id) {
    set.seed(stage_seed_fn(stage_offset, rep_id, 0L))
    flip_point <- sample_flip_point_in_cell(cell_win, treatment_time)
    set.seed(stage_seed_fn(stage_offset, rep_id, 1L))
    cat_original <- simulate_single_point_catalogue(
      params = flip_ctx$params_original,
      flip_point = flip_point,
      end_time = end_time,
      omega_win = omega_win
    )
    set.seed(stage_seed_fn(stage_offset, rep_id, 1L))
    cat_flipped <- simulate_single_point_catalogue(
      params = flip_ctx$params_flipped,
      flip_point = flip_point,
      end_time = end_time,
      omega_win = omega_win
    )
    count_original <- count_fn(cat_original, flip_point)
    count_flipped <- count_fn(cat_flipped, flip_point)
    abs_delta_n <- abs(count_flipped - count_original)
    pct_delta_n <- ifelse(
      count_original > 0,
      100 * abs_delta_n / count_original,
      NA_real_
    )
    cbind(
      data.frame(
        rep = rep_id,
        flip_x = flip_point$x[1],
        flip_y = flip_point$y[1],
        flip_t = flip_point$t[1],
        stringsAsFactors = FALSE
      ),
      bins_df,
      data.frame(
        count_original = count_original,
        count_flipped = count_flipped,
        abs_delta_n = abs_delta_n,
        pct_delta_n = pct_delta_n,
        stringsAsFactors = FALSE
      )
    )
  })
  decay_df <- do.call(rbind, rep_rows)
  log_fn(sprintf("Single-point label-flip diagnostic finished in %.1fs", proc.time()[3] - t0))
  list(
    df = decay_df,
    cell_center_x = cell_center_x,
    cell_center_y = cell_center_y
  )
}

summarize_decay_bins <- function(decay_df, bin_cols) {
  decay_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("annulus", bin_cols)))) %>%
    dplyr::summarize(
      mean_abs_delta = mean(.data$abs_delta_n),
      mean_pct_delta = mean(.data$pct_delta_n, na.rm = TRUE),
      mean_count_original = mean(.data$count_original),
      sd_abs_delta = stats::sd(.data$abs_delta_n),
      q10_abs_delta = stats::quantile(.data$abs_delta_n, 0.10, names = FALSE),
      q90_abs_delta = stats::quantile(.data$abs_delta_n, 0.90, names = FALSE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      mean_pct_delta = dplyr::if_else(
        is.finite(.data$mean_pct_delta),
        .data$mean_pct_delta,
        dplyr::if_else(
          is.finite(.data$mean_count_original) & .data$mean_count_original > 0,
          100 * .data$mean_abs_delta / .data$mean_count_original,
          NA_real_
        )
      )
    )
}

make_decay_rate_lines <- function(summary_df, control_pp, treated_pp) {
  positive <- summary_df[is.finite(summary_df$mean_abs_delta) & summary_df$mean_abs_delta > 0, , drop = FALSE]
  if (nrow(positive) < 1L) return(NULL)
  anchor <- positive[1, , drop = FALSE]
  refs <- data.frame(
    reference = c("control spatial rate", "treated spatial rate"),
    alpha = c(as.numeric(control_pp$alpha), as.numeric(treated_pp$alpha)),
    stringsAsFactors = FALSE
  )
  refs <- refs[is.finite(refs$alpha) & refs$alpha > 0, , drop = FALSE]
  if (nrow(refs) < 1L) return(NULL)
  refs <- refs[!duplicated(round(refs$alpha, 12)), , drop = FALSE]
  do.call(rbind, lapply(seq_len(nrow(refs)), function(i) {
    out <- summary_df
    out$reference <- refs$reference[i]
    out$alpha <- refs$alpha[i]
    out$rate_value <- anchor$mean_abs_delta *
      exp(-refs$alpha[i] * (out$d_mid^2 - anchor$d_mid^2))
    out
  }))
}

rebuild_partition_processes <- function(n_cells, nx, ny, treat_prop = 0.5, base_seed = 123L) {
  out <- rep("control", n_cells)
  set.seed(as.integer(base_seed + 42L))
  out[sample(seq_len(nx * ny), nx * ny * treat_prop)] <- "treated"
  out
}

run_decay_validation_scenario <- function(
    omega,
    nx,
    ny,
    treatment_time,
    end_time,
    partition_processes,
    hawkes_par_1,
    hawkes_par_2,
    sim_kernel,
    decay_reps = 200L,
    annulus_width = 1,
    flip_cell = NA_integer_,
    stage_seed_fn = function(stage_offset, sim_id = 0L, extra = 0L) {
      as.integer(123L + as.integer(stage_offset) * 100000L + as.integer(sim_id) * 1000L + as.integer(extra))
    },
    log_fn = function(...) invisible(NULL)) {
  if (!identical(normalize_hawkes_kernel(sim_kernel), "exponential") &&
      !identical(normalize_hawkes_kernel(sim_kernel), "power_law")) {
    return(NULL)
  }
  partition <- spatstat.geom::quadrats(X = omega, nx = nx, ny = ny)
  if (length(partition_processes) != partition$n) {
    partition_processes <- rebuild_partition_processes(partition$n, nx, ny)
  }
  omega_win <- as.owin(omega)
  flip_ctx <- resolve_single_point_flip_context(
    partition = partition,
    partition_processes = partition_processes,
    omega = omega,
    hawkes_par_1 = hawkes_par_1,
    hawkes_par_2 = hawkes_par_2,
    flip_cell = flip_cell
  )
  cell_win <- flip_ctx$cell_win
  max_corner_dist <- max(point_distance_to_rect(
    x = c(omega[1], omega[1], omega[2], omega[2]),
    y = c(omega[3], omega[4], omega[3], omega[4]),
    win = cell_win
  ))
  max_bin <- max(1L, ceiling(max_corner_dist / annulus_width))
  annuli <- data.frame(
    annulus = seq.int(0L, max_bin),
    d_left = seq.int(0L, max_bin) * annulus_width,
    d_mid = (seq.int(0L, max_bin) + 0.5) * annulus_width
  )
  count_fn <- function(catalogue, flip_point) {
    annulus_counts_from_point(
      catalogue = catalogue,
      x0 = flip_point$x[1],
      y0 = flip_point$y[1],
      annulus_width = annulus_width,
      max_bin = max_bin
    )
  }
  out <- run_single_point_label_flip_reps(
    flip_ctx = flip_ctx,
    omega = omega,
    omega_win = omega_win,
    treatment_time = treatment_time,
    end_time = end_time,
    decay_reps = decay_reps,
    count_fn = count_fn,
    bins_df = annuli,
    stage_seed_fn = stage_seed_fn,
    stage_offset = 13L,
    log_fn = log_fn
  )
  decay_df <- out$df
  decay_summary <- summarize_decay_bins(decay_df, c("d_left", "d_mid"))
  rate_lines <- make_decay_rate_lines(decay_summary, hawkes_par_1, hawkes_par_2)
  eps <- min(decay_summary$mean_abs_delta[decay_summary$mean_abs_delta > 0], na.rm = TRUE)
  if (!is.finite(eps)) eps <- 1e-6
  eps <- eps / 2
  decay_summary$mean_abs_delta_plot <- pmax(decay_summary$mean_abs_delta, eps)
  if (!is.null(rate_lines) && nrow(rate_lines) > 0) {
    rate_lines$rate_value_plot <- pmax(rate_lines$rate_value, eps)
  }
  list(
    specs = list(
      enabled = TRUE,
      flip_mode = "single_point_label",
      reps = decay_reps,
      annulus_width = annulus_width,
      flip_cell = flip_ctx$flip_cell,
      flip_time = treatment_time,
      original_process = flip_ctx$original_process,
      flipped_process = flip_ctx$flipped_process,
      slope_reference = "exp(-alpha * d^2), anchored to first positive empirical annulus"
    ),
    df = decay_df,
    summary = decay_summary,
    rate_lines = rate_lines
  )
}

default_decay_time_bin_width <- function(treatment_time, end_time, requested_width = NA_real_) {
  if (is.finite(requested_width) && requested_width > 0) return(requested_width)
  horizon <- end_time - treatment_time
  if (!is.finite(horizon) || horizon <= 0) return(1)
  max(1, horizon / 40)
}

run_temporal_decay_validation_scenario <- function(
    omega,
    nx,
    ny,
    treatment_time,
    end_time,
    partition_processes,
    hawkes_par_1,
    hawkes_par_2,
    sim_kernel,
    decay_reps = 200L,
    time_bin_width = NA_real_,
    flip_cell = NA_integer_,
    stage_seed_fn = function(stage_offset, sim_id = 0L, extra = 0L) {
      as.integer(123L + as.integer(stage_offset) * 100000L + as.integer(sim_id) * 1000L + as.integer(extra))
    },
    log_fn = function(...) invisible(NULL)) {
  if (!identical(normalize_hawkes_kernel(sim_kernel), "exponential") &&
      !identical(normalize_hawkes_kernel(sim_kernel), "power_law")) {
    return(NULL)
  }
  time_bin_width <- default_decay_time_bin_width(treatment_time, end_time, time_bin_width)
  partition <- spatstat.geom::quadrats(X = omega, nx = nx, ny = ny)
  if (length(partition_processes) != partition$n) {
    partition_processes <- rebuild_partition_processes(partition$n, nx, ny)
  }
  omega_win <- as.owin(omega)
  flip_ctx <- resolve_single_point_flip_context(
    partition = partition,
    partition_processes = partition_processes,
    omega = omega,
    hawkes_par_1 = hawkes_par_1,
    hawkes_par_2 = hawkes_par_2,
    flip_cell = flip_cell
  )
  post_treatment_horizon <- end_time - treatment_time
  max_bin <- max(1L, ceiling(post_treatment_horizon / time_bin_width))
  time_bins <- data.frame(
    annulus = seq.int(0L, max_bin),
    t_left = seq.int(0L, max_bin) * time_bin_width,
    t_mid = (seq.int(0L, max_bin) + 0.5) * time_bin_width
  )
  count_fn <- function(catalogue, flip_point) {
    time_bin_counts_since_event(
      catalogue = catalogue,
      event_time = flip_point$t[1],
      time_bin_width = time_bin_width,
      max_bin = max_bin
    )
  }
  log_fn(
    "Running temporal single-point label-flip diagnostic: reps=", decay_reps,
    " | flip_cell=", flip_ctx$flip_cell,
    " | time_bin_width=", time_bin_width,
    " | ", flip_ctx$original_process, "->", flip_ctx$flipped_process
  )
  out <- run_single_point_label_flip_reps(
    flip_ctx = flip_ctx,
    omega = omega,
    omega_win = omega_win,
    treatment_time = treatment_time,
    end_time = end_time,
    decay_reps = decay_reps,
    count_fn = count_fn,
    bins_df = time_bins,
    stage_seed_fn = stage_seed_fn,
    stage_offset = 14L,
    log_fn = log_fn
  )
  decay_df <- out$df
  decay_summary <- summarize_decay_bins(decay_df, c("t_left", "t_mid"))
  eps <- min(decay_summary$mean_abs_delta[decay_summary$mean_abs_delta > 0], na.rm = TRUE)
  if (!is.finite(eps)) eps <- 1e-6
  eps <- eps / 2
  decay_summary$mean_abs_delta_plot <- pmax(decay_summary$mean_abs_delta, eps)
  list(
    specs = list(
      enabled = TRUE,
      flip_mode = "single_point_label",
      axis = "lag_since_flip_event",
      reps = decay_reps,
      time_bin_width = time_bin_width,
      flip_cell = flip_ctx$flip_cell,
      flip_time = treatment_time,
      original_process = flip_ctx$original_process,
      flipped_process = flip_ctx$flipped_process
    ),
    df = decay_df,
    summary = decay_summary
  )
}

hawkes_exponential_spatial_mean <- function(alpha) {
  alpha <- as.numeric(alpha)[1]
  if (!is.finite(alpha) || alpha <= 0) return(NA_real_)
  sqrt(pi) / (2 * sqrt(alpha))
}

hawkes_power_law_spatial_mean <- function(spatial_d, q) {
  spatial_d <- as.numeric(spatial_d)[1]
  q <- as.numeric(q)[1]
  if (!is.finite(spatial_d) || spatial_d <= 0 || !is.finite(q) || q <= 1.5) return(NA_real_)
  sqrt(spatial_d) * (q - 1) * sqrt(pi) / 2 * gamma(q - 1.5) / gamma(q)
}

hawkes_power_law_spatial_d_for_mean <- function(target_mean, q) {
  target_mean <- as.numeric(target_mean)[1]
  q <- as.numeric(q)[1]
  if (!is.finite(target_mean) || target_mean <= 0 || !is.finite(q) || q <= 1.5) return(NA_real_)
  factor <- (q - 1) * sqrt(pi) / 2 * gamma(q - 1.5) / gamma(q)
  if (!is.finite(factor) || factor <= 0) return(NA_real_)
  (target_mean / factor)^2
}

hawkes_exponential_temporal_mean <- function(beta) {
  beta <- as.numeric(beta)[1]
  if (!is.finite(beta) || beta <= 0) return(NA_real_)
  1 / beta
}

hawkes_power_law_temporal_mean <- function(cc, p) {
  cc <- as.numeric(cc)[1]
  p <- as.numeric(p)[1]
  if (!is.finite(cc) || cc <= 0 || !is.finite(p) || p <= 2) return(Inf)
  cc / (p - 2)
}

hawkes_power_law_c_for_mean <- function(target_mean, p) {
  target_mean <- as.numeric(target_mean)[1]
  p <- as.numeric(p)[1]
  if (!is.finite(target_mean) || target_mean <= 0 || !is.finite(p) || p <= 2) return(NA_real_)
  target_mean * (p - 2)
}

make_hawkes_params_from_cfg <- function(
    cfg,
    k,
    sim_kernel = "exponential",
    sim_spatial_kernel = "exponential",
    spatial_q = NULL,
    spatial_d = NULL,
    power_c = NULL,
    power_p = NULL) {
  mu <- cfg$TRUE_MU
  alpha <- cfg$HAWKES_ALPHA
  as_hawkes <- getFromNamespace("as_hawkes_params", "PPDisentangle")
  normalize_spatial <- getFromNamespace("normalize_hawkes_spatial_kernel", "PPDisentangle")
  sim_kernel <- normalize_hawkes_kernel(sim_kernel)
  sim_spatial_kernel <- normalize_spatial(sim_spatial_kernel)

  if (identical(sim_kernel, "power_law")) {
    cc <- if (!is.null(power_c) && is.finite(power_c)) as.numeric(power_c) else cfg$HAWKES_POWER_C
    pp <- if (!is.null(power_p) && is.finite(power_p)) as.numeric(power_p) else cfg$HAWKES_POWER_P
    out <- list(
      mu = mu,
      alpha = alpha,
      c = cc,
      p = pp,
      beta = cc,
      K = k,
      kernel = "power_law"
    )
  } else {
    out <- list(
      mu = mu,
      alpha = alpha,
      beta = cfg$HAWKES_BETA,
      K = k,
      kernel = "exponential"
    )
  }
  out$spatial_kernel <- sim_spatial_kernel
  if (identical(sim_spatial_kernel, "power_law")) {
    out$spatial_q <- if (!is.null(spatial_q) && is.finite(spatial_q)) {
      as.numeric(spatial_q)
    } else if (!is.null(cfg$HAWKES_SPATIAL_POWER_Q)) {
      as.numeric(cfg$HAWKES_SPATIAL_POWER_Q)
    } else {
      2
    }
  }
  as_hawkes(out, out$kernel, out$spatial_kernel, out$spatial_q, spatial_d)
}

run_decay_from_cfg <- function(
    cfg,
    res,
    hawkes_par_1,
    hawkes_par_2,
    sim_kernel,
    decay_reps,
    base_seed = 123L) {
  partition <- spatstat.geom::quadrats(X = cfg$OMEGA, nx = cfg$NX, ny = cfg$NY)
  seed_base <- if (!is.null(cfg$BASE_SEED)) cfg$BASE_SEED else base_seed
  partition_processes <- rebuild_partition_processes(partition$n, cfg$NX, cfg$NY, base_seed = seed_base)
  flip_cell <- if (!is.null(res$decay_validation$specs$flip_cell)) {
    res$decay_validation$specs$flip_cell
  } else {
    cfg$DECAY_FLIP_CELL
  }
  stage_seed_fn <- function(stage_offset, sim_id = 0L, extra = 0L) {
    as.integer(seed_base + as.integer(stage_offset) * 100000L + as.integer(sim_id) * 1000L + as.integer(extra))
  }
  run_decay_validation_scenario(
    omega = cfg$OMEGA,
    nx = cfg$NX,
    ny = cfg$NY,
    treatment_time = cfg$TREATMENT_TIME,
    end_time = cfg$END_TIME,
    partition_processes = partition_processes,
    hawkes_par_1 = hawkes_par_1,
    hawkes_par_2 = hawkes_par_2,
    sim_kernel = sim_kernel,
    decay_reps = decay_reps,
    annulus_width = if (!is.null(cfg$DECAY_ANNULUS_WIDTH)) cfg$DECAY_ANNULUS_WIDTH else 1,
    flip_cell = flip_cell,
    stage_seed_fn = stage_seed_fn,
    log_fn = message
  )
}

run_temporal_decay_from_cfg <- function(
    cfg,
    res,
    hawkes_par_1,
    hawkes_par_2,
    sim_kernel,
    decay_reps,
    base_seed = 123L) {
  partition <- spatstat.geom::quadrats(X = cfg$OMEGA, nx = cfg$NX, ny = cfg$NY)
  seed_base <- if (!is.null(cfg$BASE_SEED)) cfg$BASE_SEED else base_seed
  partition_processes <- rebuild_partition_processes(partition$n, cfg$NX, cfg$NY, base_seed = seed_base)
  flip_cell <- if (!is.null(res$decay_validation$specs$flip_cell)) {
    res$decay_validation$specs$flip_cell
  } else {
    cfg$DECAY_FLIP_CELL
  }
  stage_seed_fn <- function(stage_offset, sim_id = 0L, extra = 0L) {
    as.integer(seed_base + as.integer(stage_offset) * 100000L + as.integer(sim_id) * 1000L + as.integer(extra))
  }
  time_bin_width <- if (!is.null(cfg$DECAY_TIME_BIN_WIDTH)) cfg$DECAY_TIME_BIN_WIDTH else NA_real_
  run_temporal_decay_validation_scenario(
    omega = cfg$OMEGA,
    nx = cfg$NX,
    ny = cfg$NY,
    treatment_time = cfg$TREATMENT_TIME,
    end_time = cfg$END_TIME,
    partition_processes = partition_processes,
    hawkes_par_1 = hawkes_par_1,
    hawkes_par_2 = hawkes_par_2,
    sim_kernel = sim_kernel,
    decay_reps = decay_reps,
    time_bin_width = time_bin_width,
    flip_cell = flip_cell,
    stage_seed_fn = stage_seed_fn,
    log_fn = message
  )
}

refresh_temporal_decay_for_spec <- function(
    rds_path,
    control_k,
    treated_k,
    sim_kernel,
    decay_reps,
    showcase_id,
    decay_label,
    sim_spatial_kernel = "exponential",
    spatial_q = NULL,
    spatial_d = NULL,
    power_c = NULL,
    power_p = NULL,
    base_seed = 123L) {
  if (!file.exists(rds_path)) stop(sprintf("RDS not found: %s", rds_path))
  res <- readRDS(rds_path)
  cfg <- res$config
  hawkes_par_1 <- make_hawkes_params_from_cfg(
    cfg, control_k, sim_kernel,
    sim_spatial_kernel = sim_spatial_kernel,
    spatial_q = spatial_q, spatial_d = spatial_d,
    power_c = power_c, power_p = power_p
  )
  hawkes_par_2 <- make_hawkes_params_from_cfg(
    cfg, treated_k, sim_kernel,
    sim_spatial_kernel = sim_spatial_kernel,
    spatial_q = spatial_q, spatial_d = spatial_d,
    power_c = power_c, power_p = power_p
  )
  out <- run_temporal_decay_from_cfg(
    cfg = cfg,
    res = res,
    hawkes_par_1 = hawkes_par_1,
    hawkes_par_2 = hawkes_par_2,
    sim_kernel = sim_kernel,
    decay_reps = decay_reps,
    base_seed = base_seed
  )
  if (is.null(out)) return(NULL)
  summary_df <- as.data.frame(out$summary)
  summary_df$showcase_id <- showcase_id
  summary_df$decay_label <- decay_label
  summary_df$control_k <- control_k
  summary_df$treated_k <- treated_k
  summary_df$sim_kernel <- sim_kernel
  summary_df$sim_spatial_kernel <- sim_spatial_kernel
  if (!is.null(spatial_q)) summary_df$spatial_q <- spatial_q
  if (!is.null(power_p)) summary_df$power_p <- power_p
  summary_df
}

refresh_decay_from_rds <- function(rds_path, decay_reps, base_seed = 123L) {
  if (!file.exists(rds_path)) stop(sprintf("RDS not found: %s", rds_path))
  res <- readRDS(rds_path)
  cfg <- res$config
  if (is.null(cfg$hawkes_par_1) || is.null(cfg$hawkes_par_2)) {
    stop(sprintf("RDS config missing Hawkes parameters: %s", rds_path))
  }
  run_decay_from_cfg(
    cfg = cfg,
    res = res,
    hawkes_par_1 = cfg$hawkes_par_1,
    hawkes_par_2 = cfg$hawkes_par_2,
    sim_kernel = cfg$SIM_KERNEL,
    decay_reps = decay_reps,
    base_seed = base_seed
  )
}

refresh_decay_for_spec <- function(
    rds_path,
    control_k,
    treated_k,
    sim_kernel,
    decay_reps,
    showcase_id,
    decay_label,
    sim_spatial_kernel = "exponential",
    spatial_q = NULL,
    spatial_d = NULL,
    power_c = NULL,
    power_p = NULL,
    base_seed = 123L) {
  if (!file.exists(rds_path)) stop(sprintf("RDS not found: %s", rds_path))
  res <- readRDS(rds_path)
  cfg <- res$config
  hawkes_par_1 <- make_hawkes_params_from_cfg(
    cfg, control_k, sim_kernel,
    sim_spatial_kernel = sim_spatial_kernel,
    spatial_q = spatial_q, spatial_d = spatial_d,
    power_c = power_c, power_p = power_p
  )
  hawkes_par_2 <- make_hawkes_params_from_cfg(
    cfg, treated_k, sim_kernel,
    sim_spatial_kernel = sim_spatial_kernel,
    spatial_q = spatial_q, spatial_d = spatial_d,
    power_c = power_c, power_p = power_p
  )
  out <- run_decay_from_cfg(
    cfg = cfg,
    res = res,
    hawkes_par_1 = hawkes_par_1,
    hawkes_par_2 = hawkes_par_2,
    sim_kernel = sim_kernel,
    decay_reps = decay_reps,
    base_seed = base_seed
  )
  if (is.null(out)) return(NULL)
  summary_df <- as.data.frame(out$summary)
  summary_df$showcase_id <- showcase_id
  summary_df$decay_label <- decay_label
  summary_df$control_k <- control_k
  summary_df$treated_k <- treated_k
  summary_df$sim_kernel <- sim_kernel
  summary_df$sim_spatial_kernel <- sim_spatial_kernel
  if (!is.null(spatial_q)) summary_df$spatial_q <- spatial_q
  if (!is.null(power_p)) summary_df$power_p <- power_p
  summary_df
}
