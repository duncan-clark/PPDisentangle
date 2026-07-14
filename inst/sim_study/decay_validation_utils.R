# Forward-simulation decay-validation helpers shared by sim_study.R and
# sim_study_robustness.R.

normalize_hawkes_kernel <- getFromNamespace("normalize_hawkes_kernel", "PPDisentangle")

cell_centers_matrix <- function(partition) {
  do.call(rbind, lapply(seq_len(partition$n), function(i) {
    wi <- as.owin(partition[i])
    c(x = mean(wi$xrange), y = mean(wi$yrange))
  }))
}

# Treated cells whose centroids lie in the central fraction of Omega
# (default: middle 50% box). Falls back to all treated cells if none qualify.
central_treated_cells <- function(
    partition,
    partition_processes,
    omega,
    central_frac = 0.5) {
  treated_cells <- which(partition_processes == "treated")
  if (length(treated_cells) < 1L) {
    stop("Decay flip requires at least one treated cell (treated -> control).")
  }
  frac <- as.numeric(central_frac)[1]
  if (!is.finite(frac) || frac <= 0 || frac > 1) frac <- 0.5
  centers <- cell_centers_matrix(partition)
  xr <- as.numeric(omega[1:2])
  yr <- as.numeric(omega[3:4])
  cx <- mean(xr)
  cy <- mean(yr)
  half_w <- 0.5 * frac * diff(xr)
  half_h <- 0.5 * frac * diff(yr)
  in_center <- centers[, "x"] >= (cx - half_w) & centers[, "x"] <= (cx + half_w) &
    centers[, "y"] >= (cy - half_h) & centers[, "y"] <= (cy + half_h)
  candidates <- treated_cells[in_center[treated_cells]]
  if (length(candidates) < 1L) {
    omega_center <- c(cx, cy)
    d2 <- (centers[treated_cells, "x"] - omega_center[1])^2 +
      (centers[treated_cells, "y"] - omega_center[2])^2
    return(as.integer(treated_cells[which.min(d2)]))
  }
  as.integer(candidates)
}

choose_decay_flip_cell <- function(partition, partition_processes, omega, requested_cell = NA_integer_) {
  if (is.finite(requested_cell) && !is.na(requested_cell) &&
      requested_cell >= 1L && requested_cell <= partition$n) {
    return(as.integer(requested_cell))
  }
  # Default: the unique cell whose centroid is closest to the window centre.
  # (Do not depend on which cells happen to be treated in a random allocation.)
  centers <- cell_centers_matrix(partition)
  omega_center <- c(mean(omega[1:2]), mean(omega[3:4]))
  d2 <- (centers[, "x"] - omega_center[1])^2 +
    (centers[, "y"] - omega_center[2])^2
  as.integer(which.min(d2))
}

# Ensure the chosen flip cell is treated in the baseline allocation.
force_flip_cell_treated <- function(partition_processes, flip_cell) {
  out <- as.character(partition_processes)
  out[as.integer(flip_cell)] <- "treated"
  out
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

sample_decay_flip_time <- function(treatment_time, end_time, vary_time = FALSE) {
  if (!isTRUE(vary_time)) return(as.numeric(treatment_time))
  horizon <- as.numeric(end_time) - as.numeric(treatment_time)
  if (!is.finite(horizon) || horizon <= 0) return(as.numeric(treatment_time))
  # Keep at least 25% of the post-treatment window for offspring.
  latest <- as.numeric(treatment_time) + 0.75 * horizon
  stats::runif(1L, min = as.numeric(treatment_time), max = max(as.numeric(treatment_time), latest))
}

sample_decay_flip_context <- function(
    partition,
    partition_processes,
    hawkes_par_1,
    hawkes_par_2,
    treatment_time,
    end_time,
    omega,
    vary_time = FALSE,
    central_frac = 0.5) {
  candidates <- central_treated_cells(
    partition = partition,
    partition_processes = partition_processes,
    omega = omega,
    central_frac = central_frac
  )
  flip_cell <- as.integer(candidates[[sample.int(length(candidates), 1L)]])
  cell_win <- as.owin(partition[flip_cell])
  # Always flip treated -> control so Delta N = N_control - N_treated is an expected increase.
  original_process <- "treated"
  flipped_process <- "control"
  params_original <- hawkes_par_2
  params_flipped <- hawkes_par_1
  flip_time <- sample_decay_flip_time(treatment_time, end_time, vary_time = vary_time)
  flip_point <- sample_flip_point_in_cell(cell_win, flip_time)
  list(
    flip_cell = flip_cell,
    cell_win = cell_win,
    original_process = original_process,
    flipped_process = flipped_process,
    params_original = params_original,
    params_flipped = params_flipped,
    flip_point = flip_point
  )
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
  # Always flip treated -> control so Delta N = N_control - N_treated is an expected increase.
  original_process <- "treated"
  flipped_process <- "control"
  params_original <- hawkes_par_2
  params_flipped <- hawkes_par_1
  list(
    flip_cell = flip_cell,
    cell_win = cell_win,
    original_process = original_process,
    flipped_process = flipped_process,
    params_original = params_original,
    params_flipped = params_flipped
  )
}

simulate_single_point_catalogue <- function(
    params,
    flip_point,
    end_time,
    omega_win,
    clip_spatial = TRUE,
    unclipped_pad = 500) {
  t0 <- flip_point$t[1]
  if (!is.finite(end_time) || !is.finite(t0) || end_time <= t0) {
    return(data.frame(x = numeric(0), y = numeric(0), t = numeric(0)))
  }
  win <- omega_win
  if (!isTRUE(clip_spatial)) {
    # Retain long-range offspring for decay diagnostics: simulate on a large
    # window centred at the flip point (sim_hawkes_fast always clips to windowS).
    x0 <- as.numeric(flip_point$x[1])
    y0 <- as.numeric(flip_point$y[1])
    pad <- as.numeric(unclipped_pad)[1]
    if (!is.finite(pad) || pad <= 0) pad <- 500
    if (is.finite(x0) && is.finite(y0)) {
      win <- spatstat.geom::owin(c(x0 - pad, x0 + pad), c(y0 - pad, y0 + pad))
    }
  }
  ev <- getFromNamespace("sim_hawkes_fast", "PPDisentangle")(
    params = params,
    windowT = c(t0, end_time),
    windowS = win,
    background_realization = flip_point
  )
  data.frame(x = ev$x, y = ev$y, t = ev$t)
}

# First-generation offspring only (no cascade), with no spatial window clipping.
# This isolates the temporal/spatial triggering kernel shape that the cascade
# ΔN plots otherwise obscure (and avoids domain truncation of long-range PL draws).
simulate_direct_offspring_catalogue <- function(params, flip_point, end_time) {
  empty <- data.frame(x = numeric(0), y = numeric(0), t = numeric(0))
  t0 <- as.numeric(flip_point$t[1])
  x0 <- as.numeric(flip_point$x[1])
  y0 <- as.numeric(flip_point$y[1])
  horizon <- as.numeric(end_time) - t0
  if (!is.finite(horizon) || horizon <= 0 || !is.finite(x0) || !is.finite(y0)) return(empty)

  K <- as.numeric(params$K)[1]
  if (!is.finite(K) || K <= 0) return(empty)

  `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0L && !(length(a) == 1L && is.na(a))) a else b
  kernel <- normalize_hawkes_kernel(params$kernel %||% "exponential")
  normalize_spatial <- getFromNamespace("normalize_hawkes_spatial_kernel", "PPDisentangle")
  spatial_kernel <- normalize_spatial(params$spatial_kernel %||% "exponential")

  if (identical(kernel, "power_law")) {
    cc <- as.numeric((params[["c"]] %||% params$beta))[1]
    pp <- as.numeric(params[["p"]])[1]
    if (!is.finite(cc) || cc <= 0 || !is.finite(pp) || pp <= 1) return(empty)
    cdf_h <- 1 - (1 + horizon / cc)^(1 - pp)
    if (!is.finite(cdf_h) || cdf_h <= 0) return(empty)
    n <- stats::rpois(1L, K * cdf_h)
    if (n < 1L) return(empty)
    u <- stats::runif(n, 0, cdf_h)
    dt <- cc * ((1 - u)^(-1 / (pp - 1)) - 1)
  } else {
    beta <- as.numeric(params$beta)[1]
    if (!is.finite(beta) || beta <= 0) return(empty)
    cdf_h <- 1 - exp(-beta * horizon)
    if (!is.finite(cdf_h) || cdf_h <= 0) return(empty)
    n <- stats::rpois(1L, K * cdf_h)
    if (n < 1L) return(empty)
    u <- stats::runif(n, 0, cdf_h)
    dt <- -log(pmax(1e-15, 1 - u)) / beta
  }

  if (identical(spatial_kernel, "power_law")) {
    spatial_d <- as.numeric(params$spatial_d)[1]
    spatial_q <- as.numeric(params$spatial_q)[1]
    if (!is.finite(spatial_d) || spatial_d <= 0 || !is.finite(spatial_q) || spatial_q <= 1) {
      return(empty)
    }
    u_s <- stats::runif(n)
    u_s <- pmin(pmax(u_s, 0), 1 - 1e-15)
    r2 <- spatial_d * ((1 - u_s)^(-1 / (spatial_q - 1)) - 1)
  } else {
    alpha <- as.numeric(params$alpha)[1]
    if (!is.finite(alpha) || alpha <= 0) return(empty)
    # Match C++ R::rexp(1/alpha): Exp with mean 1/alpha for r^2.
    r2 <- stats::rexp(n, rate = alpha)
  }

  ang <- stats::runif(n, 0, 2 * pi)
  rr <- sqrt(pmax(r2, 0))
  data.frame(
    x = x0 + rr * cos(ang),
    y = y0 + rr * sin(ang),
    t = t0 + dt
  )
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

# Distance-to-cell annuli: events inside the flipped cell land in bin 0.
annulus_counts_from_cell <- function(catalogue, cell_win, annulus_width, max_bin) {
  if (is.null(catalogue) || nrow(catalogue) < 1L) {
    return(as.numeric(rep(0, max_bin + 1L)))
  }
  dist <- point_distance_to_rect(catalogue$x, catalogue$y, cell_win)
  bins <- pmin(max_bin, floor(dist / annulus_width))
  as.numeric(tabulate(bins + 1L, nbins = max_bin + 1L))
}

# Area of {k w <= d_rect < (k+1) w} intersected with the study window, via MC.
# Matches the distance metric used for binning and accounts for growing ring size
# and clipping at the window boundary.
estimate_rect_annulus_areas_in_window <- function(
    cell_win,
    omega_win,
    annulus_width,
    max_bin,
    n_mc = 200000L,
    seed = 314159L) {
  xr <- as.numeric(omega_win$xrange)
  yr <- as.numeric(omega_win$yrange)
  area_omega <- diff(xr) * diff(yr)
  if (!is.finite(area_omega) || area_omega <= 0) {
    stop("Invalid study window for annulus-area estimation.")
  }
  rng_state <- if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  } else {
    NULL
  }
  on.exit({
    if (is.null(rng_state)) {
      if (exists(".Random.seed", envir = globalenv())) {
        rm(".Random.seed", envir = globalenv())
      }
    } else {
      assign(".Random.seed", rng_state, envir = globalenv())
    }
  }, add = TRUE)
  set.seed(as.integer(seed)[1])
  x <- stats::runif(n_mc, xr[[1]], xr[[2]])
  y <- stats::runif(n_mc, yr[[1]], yr[[2]])
  dist <- point_distance_to_rect(x, y, cell_win)
  bins <- pmin(max_bin, floor(dist / annulus_width))
  props <- as.numeric(tabulate(bins + 1L, nbins = max_bin + 1L)) / as.numeric(n_mc)
  pmax(props * area_omega, .Machine$double.eps)
}

# Largest d_rect before the exterior buffer hits the study-window boundary.
# Beyond this, annulus ∩ Omega is truncated and area-normalized curves are unreliable.
max_unclipped_rect_distance <- function(cell_win, omega_win) {
  gaps <- c(
    as.numeric(cell_win$xrange[1] - omega_win$xrange[1]),
    as.numeric(omega_win$xrange[2] - cell_win$xrange[2]),
    as.numeric(cell_win$yrange[1] - omega_win$yrange[1]),
    as.numeric(omega_win$yrange[2] - cell_win$yrange[2])
  )
  max(0, min(gaps, na.rm = TRUE))
}

# Largest radial distance from a point before the circle hits Omega.
max_unclipped_point_distance <- function(x0, y0, omega_win) {
  gaps <- c(
    as.numeric(x0 - omega_win$xrange[1]),
    as.numeric(omega_win$xrange[2] - x0),
    as.numeric(y0 - omega_win$yrange[1]),
    as.numeric(omega_win$yrange[2] - y0)
  )
  max(0, min(gaps, na.rm = TRUE))
}

# Unconstrained (whole-plane) radial annulus areas around a point.
point_distance_annulus_areas_unconstrained <- function(annulus_width, max_bin) {
  w <- as.numeric(annulus_width)[1]
  d_left <- seq.int(0L, max_bin) * w
  d_right <- d_left + w
  pmax(pi * (d_right^2 - d_left^2), .Machine$double.eps)
}

# Unconstrained (whole-plane) Minkowski annulus areas for a rectangle.
rect_distance_annulus_areas_unconstrained <- function(cell_win, annulus_width, max_bin) {
  A <- diff(as.numeric(cell_win$xrange)) * diff(as.numeric(cell_win$yrange))
  P <- 2 * (diff(as.numeric(cell_win$xrange)) + diff(as.numeric(cell_win$yrange)))
  w <- as.numeric(annulus_width)[1]
  d_left <- seq.int(0L, max_bin) * w
  d_right <- d_left + w
  # F(r) = area({d_rect < r}); F(0) = 0 so bin 0 includes the cell interior.
  F <- function(r) {
    r <- pmax(as.numeric(r), 0)
    ifelse(r <= 0, 0, A + P * r + pi * r^2)
  }
  pmax(F(d_right) - F(d_left), .Machine$double.eps)
}

time_bin_counts_since_event <- function(catalogue, event_time, time_bin_width, max_bin) {
  counts <- rep(0, max_bin + 1L)
  if (is.null(catalogue) || nrow(catalogue) < 1L) return(counts)
  elapsed <- catalogue$t - event_time
  keep <- is.finite(elapsed) & elapsed >= 0
  if (!any(keep)) return(counts)
  bins <- pmin(max_bin, floor(elapsed[keep] / time_bin_width))
  as.numeric(tabulate(bins + 1L, nbins = max_bin + 1L))
}

# Forward-simulate z vs z with one treated cell flipped to control (CRN), matching
# Supplement Proposition (allocation-influence decay): influence of a single
# coordinate z_j measured by |Delta N| in distance/lag bins.
run_single_cell_allocation_flip_reps <- function(
    partition,
    partition_processes,
    omega,
    omega_win,
    treatment_time,
    end_time,
    hawkes_par_1,
    hawkes_par_2,
    decay_reps,
    count_fn,
    bins_df,
    stage_seed_fn,
    stage_offset,
    flip_cell = NA_integer_,
    log_fn = function(...) invisible(NULL)) {
  flip_cell <- choose_decay_flip_cell(
    partition, partition_processes, omega, flip_cell
  )
  partition_processes <- force_flip_cell_treated(partition_processes, flip_cell)
  cell_win <- as.owin(partition[flip_cell])
  z <- partition_processes == "treated"
  z_flip <- z
  z_flip[[flip_cell]] <- FALSE
  pp <- ifelse(z, "treated", "control")
  ppf <- ifelse(z_flip, "treated", "control")
  hp <- list(control = hawkes_par_1, treated = hawkes_par_2)
  sim_hawkes <- getFromNamespace("sim_hawkes", "PPDisentangle")
  gen <- PPDisentangle::generate_inhomogeneous_hawkes

  log_fn(
    "Running single-cell allocation-flip diagnostic: reps=", decay_reps,
    " | fixed_flip_cell=", flip_cell, " | treated->control"
  )
  t0 <- proc.time()[3]
  rep_rows <- lapply(seq_len(decay_reps), function(rep_id) {
    set.seed(stage_seed_fn(stage_offset, rep_id, 0L))
    pre <- as.data.frame(sim_hawkes(
      params = hawkes_par_1,
      windowT = c(0, treatment_time),
      windowS = omega_win
    ))
    if (nrow(pre) > 0L) {
      pre$process <- "control"
      pre$location_process <- "control"
    }
    filtration <- if (nrow(pre) > 0L) pre else NULL
    seed_post <- stage_seed_fn(stage_offset, rep_id, 1L)
    set.seed(seed_post)
    cat_z <- as.data.frame(gen(
      Omega = omega,
      partition = partition,
      time_window = c(treatment_time, end_time),
      partition_processes = pp,
      hawkes_params = hp,
      filtration = filtration
    ))
    set.seed(seed_post)
    cat_flip <- as.data.frame(gen(
      Omega = omega,
      partition = partition,
      time_window = c(treatment_time, end_time),
      partition_processes = ppf,
      hawkes_params = hp,
      filtration = filtration
    ))
    count_original <- count_fn(cat_z, cell_win)
    count_flipped <- count_fn(cat_flip, cell_win)
    delta_n <- count_flipped - count_original
    abs_delta_n <- abs(delta_n)
    pct_delta_n <- ifelse(
      count_original > 0,
      100 * delta_n / count_original,
      NA_real_
    )
    cbind(
      data.frame(
        rep = rep_id,
        flip_cell = flip_cell,
        flip_x = mean(cell_win$xrange),
        flip_y = mean(cell_win$yrange),
        flip_t = treatment_time,
        original_process = "treated",
        flipped_process = "control",
        stringsAsFactors = FALSE
      ),
      bins_df,
      data.frame(
        count_original = count_original,
        count_flipped = count_flipped,
        delta_n = delta_n,
        abs_delta_n = abs_delta_n,
        pct_delta_n = pct_delta_n,
        stringsAsFactors = FALSE
      )
    )
  })
  decay_df <- do.call(rbind, rep_rows)
  log_fn(sprintf(
    "Single-cell allocation-flip diagnostic finished in %.1fs",
    proc.time()[3] - t0
  ))
  list(
    df = decay_df,
    flip_cell = flip_cell,
    cell_center_x = mean(cell_win$xrange),
    cell_center_y = mean(cell_win$yrange),
    cell_win = cell_win
  )
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
    log_fn = function(...) invisible(NULL),
    partition = NULL,
    partition_processes = NULL,
    hawkes_par_1 = NULL,
    hawkes_par_2 = NULL,
    resample_flip_each_rep = TRUE,
    vary_flip_time = FALSE,
    first_generation_only = FALSE,
    clip_spatial = FALSE) {
  # Legacy single-cell mode keeps flip_ctx fixed; preferred mode resamples a
  # treated cell and point (and optionally flip time) independently each rep.
  # Default: full Hawkes cascade. clip_spatial=FALSE expands the sim window so
  # long-range offspring are retained for decay diagnostics.
  can_resample <- isTRUE(resample_flip_each_rep) &&
    !is.null(partition) && !is.null(partition_processes) &&
    !is.null(hawkes_par_1) && !is.null(hawkes_par_2)
  cell_win <- flip_ctx$cell_win
  log_fn(
    "Running single-point label-flip diagnostic: reps=", decay_reps,
    " | resample_each_rep=", can_resample,
    " | vary_flip_time=", isTRUE(vary_flip_time),
    " | first_generation_only=", isTRUE(first_generation_only),
    " | clip_spatial=", isTRUE(clip_spatial),
    if (can_resample) "" else paste0(" | flip_cell=", flip_ctx$flip_cell),
    if (can_resample) "" else paste0(" | ", flip_ctx$original_process, "->", flip_ctx$flipped_process)
  )
  t0 <- proc.time()[3]
  rep_rows <- lapply(seq_len(decay_reps), function(rep_id) {
    set.seed(stage_seed_fn(stage_offset, rep_id, 0L))
    if (can_resample) {
      draw <- sample_decay_flip_context(
        partition = partition,
        partition_processes = partition_processes,
        hawkes_par_1 = hawkes_par_1,
        hawkes_par_2 = hawkes_par_2,
        treatment_time = treatment_time,
        end_time = end_time,
        omega = omega,
        vary_time = vary_flip_time
      )
      flip_point <- draw$flip_point
      params_original <- draw$params_original
      params_flipped <- draw$params_flipped
      flip_cell <- draw$flip_cell
      original_process <- draw$original_process
      flipped_process <- draw$flipped_process
    } else {
      flip_point <- sample_flip_point_in_cell(cell_win, treatment_time)
      params_original <- flip_ctx$params_original
      params_flipped <- flip_ctx$params_flipped
      flip_cell <- flip_ctx$flip_cell
      original_process <- flip_ctx$original_process
      flipped_process <- flip_ctx$flipped_process
    }
    set.seed(stage_seed_fn(stage_offset, rep_id, 1L))
    if (isTRUE(first_generation_only)) {
      cat_original <- simulate_direct_offspring_catalogue(
        params = params_original,
        flip_point = flip_point,
        end_time = end_time
      )
      set.seed(stage_seed_fn(stage_offset, rep_id, 1L))
      cat_flipped <- simulate_direct_offspring_catalogue(
        params = params_flipped,
        flip_point = flip_point,
        end_time = end_time
      )
    } else {
      cat_original <- simulate_single_point_catalogue(
        params = params_original,
        flip_point = flip_point,
        end_time = end_time,
        omega_win = omega_win,
        clip_spatial = clip_spatial
      )
      set.seed(stage_seed_fn(stage_offset, rep_id, 1L))
      cat_flipped <- simulate_single_point_catalogue(
        params = params_flipped,
        flip_point = flip_point,
        end_time = end_time,
        omega_win = omega_win,
        clip_spatial = clip_spatial
      )
    }
    count_original <- count_fn(cat_original, flip_point)
    count_flipped <- count_fn(cat_flipped, flip_point)
    # Signed: flipped (control) minus original (treated); expected positive under K_0 > K_1.
    delta_n <- count_flipped - count_original
    abs_delta_n <- abs(delta_n)
    pct_delta_n <- ifelse(
      count_original > 0,
      100 * delta_n / count_original,
      NA_real_
    )
    cbind(
      data.frame(
        rep = rep_id,
        flip_cell = flip_cell,
        flip_x = flip_point$x[1],
        flip_y = flip_point$y[1],
        flip_t = flip_point$t[1],
        original_process = original_process,
        flipped_process = flipped_process,
        stringsAsFactors = FALSE
      ),
      bins_df,
      data.frame(
        count_original = count_original,
        count_flipped = count_flipped,
        delta_n = delta_n,
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
    cell_center_x = mean(cell_win$xrange),
    cell_center_y = mean(cell_win$yrange)
  )
}

summarize_decay_bins <- function(decay_df, bin_cols) {
  decay_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("annulus", bin_cols)))) %>%
    dplyr::summarize(
      mean_delta = mean(.data$delta_n),
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
          100 * .data$mean_delta / .data$mean_count_original,
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
  # Cover the full window diameter so far-field annuli remain defined.
  max_corner_dist <- sqrt((omega[2] - omega[1])^2 + (omega[4] - omega[3])^2)
  max_bin <- max(1L, ceiling(max_corner_dist / annulus_width))
  # Resolve the fixed flip cell early so annulus areas match that cell.
  flip_cell_fixed <- choose_decay_flip_cell(
    partition, partition_processes, omega, flip_cell
  )
  cell_win_areas <- as.owin(partition[flip_cell_fixed])
  annulus_areas <- estimate_rect_annulus_areas_in_window(
    cell_win = cell_win_areas,
    omega_win = omega_win,
    annulus_width = annulus_width,
    max_bin = max_bin
  )
  area_free <- rect_distance_annulus_areas_unconstrained(
    cell_win_areas, annulus_width, max_bin
  )
  max_unclipped <- max_unclipped_rect_distance(cell_win_areas, omega_win)
  annuli <- data.frame(
    annulus = seq.int(0L, max_bin),
    d_left = seq.int(0L, max_bin) * annulus_width,
    d_mid = (seq.int(0L, max_bin) + 0.5) * annulus_width,
    annulus_area = annulus_areas,
    annulus_area_unconstrained = area_free,
    area_keep_frac = annulus_areas / area_free,
    max_unclipped_distance = max_unclipped
  )
  count_fn <- function(catalogue, cell_win) {
    counts <- annulus_counts_from_cell(
      catalogue = catalogue,
      cell_win = cell_win,
      annulus_width = annulus_width,
      max_bin = max_bin
    )
    # Density: count per unit area in the distance ring ∩ study window.
    counts / annulus_areas
  }
  out <- run_single_cell_allocation_flip_reps(
    partition = partition,
    partition_processes = partition_processes,
    omega = omega,
    omega_win = omega_win,
    treatment_time = treatment_time,
    end_time = end_time,
    hawkes_par_1 = hawkes_par_1,
    hawkes_par_2 = hawkes_par_2,
    decay_reps = decay_reps,
    count_fn = count_fn,
    bins_df = annuli,
    stage_seed_fn = stage_seed_fn,
    stage_offset = 13L,
    flip_cell = flip_cell_fixed,
    log_fn = log_fn
  )
  decay_df <- out$df
  decay_summary <- summarize_decay_bins(
    decay_df,
    c("d_left", "d_mid", "annulus_area", "annulus_area_unconstrained",
      "area_keep_frac", "max_unclipped_distance")
  )
  # Keep only annuli whose exterior buffer has not yet hit the study-window
  # boundary (spatial clipping). Beyond this, Ω-truncation dominates.
  keep <- is.finite(decay_summary$d_mid) &
    decay_summary$d_mid <= decay_summary$max_unclipped_distance &
    (!("area_keep_frac" %in% names(decay_summary)) |
       is.na(decay_summary$area_keep_frac) |
       decay_summary$area_keep_frac >= 0.95)
  decay_summary <- decay_summary[keep, , drop = FALSE]
  rate_lines <- make_decay_rate_lines(decay_summary, hawkes_par_1, hawkes_par_2)
  pos <- decay_summary$mean_abs_delta[
    is.finite(decay_summary$mean_abs_delta) & decay_summary$mean_abs_delta > 0
  ]
  eps <- if (length(pos) > 0L) min(pos, na.rm = TRUE) / 2 else 1e-6
  if (!is.finite(eps)) eps <- 1e-6
  decay_summary$mean_delta_plot <- pmax(decay_summary$mean_abs_delta, eps)
  decay_summary$mean_abs_delta_plot <- decay_summary$mean_delta_plot
  if (!is.null(rate_lines) && nrow(rate_lines) > 0) {
    rate_lines$rate_value_plot <- pmax(rate_lines$rate_value, eps)
  }
  list(
    specs = list(
      enabled = TRUE,
      flip_mode = "single_cell_allocation_fixed_central_treated_to_control",
      flip_direction = "treated_to_control",
      reps = decay_reps,
      annulus_width = annulus_width,
      normalize_by_annulus_area = TRUE,
      spatial_clip_to_untruncated_annuli = TRUE,
      max_unclipped_distance = max_unclipped,
      flip_cell = out$flip_cell,
      flip_cell_fixed = TRUE,
      flip_time = treatment_time,
      original_process = "treated",
      flipped_process = "control",
      slope_reference = "area-normalized |Delta N| per untruncated annulus under CRN allocation flip"
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
  post_treatment_horizon <- end_time - treatment_time
  max_bin <- max(1L, ceiling(post_treatment_horizon / time_bin_width))
  time_bins <- data.frame(
    annulus = seq.int(0L, max_bin),
    t_left = seq.int(0L, max_bin) * time_bin_width,
    t_mid = (seq.int(0L, max_bin) + 0.5) * time_bin_width
  )
  count_fn <- function(catalogue, cell_win) {
    time_bin_counts_since_event(
      catalogue = catalogue,
      event_time = treatment_time,
      time_bin_width = time_bin_width,
      max_bin = max_bin
    )
  }
  out <- run_single_cell_allocation_flip_reps(
    partition = partition,
    partition_processes = partition_processes,
    omega = omega,
    omega_win = omega_win,
    treatment_time = treatment_time,
    end_time = end_time,
    hawkes_par_1 = hawkes_par_1,
    hawkes_par_2 = hawkes_par_2,
    decay_reps = decay_reps,
    count_fn = count_fn,
    bins_df = time_bins,
    stage_seed_fn = stage_seed_fn,
    stage_offset = 14L,
    flip_cell = flip_cell,
    log_fn = log_fn
  )
  decay_df <- out$df
  decay_summary <- summarize_decay_bins(decay_df, c("t_left", "t_mid"))
  pos <- decay_summary$mean_abs_delta[
    is.finite(decay_summary$mean_abs_delta) & decay_summary$mean_abs_delta > 0
  ]
  eps <- if (length(pos) > 0L) min(pos, na.rm = TRUE) / 2 else 1e-6
  if (!is.finite(eps)) eps <- 1e-6
  decay_summary$mean_delta_plot <- pmax(decay_summary$mean_abs_delta, eps)
  decay_summary$mean_abs_delta_plot <- decay_summary$mean_delta_plot
  list(
    specs = list(
      enabled = TRUE,
      flip_mode = "single_cell_allocation_fixed_central_treated_to_control",
      flip_direction = "treated_to_control",
      axis = "lag_since_treatment_time",
      reps = decay_reps,
      time_bin_width = time_bin_width,
      flip_cell = out$flip_cell,
      flip_cell_fixed = TRUE,
      flip_time = treatment_time,
      original_process = "treated",
      flipped_process = "control"
    ),
    df = decay_df,
    summary = decay_summary
  )
}

run_label_flip_decay_validation_scenario <- function(
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
  max_corner_dist <- sqrt((omega[2] - omega[1])^2 + (omega[4] - omega[3])^2)
  max_bin <- max(1L, ceiling(max_corner_dist / annulus_width))
  flip_ctx <- resolve_single_point_flip_context(
    partition = partition,
    partition_processes = partition_processes,
    omega = omega,
    hawkes_par_1 = hawkes_par_1,
    hawkes_par_2 = hawkes_par_2,
    flip_cell = flip_cell
  )
  # Clip using the fixed cell centre (flip points are drawn inside this cell).
  x0 <- mean(as.numeric(flip_ctx$cell_win$xrange))
  y0 <- mean(as.numeric(flip_ctx$cell_win$yrange))
  max_unclipped <- max_unclipped_point_distance(x0, y0, omega_win)
  area_free <- point_distance_annulus_areas_unconstrained(annulus_width, max_bin)
  # Approximate in-window ring areas by MC distance-to-point from the cell centre.
  annulus_areas <- local({
    xr <- as.numeric(omega_win$xrange)
    yr <- as.numeric(omega_win$yrange)
    area_omega <- diff(xr) * diff(yr)
    rng_state <- if (exists(".Random.seed", envir = globalenv())) {
      get(".Random.seed", envir = globalenv())
    } else {
      NULL
    }
    on.exit({
      if (is.null(rng_state)) {
        if (exists(".Random.seed", envir = globalenv())) {
          rm(".Random.seed", envir = globalenv())
        }
      } else {
        assign(".Random.seed", rng_state, envir = globalenv())
      }
    }, add = TRUE)
    set.seed(314159L)
    n_mc <- 200000L
    xs <- stats::runif(n_mc, xr[[1]], xr[[2]])
    ys <- stats::runif(n_mc, yr[[1]], yr[[2]])
    dist <- point_distance_to_point(xs, ys, x0, y0)
    bins <- pmin(max_bin, floor(dist / annulus_width))
    props <- as.numeric(tabulate(bins + 1L, nbins = max_bin + 1L)) / as.numeric(n_mc)
    pmax(props * area_omega, .Machine$double.eps)
  })
  annuli <- data.frame(
    annulus = seq.int(0L, max_bin),
    d_left = seq.int(0L, max_bin) * annulus_width,
    d_mid = (seq.int(0L, max_bin) + 0.5) * annulus_width,
    annulus_area = annulus_areas,
    annulus_area_unconstrained = area_free,
    area_keep_frac = annulus_areas / area_free,
    max_unclipped_distance = max_unclipped
  )
  count_fn <- function(catalogue, flip_point) {
    counts <- annulus_counts_from_point(
      catalogue = catalogue,
      x0 = flip_point$x[1],
      y0 = flip_point$y[1],
      annulus_width = annulus_width,
      max_bin = max_bin
    )
    counts / annulus_areas
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
    stage_offset = 15L,
    log_fn = log_fn,
    partition = partition,
    partition_processes = force_flip_cell_treated(partition_processes, flip_ctx$flip_cell),
    hawkes_par_1 = hawkes_par_1,
    hawkes_par_2 = hawkes_par_2,
    resample_flip_each_rep = FALSE,
    vary_flip_time = FALSE,
    first_generation_only = FALSE,
    clip_spatial = TRUE
  )
  decay_df <- out$df
  decay_summary <- summarize_decay_bins(
    decay_df,
    c("d_left", "d_mid", "annulus_area", "annulus_area_unconstrained",
      "area_keep_frac", "max_unclipped_distance")
  )
  keep <- is.finite(decay_summary$d_mid) &
    decay_summary$d_mid <= decay_summary$max_unclipped_distance &
    (!("area_keep_frac" %in% names(decay_summary)) |
       is.na(decay_summary$area_keep_frac) |
       decay_summary$area_keep_frac >= 0.95)
  decay_summary <- decay_summary[keep, , drop = FALSE]
  pos <- decay_summary$mean_abs_delta[
    is.finite(decay_summary$mean_abs_delta) & decay_summary$mean_abs_delta > 0
  ]
  eps <- if (length(pos) > 0L) min(pos, na.rm = TRUE) / 2 else 1e-6
  if (!is.finite(eps)) eps <- 1e-6
  decay_summary$mean_delta_plot <- pmax(decay_summary$mean_abs_delta, eps)
  decay_summary$mean_abs_delta_plot <- decay_summary$mean_delta_plot
  list(
    specs = list(
      enabled = TRUE,
      flip_mode = "single_point_label_fixed_central_treated_to_control",
      flip_direction = "treated_to_control",
      reps = decay_reps,
      annulus_width = annulus_width,
      normalize_by_annulus_area = TRUE,
      spatial_clip_to_untruncated_annuli = TRUE,
      max_unclipped_distance = max_unclipped,
      flip_cell = flip_ctx$flip_cell,
      flip_cell_fixed = TRUE,
      flip_time = treatment_time,
      original_process = "treated",
      flipped_process = "control"
    ),
    df = decay_df,
    summary = decay_summary
  )
}

run_label_flip_temporal_decay_validation_scenario <- function(
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
  post_treatment_horizon <- end_time - treatment_time
  max_bin <- max(1L, ceiling(post_treatment_horizon / time_bin_width))
  time_bins <- data.frame(
    annulus = seq.int(0L, max_bin),
    t_left = seq.int(0L, max_bin) * time_bin_width,
    t_mid = (seq.int(0L, max_bin) + 0.5) * time_bin_width
  )
  flip_ctx <- resolve_single_point_flip_context(
    partition = partition,
    partition_processes = partition_processes,
    omega = omega,
    hawkes_par_1 = hawkes_par_1,
    hawkes_par_2 = hawkes_par_2,
    flip_cell = flip_cell
  )
  count_fn <- function(catalogue, flip_point) {
    time_bin_counts_since_event(
      catalogue = catalogue,
      event_time = flip_point$t[1],
      time_bin_width = time_bin_width,
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
    bins_df = time_bins,
    stage_seed_fn = stage_seed_fn,
    stage_offset = 16L,
    log_fn = log_fn,
    partition = partition,
    partition_processes = force_flip_cell_treated(partition_processes, flip_ctx$flip_cell),
    hawkes_par_1 = hawkes_par_1,
    hawkes_par_2 = hawkes_par_2,
    resample_flip_each_rep = FALSE,
    vary_flip_time = FALSE,
    first_generation_only = FALSE,
    clip_spatial = TRUE
  )
  decay_df <- out$df
  decay_summary <- summarize_decay_bins(decay_df, c("t_left", "t_mid"))
  pos <- decay_summary$mean_abs_delta[
    is.finite(decay_summary$mean_abs_delta) & decay_summary$mean_abs_delta > 0
  ]
  eps <- if (length(pos) > 0L) min(pos, na.rm = TRUE) / 2 else 1e-6
  if (!is.finite(eps)) eps <- 1e-6
  decay_summary$mean_delta_plot <- pmax(decay_summary$mean_abs_delta, eps)
  decay_summary$mean_abs_delta_plot <- decay_summary$mean_delta_plot
  list(
    specs = list(
      enabled = TRUE,
      flip_mode = "single_point_label_fixed_central_treated_to_control",
      flip_direction = "treated_to_control",
      axis = "lag_since_flip_time",
      reps = decay_reps,
      time_bin_width = time_bin_width,
      flip_cell = flip_ctx$flip_cell,
      flip_cell_fixed = TRUE,
      flip_time = treatment_time,
      original_process = "treated",
      flipped_process = "control"
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
  # Prefer an explicit config cell; otherwise always the fixed geometric centre cell.
  flip_cell <- if (!is.null(cfg$DECAY_FLIP_CELL) &&
                   is.finite(suppressWarnings(as.integer(cfg$DECAY_FLIP_CELL)[1]))) {
    as.integer(cfg$DECAY_FLIP_CELL)[1]
  } else {
    NA_integer_
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
    base_seed = 123L,
    time_bin_width = NA_real_) {
  partition <- spatstat.geom::quadrats(X = cfg$OMEGA, nx = cfg$NX, ny = cfg$NY)
  seed_base <- if (!is.null(cfg$BASE_SEED)) cfg$BASE_SEED else base_seed
  partition_processes <- rebuild_partition_processes(partition$n, cfg$NX, cfg$NY, base_seed = seed_base)
  # Prefer an explicit config cell; otherwise always the fixed geometric centre cell.
  flip_cell <- if (!is.null(cfg$DECAY_FLIP_CELL) &&
                   is.finite(suppressWarnings(as.integer(cfg$DECAY_FLIP_CELL)[1]))) {
    as.integer(cfg$DECAY_FLIP_CELL)[1]
  } else {
    NA_integer_
  }
  stage_seed_fn <- function(stage_offset, sim_id = 0L, extra = 0L) {
    as.integer(seed_base + as.integer(stage_offset) * 100000L + as.integer(sim_id) * 1000L + as.integer(extra))
  }
  if (!is.finite(time_bin_width) || time_bin_width <= 0) {
    time_bin_width <- if (!is.null(cfg$DECAY_TIME_BIN_WIDTH)) cfg$DECAY_TIME_BIN_WIDTH else NA_real_
  }
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

run_label_flip_decay_from_cfg <- function(
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
  flip_cell <- if (!is.null(cfg$DECAY_FLIP_CELL) &&
                   is.finite(suppressWarnings(as.integer(cfg$DECAY_FLIP_CELL)[1]))) {
    as.integer(cfg$DECAY_FLIP_CELL)[1]
  } else {
    NA_integer_
  }
  stage_seed_fn <- function(stage_offset, sim_id = 0L, extra = 0L) {
    as.integer(seed_base + as.integer(stage_offset) * 100000L + as.integer(sim_id) * 1000L + as.integer(extra))
  }
  run_label_flip_decay_validation_scenario(
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

run_label_flip_temporal_decay_from_cfg <- function(
    cfg,
    res,
    hawkes_par_1,
    hawkes_par_2,
    sim_kernel,
    decay_reps,
    base_seed = 123L,
    time_bin_width = NA_real_) {
  partition <- spatstat.geom::quadrats(X = cfg$OMEGA, nx = cfg$NX, ny = cfg$NY)
  seed_base <- if (!is.null(cfg$BASE_SEED)) cfg$BASE_SEED else base_seed
  partition_processes <- rebuild_partition_processes(partition$n, cfg$NX, cfg$NY, base_seed = seed_base)
  flip_cell <- if (!is.null(cfg$DECAY_FLIP_CELL) &&
                   is.finite(suppressWarnings(as.integer(cfg$DECAY_FLIP_CELL)[1]))) {
    as.integer(cfg$DECAY_FLIP_CELL)[1]
  } else {
    NA_integer_
  }
  stage_seed_fn <- function(stage_offset, sim_id = 0L, extra = 0L) {
    as.integer(seed_base + as.integer(stage_offset) * 100000L + as.integer(sim_id) * 1000L + as.integer(extra))
  }
  if (!is.finite(time_bin_width) || time_bin_width <= 0) {
    time_bin_width <- if (!is.null(cfg$DECAY_TIME_BIN_WIDTH)) cfg$DECAY_TIME_BIN_WIDTH else NA_real_
  }
  run_label_flip_temporal_decay_validation_scenario(
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
    base_seed = 123L,
    time_bin_width = 0.2,
    exponential_beta = NULL) {
  if (!file.exists(rds_path)) stop(sprintf("RDS not found: %s", rds_path))
  res <- readRDS(rds_path)
  cfg <- res$config
  if (!is.null(exponential_beta) && is.finite(exponential_beta) && exponential_beta > 0) {
    cfg$HAWKES_BETA <- as.numeric(exponential_beta)[1]
  }
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
    base_seed = base_seed,
    time_bin_width = time_bin_width
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

refresh_label_flip_decay_for_spec <- function(
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
  out <- run_label_flip_decay_from_cfg(
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

refresh_label_flip_temporal_decay_for_spec <- function(
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
    base_seed = 123L,
    time_bin_width = 0.2,
    exponential_beta = NULL) {
  if (!file.exists(rds_path)) stop(sprintf("RDS not found: %s", rds_path))
  res <- readRDS(rds_path)
  cfg <- res$config
  if (!is.null(exponential_beta) && is.finite(exponential_beta) && exponential_beta > 0) {
    cfg$HAWKES_BETA <- as.numeric(exponential_beta)[1]
  }
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
  out <- run_label_flip_temporal_decay_from_cfg(
    cfg = cfg,
    res = res,
    hawkes_par_1 = hawkes_par_1,
    hawkes_par_2 = hawkes_par_2,
    sim_kernel = sim_kernel,
    decay_reps = decay_reps,
    base_seed = base_seed,
    time_bin_width = time_bin_width
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
