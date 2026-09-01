# Support routines for oklahoma_report.qmd when results RDS lacks pp_data / counties
# (e.g. older pre-bootstrap checkpoints). Mirrors oklahoma_analysis.R Steps 1–2 + pre split.

local({
  # Source shared geometry helpers from the same directory when available.
  this_file <- tryCatch(
    normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE),
    error = function(e) NA_character_
  )
  candidates <- c(
    if (is.character(this_file) && nzchar(this_file)) {
      file.path(dirname(this_file), "oklahoma_geometry.R")
    } else {
      character(0)
    },
    file.path(getwd(), "inst/oklahoma/oklahoma_geometry.R"),
    file.path(getwd(), "oklahoma_geometry.R")
  )
  if (!exists("oklahoma_sf_to_owin_km", mode = "function", inherits = TRUE)) {
    found <- candidates[file.exists(candidates)]
    if (!length(found)) {
      stop("Could not locate oklahoma_geometry.R; source it before oklahoma_report_support.R.")
    }
    source(found[[1]], local = FALSE)
  }
})

oklahoma_report_rebuild_pp_and_counties <- function(data_dir, etas_m0, crs_proj = 5070L) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table required")
  if (!requireNamespace("sf", quietly = TRUE)) stop("sf required")
  if (!requireNamespace("tigris", quietly = TRUE)) stop("tigris required")
  if (!requireNamespace("spatstat.geom", quietly = TRUE)) stop("spatstat required")
  if (!exists("oklahoma_sf_to_owin_km", mode = "function", inherits = TRUE)) {
    stop("oklahoma_geometry.R helpers are not loaded.")
  }

  meta <- jsonlite::fromJSON(readLines(file.path(data_dir, "metadata.json")))
  ev_all <- data.table::fread(file.path(data_dir, "events_all.csv"))
  t_star_utc <- as.POSIXct(meta$design$t_star_utc, tz = "UTC")
  ev_all[, time_utc := as.POSIXct(time_utc, tz = "UTC")]
  ev_all[, t_days := as.numeric(difftime(time_utc, t_star_utc, units = "days"))]
  ev_all[, x_km := x_m / 1000]
  ev_all[, y_km := y_m / 1000]
  ev_all <- ev_all[mag >= etas_m0]

  post_end_days <- as.numeric(difftime(
    as.POSIXct(meta$design$post_end_utc, tz = "UTC"), t_star_utc, units = "days"))

  options(tigris_use_cache = TRUE)
  counties_sf <- tigris::counties(state = "OK", cb = TRUE, year = 2022)
  counties_sf <- sf::st_transform(counties_sf, crs_proj)
  counties_sf <- sf::st_make_valid(counties_sf)

  ok_boundary <- sf::st_union(counties_sf)
  ok_boundary <- sf::st_make_valid(ok_boundary)
  win_km <- oklahoma_sf_to_owin_km(ok_boundary)

  county_owins <- oklahoma_sf_features_to_owins_km(counties_sf, name_col = "NAME")
  valid_idx <- !vapply(county_owins, is.null, logical(1))
  county_owins_valid <- county_owins[valid_idx]
  counties_sf_valid <- counties_sf[valid_idx, ]
  partition <- spatstat.geom::tess(tiles = county_owins_valid, window = win_km)

  aoi_sf <- sf::st_read(file.path(data_dir, "occ_aoi_layer_2.geojson"), quiet = TRUE)
  aoi_sf <- sf::st_transform(aoi_sf, crs_proj)
  aoi_sf <- sf::st_make_valid(aoi_sf)
  aoi_union <- sf::st_union(aoi_sf)
  county_centroids <- sf::st_centroid(counties_sf_valid)
  inside_aoi <- lengths(sf::st_within(county_centroids, aoi_union)) > 0
  partition_processes <- ifelse(inside_aoi, "treated", "control")
  names(partition_processes) <- counties_sf_valid$NAME
  treated_idx <- partition_processes == "treated"
  treated_names <- names(partition_processes)[treated_idx]
  control_ss <- spatstat.geom::as.owin(partition[!treated_idx])
  treated_ss <- spatstat.geom::as.owin(partition[treated_idx])

  assign_county <- function(df) {
    df$location_process <- oklahoma_assign_partition_process(
      df$x, df$y, partition, partition_processes
    )
    df$W <- 1.0
    df$n <- nrow(df)
    df$background <- TRUE
    df
  }

  pp_pre <- assign_county(as.data.frame(ev_all[t_days < 0,
    list(x = x_km, y = y_km, t = t_days, mag = mag)]))
  pp_post <- assign_county(as.data.frame(ev_all[t_days >= 0 & t_days <= post_end_days,
    list(x = x_km, y = y_km, t = t_days, mag = mag)]))
  pp_pre <- pp_pre[!is.na(pp_pre$location_process), , drop = FALSE]
  pp_post <- pp_post[!is.na(pp_post$location_process), , drop = FALSE]
  oklahoma_assert_label_support(pp_pre, control_ss, treated_ss, context = "pp_pre")
  oklahoma_assert_label_support(pp_post, control_ss, treated_ss, context = "pp_post")

  pp_pre_all <- pp_pre[order(pp_pre$t), , drop = FALSE]
  n_pre_total <- nrow(pp_pre_all)
  n_pre_holdout <- floor(n_pre_total * 0.5)
  if (n_pre_total > 0 && n_pre_holdout < 1) n_pre_holdout <- 1
  holdout_idx <- if (n_pre_holdout > 0) seq_len(n_pre_holdout) else integer(0)
  keep_idx <- if (n_pre_total > n_pre_holdout) (n_pre_holdout + 1):n_pre_total else integer(0)
  if (length(keep_idx) < 1) {
    stop("Pre-treatment split left zero estimation events; cannot rebuild pp_data for report.")
  }
  pp_pre_holdout <- pp_pre_all[holdout_idx, , drop = FALSE]
  pp_pre <- pp_pre_all[keep_idx, , drop = FALSE]

  pp_pre$process <- "control"
  pp_post$process <- pp_post$location_process
  pp_pre$inferred_process <- "control"
  pp_post$inferred_process <- pp_post$location_process

  counties <- list(
    names = counties_sf_valid$NAME,
    treated_names = treated_names,
    n_counties = partition$n,
    n_treated = sum(treated_idx)
  )
  pp_data <- list(
    pp_pre = pp_pre,
    pp_pre_holdout = pp_pre_holdout,
    pp_post = pp_post
  )
  list(counties = counties, pp_data = pp_data, meta = meta, post_end_days = post_end_days)
}

oklahoma_report_grid_partition <- function(diameter, win, aoi_owin, label, max_tiles = 5000L) {
  if (is.null(diameter) || length(diameter) < 1L || !is.finite(diameter[[1]])) {
    target_tiles <- 100L
    aspect <- diff(win$xrange) / diff(win$yrange)
    nx <- max(2L, round(sqrt(target_tiles * aspect)))
    ny <- max(2L, round(target_tiles / nx))
  } else {
    nx <- max(2L, ceiling(diff(win$xrange) / diameter))
    ny <- max(2L, ceiling(diff(win$yrange) / diameter))
    if ((nx * ny) > max_tiles) {
      shrink <- sqrt((nx * ny) / max_tiles)
      nx <- max(2L, ceiling(nx / shrink))
      ny <- max(2L, ceiling(ny / shrink))
    }
  }
  grid_part <- spatstat.geom::quadrats(win, nx = nx, ny = ny)
  n_tiles <- grid_part$n
  grid_procs <- character(n_tiles)
  tile_list <- spatstat.geom::tiles(grid_part)
  for (i in seq_len(n_tiles)) {
    tile_cent <- c(mean(tile_list[[i]]$xrange), mean(tile_list[[i]]$yrange))
    grid_procs[i] <- if (spatstat.geom::inside.owin(tile_cent[1], tile_cent[2], aoi_owin)) {
      "treated"
    } else {
      "control"
    }
  }
  names(grid_procs) <- spatstat.geom::tilenames(grid_part)
  list(partition = grid_part, processes = grid_procs, label = label)
}

oklahoma_report_tess_rects <- function(tess_obj, processes) {
  tile_list <- spatstat.geom::tiles(tess_obj)
  nms <- names(tile_list)
  do.call(rbind, lapply(seq_along(tile_list), function(i) {
    ow <- tile_list[[i]]
    nm <- if (!is.null(nms) && length(nms) >= i) nms[[i]] else as.character(i)
    data.frame(
      xmin = ow$xrange[1], xmax = ow$xrange[2],
      ymin = ow$yrange[1], ymax = ow$yrange[2],
      treatment = unname(processes[[nm]]),
      stringsAsFactors = FALSE
    )
  }))
}

oklahoma_report_plot_partition_pp <- function(base_map, pp_pre, pp_post) {
  p1 <- base_map +
    ggplot2::geom_point(
      data = pp_pre, ggplot2::aes(x = x, y = y, alpha = t),
      colour = "#2166ac", size = 0.8
    ) +
    ggplot2::scale_alpha_continuous(range = c(0.15, 0.8), guide = "none") +
    ggplot2::labs(title = "Pre-treatment", x = "X (km)", y = "Y (km)")
  p2 <- base_map +
    ggplot2::geom_point(
      data = pp_post,
      ggplot2::aes(x = x, y = y, colour = Process, alpha = t),
      size = 0.8
    ) +
    ggplot2::scale_colour_manual(
      values = c(control = "#2166ac", treated = "#b2182b"),
      guide = "none"
    ) +
    ggplot2::scale_alpha_continuous(range = c(0.3, 1), guide = "none") +
    ggplot2::labs(title = "Post-treatment", x = "X (km)", y = "Y (km)")
  gridExtra::grid.arrange(p1, p2, ncol = 2)
}

# ---------------------------------------------------------------------------
# SEM D log-likelihood on a common AOI geometry
#
# Last-labelling vs AOI-boundary process labels are scored with the same
# fitted SEM-D parameters, AOI control/treated regions, and KDE field.
# Only inferred_process changes. Pre-treatment events stay control.
# ---------------------------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a)) a else b

oklahoma_report_cfg_scalar <- function(cfg, name, default) {
  if (is.null(cfg) || is.null(cfg[[name]])) return(default)
  v <- suppressWarnings(as.numeric(unlist(cfg[[name]], use.names = FALSE)))
  v <- v[is.finite(v)]
  if (length(v) < 1L) return(default)
  v[[1]]
}

oklahoma_report_event_key <- function(df) {
  paste(
    sprintf("%.8f", as.numeric(df$t)),
    sprintf("%.6f", as.numeric(df$x)),
    sprintf("%.6f", as.numeric(df$y)),
    sep = "|"
  )
}

oklahoma_report_aoi_state_spaces <- function(win_km, aoi_owin) {
  treat <- tryCatch(
    spatstat.geom::intersect.owin(aoi_owin, win_km),
    error = function(e) aoi_owin
  )
  ctrl <- tryCatch(
    spatstat.geom::setminus.owin(win_km, treat),
    error = function(e) NULL
  )
  if (is.null(ctrl) || is.null(treat)) return(NULL)
  list(control = ctrl, treated = treat)
}

oklahoma_report_normalize_bg_weights <- function(df_sub, win_sub, covariate_im,
                                                 mark_name = "W") {
  if (is.null(df_sub) || nrow(df_sub) == 0) {
    return(list(new_df = df_sub, norm = 0))
  }
  cov_in_window <- covariate_im[win_sub, drop = FALSE]
  total_mass_raw <- spatstat.geom::integral.im(cov_in_window)
  target_area <- spatstat.geom::area(win_sub)
  if (!is.finite(total_mass_raw) || total_mass_raw <= 0 ||
      !is.finite(target_area) || target_area <= 0) {
    stop("Cannot normalize KDE background weights on a state space.")
  }
  norm_factor <- target_area / total_mass_raw
  vals_raw <- spatstat.geom::interp.im(covariate_im, df_sub$x, df_sub$y)
  vals_raw[!is.finite(vals_raw)] <- 0
  df_sub[[mark_name]] <- vals_raw * norm_factor
  pos <- df_sub[[mark_name]][df_sub[[mark_name]] > 0]
  min_val <- if (length(pos) > 0L) min(pos, na.rm = TRUE) else 1e-9
  if (!is.finite(min_val) || min_val <= 0) min_val <- 1e-9
  df_sub[[mark_name]][df_sub[[mark_name]] <= 0] <- min_val
  list(new_df = df_sub, norm = norm_factor)
}

oklahoma_report_kde_pre_mass_ratio <- function(covariate_im, ctrl_win) {
  total_mass <- spatstat.geom::integral.im(covariate_im)
  ctrl_mass <- spatstat.geom::integral.im(covariate_im[ctrl_win, drop = FALSE])
  if (!is.finite(total_mass) || !is.finite(ctrl_mass) || ctrl_mass <= 0) {
    stop("Cannot compute control-background pre-treatment mass ratio.")
  }
  max(1, total_mass / ctrl_mass)
}

oklahoma_report_rebuild_lambda_im <- function(holdout, win_km, kde_info) {
  if (is.null(holdout) || nrow(as.data.frame(holdout)) < 2L || is.null(win_km)) {
    return(NULL)
  }
  holdout <- as.data.frame(holdout)
  xy <- if (!is.null(kde_info$bw_sigma_xy)) {
    suppressWarnings(as.numeric(kde_info$bw_sigma_xy))
  } else {
    numeric(0)
  }
  xy <- xy[is.finite(xy) & xy > 0]
  if (length(xy) < 1L && !is.null(kde_info$bw_sigma)) {
    xy <- suppressWarnings(as.numeric(kde_info$bw_sigma))
    xy <- xy[is.finite(xy) & xy > 0]
  }
  if (length(xy) < 1L) return(NULL)
  sigma <- if (length(xy) >= 2L) xy[seq_len(2L)] else xy[[1]]
  X <- spatstat.geom::ppp(holdout$x, holdout$y, window = win_km)
  dens_fn <- if (is.function(get0("density.ppp", mode = "function", inherits = TRUE))) {
    get("density.ppp", mode = "function")
  } else if (requireNamespace("spatstat.explore", quietly = TRUE) &&
             is.function(spatstat.explore::density.ppp)) {
    spatstat.explore::density.ppp
  } else {
    stats::density
  }
  lambda_im <- dens_fn(X, sigma = sigma, edge = TRUE, at = "pixels")
  pos <- lambda_im$v[lambda_im$v > 0]
  if (length(pos) < 1L) return(NULL)
  min_nz <- min(pos, na.rm = TRUE)
  lambda_im$v[lambda_im$v <= 0] <- min_nz
  lambda_im
}

oklahoma_report_sem_last_labelling <- function(info, sem_obj = NULL) {
  pull_lab <- function(obj) {
    if (is.null(obj)) return(NULL)
    if (!is.null(obj$adaptive$adaptive_labelling)) {
      return(as.data.frame(obj$adaptive$adaptive_labelling))
    }
    if (!is.null(obj$fit) && !is.null(obj$fit$adaptive$adaptive_labelling)) {
      return(as.data.frame(obj$fit$adaptive$adaptive_labelling))
    }
    if (!is.null(obj$sem) && !is.null(obj$sem$adaptive$adaptive_labelling)) {
      return(as.data.frame(obj$sem$adaptive$adaptive_labelling))
    }
    NULL
  }
  lab <- pull_lab(sem_obj)
  if (is.null(lab)) lab <- pull_lab(info$semF)
  if (is.null(lab)) lab <- pull_lab(info$sem)
  if (is.null(lab)) lab <- pull_lab(info$fit)
  if (is.null(lab) && !is.null(info$pp_post_sem)) {
    lab <- as.data.frame(info$pp_post_sem)
  }
  lab
}

oklahoma_report_apply_last_labels <- function(catalog, last_lab) {
  out <- as.data.frame(catalog)
  loc <- if ("location_process" %in% names(out)) {
    as.character(out$location_process)
  } else {
    rep("control", nrow(out))
  }
  out$inferred_process <- ifelse(out$t < 0, "control", loc)
  if (is.null(last_lab) || nrow(as.data.frame(last_lab)) < 1L) {
    out$process <- out$inferred_process
    return(out)
  }
  last_lab <- as.data.frame(last_lab)
  lab_col <- if ("inferred_process" %in% names(last_lab)) {
    "inferred_process"
  } else if ("process" %in% names(last_lab)) {
    "process"
  } else {
    NULL
  }
  if (is.null(lab_col)) {
    out$process <- out$inferred_process
    return(out)
  }
  keys_cat <- oklahoma_report_event_key(out)
  keys_lab <- oklahoma_report_event_key(last_lab)
  map <- as.character(last_lab[[lab_col]])
  names(map) <- keys_lab
  hit <- keys_cat %in% names(map)
  n_post <- sum(out$t >= 0, na.rm = TRUE)
  n_hit_post <- sum(hit & out$t >= 0, na.rm = TRUE)
  if (n_post > 0L && n_hit_post < floor(0.9 * n_post)) {
    post_cat <- which(out$t >= 0)
    post_lab <- which(last_lab$t >= 0)
    if (length(post_cat) == length(post_lab) && length(post_cat) > 0L) {
      ord_c <- post_cat[order(out$t[post_cat], out$x[post_cat], out$y[post_cat])]
      ord_l <- post_lab[order(last_lab$t[post_lab], last_lab$x[post_lab], last_lab$y[post_lab])]
      out$inferred_process[ord_c] <- as.character(last_lab[[lab_col]][ord_l])
      hit <- rep(TRUE, nrow(out))
    }
  } else if (any(hit)) {
    out$inferred_process[hit] <- unname(map[keys_cat[hit]])
  }
  out$inferred_process[out$t < 0] <- "control"
  out$inferred_process[!out$inferred_process %in% c("control", "treated")] <- "control"
  out$process <- out$inferred_process
  out
}

oklahoma_report_aoi_boundary_labels <- function(catalog, aoi_treat_ss) {
  out <- as.data.frame(catalog)
  in_aoi <- spatstat.geom::inside.owin(out$x, out$y, aoi_treat_ss)
  out$inferred_process <- ifelse(out$t >= 0 & in_aoi, "treated", "control")
  out$process <- out$inferred_process
  out
}

oklahoma_report_aoi_weighted_catalog <- function(pp_pre, pp_post, aoi_ss, lambda_im) {
  pre <- as.data.frame(pp_pre)
  post <- as.data.frame(pp_post)
  keep <- c("x", "y", "t", "mag")
  miss <- setdiff(keep, names(pre))
  if (length(miss) > 0L) stop("pp_pre missing columns: ", paste(miss, collapse = ", "))
  miss <- setdiff(keep, names(post))
  if (length(miss) > 0L) stop("pp_post missing columns: ", paste(miss, collapse = ", "))
  pre <- pre[, keep, drop = FALSE]
  post <- post[, keep, drop = FALSE]
  in_aoi <- spatstat.geom::inside.owin(post$x, post$y, aoi_ss$treated)
  pre_w <- oklahoma_report_normalize_bg_weights(pre, aoi_ss$control, lambda_im)$new_df
  post_c <- oklahoma_report_normalize_bg_weights(
    post[!in_aoi, , drop = FALSE], aoi_ss$control, lambda_im
  )$new_df
  post_t <- oklahoma_report_normalize_bg_weights(
    post[in_aoi, , drop = FALSE], aoi_ss$treated, lambda_im
  )$new_df
  post_w <- rbind(post_c, post_t)
  all_df <- rbind(pre_w, post_w)
  all_df[order(all_df$t), , drop = FALSE]
}

oklahoma_report_eval_biv_kde_loglik <- function(params,
                                                realiz,
                                                windowT,
                                                windowS,
                                                aoi_ss,
                                                mass_ratio,
                                                cfg) {
  if (is.null(params) || is.null(realiz) || nrow(realiz) < 1L) return(NA_real_)
  pv <- if (is.list(params)) unlist(params) else params
  pv <- suppressWarnings(as.numeric(pv))
  if (!is.null(names(params))) names(pv) <- names(params)
  if (is.list(params) && is.null(names(pv))) names(pv) <- names(unlist(params))
  enforce <- isTRUE(cfg$ETAS_ENFORCE_FINITE_TRIGGER_MOMENTS)
  p_lo <- oklahoma_report_cfg_scalar(cfg, "ETAS_P_LOWER_BOUND", 0)
  q_lo <- oklahoma_report_cfg_scalar(cfg, "ETAS_Q_LOWER_BOUND", 1)
  rho_max <- oklahoma_report_cfg_scalar(cfg, "ETAS_BRANCHING_MAX", 0.98)
  t_trunc <- oklahoma_report_cfg_scalar(cfg, "SEM_T_TRUNC_DAYS", NA_real_)
  if (!is.finite(t_trunc) || t_trunc <= 0) t_trunc <- NULL
  m0 <- oklahoma_report_cfg_scalar(cfg, "ETAS_M0", min(realiz$mag, na.rm = TRUE))
  beta_gr <- oklahoma_report_cfg_scalar(cfg, "BETA_GR", NA_real_)
  ll <- tryCatch(
    loglik_etas_bivariate(
      params = pv,
      realiz = realiz,
      windowT = windowT,
      windowS = windowS,
      m0 = m0,
      control_state_space = aoi_ss$control,
      treated_state_space = aoi_ss$treated,
      background_rate_var = "W",
      treated_background_zero_before = 0,
      control_background_everywhere_before = 0,
      control_background_pre_mass_ratio = mass_ratio,
      beta_gr = beta_gr,
      enforce_finite_trigger_moments = enforce,
      p_lower_bound = p_lo,
      q_lower_bound = q_lo,
      max_branching_radius = rho_max,
      t_trunc = t_trunc
    ),
    error = function(e) NA_real_
  )
  as.numeric(ll)
}

oklahoma_report_fit_d_params <- function(info) {
  if (!is.null(info$F_params)) return(info$F_params)
  if (!is.null(info$params)) return(info$params)
  if (!is.null(info$semF) && !is.null(info$semF$etas_bivariate_params)) {
    return(info$semF$etas_bivariate_params)
  }
  if (!is.null(info$fit) && !is.null(info$fit$etas_bivariate_params)) {
    return(info$fit$etas_bivariate_params)
  }
  NULL
}

oklahoma_report_partition_sem_aoi_loglik <- function(partition_results,
                                                     primary_id,
                                                     primary_info,
                                                     pp_pre,
                                                     pp_post,
                                                     pp_pre_holdout,
                                                     win_km,
                                                     aoi_owin,
                                                     kde_info,
                                                     cfg,
                                                     format_label = identity) {
  empty <- data.frame(
    PartitionID = character(),
    Partition = character(),
    Last_labelling_LL = numeric(),
    AOI_boundary_LL = numeric(),
    Delta_last_minus_AOI = numeric(),
    Post_last_ne_AOI = integer(),
    N_post = integer(),
    stringsAsFactors = FALSE
  )
  if (is.null(win_km) || is.null(aoi_owin) || is.null(pp_pre) || is.null(pp_post)) {
    return(empty)
  }
  aoi_ss <- oklahoma_report_aoi_state_spaces(win_km, aoi_owin)
  if (is.null(aoi_ss)) return(empty)
  lambda_im <- oklahoma_report_rebuild_lambda_im(pp_pre_holdout, win_km, kde_info)
  if (is.null(lambda_im)) return(empty)

  catalog <- oklahoma_report_aoi_weighted_catalog(pp_pre, pp_post, aoi_ss, lambda_im)
  mass_ratio <- oklahoma_report_kde_pre_mass_ratio(lambda_im, aoi_ss$control)
  post_end <- oklahoma_report_cfg_scalar(cfg, "windowT_post", max(pp_post$t, na.rm = TRUE))
  if (!is.null(cfg$windowT_post)) {
    wt <- suppressWarnings(as.numeric(unlist(cfg$windowT_post, use.names = FALSE)))
    wt <- wt[is.finite(wt)]
    if (length(wt) >= 2L) post_end <- wt[[2]]
  }
  windowT <- c(min(pp_pre$t, na.rm = TRUE), post_end)
  aoi_realiz <- oklahoma_report_aoi_boundary_labels(catalog, aoi_ss$treated)

  pr <- if (is.null(partition_results)) list() else partition_results
  if (!is.null(primary_info) && !primary_id %in% names(pr)) {
    pr <- c(setNames(list(primary_info), primary_id), pr)
  } else if (!is.null(primary_info) && primary_id %in% names(pr)) {
    if (is.null(pr[[primary_id]]$F_params) && !is.null(primary_info$F_params)) {
      pr[[primary_id]]$F_params <- primary_info$F_params
    }
    if (is.null(pr[[primary_id]]$semF) && !is.null(primary_info$semF)) {
      pr[[primary_id]]$semF <- primary_info$semF
    }
    if (is.null(pr[[primary_id]]$pp_post_sem) && !is.null(primary_info$pp_post_sem)) {
      pr[[primary_id]]$pp_post_sem <- primary_info$pp_post_sem
    }
  }

  rows <- list()
  for (pname in names(pr)) {
    info <- pr[[pname]]
    if (is.null(info)) next
    params <- oklahoma_report_fit_d_params(info)
    if (is.null(params)) next
    pid <- as.character(info$label %||% pname)
    sem_for_row <- if (identical(pid, as.character(primary_id))) {
      info$semF %||% (if (!is.null(primary_info)) primary_info$semF else NULL)
    } else {
      info$semF
    }
    last_lab <- oklahoma_report_sem_last_labelling(info, sem_obj = sem_for_row)
    last_realiz <- oklahoma_report_apply_last_labels(catalog, last_lab)
    ll_last <- oklahoma_report_eval_biv_kde_loglik(
      params, last_realiz, windowT, win_km, aoi_ss, mass_ratio, cfg
    )
    ll_aoi <- oklahoma_report_eval_biv_kde_loglik(
      params, aoi_realiz, windowT, win_km, aoi_ss, mass_ratio, cfg
    )
    post_last <- last_realiz$inferred_process[last_realiz$t >= 0]
    post_aoi <- aoi_realiz$inferred_process[aoi_realiz$t >= 0]
    n_post <- length(post_aoi)
    n_diff <- if (length(post_last) == length(post_aoi)) {
      sum(post_last != post_aoi, na.rm = TRUE)
    } else {
      NA_integer_
    }
    rows[[length(rows) + 1L]] <- data.frame(
      PartitionID = pid,
      Partition = format_label(pid),
      Last_labelling_LL = ll_last,
      AOI_boundary_LL = ll_aoi,
      Delta_last_minus_AOI = ll_last - ll_aoi,
      Post_last_ne_AOI = as.integer(n_diff),
      N_post = as.integer(n_post),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) < 1L) return(empty)
  out <- do.call(rbind, rows)
  ids <- unique(as.character(out$PartitionID))
  out$PartitionID <- factor(
    out$PartitionID,
    levels = c(primary_id, sort(setdiff(ids, primary_id)))
  )
  out[order(out$PartitionID), , drop = FALSE]
}
