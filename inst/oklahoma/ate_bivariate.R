# Bivariate ATE via forward simulation of the full 2x2 ETAS law (cross-excitation on).
#
# contrast:
#   "observed"       — control-everywhere vs observed mixed allocation
#   "all_or_nothing" — control-everywhere vs treated-everywhere
#                      (same AoN estimand as the original marginal code, but bivariate)

ate_estim_bivariate <- function(
    biv_params,
    windowT,
    windowS,
    state_spaces_obs,
    label = "bivariate",
    n_sims = 500L,
    n_cores = 1L,
    m0 = 2.5,
    beta_gr = 2.3,
    filtration_history = NULL,
    t_trunc = NULL,
    n_tiles = 1L,
    crn_base_seed = 100000L,
    use_crn = TRUE,
    crn_pair = TRUE,
    quiet = FALSE,
    contrast = c("all_or_nothing", "observed"),
    # Optional named list of background lookups (control/treated), each
    # normalized to spatial mean one on its observed region but evaluable on
    # the whole domain (e.g. the fitted KDE).
    covariate_lookup = NULL,
    # Density-continuous everywhere-worlds (matches the fitted law where the
    # control background covers the whole domain pre-treatment at density
    # mu_0/|S0_obs|): each counterfactual background keeps its fitted density
    # and extends its support, so its total rate scales with the support
    # mass. FALSE reproduces the legacy total-rate-preserving worlds.
    density_reference = TRUE
) {
  contrast <- match.arg(contrast)
  if (is.null(biv_params)) stop("biv_params is NULL")
  pv <- unlist(biv_params)
  needed <- c(
    "mu_0", "mu_1", "A_00", "A_11", "A_01", "A_10",
    "alpha_m_00", "alpha_m_11", "alpha_m_01", "alpha_m_10",
    "c", "p", "D", "gamma", "q"
  )
  miss <- setdiff(needed, names(pv))
  if (length(miss) > 0L) {
    stop("Missing bivariate params: ", paste(miss, collapse = ", "))
  }
  params <- as.list(pv[needed])
  eta_channel <- function(A_name, alpha_name) {
    gap <- beta_gr - as.numeric(params[[alpha_name]])
    A <- as.numeric(params[[A_name]])
    if (!is.finite(A) || A < 0 || !is.finite(gap) || gap <= 1e-8) return(Inf)
    A * beta_gr / gap
  }
  offspring <- matrix(
    c(
      eta_channel("A_00", "alpha_m_00"),
      eta_channel("A_01", "alpha_m_01"),
      eta_channel("A_10", "alpha_m_10"),
      eta_channel("A_11", "alpha_m_11")
    ),
    nrow = 2, byrow = TRUE
  )
  rho <- if (all(is.finite(offspring))) {
    max(Re(eigen(offspring, only.values = TRUE)$values))
  } else {
    Inf
  }
  if (!is.finite(rho) || rho >= 1) {
    stop(sprintf(
      "%s has an explosive bivariate ATE law under beta_gr=%.4g (rho=%s).",
      label, beta_gr, format(rho, digits = 6)
    ))
  }

  if (is.null(filtration_history)) {
    pre_history <- data.frame(
      x = numeric(0), y = numeric(0), t = numeric(0),
      mag = numeric(0), process_id = integer(0)
    )
  } else {
    fh <- as.data.frame(filtration_history)
    pre_history <- data.frame(
      x = fh$x, y = fh$y, t = fh$t,
      mag = if ("mag" %in% names(fh)) fh$mag else rep(m0, nrow(fh)),
      process_id = if ("process_id" %in% names(fh)) {
        as.integer(fh$process_id)
      } else if ("inferred_process" %in% names(fh)) {
        as.integer(fh$inferred_process == "treated")
      } else {
        # Pre-treatment history is control-component by convention.
        rep(0L, nrow(fh))
      }
    )
  }

  # Control everywhere: immigrants only on full window; treated support empty.
  # Right-hand regime: observed mixed supports, or treated-everywhere.
  # Triggering (incl. cross) unrestricted in both cases.
  obs_ctrl <- if (!is.null(state_spaces_obs)) state_spaces_obs$control else NULL
  obs_treat <- if (!is.null(state_spaces_obs)) state_spaces_obs$treated else NULL
  if (isTRUE(density_reference) && (is.null(obs_ctrl) || is.null(obs_treat))) {
    stop("density_reference=TRUE requires state_spaces_obs$control/treated.")
  }
  everywhere_win <- if (isTRUE(density_reference) && is.null(covariate_lookup)) {
    # Flat background: extend over the modelled domain (union of observed
    # supports), matching the fit's (|S0|+|S1|)/|S0| mass-ratio convention.
    spatstat.geom::union.owin(obs_ctrl, obs_treat)
  } else {
    # KDE background: the fitted field lives on the full window, so extend
    # there (consistent with a full-domain KDE mass ratio in fitting).
    windowS
  }
  ss_all_control <- list(control = everywhere_win, treated = NULL)
  ref_all_control <- if (isTRUE(density_reference)) {
    list(control = spatstat.geom::area(obs_ctrl))
  } else {
    NULL
  }
  if (identical(contrast, "all_or_nothing")) {
    ss_right <- list(control = NULL, treated = everywhere_win)
    ref_right <- if (isTRUE(density_reference)) {
      list(treated = spatstat.geom::area(obs_treat))
    } else {
      NULL
    }
  } else {
    if (is.null(state_spaces_obs)) {
      stop("state_spaces_obs required for contrast='observed'")
    }
    ss_right <- list(control = obs_ctrl, treated = obs_treat)
    ref_right <- NULL
  }
  right_lab <- if (identical(contrast, "all_or_nothing")) "all-treated" else "obs-regime"

  sim_one <- function(ss, ref_areas) {
    sim_etas_bivariate(
      params = params,
      windowT = windowT,
      windowS = windowS,
      state_spaces = ss,
      m0 = m0,
      beta_gr = beta_gr,
      filtration = pre_history,
      covariate_lookup = covariate_lookup,
      bg_ref_areas = ref_areas,
      t_trunc = t_trunc
    )
  }

  n_sims_i <- max(1L, as.integer(n_sims))
  n_cores_use <- suppressWarnings(as.integer(n_cores))
  if (!is.finite(n_cores_use) || is.na(n_cores_use) || n_cores_use < 1L) n_cores_use <- 1L
  n_cores_use <- max(1L, min(n_cores_use, n_sims_i))

  if (!isTRUE(quiet)) {
    cat(sprintf(
      "  [ATE:bivariate/%s] %s: sims=%d cores=%d window=[%.1f,%.1f] t_trunc=%s\n",
      contrast, label, n_sims_i, n_cores_use, windowT[1], windowT[2],
      if (is.null(t_trunc)) "NULL" else format(t_trunc, digits = 4)
    ))
  }

  run_one <- function(s) {
    s_int <- as.integer(s)
    if (isTRUE(use_crn)) {
      seed_s <- as.integer(crn_base_seed + s_int)
      if (isTRUE(crn_pair)) {
        set.seed(seed_s)
        c_sim <- sim_one(ss_all_control, ref_all_control)
        set.seed(seed_s)
        t_sim <- sim_one(ss_right, ref_right)
      } else {
        set.seed(seed_s)
        c_sim <- sim_one(ss_all_control, ref_all_control)
        set.seed(as.integer(seed_s + 1000000L))
        t_sim <- sim_one(ss_right, ref_right)
      }
    } else {
      c_sim <- sim_one(ss_all_control, ref_all_control)
      t_sim <- sim_one(ss_right, ref_right)
    }
    c(c_count = nrow(c_sim), t_count = nrow(t_sim))
  }

  sim_ids <- as.list(seq_len(n_sims_i))
  ate_label_slug <- gsub("[^A-Za-z0-9]+", "_", label)
  ate_label_slug <- gsub("^_+|_+$", "", ate_label_slug)
  if (!nzchar(ate_label_slug)) ate_label_slug <- "model"
  parallel_label <- sprintf(
    "ate-biv-%s-%s", contrast, tolower(substr(ate_label_slug, 1L, 32L))
  )
  # Prefer the Oklahoma analysis runner when available (PSOCK/fork + progress).
  # Fall back to mclapply, then sequential lapply.
  sim_results <- if (n_cores_use > 1L && exists("run_parallel", mode = "function")) {
    run_parallel(sim_ids, run_one, cores = n_cores_use, label = parallel_label)
  } else if (n_cores_use > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(sim_ids, run_one, mc.cores = n_cores_use)
  } else {
    lapply(sim_ids, run_one)
  }
  c_counts <- vapply(sim_results, function(z) as.numeric(z[["c_count"]]), numeric(1))
  t_counts <- vapply(sim_results, function(z) as.numeric(z[["t_count"]]), numeric(1))
  total_saved <- c_counts - t_counts

  if (!isTRUE(quiet)) {
    cat(sprintf(
      "    %s: ctrl mean=%.1f | %s mean=%.1f | saved mean=%.1f\n",
      label,
      mean(c_counts, na.rm = TRUE),
      right_lab,
      mean(t_counts, na.rm = TRUE),
      mean(total_saved, na.rm = TRUE)
    ))
  }

  n_tiles_used <- max(1L, as.integer(n_tiles))
  all_nothing_sim <- data.frame(
    c_total = c_counts,
    t_total = t_counts,
    total_saved = total_saved,
    total_effect = total_saved,
    c_mean = c_counts / n_tiles_used,
    t_mean = t_counts / n_tiles_used,
    saved_per_tile = total_saved / n_tiles_used,
    ATE = total_saved / n_tiles_used
  )

  list(
    all_nothing_sim = all_nothing_sim,
    saved_naive = NA_real_,
    saved_spillover = NA_real_,
    total_saved_naive = NA_real_,
    total_saved_spillover = NA_real_,
    ATE_naive = NA_real_,
    ATE_spillover = NA_real_,
    total_naive = NA_real_,
    total_spillover = NA_real_,
    n_tiles_used = n_tiles_used,
    treated_pp = NULL,
    control_pp = NULL,
    analytic = NULL,
    analytic_saved = c(
      eta_ctrl_minus_treat = NA_real_,
      total_ctrl_minus_treat = NA_real_
    ),
    ate_method = paste0("bivariate_", contrast),
    contrast = contrast,
    branching_radius = rho,
    biv_params = params
  )
}

# Rebuild spatstat window + observed control/treated supports (km).
oklahoma_build_state_spaces <- function(data_dir, crs_proj = 5070L) {
  if (!requireNamespace("sf", quietly = TRUE)) stop("sf required")
  if (!requireNamespace("tigris", quietly = TRUE)) stop("tigris required")
  if (!requireNamespace("spatstat.geom", quietly = TRUE)) stop("spatstat required")
  if (!exists("oklahoma_sf_features_to_owins_km", mode = "function", inherits = TRUE)) {
    # When sourced from oklahoma_analysis.R, oklahoma_geometry.R is already loaded.
    # Standalone use: try adjacent helper file.
    geom_path <- file.path(dirname(data_dir), "oklahoma_geometry.R")
    if (!file.exists(geom_path)) {
      geom_path <- file.path(getwd(), "inst/oklahoma/oklahoma_geometry.R")
    }
    if (!file.exists(geom_path)) {
      stop("oklahoma_geometry.R helpers are not loaded.")
    }
    source(geom_path, local = FALSE)
  }

  options(tigris_use_cache = TRUE)
  counties_sf <- tigris::counties(state = "OK", cb = TRUE, year = 2022)
  counties_sf <- sf::st_transform(counties_sf, crs_proj)
  counties_sf <- sf::st_make_valid(counties_sf)
  ok_boundary <- sf::st_make_valid(sf::st_union(counties_sf))
  bb <- sf::st_bbox(ok_boundary)
  win_km <- spatstat.geom::owin(
    xrange = c(bb["xmin"], bb["xmax"]) / 1000,
    yrange = c(bb["ymin"], bb["ymax"]) / 1000
  )

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

  list(
    win_km = win_km,
    partition = partition,
    partition_processes = partition_processes,
    treated_idx = treated_idx,
    state_spaces = list(
      control = spatstat.geom::as.owin(partition[!treated_idx]),
      treated = spatstat.geom::as.owin(partition[treated_idx])
    ),
    n_tiles = partition$n
  )
}
