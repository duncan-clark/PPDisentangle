# Gambia Radii Study: Visualizations and CIs
#
# Run after gambia_radii_study.R (IPD) and/or gambia_radii_study_nonIPD.R (non-IPD).
# Produces:
# 1. Cumulative label changes over SEM iterations
# 2. Radius vs ATE (naive vs SEM)
# 3. Radius vs effect with 95% CIs (vanilla = bootstrap; SEM = Louis)
#
# DATA_TYPE: "IPD" | "nonIPD" | "both"
#   - "IPD": plot IPD results only
#  - "nonIPD": plot non-IPD (placebo) results only
#  - "both": plot both and add combined comparison
#
# Usage: source("gambia_radii_plots.R") from inst/gambia/ with package loaded

library(PPDisentangle)
library(ggplot2)
library(dplyr)
library(tidyr)
library(spatstat)
library(terra)
library(raster)
library(lubridate)

# =========================================================
# Configuration
# =========================================================
DATA_TYPE <- "IPD"   # "IPD" | "nonIPD" | "both"

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
OUT_DIR <- file.path(SCRIPT_DIR, "output")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

RESULTS_FILE_IPD    <- file.path(OUT_DIR, "gambia_radii_results.rds")
RESULTS_FILE_NONIPD <- file.path(OUT_DIR, "gambia_radii_results_nonIPD.rds")
BOOTSTRAP_B  <- 50   # set higher (e.g. 200) for production
BOOTSTRAP_SEED <- 12345

# =========================================================
# Load Results
# =========================================================
load_results <- function(path, label) {
  if (!file.exists(path)) {
    stop("Run the ", label, " radii study first to create ", path)
  }
  res <- readRDS(path)
  if (nrow(res$results) == 0) stop("No results in ", path)
  res$data_label <- label
  res
}

results_ipd <- NULL
results_nonipd <- NULL

if (DATA_TYPE %in% c("IPD", "both")) {
  results_ipd <- load_results(RESULTS_FILE_IPD, "IPD")
}
if (DATA_TYPE %in% c("nonIPD", "both")) {
  results_nonipd <- load_results(RESULTS_FILE_NONIPD, "non-IPD")
}
if (DATA_TYPE == "both" && (is.null(results_ipd) || is.null(results_nonipd))) {
  stop("DATA_TYPE='both' requires both IPD and non-IPD result files")
}

# For single-type mode, use the appropriate results as primary
if (DATA_TYPE == "IPD") {
  results_full <- results_ipd
} else if (DATA_TYPE == "nonIPD") {
  results_full <- results_nonipd
} else {
  results_full <- results_ipd  # primary for shared setup
}

df <- results_full$results
adaptive_by_radius <- results_full$adaptive_by_radius
louis_by_radius <- if (!is.null(results_full$louis_by_radius)) results_full$louis_by_radius else list()
config <- results_full$config
data_label <- results_full$data_label
file_suffix <- if (DATA_TYPE == "both") "" else paste0("_", tolower(gsub("-", "", data_label)))

# =========================================================
# Data Loading (same as radii study, for bootstrap)
# =========================================================
gambia_all <- dget(system.file("extdata/gambia/DISCLEAN.txt", package = "PPDisentangle"))
gambia_ipd <- gambia_all %>% filter(IPDEvent == "Y")
gambia_nonipd <- gambia_all %>% filter(IPDEvent == "N")
v_all <- vect(gambia_all, geom = c("long","lat"), crs = "EPSG:4326")
v_all_utm <- project(v_all, "EPSG:32628")
coords_all <- crds(v_all_utm)
gambia_all$easting <- coords_all[,1]
gambia_all$northing <- coords_all[,2]

TREATMENT_DATE <- as_datetime("2011-05-01")
t0 <- min(as.Date(gambia_ipd$date_IPD_Pne_Event))
TREATMENT_TIME <- as.numeric(as.Date(TREATMENT_DATE) - t0)
STUDY_END <- as.numeric(max(as.Date(gambia_ipd$date_IPD_Pne_Event)) - t0)
duration_post <- STUDY_END - TREATMENT_TIME

win_orig <- dget(system.file("extdata/gambia/GAMWIN_PROJ.txt", package = "PPDisentangle"))
S <- diag(c(1000, 1000)); t <- c(0, 0)
win_full <- affine.owin(win_orig, mat = S, vec = t)
hc <- dget(system.file("extdata/gambia/hc.txt", package = "PPDisentangle"))
hc <- affine(hc, mat = S, vec = t)

# IPD point data (for IPD bootstrap)
ipd_all <- gambia_all %>% filter(IPDEvent == "Y")
ipd_all <- ipd_all[inside.owin(ipd_all$easting, ipd_all$northing, win_full), ]
pp_data_ipd <- data.frame(
  x = ipd_all$easting, y = ipd_all$northing,
  t = as.numeric(as.Date(ipd_all$date_IPD_Pne_Event)) - as.numeric(t0),
  IPD = TRUE
)
pp_data_ipd <- pp_data_ipd[inside.owin(pp_data_ipd$x, pp_data_ipd$y, win_full), ]
pp_data_ipd$t <- pp_data_ipd$t - min(pp_data_ipd$t)
pp_data_ipd <- pp_data_ipd %>% filter(t <= STUDY_END)

# Non-IPD point data: full pool or reconstructed sample (for non-IPD bootstrap)
non_ipd_all <- gambia_all %>% filter(IPDEvent == "N")
non_ipd_all <- non_ipd_all[inside.owin(non_ipd_all$easting, non_ipd_all$northing, win_full), ]
pp_nonipd_full <- data.frame(
  x = non_ipd_all$easting, y = non_ipd_all$northing,
  t = as.numeric(as.Date(non_ipd_all$date_IPD_Pne_Event)) - as.numeric(t0),
  IPD = FALSE
)
pp_nonipd_full <- pp_nonipd_full[inside.owin(pp_nonipd_full$x, pp_nonipd_full$y, win_full), ]
pp_nonipd_full$t <- pp_nonipd_full$t - min(pp_nonipd_full$t)
pp_nonipd_full <- pp_nonipd_full %>% filter(t <= STUDY_END)

N_SAMPLE <- 1000
SAMPLE_SEED <- 2024
if (!is.null(results_nonipd) && !is.null(results_nonipd$config$N_SAMPLE)) {
  N_SAMPLE <- results_nonipd$config$N_SAMPLE
  SAMPLE_SEED <- if (!is.null(results_nonipd$config$SAMPLE_SEED)) results_nonipd$config$SAMPLE_SEED else 2024
}
reconstruct_nonipd_sample <- function() {
  pre_idx  <- which(pp_nonipd_full$t <  TREATMENT_TIME)
  post_idx <- which(pp_nonipd_full$t >= TREATMENT_TIME)
  set.seed(SAMPLE_SEED)
  n_pre  <- min(N_SAMPLE, length(pre_idx))
  n_post <- min(N_SAMPLE, length(post_idx))
  sampled <- c(sample(pre_idx, n_pre), sample(post_idx, n_post))
  pp_nonipd_full[sort(sampled), ]
}
pp_data_nonipd <- reconstruct_nonipd_sample()

# Background KDE: IPD uses non-IPD; non-IPD uses IPD
X_bg_ipd <- ppp(x = non_ipd_all$easting, y = non_ipd_all$northing, window = win_full)
X_bg_nonipd <- ppp(x = ipd_all$easting, y = ipd_all$northing, window = win_full)
bw_ipd <- suppressWarnings(bw.diggle(X_bg_ipd))
bw_nonipd <- suppressWarnings(bw.diggle(X_bg_nonipd))
lambda_im_ipd <- density(X_bg_ipd, sigma = bw_ipd, edge = TRUE, at = "pixels")
lambda_im_nonipd <- density(X_bg_nonipd, sigma = bw_nonipd, edge = TRUE, at = "pixels")
min_nz_ipd <- min(lambda_im_ipd$v[lambda_im_ipd$v > 0], na.rm = TRUE)
min_nz_nonipd <- min(lambda_im_nonipd$v[lambda_im_nonipd$v > 0], na.rm = TRUE)
lambda_im_ipd$v[lambda_im_ipd$v <= 0] <- min_nz_ipd
lambda_im_nonipd$v[lambda_im_nonipd$v <= 0] <- min_nz_nonipd

normalize_marks_gambia <- function(df_sub, win_sub, covariate_im, mark_name = "W") {
  if (nrow(df_sub) == 0) return(list(new_df = df_sub, mass = 0, norm = 0))
  cov_in_window <- covariate_im[win_sub, drop = FALSE]
  total_mass_raw <- integral.im(cov_in_window)
  target_area <- spatstat.geom::area(win_sub)
  norm_factor <- target_area / total_mass_raw
  vals_raw <- raster::extract(raster::raster(covariate_im), df_sub[, c("x", "y")])
  vals_raw[is.na(vals_raw)] <- 0
  df_sub[[mark_name]] <- vals_raw * norm_factor
  min_val <- min(df_sub[[mark_name]][df_sub[[mark_name]] > 0], na.rm = TRUE)
  if (is.infinite(min_val) | is.na(min_val)) min_val <- 1e-9
  df_sub[[mark_name]][df_sub[[mark_name]] <= 0] <- min_val
  list(new_df = df_sub, mass = total_mass_raw, norm = norm_factor)
}

# =========================================================
# Helper: estimate_savings_simple (match radii study)
# =========================================================
estimate_savings_simple <- function(p_ctrl, p_trtd, m_ctrl, m_trtd, m_target, duration) {
  denom_c <- 1 - p_ctrl[["K"]]
  denom_t <- 1 - p_trtd[["K"]]
  if (denom_c < 1e-3 || denom_t < 1e-3 || m_ctrl < 1e-10 || m_trtd < 1e-10 ||
      p_ctrl[["mu"]] < 1e-12 || p_trtd[["mu"]] < 1e-12)
    return(list(pct = NA))
  rate_c <- p_ctrl[["mu"]] / m_ctrl
  rate_t <- p_trtd[["mu"]] / m_trtd
  bg_c <- rate_c * m_target * duration
  bg_t <- rate_t * m_target * duration
  exp_c <- bg_c / denom_c
  exp_t <- bg_t / denom_t
  pct <- if (exp_c > 0) (exp_c - exp_t) / exp_c else NA
  if (!is.na(pct) && (pct < -10 || pct > 10)) return(list(pct = NA))
  return(list(pct = pct))
}

# =========================================================
# Viz 1: Cumulative points changed over SEM adaptive steps (per radius)
# =========================================================
build_flips_df <- function(adaptive_by_radius, data_label = "") {
  rows <- list()
  for (r_km in names(adaptive_by_radius)) {
    info <- adaptive_by_radius[[r_km]]
    n_post <- info$n_post
    global_iter <- 0L
    cum_flips <- 0
    for (a in seq_along(info$history)) {
      adapt <- info$history[[a]]
      if (is.null(adapt$max_metric_flips)) next
      for (i in seq_along(adapt$max_metric_flips)) {
        global_iter <- global_iter + 1L
        cum_flips <- cum_flips + adapt$max_metric_flips[i]
        rows[[length(rows) + 1]] <- data.frame(
          radius_km = as.numeric(r_km),
          iteration = global_iter,
          cum_changed = cum_flips,
          cum_pct = 100 * cum_flips / n_post,
          n_post = n_post,
          data_type = data_label
        )
      }
    }
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

plot_viz1 <- function(adaptive_by_radius, data_label, suffix) {
  flips_df <- build_flips_df(adaptive_by_radius, data_label)
  if (is.null(flips_df)) return(invisible(NULL))
  flips_df$radius_label <- paste0(flips_df$radius_km, " km")
  flips_df$radius_label <- factor(flips_df$radius_label,
    levels = paste0(sort(unique(flips_df$radius_km)), " km"))
  p1 <- ggplot(flips_df, aes(x = iteration, y = cum_pct)) +
    geom_line(linewidth = 0.8, color = "#377EB8") +
    geom_point(size = 1.2, color = "#377EB8") +
    facet_wrap(~radius_label, scales = "free_y", ncol = 4) +
    labs(
      title = paste0("Cumulative label changes over SEM iterations (", data_label, ")"),
      subtitle = "Flattening indicates convergence of the labelling",
      x = "Iteration", y = "Cumulative % points relabelled"
    ) +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold"))
  ggsave(file.path(OUT_DIR, paste0("gambia_viz1_pct_changed", suffix, ".png")), p1, width = 10, height = 8, dpi = 150)
  cat("Saved", file.path(OUT_DIR, paste0("gambia_viz1_pct_changed", suffix, ".png")), "\n")
}

if (DATA_TYPE == "both") {
  if (!is.null(results_ipd$adaptive_by_radius))
    plot_viz1(results_ipd$adaptive_by_radius, "IPD", "_IPD")
  if (!is.null(results_nonipd$adaptive_by_radius))
    plot_viz1(results_nonipd$adaptive_by_radius, "non-IPD (placebo)", "_nonIPD")
} else {
  if (!is.null(adaptive_by_radius)) {
    plot_viz1(adaptive_by_radius, data_label, file_suffix)
  } else {
    cat("No adaptive_history data for viz 1. Re-run radii study with updated package.\n")
  }
}

# =========================================================
# Viz 2: Radius vs ATE (naive vs SEM)
# =========================================================
plot_viz2 <- function(df_use, data_label, suffix) {
  df_plot <- df_use %>%
    filter(!sem_degenerate) %>%
    dplyr::select(radius_km, vanilla_savings_pct, sem_savings_pct) %>%
    tidyr::pivot_longer(cols = c(vanilla_savings_pct, sem_savings_pct),
                        names_to = "method", values_to = "savings_pct") %>%
    mutate(method = ifelse(method == "vanilla_savings_pct", "Naive", "SEM"))
  df_plot <- df_plot %>% filter(!is.na(savings_pct))
  if (nrow(df_plot) == 0) return(invisible(NULL))
  p2 <- ggplot(df_plot, aes(x = radius_km, y = savings_pct * 100, color = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = c(Naive = "#E41A1C", SEM = "#377EB8")) +
    labs(
      title = paste0("Effect (savings %) by radius — ", data_label),
      x = "Radius (km)", y = "Savings (%)",
      color = "Method"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  ggsave(file.path(OUT_DIR, paste0("gambia_viz2_radius_vs_ate", suffix, ".png")), p2, width = 7, height = 5, dpi = 150)
  cat("Saved", file.path(OUT_DIR, paste0("gambia_viz2_radius_vs_ate", suffix, ".png")), "\n")
}

if (DATA_TYPE == "both") {
  df_combined <- bind_rows(
    results_ipd$results %>% filter(!sem_degenerate) %>% mutate(data_type = "IPD"),
    results_nonipd$results %>% filter(!sem_degenerate) %>% mutate(data_type = "non-IPD (placebo)")
  )
  df_plot2 <- df_combined %>%
    dplyr::select(radius_km, vanilla_savings_pct, sem_savings_pct, data_type) %>%
    tidyr::pivot_longer(cols = c(vanilla_savings_pct, sem_savings_pct),
                        names_to = "method", values_to = "savings_pct") %>%
    mutate(method = ifelse(method == "vanilla_savings_pct", "Naive", "SEM"))
  df_plot2 <- df_plot2 %>% filter(!is.na(savings_pct))
  if (nrow(df_plot2) > 0) {
    p2b <- ggplot(df_plot2, aes(x = radius_km, y = savings_pct * 100, color = method, linetype = data_type)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2.5) +
      scale_color_manual(values = c(Naive = "#E41A1C", SEM = "#377EB8")) +
      scale_linetype_manual(values = c("IPD" = "solid", "non-IPD (placebo)" = "dashed")) +
      labs(title = "IPD vs non-IPD (placebo): savings % by radius", x = "Radius (km)", y = "Savings (%)") +
      theme_minimal() + theme(legend.position = "bottom")
    ggsave(file.path(OUT_DIR, "gambia_viz2_radius_vs_ate_combined.png"), p2b, width = 8, height = 5, dpi = 150)
    cat("Saved", file.path(OUT_DIR, "gambia_viz2_radius_vs_ate_combined.png"), "\n")
  }
  plot_viz2(results_ipd$results, "IPD", "_IPD")
  plot_viz2(results_nonipd$results, "non-IPD (placebo)", "_nonIPD")
} else {
  plot_viz2(df, data_label, file_suffix)
}

# =========================================================
# Viz 3: Parametric bootstrap for ATE + error bars
# =========================================================
# Bootstrap: for each radius, simulate from fitted vanilla Hawkes, refit, recompute savings.
# We need partition (treated_ss, control_ss), covariate_lookup, masses per radius.
# The radii study doesn't save these - we'd need to recompute or save them.
# Simplified: run bootstrap within a modified flow that has access to partition.
# For now, implement bootstrap that re-runs the data prep per radius (lightweight)
# and simulates from vanilla params.

run_bootstrap_replicate <- function(treated_ss, control_ss, lambda_im,
                                    res_pre_trtd, res_pre_ctrl, total_m, duration_post,
                                    cov_lookup_fn, win_full, params_ctrl, params_trtd, seed, cfg) {
  partition <- tess(tiles = list("control" = control_ss, "treated" = treated_ss))
  state_spaces <- list(control_ss, treated_ss)
  set.seed(seed)
  sim_df <- tryCatch({
    sim <- generate_inhomogeneous_hawkes(
      Omega = win_full,
      partition = partition,
      time_window = c(TREATMENT_TIME, STUDY_END),
      partition_processes = c("control", "treated"),
      hawkes_params = list(control = params_ctrl, treated = params_trtd),
      state_spaces = state_spaces,
      covariate_lookup = cov_lookup_fn,
      t_trunc = cfg$T_TRUNC
    )
    as.data.frame(sim)
  }, error = function(e) NULL)
  if (is.null(sim_df) || nrow(sim_df) < 10) return(NA)

  sim_df$location_process <- ifelse(
    inside.owin(sim_df$x, sim_df$y, treated_ss), "treated", "control")
  sim_post_trtd <- sim_df %>% filter(t >= TREATMENT_TIME, location_process == "treated")
  sim_post_ctrl <- sim_df %>% filter(t >= TREATMENT_TIME, location_process == "control")
  if (nrow(sim_post_trtd) < 5 || nrow(sim_post_ctrl) < 5) return(NA)

  norm_marks <- function(df_sub, win_sub) {
    cov_in_window <- lambda_im[win_sub, drop = FALSE]
    total_mass_raw <- integral.im(cov_in_window)
    target_area <- spatstat.geom::area(win_sub)
    norm_factor <- target_area / total_mass_raw
    vals_raw <- raster::extract(raster::raster(lambda_im), df_sub[, c("x", "y")])
    vals_raw[is.na(vals_raw)] <- 0
    df_sub$W <- vals_raw * norm_factor
    min_val <- min(df_sub$W[df_sub$W > 0], na.rm = TRUE)
    if (is.infinite(min_val) | is.na(min_val)) min_val <- 1e-9
    df_sub$W[df_sub$W <= 0] <- min_val
    list(new_df = df_sub, mass = total_mass_raw)
  }
  r_trtd <- norm_marks(sim_post_trtd, treated_ss)
  r_ctrl <- norm_marks(sim_post_ctrl, control_ss)

  fp <- if (!is.null(cfg$FIXED_BETA)) list(beta = cfg$FIXED_BETA) else NULL
  fit_trtd <- tryCatch(
    suppressWarnings(fit_hawkes(
      params_init = list(mu = 0.001, alpha = 1e-6, beta = 0.05, K = 0.2),
      realiz = r_trtd$new_df, windowT = range(r_trtd$new_df$t),
      windowS = treated_ss, background_rate_var = "W", maxit = 500,
      use_fast = TRUE, method = "Nelder-Mead", fixed_params = fp,
      t_trunc = cfg$T_TRUNC)),
    error = function(e) NULL
  )
  fit_ctrl <- tryCatch(
    suppressWarnings(fit_hawkes(
      params_init = list(mu = 0.001, alpha = 1e-6, beta = 0.05, K = 0.2),
      realiz = r_ctrl$new_df, windowT = range(r_ctrl$new_df$t),
      windowS = control_ss, background_rate_var = "W", maxit = 500,
      use_fast = TRUE, method = "Nelder-Mead", fixed_params = fp,
      t_trunc = cfg$T_TRUNC)),
    error = function(e) NULL
  )
  if (is.null(fit_trtd) || is.null(fit_ctrl)) return(NA)
  par_trtd <- as.list(as.numeric(fit_trtd$par))
  names(par_trtd) <- c("mu", "alpha", "beta", "K")
  par_ctrl <- as.list(as.numeric(fit_ctrl$par))
  names(par_ctrl) <- c("mu", "alpha", "beta", "K")
  if (par_ctrl[["K"]] >= 0.999 || par_trtd[["K"]] >= 0.999) return(NA)

  sav <- estimate_savings_simple(par_ctrl, par_trtd,
    r_ctrl$mass, r_trtd$mass, total_m, duration_post)
  sav$pct
}

run_bootstrap_for_dataset <- function(df_use, pp_data, lambda_im, louis_by_radius_use, cfg, data_label) {
  boot_results <- list()
  louis_sem_results <- list()
  cat("Running parametric bootstrap (B=", BOOTSTRAP_B, ") for ", data_label, " vanilla ATE...\n", sep = "")
  for (i in seq_len(nrow(df_use))) {
    r_km <- df_use$radius_km[i]
    r_m <- r_km * 1000
    discs <- lapply(seq_len(npoints(hc)), function(j)
      disc(radius = r_m, centre = c(hc$x[j], hc$y[j]), npoly = 128))
    treated_ss <- Reduce(union.owin, discs)
    control_ss <- setminus.owin(win_full, treated_ss)

    pp_study <- pp_data
    pp_study$location_process <- ifelse(inside.owin(pp_study$x, pp_study$y, treated_ss), "treated", "control")
    res_pre_trtd <- normalize_marks_gambia(
      pp_study %>% filter(t < TREATMENT_TIME, location_process == "treated"), treated_ss, lambda_im)
    res_pre_ctrl <- normalize_marks_gambia(
      pp_study %>% filter(t < TREATMENT_TIME, location_process == "control"), control_ss, lambda_im)
    total_m <- res_pre_trtd$mass + res_pre_ctrl$mass

    cov_lookup <- function(x, y) {
      rv <- as.numeric(spatstat.geom::interp.im(lambda_im, list(x = x, y = y)))
      rv[is.na(rv)] <- 0
      is_t <- inside.owin(x, y, treated_ss)
      rv[is_t] <- rv[is_t] * res_pre_trtd$norm
      rv[!is_t] <- rv[!is_t] * res_pre_ctrl$norm
      rv
    }

    params_ctrl <- list(mu = df_use$vanilla_mu_ctrl[i], alpha = df_use$vanilla_alpha_ctrl[i],
                        beta = df_use$vanilla_beta_ctrl[i], K = df_use$vanilla_K_ctrl[i])
    params_trtd <- list(mu = df_use$vanilla_mu_trtd[i], alpha = df_use$vanilla_alpha_trtd[i],
                        beta = df_use$vanilla_beta_trtd[i], K = df_use$vanilla_K_trtd[i])

    reps <- sapply(seq_len(BOOTSTRAP_B), function(b) {
      run_bootstrap_replicate(treated_ss, control_ss, lambda_im,
        res_pre_trtd, res_pre_ctrl, total_m, duration_post,
        cov_lookup, win_full, params_ctrl, params_trtd,
        BOOTSTRAP_SEED + b * 1000 + i, cfg)
    })
    valid <- reps[!is.na(reps) & is.finite(reps) & reps > -10 & reps < 10]
    boot_results[[as.character(r_km)]] <- list(
      mean = mean(valid),
      lwr = if (length(valid) >= 10) quantile(valid, 0.025) else NA,
      upr = if (length(valid) >= 10) quantile(valid, 0.975) else NA,
      n_valid = length(valid)
    )
    cat(sprintf("  %.1f km vanilla: %d/%d valid, mean=%.1f%% [%.1f, %.1f]\n",
        r_km, length(valid), BOOTSTRAP_B,
        if (length(valid) > 0) mean(valid) * 100 else NA,
        if (length(valid) >= 10) quantile(valid, 0.025) * 100 else NA,
        if (length(valid) >= 10) quantile(valid, 0.975) * 100 else NA))

    lr <- louis_by_radius_use[[as.character(r_km)]]
    if (!is.null(lr) && !is.null(lr$vcov) && nrow(lr$vcov) == 8) {
      set.seed(BOOTSTRAP_SEED + i)
      n_samp <- BOOTSTRAP_B * 3
      theta_samp <- tryCatch(
        MASS::mvrnorm(n_samp, lr$theta_hat, lr$vcov),
        error = function(e) NULL
      )
      if (!is.null(theta_samp)) {
        ate_samp <- apply(theta_samp, 1, function(th) {
          pc <- list(mu = th[1], alpha = th[2], beta = th[3], K = th[4])
          pt <- list(mu = th[5], alpha = th[6], beta = th[7], K = th[8])
          sav <- estimate_savings_simple(pc, pt, res_pre_ctrl$mass, res_pre_trtd$mass, total_m, duration_post)
          sav$pct
        })
        ate_valid <- ate_samp[!is.na(ate_samp) & is.finite(ate_samp) & ate_samp > -10 & ate_samp < 10]
        if (length(ate_valid) >= 10) {
          louis_sem_results[[as.character(r_km)]] <- list(
            lwr = quantile(ate_valid, 0.025),
            upr = quantile(ate_valid, 0.975),
            n_valid = length(ate_valid)
          )
          cat(sprintf("       SEM (Louis): %d valid samples, 95%% CI [%.1f, %.1f]\n",
              length(ate_valid), quantile(ate_valid, 0.025) * 100, quantile(ate_valid, 0.975) * 100))
        }
      }
    }
  }
  list(boot_results = boot_results, louis_sem_results = louis_sem_results)
}

plot_viz3 <- function(df_use, boot_results, louis_sem_results, data_label, suffix) {
  boot_df <- do.call(rbind, lapply(names(boot_results), function(r) {
    br <- boot_results[[r]]
    data.frame(radius_km = as.numeric(r), boot_mean = br$mean, boot_lwr = br$lwr, boot_upr = br$upr)
  }))
  louis_df <- do.call(rbind, lapply(names(louis_sem_results), function(r) {
    lr <- louis_sem_results[[r]]
    data.frame(radius_km = as.numeric(r), sem_lwr = lr$lwr, sem_upr = lr$upr)
  }))
  if (is.null(louis_df) || nrow(louis_df) == 0)
    louis_df <- data.frame(radius_km = numeric(), sem_lwr = numeric(), sem_upr = numeric())

  df_plot3 <- df_use %>%
    filter(!sem_degenerate) %>%
    dplyr::select(radius_km, vanilla_savings_pct, sem_savings_pct) %>%
    left_join(boot_df, by = "radius_km") %>%
    left_join(louis_df, by = "radius_km")

  df_long <- df_plot3 %>%
    filter(!is.na(vanilla_savings_pct)) %>%
    tidyr::pivot_longer(
      cols = c(vanilla_savings_pct, sem_savings_pct),
      names_to = "method", values_to = "savings_pct"
    ) %>%
    mutate(
      method = ifelse(method == "vanilla_savings_pct", "Vanilla", "SEM"),
      lwr = ifelse(method == "Vanilla", boot_lwr, sem_lwr),
      upr = ifelse(method == "Vanilla", boot_upr, sem_upr),
      x_dodge = radius_km + ifelse(method == "Vanilla", -0.15, 0.15)
    ) %>%
    arrange(method, radius_km)

  if (nrow(df_long) == 0) return(invisible(NULL))
  p3 <- ggplot(df_long, aes(x = x_dodge, y = savings_pct * 100, color = method)) +
    geom_errorbar(aes(ymin = lwr * 100, ymax = upr * 100), linewidth = 0.8, width = 0.2, na.rm = TRUE) +
    geom_point(size = 3, na.rm = TRUE) +
    geom_line(aes(group = method), linetype = "dashed", linewidth = 0.5, alpha = 0.7, na.rm = TRUE) +
    scale_x_continuous(breaks = unique(df_plot3$radius_km), minor_breaks = NULL) +
    scale_color_manual(values = c(Vanilla = "#E41A1C", SEM = "#377EB8")) +
    labs(
      title = paste0("Savings (%) by radius: Vanilla vs SEM with 95% CIs (", data_label, ")"),
      x = "Radius (km)", y = "Savings (%)",
      color = "Method"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  ggsave(file.path(OUT_DIR, paste0("gambia_viz3_radius_vs_ate_bootstrap", suffix, ".png")), p3, width = 7, height = 5, dpi = 150)
  cat("Saved", file.path(OUT_DIR, paste0("gambia_viz3_radius_vs_ate_bootstrap", suffix, ".png")), "\n")
}

# Run bootstrap and plot
`%||%` <- function(x, y) if (is.null(x)) y else x
if (DATA_TYPE == "both") {
  bl_ipd <- run_bootstrap_for_dataset(
    results_ipd$results, pp_data_ipd, lambda_im_ipd,
    results_ipd$louis_by_radius %||% list(), results_ipd$config, "IPD")
  bl_nonipd <- run_bootstrap_for_dataset(
    results_nonipd$results, pp_data_nonipd, lambda_im_nonipd,
    results_nonipd$louis_by_radius %||% list(), results_nonipd$config, "non-IPD")
  plot_viz3(results_ipd$results, bl_ipd$boot_results, bl_ipd$louis_sem_results, "IPD", "_IPD")
  plot_viz3(results_nonipd$results, bl_nonipd$boot_results, bl_nonipd$louis_sem_results, "non-IPD (placebo)", "_nonIPD")
  saveRDS(list(IPD = list(boot_results = bl_ipd$boot_results, louis_sem_results = bl_ipd$louis_sem_results),
               nonIPD = list(boot_results = bl_nonipd$boot_results, louis_sem_results = bl_nonipd$louis_sem_results)),
          file.path(OUT_DIR, "gambia_bootstrap_results.rds"))
} else {
  louis_use <- louis_by_radius %||% list()
  pp_use <- if (data_label == "IPD") pp_data_ipd else pp_data_nonipd
  lambda_use <- if (data_label == "IPD") lambda_im_ipd else lambda_im_nonipd
  bl <- run_bootstrap_for_dataset(df, pp_use, lambda_use, louis_use, config, data_label)
  plot_viz3(df, bl$boot_results, bl$louis_sem_results, data_label, file_suffix)
  saveRDS(list(boot_results = bl$boot_results, louis_sem_results = bl$louis_sem_results),
          file.path(OUT_DIR, paste0("gambia_bootstrap_results", file_suffix, ".rds")))
}
cat("Saved bootstrap results.\n")
