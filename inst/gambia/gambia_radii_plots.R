# Gambia Radii Study: Visualizations and CIs
#
# Run after gambia_radii_study.R. Loads results and produces:
# 1. % points changed over SEM iterations (all radii combined)
# 2. Radius vs ATE (naive vs SEM)
# 3. Radius vs effect: vanilla = parametric bootstrap CI; SEM = Louis-method CI
#    (sample params from N(MLE, vcov), compute ATE for each)
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
SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
OUT_DIR <- file.path(SCRIPT_DIR, "output")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

RESULTS_FILE <- file.path(OUT_DIR, "gambia_radii_results.rds")
BOOTSTRAP_B  <- 50   # set higher (e.g. 200) for production
BOOTSTRAP_SEED <- 12345

# =========================================================
# Load Results
# =========================================================
if (!file.exists(RESULTS_FILE)) {
  stop("Run gambia_radii_study.R first to create ", RESULTS_FILE)
}
results_full <- readRDS(RESULTS_FILE)
df <- results_full$results
adaptive_by_radius <- results_full$adaptive_by_radius
louis_by_radius <- if (!is.null(results_full$louis_by_radius)) results_full$louis_by_radius else list()
config <- results_full$config

if (nrow(df) == 0) {
  stop("No results in ", RESULTS_FILE)
}

# =========================================================
# Data Loading (same as radii study, for bootstrap)
# =========================================================
gambia_all <- dget(system.file("extdata/gambia/DISCLEAN.txt", package = "PPDisentangle"))
gambia_data <- gambia_all %>% filter(IPDEvent == "Y")
v_all <- vect(gambia_all, geom = c("long","lat"), crs = "EPSG:4326")
v_all_utm <- project(v_all, "EPSG:32628")
coords_all <- crds(v_all_utm)
gambia_all$easting <- coords_all[,1]
gambia_all$northing <- coords_all[,2]

TREATMENT_DATE <- as_datetime("2011-05-01")
TREATMENT_TIME <- as.numeric(as.Date(TREATMENT_DATE) - min(as.Date(gambia_data$date_IPD_Pne_Event)))
STUDY_END <- as.numeric(max(as.Date(gambia_data$date_IPD_Pne_Event)) - min(as.Date(gambia_data$date_IPD_Pne_Event)))
duration_post <- STUDY_END - TREATMENT_TIME

win_orig <- dget(system.file("extdata/gambia/GAMWIN_PROJ.txt", package = "PPDisentangle"))
S <- diag(c(1000, 1000)); t <- c(0, 0)
win_full <- affine.owin(win_orig, mat = S, vec = t)
hc <- dget(system.file("extdata/gambia/hc.txt", package = "PPDisentangle"))
hc <- affine(hc, mat = S, vec = t)

pp_data_full <- data.frame(
  x = gambia_all$easting[gambia_all$IPDEvent == "Y"],
  y = gambia_all$northing[gambia_all$IPDEvent == "Y"],
  t = as.numeric(as.Date(gambia_all$date_IPD_Pne_Event[gambia_all$IPDEvent == "Y"])) -
      min(as.numeric(as.Date(gambia_data$date_IPD_Pne_Event))),
  IPD = TRUE
)
pp_data_full <- pp_data_full[inside.owin(pp_data_full$x, pp_data_full$y, win_full),]
pp_data_full$t <- pp_data_full$t - min(pp_data_full$t)
pp_data_full <- pp_data_full %>% filter(t <= STUDY_END)

non_ipd_all <- gambia_all %>% filter(IPDEvent == "N")
non_ipd_all <- non_ipd_all[inside.owin(non_ipd_all$easting, non_ipd_all$northing, win_full), ]
X_bg <- ppp(x = non_ipd_all$easting, y = non_ipd_all$northing, window = win_full)
bw_sigma <- suppressWarnings(bw.diggle(X_bg))
lambda_im <- density(X_bg, sigma = bw_sigma, edge = TRUE, at = "pixels")
min_nz <- min(lambda_im$v[lambda_im$v > 0], na.rm = TRUE)
lambda_im$v[lambda_im$v <= 0] <- min_nz

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
build_flips_df <- function(adaptive_by_radius) {
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
          n_post = n_post
        )
      }
    }
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

flips_df <- build_flips_df(adaptive_by_radius)
if (!is.null(flips_df)) {
  flips_df$radius_label <- paste0(flips_df$radius_km, " km")
  flips_df$radius_label <- factor(flips_df$radius_label,
    levels = paste0(sort(unique(flips_df$radius_km)), " km"))

  p1 <- ggplot(flips_df, aes(x = iteration, y = cum_pct)) +
    geom_line(linewidth = 0.8, color = "#377EB8") +
    geom_point(size = 1.2, color = "#377EB8") +
    facet_wrap(~radius_label, scales = "free_y", ncol = 4) +
    labs(
      title = "Cumulative label changes over SEM iterations",
      subtitle = "Flattening indicates convergence of the labelling",
      x = "Iteration", y = "Cumulative % points relabelled"
    ) +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold"))
  ggsave(file.path(OUT_DIR, "gambia_viz1_pct_changed.png"), p1, width = 10, height = 8, dpi = 150)
  cat("Saved", file.path(OUT_DIR, "gambia_viz1_pct_changed.png"), "\n")
} else {
  cat("No adaptive_history data for viz 1. Re-run radii study with updated package.\n")
}

# =========================================================
# Viz 2: Radius vs ATE (naive vs SEM)
# =========================================================
df_plot <- df %>%
  filter(!sem_degenerate) %>%
  dplyr::select(radius_km, vanilla_savings_pct, sem_savings_pct) %>%
  tidyr::pivot_longer(cols = c(vanilla_savings_pct, sem_savings_pct),
                      names_to = "method", values_to = "savings_pct") %>%
  mutate(method = ifelse(method == "vanilla_savings_pct", "Naive", "SEM"))

p2 <- ggplot(df_plot %>% filter(!is.na(savings_pct)), aes(x = radius_km, y = savings_pct * 100, color = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c(Naive = "#E41A1C", SEM = "#377EB8")) +
  labs(
    title = "Vaccine effect (savings %) by radius",
    x = "Radius (km)", y = "Savings (%)",
    color = "Method"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "gambia_viz2_radius_vs_ate.png"), p2, width = 7, height = 5, dpi = 150)
cat("Saved", file.path(OUT_DIR, "gambia_viz2_radius_vs_ate.png"), "\n")

# =========================================================
# Viz 3: Parametric bootstrap for ATE + error bars
# =========================================================
# Bootstrap: for each radius, simulate from fitted vanilla Hawkes, refit, recompute savings.
# We need partition (treated_ss, control_ss), covariate_lookup, masses per radius.
# The radii study doesn't save these - we'd need to recompute or save them.
# Simplified: run bootstrap within a modified flow that has access to partition.
# For now, implement bootstrap that re-runs the data prep per radius (lightweight)
# and simulates from vanilla params.

run_bootstrap_replicate <- function(r_km, params_ctrl, params_trtd,
                                    treated_ss, control_ss, lambda_im,
                                    res_pre_trtd, res_pre_ctrl, total_m, duration_post,
                                    cov_lookup_fn, win_full, seed) {
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
      t_trunc = config$T_TRUNC
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

  fp <- if (!is.null(config$FIXED_BETA)) list(beta = config$FIXED_BETA) else NULL
  fit_trtd <- tryCatch(
    suppressWarnings(fit_hawkes(
      params_init = list(mu = 0.001, alpha = 1e-6, beta = 0.05, K = 0.2),
      realiz = r_trtd$new_df, windowT = range(r_trtd$new_df$t),
      windowS = treated_ss, background_rate_var = "W", maxit = 500,
      use_fast = TRUE, method = "Nelder-Mead", fixed_params = fp,
      t_trunc = config$T_TRUNC)),
    error = function(e) NULL
  )
  fit_ctrl <- tryCatch(
    suppressWarnings(fit_hawkes(
      params_init = list(mu = 0.001, alpha = 1e-6, beta = 0.05, K = 0.2),
      realiz = r_ctrl$new_df, windowT = range(r_ctrl$new_df$t),
      windowS = control_ss, background_rate_var = "W", maxit = 500,
      use_fast = TRUE, method = "Nelder-Mead", fixed_params = fp,
      t_trunc = config$T_TRUNC)),
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

# We need to rebuild geometry and lambda_im per radius. Load full data.
non_ipd_all <- gambia_all %>% filter(IPDEvent == "N")
non_ipd_all <- non_ipd_all[inside.owin(non_ipd_all$easting, non_ipd_all$northing, win_full), ]
X_bg <- ppp(x = non_ipd_all$easting, y = non_ipd_all$northing, window = win_full)
bw_sigma <- suppressWarnings(bw.diggle(X_bg))
lambda_im <- density(X_bg, sigma = bw_sigma, edge = TRUE, at = "pixels")
min_nz <- min(lambda_im$v[lambda_im$v > 0], na.rm = TRUE)
lambda_im$v[lambda_im$v <= 0] <- min_nz

pp_data_full <- data.frame(
  x = gambia_all$easting[gambia_all$IPDEvent == "Y"],
  y = gambia_all$northing[gambia_all$IPDEvent == "Y"],
  t = as.numeric(as.Date(gambia_all$date_IPD_Pne_Event[gambia_all$IPDEvent == "Y"])) -
      min(as.numeric(as.Date(gambia_data$date_IPD_Pne_Event))),
  IPD = TRUE
)
pp_data_full <- pp_data_full[inside.owin(pp_data_full$x, pp_data_full$y, win_full), ]
pp_data_full$t <- pp_data_full$t - min(pp_data_full$t)
pp_data_full <- pp_data_full %>% filter(t <= STUDY_END)

cat("Running parametric bootstrap (B=", BOOTSTRAP_B, ") for vanilla ATE...\n", sep = "")
boot_results <- list()
louis_sem_results <- list()
for (i in seq_len(nrow(df))) {
  r_km <- df$radius_km[i]
  r_m <- r_km * 1000
  discs <- lapply(seq_len(npoints(hc)), function(j)
    disc(radius = r_m, centre = c(hc$x[j], hc$y[j]), npoly = 128))
  treated_ss <- Reduce(union.owin, discs)
  control_ss <- setminus.owin(win_full, treated_ss)

  pp_study <- pp_data_full
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

  params_ctrl <- list(mu = df$vanilla_mu_ctrl[i], alpha = df$vanilla_alpha_ctrl[i],
                      beta = df$vanilla_beta_ctrl[i], K = df$vanilla_K_ctrl[i])
  params_trtd <- list(mu = df$vanilla_mu_trtd[i], alpha = df$vanilla_alpha_trtd[i],
                      beta = df$vanilla_beta_trtd[i], K = df$vanilla_K_trtd[i])

  reps <- sapply(seq_len(BOOTSTRAP_B), function(b) {
    run_bootstrap_replicate(r_km, params_ctrl, params_trtd,
      treated_ss, control_ss, lambda_im,
      res_pre_trtd, res_pre_ctrl, total_m, duration_post,
      cov_lookup, win_full, BOOTSTRAP_SEED + b * 1000 + i)
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

  # Louis-based SEM: sample params from N(theta_hat, vcov), compute ATE
  lr <- louis_by_radius[[as.character(r_km)]]
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

# Build bootstrap df for plotting
boot_df <- do.call(rbind, lapply(names(boot_results), function(r) {
  br <- boot_results[[r]]
  data.frame(radius_km = as.numeric(r), boot_mean = br$mean, boot_lwr = br$lwr, boot_upr = br$upr)
}))

louis_df <- do.call(rbind, lapply(names(louis_sem_results), function(r) {
  lr <- louis_sem_results[[r]]
  data.frame(radius_km = as.numeric(r), sem_lwr = lr$lwr, sem_upr = lr$upr)
}))
if (is.null(louis_df)) louis_df <- data.frame(radius_km = numeric(), sem_lwr = numeric(), sem_upr = numeric())

df_plot3 <- df %>%
  filter(!sem_degenerate) %>%
  dplyr::select(radius_km, vanilla_savings_pct, sem_savings_pct) %>%
  left_join(boot_df, by = "radius_km") %>%
  left_join(louis_df, by = "radius_km")

# Long format for cleaner plotting: vanilla + SEM with error bars side-by-side
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

p3 <- ggplot(df_long, aes(x = x_dodge, y = savings_pct * 100, color = method)) +
  geom_errorbar(aes(ymin = lwr * 100, ymax = upr * 100), linewidth = 0.8, width = 0.2, na.rm = TRUE) +
  geom_point(size = 3, na.rm = TRUE) +
  geom_line(aes(group = method), linetype = "dashed", linewidth = 0.5, alpha = 0.7, na.rm = TRUE) +
  scale_x_continuous(breaks = unique(df_plot3$radius_km), minor_breaks = NULL) +
  scale_color_manual(values = c(Vanilla = "#E41A1C", SEM = "#377EB8")) +
  labs(
    title = "Savings (%) by radius: Vanilla vs SEM with 95% CIs",
    x = "Radius (km)", y = "Savings (%)",
    color = "Method"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "gambia_viz3_radius_vs_ate_bootstrap.png"), p3, width = 7, height = 5, dpi = 150)
cat("Saved", file.path(OUT_DIR, "gambia_viz3_radius_vs_ate_bootstrap.png"), "\n")

saveRDS(list(boot_results = boot_results, boot_df = boot_df,
             louis_sem_results = louis_sem_results, louis_df = louis_df),
        file.path(OUT_DIR, "gambia_bootstrap_results.rds"))
cat("Saved", file.path(OUT_DIR, "gambia_bootstrap_results.rds"), "\n")
