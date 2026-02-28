library(PPDisentangle)
library(spatstat)
library(ggplot2)
library(dplyr)
library(data.table)
library(pbapply)
library(parallel)
library(doParallel)
library(lubridate)
library(terra)
library(raster)

# =========================================================
# Configuration & Tuning Parameters
# =========================================================

TREATMENT_DATE  <- as_datetime("2011-05-01")
TRUNC_DAYS      <- 30            # temporal kernel truncation (days)
FIXED_BETA      <- NULL          # beta is FREE with truncation

RADII_KM <- c(
  seq(0.5, 3.0, by = 0.5),
  seq(4, 10, by = 1),
  seq(12, 20, by = 2)
)

# --- Vanilla Hawkes fit ---
VANILLA_MAXIT   <- 5000
VANILLA_STARTS  <- list(
  list(mu = 0.0005, alpha = 1e-7, beta = 0.05, K = 0.3),
  list(mu = 0.001,  alpha = 1e-6, beta = 0.03, K = 0.2),
  list(mu = 0.0005, alpha = 1e-8, beta = 0.1,  K = 0.4)
)

# --- SEM outer loop ---
SEM_N_LABELLINGS <- 20
SEM_N_ITER       <- 2

# --- SEM adaptive (inner) step ---
SEM_INNER_ITER          <- 50
SEM_INNER_N_PROPS       <- 20
SEM_PARAM_UPDATE_CADENCE <- 10
SEM_CHANGE_FACTOR       <- 0.05
SEM_INCLUDE_STARTING    <- TRUE
SEM_UPDATE_STARTING     <- TRUE
SEM_VERBOSE             <- TRUE

# =========================================================
# Data Loading & Preprocessing
# =========================================================
gambia_all <- dget(system.file("extdata/gambia/DISCLEAN.txt", package = "PPDisentangle"))
gambia_data <- gambia_all %>% filter(IPDEvent == "Y")

TREATMENT_TIME <- as.numeric(as.Date(TREATMENT_DATE) -
                               min(as.Date(gambia_data$date_IPD_Pne_Event)))
STUDY_END <- as.numeric(max(as.Date(gambia_data$date_IPD_Pne_Event)) -
                          min(as.Date(gambia_data$date_IPD_Pne_Event)))

win_orig <- dget(system.file("extdata/gambia/GAMWIN_PROJ.txt", package = "PPDisentangle"))
S  <- diag(c(1000, 1000))
t  <- c(0, 0)
win_full <- affine.owin(win_orig, mat = S, vec = t)
hc <- dget(system.file("extdata/gambia/hc.txt", package = "PPDisentangle"))
hc <- affine(hc, mat = S, vec = t)

v_all <- vect(gambia_all, geom = c("long","lat"), crs = "EPSG:4326")
v_all_utm <- project(v_all, "EPSG:32628")
coords_all <- crds(v_all_utm)
gambia_all$easting <- coords_all[,1]
gambia_all$northing <- coords_all[,2]

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

# =========================================================
# Helper Functions
# =========================================================

normalize_marks_gambia <- function(df_sub, win_sub, covariate_im, mark_name = "W") {
  cov_in_window <- covariate_im[win_sub, drop = FALSE]
  total_mass_raw <- integral.im(cov_in_window)
  target_area <- spatstat.geom::area(win_sub)
  norm_factor <- target_area / total_mass_raw
  vals_raw <- raster::extract(raster(covariate_im), df_sub[, c("x", "y")])
  vals_raw[is.na(vals_raw)] <- 0
  df_sub[[mark_name]] <- vals_raw * norm_factor
  min_val <- min(df_sub[[mark_name]][df_sub[[mark_name]] > 0], na.rm = TRUE)
  if(is.infinite(min_val) | is.na(min_val)) min_val <- 1e-9
  df_sub[[mark_name]][df_sub[[mark_name]] <= 0] <- min_val
  return(list(new_df = df_sub, mass = total_mass_raw, norm = norm_factor))
}

run_full_fit <- function(df, win, label, fixed_beta = FIXED_BETA,
                         starts = VANILLA_STARTS, maxit = VANILLA_MAXIT,
                         t_trunc = TRUNC_DAYS) {
  fp <- if (!is.null(fixed_beta)) list(beta = fixed_beta) else NULL
  best_fit <- NULL
  best_val <- -Inf
  for(i in seq_along(starts)) {
    fit <- tryCatch(
      fit_hawkes(params_init = starts[[i]], realiz = df, windowT = range(df$t),
                 windowS = win, background_rate_var = 'W', maxit = maxit,
                 use_fast = TRUE, method = "Nelder-Mead", fixed_params = fp,
                 t_trunc = t_trunc),
      error = function(e) NULL
    )
    if(!is.null(fit) && fit$value > best_val) {
      best_val <- fit$value
      best_fit <- fit
    }
  }
  pars <- as.numeric(best_fit$par)
  names(pars) <- c("mu", "alpha", "beta", "K")
  return(as.list(pars))
}

as_parlist <- function(p) {
  v <- as.numeric(p); names(v) <- c("mu", "alpha", "beta", "K"); as.list(v)
}

estimate_savings_simple <- function(p_ctrl, p_trtd, m_ctrl, m_trtd, m_target, duration) {
  rate_c <- p_ctrl[["mu"]] / m_ctrl
  rate_t <- p_trtd[["mu"]] / m_trtd
  bg_c <- rate_c * m_target * duration
  bg_t <- rate_t * m_target * duration
  exp_c <- bg_c / (1 - p_ctrl[["K"]])
  exp_t <- bg_t / (1 - p_trtd[["K"]])
  return(list(exp_c = exp_c, exp_t = exp_t, savings = exp_c - exp_t, pct = (exp_c - exp_t)/exp_c))
}

# =========================================================
# Radii Study Loop (TRUNCATED KERNEL)
# =========================================================
results_list <- list()
t_study_start <- proc.time()[3]

cat(sprintf("Starting radii study (TRUNCATED KERNEL, T_max=%d days): %d radii, SEM inner=%d iters x %d props, outer=%d iters x %d labellings\n",
            TRUNC_DAYS, length(RADII_KM), SEM_INNER_ITER, SEM_INNER_N_PROPS, SEM_N_ITER, SEM_N_LABELLINGS))
cat(sprintf("Beta is %s\n", if (is.null(FIXED_BETA)) "FREE" else sprintf("FIXED at %g", FIXED_BETA)))

for (r_km in RADII_KM) {
  r_m <- r_km * 1000
  t_radius_start <- proc.time()[3]
  cat(sprintf("\n=== Radius: %.1f km ===\n", r_km))

  # 1. Define partitions
  discs <- lapply(seq_len(npoints(hc)), function(i) disc(radius = r_m, centre = c(hc$x[i], hc$y[i]), npoly = 128))
  treated_ss <- Reduce(union.owin, discs)
  win_study <- win_full
  control_ss <- setminus.owin(win_study, treated_ss)

  # 2. Filter data
  pp_study <- pp_data_full[inside.owin(pp_data_full$x, pp_data_full$y, win_study),]
  dups <- duplicated(pp_study[, c("x", "y", "t")])
  if (any(dups)) {
    cat(sprintf("  Removing %d duplicated points...\n", sum(dups)))
    pp_study <- pp_study[!dups, ]
  }
  pp_study$location_process <- ifelse(inside.owin(pp_study$x, pp_study$y, treated_ss), "treated", "control")
  pp_study$process <- pp_study$location_process
  pp_study$background <- TRUE

  n_pre_t  <- sum(pp_study$t <  TREATMENT_TIME & pp_study$location_process == "treated")
  n_pre_c  <- sum(pp_study$t <  TREATMENT_TIME & pp_study$location_process == "control")
  n_post_t <- sum(pp_study$t >= TREATMENT_TIME & pp_study$location_process == "treated")
  n_post_c <- sum(pp_study$t >= TREATMENT_TIME & pp_study$location_process == "control")
  cat(sprintf("  Points: pre-trt=%d pre-ctrl=%d post-trt=%d post-ctrl=%d (total=%d)\n",
              n_pre_t, n_pre_c, n_post_t, n_post_c, nrow(pp_study)))

  # 3. Background rate from ALL non-IPD cases
  non_ipd_study <- non_ipd_all[inside.owin(non_ipd_all$easting, non_ipd_all$northing, win_study), ]
  X_bg <- ppp(x = non_ipd_study$easting, y = non_ipd_study$northing, window = win_study)
  bw_sigma <- bw.diggle(X_bg)
  lambda_im <- density(X_bg, sigma = bw_sigma, edge = TRUE, at = "pixels")
  min_nz <- min(lambda_im$v[lambda_im$v > 0], na.rm = TRUE)
  lambda_im$v[lambda_im$v <= 0] <- min_nz

  # 4. Normalize and split
  res_pre_trtd  <- normalize_marks_gambia(pp_study %>% filter(t <  TREATMENT_TIME, location_process == "treated"), treated_ss, lambda_im)
  res_pre_ctrl  <- normalize_marks_gambia(pp_study %>% filter(t <  TREATMENT_TIME, location_process == "control"), control_ss, lambda_im)
  res_post_trtd <- normalize_marks_gambia(pp_study %>% filter(t >= TREATMENT_TIME, location_process == "treated"), treated_ss, lambda_im)
  res_post_ctrl <- normalize_marks_gambia(pp_study %>% filter(t >= TREATMENT_TIME, location_process == "control"), control_ss, lambda_im)

  # 5. Vanilla Hawkes Fits (Post-treatment) with truncation
  cat("  Vanilla fits (truncated)...")
  vanilla_trtd <- run_full_fit(res_post_trtd$new_df, treated_ss, "Vanilla-Treated")
  vanilla_ctrl <- run_full_fit(res_post_ctrl$new_df, control_ss, "Vanilla-Control")
  cat(sprintf(" done (K_trt=%.3f K_ctrl=%.3f beta_trt=%.4f beta_ctrl=%.4f)\n",
              vanilla_trtd[["K"]], vanilla_ctrl[["K"]], vanilla_trtd[["beta"]], vanilla_ctrl[["beta"]]))

  # 6. SEM Fit with truncation
  cat("  SEM fit (truncated)...")
  pp_sem <- rbind(res_pre_trtd$new_df, res_pre_ctrl$new_df, res_post_trtd$new_df, res_post_ctrl$new_df)
  cov_lookup <- function(x, y) {
    rv <- as.numeric(spatstat.geom::interp.im(lambda_im, list(x = x, y = y)))
    rv[is.na(rv)] <- 0
    is_t <- inside.owin(x, y, treated_ss)
    rv[is_t] <- rv[is_t] * res_pre_trtd$norm
    rv[!is_t] <- rv[!is_t] * res_pre_ctrl$norm
    return(rv)
  }

  sem_res <- adaptive_SEM(
    pp_data = pp_sem,
    partition = tess(tiles = list("control" = control_ss, "treated" = treated_ss)),
    partition_processes = c("control", "treated"),
    statespace = win_study,
    time_window = c(0, STUDY_END),
    treatment_time = TREATMENT_TIME,
    hawkes_params_control = vanilla_ctrl,
    hawkes_params_treated = vanilla_trtd,
    N_labellings = SEM_N_LABELLINGS,
    N_iter = SEM_N_ITER,
    covariate_lookup = cov_lookup,
    background_rate_var = 'W',
    t_trunc = TRUNC_DAYS,
    adaptive_control = list(
      update_control_params    = TRUE,
      param_update_cadence     = SEM_PARAM_UPDATE_CADENCE,
      proposal_method          = "simulation",
      fixed_params             = if (!is.null(FIXED_BETA)) list(beta = FIXED_BETA) else NULL,
      state_spaces             = list(control_ss, treated_ss),
      iter                     = SEM_INNER_ITER,
      n_props                  = SEM_INNER_N_PROPS,
      change_factor            = SEM_CHANGE_FACTOR,
      include_starting_data    = SEM_INCLUDE_STARTING,
      update_starting_data     = SEM_UPDATE_STARTING,
      verbose                  = SEM_VERBOSE
    )
  )

  sem_ctrl <- as_parlist(sem_res$hawkes_params_control)
  sem_trtd <- as_parlist(sem_res$hawkes_params_treated)
  cat(sprintf(" done (K_trt=%.3f K_ctrl=%.3f beta_trt=%.4f beta_ctrl=%.4f)\n",
              sem_trtd[["K"]], sem_ctrl[["K"]], sem_trtd[["beta"]], sem_ctrl[["beta"]]))

  # 7. Calculate Effects
  duration_post <- STUDY_END - TREATMENT_TIME
  total_m <- res_pre_trtd$mass + res_pre_ctrl$mass

  sav_vanilla <- estimate_savings_simple(vanilla_ctrl, vanilla_trtd, res_pre_ctrl$mass, res_pre_trtd$mass, total_m, duration_post)
  sav_sem     <- estimate_savings_simple(sem_ctrl, sem_trtd, res_pre_ctrl$mass, res_pre_trtd$mass, total_m, duration_post)

  # 8. Store Results
  elapsed_radius <- round((proc.time()[3] - t_radius_start) / 60, 1)
  cat(sprintf("  Savings: vanilla=%.1f%% sem=%.1f%% (%.1f min)\n",
              sav_vanilla$pct * 100, sav_sem$pct * 100, elapsed_radius))

  results_list[[as.character(r_km)]] <- data.frame(
    radius_km = r_km,
    n_post_trt = n_post_t, n_post_ctrl = n_post_c,
    vanilla_K_ctrl = vanilla_ctrl[["K"]],
    vanilla_K_trtd = vanilla_trtd[["K"]],
    vanilla_mu_ctrl = vanilla_ctrl[["mu"]],
    vanilla_mu_trtd = vanilla_trtd[["mu"]],
    vanilla_alpha_ctrl = vanilla_ctrl[["alpha"]],
    vanilla_alpha_trtd = vanilla_trtd[["alpha"]],
    vanilla_beta_ctrl = vanilla_ctrl[["beta"]],
    vanilla_beta_trtd = vanilla_trtd[["beta"]],
    sem_K_ctrl = sem_ctrl[["K"]],
    sem_K_trtd = sem_trtd[["K"]],
    sem_mu_ctrl = sem_ctrl[["mu"]],
    sem_mu_trtd = sem_trtd[["mu"]],
    sem_alpha_ctrl = sem_ctrl[["alpha"]],
    sem_alpha_trtd = sem_trtd[["alpha"]],
    sem_beta_ctrl = sem_ctrl[["beta"]],
    sem_beta_trtd = sem_trtd[["beta"]],
    vanilla_savings_pct = sav_vanilla$pct,
    sem_savings_pct = sav_sem$pct
  )
}

elapsed_total <- round((proc.time()[3] - t_study_start) / 60, 1)

# =========================================================
# Final Comparison Table
# =========================================================
final_results <- do.call(rbind, results_list)
final_results$sem_minus_vanilla <- final_results$sem_savings_pct - final_results$vanilla_savings_pct

cat("\n\n================================================================\n")
cat(sprintf("  RADII STUDY RESULTS: IPD (TRUNCATED KERNEL, T_max=%d days)\n", TRUNC_DAYS))
cat(sprintf("  Beta: %s\n", if (is.null(FIXED_BETA)) "FREE" else sprintf("fixed at %g", FIXED_BETA)))
cat(sprintf("  Total runtime: %.1f minutes\n", elapsed_total))
cat("================================================================\n")
print(final_results, row.names = FALSE)
cat("================================================================\n")

saveRDS(final_results, "gambia_radii_results_trunc30.rds")
cat("\nResults saved to gambia_radii_results_trunc30.rds\n")
