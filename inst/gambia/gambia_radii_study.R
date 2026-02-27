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
# Configuration & Constants
# =========================================================
TREATMENT_DATE <- as_datetime("2011-05-01")
FIXED_BETA <- 0.05  # ~20 day mean trigger time
RADII_KM <- seq(0.5, 20, length.out = 10) # 10 radii from 500m to 20km
N_BG_SAMPLE <- 1000 # Sample size for background KDE if needed, but we'll use all non-IPD

# =========================================================
# Data Loading & Preprocessing
# =========================================================
gambia_all <- dget(system.file("extdata/gambia/DISCLEAN.txt", package = "PPDisentangle"))
gambia_data <- gambia_all %>% filter(IPDEvent == "Y")

TREATMENT_TIME <- as.numeric(as.Date(TREATMENT_DATE) -
                               min(as.Date(gambia_data$date_IPD_Pne_Event)))
STUDY_END <- as.numeric(max(as.Date(gambia_data$date_IPD_Pne_Event)) -
                          min(as.Date(gambia_data$date_IPD_Pne_Event)))

# Window and Health Centers
win_orig <- dget(system.file("extdata/gambia/GAMWIN_PROJ.txt", package = "PPDisentangle"))
S  <- diag(c(1000, 1000))
t  <- c(0, 0)
win_full <- affine.owin(win_orig, mat = S, vec = t)
hc <- dget(system.file("extdata/gambia/hc.txt", package = "PPDisentangle"))
hc <- affine(hc, mat = S, vec = t)

# Project to UTM
v_all <- vect(gambia_all, geom = c("long","lat"), crs = "EPSG:4326")
v_all_utm <- project(v_all, "EPSG:32628")
coords_all <- crds(v_all_utm)
gambia_all$easting <- coords_all[,1]
gambia_all$northing <- coords_all[,2]

# IPD data frame
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

# Non-IPD for background estimation
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

run_full_fit <- function(df, win, label, fixed_beta = FIXED_BETA) {
  fp <- if (!is.null(fixed_beta)) list(beta = fixed_beta) else NULL
  starts <- list(list(mu = 0.0005, alpha = 1e-7, beta = 0.05, K = 0.3))
  best_fit <- NULL
  best_val <- -Inf
  for(i in seq_along(starts)) {
    fit <- tryCatch(
      fit_hawkes(params_init = starts[[i]], realiz = df, windowT = range(df$t),
                 windowS = win, background_rate_var = 'W', maxit = 2000,
                 use_fast = TRUE, method = "Nelder-Mead", fixed_params = fp),
      error = function(e) NULL
    )
    if(!is.null(fit) && fit$value > best_val) {
      best_val <- fit$value
      best_fit <- fit
    }
  }
  pars <- best_fit$par
  return(list(mu = pars['mu'], alpha = pars['alpha'], beta = pars['beta'], K = pars['K']))
}

estimate_savings_simple <- function(p_ctrl, p_trtd, m_ctrl, m_trtd, m_target, duration) {
  rate_c <- p_ctrl$mu / m_ctrl
  rate_t <- p_trtd$mu / m_trtd
  bg_c <- rate_c * m_target * duration
  bg_t <- rate_t * m_target * duration
  exp_c <- bg_c / (1 - p_ctrl$K)
  exp_t <- bg_t / (1 - p_trtd$K)
  return(list(exp_c = exp_c, exp_t = exp_t, savings = exp_c - exp_t, pct = (exp_c - exp_t)/exp_c))
}

# =========================================================
# Radii Study Loop
# =========================================================
results_list <- list()

for (r_km in RADII_KM) {
  r_m <- r_km * 1000
  cat("\n--- Processing Radius:", r_km, "km ---\n")
  
  # 1. Define partitions
  discs <- lapply(seq_len(npoints(hc)), function(i) disc(radius = r_m, centre = c(hc$x[i], hc$y[i]), npoly = 128))
  treated_ss <- Reduce(union.owin, discs)
  # Sample space: the full window
  win_study <- win_full
  control_ss <- setminus.owin(win_study, treated_ss)
  
  # 2. Filter data for this study window
  pp_study <- pp_data_full[inside.owin(pp_data_full$x, pp_data_full$y, win_study),]
  pp_study$location_process <- ifelse(inside.owin(pp_study$x, pp_study$y, treated_ss), "treated", "control")
  
  # 3. Background rate from ALL non-IPD cases in this window
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
  
  # 5. Vanilla Hawkes Fits (Post-treatment)
  cat("  Running Vanilla Fits...\n")
  vanilla_trtd <- run_full_fit(res_post_trtd$new_df, treated_ss, "Vanilla-Treated")
  vanilla_ctrl <- run_full_fit(res_post_ctrl$new_df, control_ss, "Vanilla-Control")
  
  # 6. SEM Fit
  cat("  Running SEM (simulation)...\n")
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
    hawkes_params_control = vanilla_ctrl, # Start with vanilla
    hawkes_params_treated = vanilla_trtd,
    N_labellings = 10,
    N_iter = 1,
    covariate_lookup = cov_lookup,
    background_rate_var = 'W',
    adaptive_control = list(
      update_control_params = TRUE, param_update_cadence = 5,
      proposal_method = "simulation", fixed_params = list(beta = FIXED_BETA),
      state_spaces = list(control_ss, treated_ss),
      iter = 20, n_props = 10, change_factor = 0.05, verbose = FALSE
    )
  )
  
  # 7. Calculate Effects
  duration_post <- STUDY_END - TREATMENT_TIME
  total_m <- res_pre_trtd$mass + res_pre_ctrl$mass
  
  sav_vanilla <- estimate_savings_simple(vanilla_ctrl, vanilla_trtd, res_pre_ctrl$mass, res_pre_trtd$mass, total_m, duration_post)
  sav_sem     <- estimate_savings_simple(sem_res$hawkes_params_control, sem_res$hawkes_params_treated, res_pre_ctrl$mass, res_pre_trtd$mass, total_m, duration_post)
  
  # 8. Store Results
  results_list[[as.character(r_km)]] <- data.frame(
    radius_km = r_km,
    vanilla_K_ctrl = vanilla_ctrl$K,
    vanilla_K_trtd = vanilla_trtd$K,
    vanilla_mu_ctrl = vanilla_ctrl$mu,
    vanilla_mu_trtd = vanilla_trtd$mu,
    sem_K_ctrl = sem_res$hawkes_params_control$K,
    sem_K_trtd = sem_res$hawkes_params_treated$K,
    sem_mu_ctrl = sem_res$hawkes_params_control$mu,
    sem_mu_trtd = sem_res$hawkes_params_treated$mu,
    vanilla_savings_pct = sav_vanilla$pct,
    sem_savings_pct = sav_sem$pct
  )
}

# =========================================================
# Final Comparison Table
# =========================================================
final_results <- do.call(rbind, results_list)
cat("\n\n================================================================\n")
cat("        RADII STUDY RESULTS: IPD ANALYSIS\n")
cat("================================================================\n")
print(final_results, row.names = FALSE)
cat("================================================================\n")

# Save results
saveRDS(final_results, "gambia_radii_results.rds")
cat("\nResults saved to gambia_radii_results.rds\n")
