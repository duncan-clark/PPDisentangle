library(PPDisentangle)
library(spatstat)
library(ggplot2)
library(dplyr)
library(data.table)
library(pbapply)
library(parallel)
library(doParallel)
library(R.utils)
library(caret)
library(terra)
library(raster)
library(lubridate)

# TODO
# - All time just IPD,
# - Correct treatment time
# - Grid over diameters (500m -> 20km)

TREATMENT_DATE <- as_datetime("2011-05-01")


# 2 years of data # to make smaller
# Define treatment as "close to health center" distance in km
HEALTH_CENTER_DIST <- 10000 # 10 km

# Get data:
gambia_data <- dget(system.file("extdata/gambia/DISCLEAN.txt", package = "PPDisentangle"))
gambia_data <- gambia_data %>% filter(IPDEvent == "Y")
TREATMENT_TIME <- as.numeric(as.Date(TREATMENT_DATE) -
                               min(as.Date(gambia_data$date_IPD_Pne_Event)))
STUDY_END <- as.numeric(max(as.Date(gambia_data$date_IPD_Pne_Event)) -
                          min(as.Date(gambia_data$date_IPD_Pne_Event)))
# window
win <- dget(system.file("extdata/gambia/GAMWIN_PROJ.txt", package = "PPDisentangle"))
S  <- diag(c(1000, 1000))
t  <- c(0, 0)
win <- affine.owin(win, mat = S, vec = t)
hc <- dget(system.file("extdata/gambia/hc.txt", package = "PPDisentangle"))
hc <- affine(hc, mat = S, vec = t)

# make the northing and easting a projection with terra - so we can recreate later:
# 1. Create a SpatVector of points, original CRS WGS84:
v <- vect(gambia_data, geom = c("long","lat"), crs = "EPSG:4326")
# 2. Project to UTM Zone 28N:
v_utm <- project(v, "EPSG:32628")
# 3. Extract coordinates: in km
coords <- crds(v_utm)
gambia_data$easting <- coords[,1]
gambia_data$northing <- coords[,2]

pp_data <- data.frame(
  x = gambia_data$easting,
  y = gambia_data$northing,
  lat = gambia_data$lat,
  long = gambia_data$long,
  t = as.numeric(as.Date(gambia_data$date_IPD_Pne_Event)) - min(as.numeric(as.Date(gambia_data$date_IPD_Pne_Event))),
  IPD = gambia_data$IPDEvent == "Y"
)
# drop data not in the window
pp_data <- pp_data[inside.owin(pp_data$x, pp_data$y, win),]
# make time into days since data started
pp_data$t <- pp_data$t - min(pp_data$t)
# Drops from 15000 -> 3000 cases
print(paste("Dropping ", sum(pp_data$t > STUDY_END), " cases past end of study"))
pp_data <- pp_data %>% filter(t <= STUDY_END)

# count IPD non/IPD cases:
table(pp_data$IPD)

# =============================================
# Partition the space
# =============================================
# 1) Distance to nearest X
dmap <- distmap(hc)
bands <- cut(as.im(dmap),
             breaks = c(0,HEALTH_CENTER_DIST,Inf),
             labels = c("close","far"))
band_names <- c("close","far")
tiles_list <- lapply(band_names, function(lbl) {
  as.owin(bands == lbl)
})
plot(bands, main="Distance bands from nearest healthcare center",
     xlab="Easting", ylab="Northing")
contour(dmap, add=TRUE, col="white")

# Make actual state-spaces:
discs_list <- lapply(seq_len(npoints(hc)), function(i) {
  disc(radius = HEALTH_CENTER_DIST,
       centre = c(hc$x[i], hc$y[i]),
       npoly = 200)    # increase npoly for smoother circles
})

disc_list_double <- lapply(seq_len(npoints(hc)), function(i) {
  disc(radius = HEALTH_CENTER_DIST * 2,
       centre = c(hc$x[i], hc$y[i]),
       npoly = 200)    # increase npoly for smoother circles
})

# if using whole win
treated_state_space <- Reduce(union.owin, discs_list)
control_state_space <- setminus.owin(win, treated_state_space)

# if do far near from health center
control_state_space <- setminus.owin(Reduce(union.owin, disc_list_double), treated_state_space)
win <- Reduce(union.owin, disc_list_double)
pp_data <- pp_data[inside.owin(pp_data$x, pp_data$y, win),]

# make this into a tessalation object
tiles <- discs_list
tiles[[length(tiles)+1]] <- control_state_space
partition <- tess(tiles = tiles)
partition_processes <- c(rep("treated",length(tiles)-1),"control")

# examine that its working out:
plot(partition,
     "treatment and control partition",
     xlab = "Easting",
     ylab = "Northing",
     lwd = 2
     )

pp_data$location_process <- sapply(1:dim(pp_data)[1],function(i){
  x <- pp_data$x[i]
  y <- pp_data$y[i]
  # Check if point is in treated state space
  if(inside.owin(x,y,treated_state_space)){
    return("treated")
  } else {
    return("control")
  }
})
pp_data$process[pp_data$t < TREATMENT_TIME] <- "control"
pp_data$background <- TRUE
pp_data$process <- pp_data$location_process
plot_pp(pp_data,partition)

# Jitter to break village-level geocoding ties. sd=500m is the approximate
# radius of a Gambian village catchment area; 1m was far too small and
# caused alpha to collapse onto near-duplicate locations.
pp_data$x <- pp_data$x + rnorm(nrow(pp_data), mean = 0, sd = 500)
pp_data$y <- pp_data$y + rnorm(nrow(pp_data), mean = 0, sd = 500)
plot_pp(pp_data,partition)
pp_data <- pp_data[inside.owin(pp_data$x, pp_data$y, win),]

# =========================================================
# 1. Background Rate Estimation (Pre-Treatment Kernel)
# =========================================================
# Establish "historical" risk using pre-treatment points across the whole window
X_pre <- ppp(
  x = pp_data$x[pp_data$t < TREATMENT_TIME],
  y = pp_data$y[pp_data$t < TREATMENT_TIME],
  window = win
)

# Use automated bandwidth selection (Diggle's method)
bw_sigma <- bw.diggle(X_pre)*20
lambda_space_gambia <- density(X_pre, sigma = bw_sigma, edge = TRUE, at = "pixels")

# Ensure positivity for numerical stability
min_non_zero <- min(lambda_space_gambia$v[lambda_space_gambia$v > 0], na.rm = TRUE)
lambda_space_gambia$v[lambda_space_gambia$v <= 0] <- min_non_zero

plot(lambda_space_gambia, main = "Gambia: Kernel-smoothed Background Rate (Pre-Treatment)")

# =========================================================
# 2. Regional Normalization (The NYC Logic)
# =========================================================
# This ensures that the background rate W integrates to the area of the window.
# This makes 'mu' comparable between regions with different base densities.

normalize_marks_gambia <- function(df_sub, win_sub, covariate_im, mark_name = "W") {
  # 1. Calculate the Raw Mass in this specific window
  cov_in_window <- covariate_im[win_sub, drop = FALSE]
  total_mass_raw <- integral.im(cov_in_window)

  # 2. Get the target Area
  target_area <- spatstat.geom::area(win_sub)

  # 3. Calculate Normalization Factor (Integral(W) == Area)
  norm_factor <- target_area / total_mass_raw

  # 4. Extract Raw values at the points
  vals_raw <- raster::extract(raster(covariate_im), df_sub[, c("x", "y")])
  vals_raw[is.na(vals_raw)] <- 0

  # 5. Apply Normalization
  df_sub[[mark_name]] <- vals_raw * norm_factor

  # Handle zeros/small values for numerical stability
  min_val <- min(df_sub[[mark_name]][df_sub[[mark_name]] > 0], na.rm = TRUE)
  if(is.infinite(min_val) | is.na(min_val)) min_val <- 1e-9
  df_sub[[mark_name]][df_sub[[mark_name]] <= 0] <- min_val

  return(list(new_df = df_sub, mass = total_mass_raw, norm = norm_factor))
}

# Split and Normalize into 4 groups (Pre/Post x Treated/Control)
res_pre_trtd  <- normalize_marks_gambia(pp_data %>% filter(t <  TREATMENT_TIME, location_process == "treated"), treated_state_space, lambda_space_gambia)
res_pre_ctrl  <- normalize_marks_gambia(pp_data %>% filter(t <  TREATMENT_TIME, location_process == "control"), control_state_space, lambda_space_gambia)
res_post_trtd <- normalize_marks_gambia(pp_data %>% filter(t >= TREATMENT_TIME, location_process == "treated"), treated_state_space, lambda_space_gambia)
res_post_ctrl <- normalize_marks_gambia(pp_data %>% filter(t >= TREATMENT_TIME, location_process == "control"), control_state_space, lambda_space_gambia)

# Combine for the SEM
pp_final <- rbind(res_pre_trtd$new_df,
                  res_pre_ctrl$new_df,
                  res_post_trtd$new_df,
                  res_post_ctrl$new_df)

# =========================================================
# 3. Initial Hawkes Fits (Pre-Treatment Baseline)
# =========================================================
# =========================================================
# 2. Perform the 4 Segmented Fits
# =========================================================

# Helper to run fit and return a named vector + diagnostics
# alpha ~ 1e-7 gives mean trigger dist ~ sqrt(1/alpha) ~ 3.2 km (meter-scale coords)
# beta ~ 0.05 gives mean trigger time ~ 20 days
run_full_fit <- function(df, win, label,
                         starts = list(
                           list(mu = 0.0001, alpha = 1e-7, beta = 0.05, K = 0.3),
                           list(mu = 0.001,  alpha = 1e-6, beta = 0.03, K = 0.2),
                           list(mu = 0.0005, alpha = 1e-8, beta = 0.1,  K = 0.4)
                         )) {
  cat("\nFitting:", label, "with", length(starts), "starting points...\n")
  best_fit <- NULL
  best_val <- -Inf
  for(i in seq_along(starts)) {
    fit <- tryCatch(
      fit_hawkes(
        params_init = starts[[i]],
        realiz = df,
        windowT = range(df$t),
        windowS = win,
        background_rate_var = 'W',
        maxit = 2000,
        use_fast = TRUE,
        method = "Nelder-Mead"
      ),
      error = function(e) { cat("  Start", i, "failed:", e$message, "\n"); NULL }
    )
    if(!is.null(fit) && fit$value > best_val) {
      best_val <- fit$value
      best_fit <- fit
      cat("  Start", i, "=> loglik =", fit$value,
          " K =", fit$par[4], " beta =", fit$par[3],
          " mean_time =", round(1/fit$par[3], 1), "days\n")
    }
  }
  pars <- best_fit$par
  cat("  BEST:", label, "=> K =", round(pars[4], 4),
      " beta =", round(pars[3], 6), " mean_time =", round(1/pars[3], 1),
      " alpha =", signif(pars[2], 3), " mean_dist =", round(sqrt(1/pars[2])), "m\n")
  return(data.frame(
    Segment = label,
    mu = pars['mu'],
    alpha = pars['alpha'],
    beta = pars['beta'],
    K = pars['K'],
    mean_dist_m = sqrt(1/pars['alpha']),
    mean_time_days = 1/pars['beta'],
    row.names = NULL
  ))
}

# Execute the fits
fit_pre_ctrl  <- run_full_fit(res_pre_ctrl$new_df,  control_state_space, "Pre-Control")
fit_pre_trtd  <- run_full_fit(res_pre_trtd$new_df,  treated_state_space, "Pre-Treated")
fit_post_ctrl <- run_full_fit(res_post_ctrl$new_df, control_state_space, "Post-Control")
fit_post_trtd <- run_full_fit(res_post_trtd$new_df, treated_state_space, "Post-Treated")

# =========================================================
# 3. Rbind the Results Dataframe for Inspection
# =========================================================

results_df <- rbind(
  fit_pre_ctrl,
  fit_pre_trtd,
  fit_post_ctrl,
  fit_post_trtd
)
results_df

# =========================================================
# 4. Adaptive SEM (Spillover Correction)
# =========================================================

# Define piecewise covariate lookup for the SEM
# It must apply the correct regional normalization factor based on location
covariate_lookup_func <- function(x, y) {
  raw_vals <- spatstat.geom::interp.im(lambda_space_gambia, list(x = x, y = y))
  raw_vals <- as.numeric(raw_vals)
  raw_vals[is.na(raw_vals)] <- 0

  is_trtd <- inside.owin(x, y, treated_state_space)
  w_vals <- raw_vals
  w_vals[is_trtd]  <- raw_vals[is_trtd] * res_pre_trtd$norm
  w_vals[!is_trtd] <- raw_vals[!is_trtd] * res_pre_ctrl$norm
  return(w_vals)
}

result_sem <- adaptive_SEM(
  pp_data = pp_final,
  partition = tess(tiles = list("control" = control_state_space,
                                "treated" = treated_state_space)),
  partition_processes = c("control", "treated"),
  statespace = win,
  time_window = c(0, STUDY_END),
  treatment_time = TREATMENT_TIME,
  hawkes_params_control = as.list(fit_pre_ctrl$par),
  hawkes_params_treated = as.list(fit_pre_ctrl$par),
  N_labellings = 10,
  N_iter = 1,
  covariate_lookup = covariate_lookup_func,
  background_rate_var = 'W',
  adaptive_control = list(
    update_control_params = TRUE,
    param_update_cadence = 10,
    proposal_update_cadence = 1,
    state_spaces = list(control_state_space,treated_state_space),
    iter  = 100,
    n_props = 100,
    change_factor = 0.1,
    include_starting_data = F,
    update_starting_data = T,
    verbose = T
  )
)

sem_pars_control <- result_sem$hawkes_params_control
sem_pars_treated <- result_sem$hawkes_params_treated

tmp <- data.frame(
  Segment = c("SEM-Adaptive-Control", "SEM-Adaptive-Treated"),
  mu    = c(sem_pars_control[["mu"]],    sem_pars_treated[["mu"]]),
  alpha = c(sem_pars_control[["alpha"]], sem_pars_treated[["alpha"]]),
  beta  = c(sem_pars_control[["beta"]],  sem_pars_treated[["beta"]]),
  K     = c(sem_pars_control[["K"]],     sem_pars_treated[["K"]])
)

# 3. Combine with your existing results_df (taking only the first 5 columns)
results_all <- rbind(results_df[,1:5], tmp)

# 4. Clean and Round: Convert columns 2:5 to numeric then round
# Using mutate(across(...)) is much safer than as.numeric on a whole block
results_all <- results_all %>%
  mutate(across(2:5, ~ round(as.numeric(.), 4)))

# 5. Add interpretation columns (Radiuses) for the SEM rows to match your inspection_results
results_all <- results_all %>%
  mutate(
    mean_dist_m = round(sqrt(1/alpha), 0),
    mean_time_days = round(1/beta, 1)
  )
print(results_all)



# =========================================================
# 5. Counterfactual Savings Calculation
# =========================================================
# Project the control model into the treated region's background risk
duration_post <- STUDY_END - TREATMENT_TIME
total_mass_gambia <- res_pre_trtd$mass + res_pre_ctrl$mass


estimate_hawkes_savings <- function(par_control,
                                    par_treated,
                                    mass_origin_control,
                                    mass_origin_treated,
                                    mass_target,
                                    duration,
                                    verbose = TRUE) {

  # Helper to safely extract K and mu
  get_vals <- function(p) {
    if(is.list(p)) return(list(mu=p$mu, K=p$K))
    if(!is.null(names(p))) return(list(mu=p[["mu"]], K=p[["K"]]))
    return(list(mu=p[1], K=p[4])) # Assume numeric vector order: mu, alpha, beta, K
  }

  p_ctrl <- get_vals(par_control)
  p_trtd <- get_vals(par_treated)

  # 1. Normalize mu to get the "per unit mass" background rate
  # rate = mu_fitted / mass_fitted_region
  rate_control_base <- p_ctrl$mu / mass_origin_control
  rate_treated_base <- p_trtd$mu / mass_origin_treated

  # 2. Project Background to the Target Area (Whole City)
  # mu_projected = rate * total_mass * duration
  bg_projected_control <- rate_control_base * mass_target * duration
  bg_projected_treated <- rate_treated_base * mass_target * duration

  # 3. Calculate Total Expected Events (Branching Process)
  # Total = Background / (1 - K)

  # Safety check for critical/super-critical regimes
  calc_total <- function(bg, K) {
    if(K >= 1) return(Inf)
    return(bg * (1 / (1 - K)))
  }

  total_exp_control <- calc_total(bg_projected_control, p_ctrl$K)
  total_exp_treated <- calc_total(bg_projected_treated, p_trtd$K)

  # 4. Calculate Savings
  savings <- total_exp_control - total_exp_treated
  pct_reduction <- savings / total_exp_control

  # 5. Output
  if(verbose) {
    cat("--- Counterfactual Annual Projection---\n")
    cat(sprintf("Projection Duration:                %.1f days\n", duration))
    cat(sprintf("Expected Sightings (No Intervention): %s\n", format(round(total_exp_control), big.mark=",")))
    cat(sprintf("Expected Sightings (With Intervention): %s\n", format(round(total_exp_treated), big.mark=",")))
    cat(sprintf("Total Sightings Prevented:              %s\n", format(round(savings), big.mark=",")))
    cat(sprintf("Percentage Reduction:                   %.2f%%\n", pct_reduction * 100))
    cat("----------------------------------------------\n")
  }

  return(list(
    expected_control = total_exp_control,
    expected_treated = total_exp_treated,
    savings = savings,
    pct_reduction = pct_reduction,
    rates = c(control_rate = rate_control_base, treated_rate = rate_treated_base)
  ))
}
# ==========
# SEM Report
# ==========

savings_report <- estimate_hawkes_savings(
  par_control = result_sem$adaptive$control_par[[length(result_sem$adaptive$control_par)]],
  par_treated = result_sem$hawkes_params_treated,
  mass_origin_control = res_pre_ctrl$mass,
  mass_origin_treated = res_pre_trtd$mass,
  mass_target = total_mass_gambia,
  duration = duration_post,
  verbose = T
)

# ============
# Naive Report
# ============

savings_report_naive <- estimate_hawkes_savings(
  par_control = fit_post_ctrl %>% dplyr::select(mu, K) %>% as.list(),
  par_treated = fit_post_trtd %>% dplyr::select(mu, K) %>% as.list(),
  mass_origin_control = res_pre_ctrl$mass,
  mass_origin_treated = res_pre_trtd$mass,
  mass_target = total_mass_gambia,
  duration = duration_post,
  verbose = T
)
