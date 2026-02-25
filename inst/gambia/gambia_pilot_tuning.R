library(PPDisentangle)
library(spatstat)
library(ggplot2)
library(dplyr)
library(data.table)
library(terra)
library(raster)
library(lubridate)

# =========================================================
# 1. Data Setup (Copied from gambia_fit_inhom.R)
# =========================================================
TREATMENT_DATE <- as_datetime("2011-05-01")
HEALTH_CENTER_DIST <- 10000 # 10 km

gambia_data <- dget(system.file("extdata/gambia/DISCLEAN.txt", package = "PPDisentangle"))
gambia_data <- gambia_data %>% filter(IPDEvent == "Y")
TREATMENT_TIME <- as.numeric(as.Date(TREATMENT_DATE) - min(as.Date(gambia_data$date_IPD_Pne_Event)))
STUDY_END <- as.numeric(max(as.Date(gambia_data$date_IPD_Pne_Event)) - min(as.Date(gambia_data$date_IPD_Pne_Event)))

win <- dget(system.file("extdata/gambia/GAMWIN_PROJ.txt", package = "PPDisentangle"))
S  <- diag(c(1000, 1000))
t  <- c(0, 0)
win <- affine.owin(win, mat = S, vec = t)
hc <- dget(system.file("extdata/gambia/hc.txt", package = "PPDisentangle"))
hc <- affine(hc, mat = S, vec = t)

v <- vect(gambia_data, geom = c("long","lat"), crs = "EPSG:4326")
v_utm <- project(v, "EPSG:32628")
coords <- crds(v_utm)
gambia_data$easting <- coords[,1]
gambia_data$northing <- coords[,2]

pp_data <- data.frame(
  x = gambia_data$easting,
  y = gambia_data$northing,
  t = as.numeric(as.Date(gambia_data$date_IPD_Pne_Event)) - min(as.numeric(as.Date(gambia_data$date_IPD_Pne_Event))),
  IPD = gambia_data$IPDEvent == "Y"
)
pp_data <- pp_data[inside.owin(pp_data$x, pp_data$y, win),]
pp_data$t <- pp_data$t - min(pp_data$t)
pp_data <- pp_data %>% filter(t <= STUDY_END)

dmap <- distmap(hc)
bands <- cut(as.im(dmap), breaks = c(0,HEALTH_CENTER_DIST,Inf), labels = c("close","far"))
discs_list <- lapply(seq_len(npoints(hc)), function(i) {
  disc(radius = HEALTH_CENTER_DIST, centre = c(hc$x[i], hc$y[i]), npoly = 200)
})
disc_list_double <- lapply(seq_len(npoints(hc)), function(i) {
  disc(radius = HEALTH_CENTER_DIST * 2, centre = c(hc$x[i], hc$y[i]), npoly = 200)
})
treated_state_space <- Reduce(union.owin, discs_list)
win_expanded <- Reduce(union.owin, disc_list_double)
control_state_space <- setminus.owin(win_expanded, treated_state_space)
pp_data <- pp_data[inside.owin(pp_data$x, pp_data$y, win_expanded),]

partition <- tess(tiles = list("control" = control_state_space, "treated" = treated_state_space))
partition_processes <- c("control", "treated")

pp_data$location_process <- ifelse(inside.owin(pp_data$x, pp_data$y, treated_state_space), "treated", "control")
pp_data$process <- pp_data$location_process
pp_data$inferred_process <- pp_data$location_process
pp_data$background <- TRUE

# Jitter
pp_data$x <- pp_data$x + rnorm(nrow(pp_data), mean = 0, sd = 1)
pp_data$y <- pp_data$y + rnorm(nrow(pp_data), mean = 0, sd = 1)
pp_data <- pp_data[inside.owin(pp_data$x, pp_data$y, win_expanded),]

# Background rate
X_pre <- ppp(x = pp_data$x[pp_data$t < TREATMENT_TIME], y = pp_data$y[pp_data$t < TREATMENT_TIME], window = win_expanded)
bw_sigma <- bw.diggle(X_pre)*20
lambda_space_gambia <- density(X_pre, sigma = bw_sigma, edge = TRUE, at = "pixels")

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

res_pre_trtd  <- normalize_marks_gambia(pp_data %>% filter(t <  TREATMENT_TIME, location_process == "treated"), treated_state_space, lambda_space_gambia)
res_pre_ctrl  <- normalize_marks_gambia(pp_data %>% filter(t <  TREATMENT_TIME, location_process == "control"), control_state_space, lambda_space_gambia)
res_post_trtd <- normalize_marks_gambia(pp_data %>% filter(t >= TREATMENT_TIME, location_process == "treated"), treated_state_space, lambda_space_gambia)
res_post_ctrl <- normalize_marks_gambia(pp_data %>% filter(t >= TREATMENT_TIME, location_process == "control"), control_state_space, lambda_space_gambia)

pp_final <- rbind(res_pre_trtd$new_df, res_pre_ctrl$new_df, res_post_trtd$new_df, res_post_ctrl$new_df)

# Initial fits
fit_init_ctrl <- fit_hawkes(params_init = list(mu = 0.0001, alpha = 0.01, beta = 0.1, K = 0.2),
                            realiz = res_pre_ctrl$new_df, windowT = range(res_pre_ctrl$new_df$t),
                            windowS = control_state_space, background_rate_var = 'W')

# =========================================================
# 2. Pilot Run Procedure
# =========================================================
cat("Starting pilot runs to tune change_factor...\n")

change_factors <- c(0.01, 0.05, 0.1, 0.2)
n_pilot_iter <- 50
pilot_results <- list()

for (cf in change_factors) {
  cat(paste0("Testing change_factor = ", cf, "...\n"))
  
  # Run em_style_labelling
  res <- em_style_labelling(
    pp_data = pp_final,
    partition = partition,
    partition_processes = partition_processes,
    statespace = win_expanded,
    time_window = c(0, STUDY_END),
    treatment_time = TREATMENT_TIME,
    hawkes_params_control = as.list(fit_init_ctrl$par),
    hawkes_params_treated = as.list(fit_init_ctrl$par),
    update_control_params = TRUE,
    param_update_cadence = 10,
    proposal_update_cadence = 1,
    update_starting_data = TRUE,
    include_starting_data = FALSE,
    metric_name = "post_likelihood",
    iter = n_pilot_iter,
    n_props = 20, # Smaller for speed in pilot
    change_factor = cf,
    verbose = FALSE,
    background_rate_var = 'W'
  )
  
  pilot_results[[as.character(cf)]] <- data.frame(
    iteration = 1:n_pilot_iter,
    average_flips = res$average_flips,
    max_metric_flips = res$max_metric_flips,
    metric = res$metrics,
    change_factor = as.factor(cf)
  )
}

pilot_df <- do.call(rbind, pilot_results)

# Plot Average Flips
p1 <- ggplot(pilot_df, aes(x = iteration, y = average_flips, color = change_factor)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  labs(title = "Average Flips per Iteration (Pilot Run)",
       subtitle = "Goal: ~1 flip per iteration for long-term stability",
       x = "Iteration", y = "Average Flips") +
  theme_minimal()

# Plot Max Metric Flips (the flips in the chosen labelling)
p2 <- ggplot(pilot_df, aes(x = iteration, y = max_metric_flips, color = change_factor)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  labs(title = "Max Metric Flips per Iteration (Pilot Run)",
       subtitle = "Flips in the best proposal chosen for next iteration",
       x = "Iteration", y = "Max Metric Flips") +
  theme_minimal()

# Plot Metric (Log-Likelihood)
p3 <- ggplot(pilot_df, aes(x = iteration, y = metric, color = change_factor)) +
  geom_line() +
  labs(title = "Log-Likelihood (Metric) over Iterations",
       subtitle = "Higher is better; check for convergence",
       x = "Iteration", y = "Metric") +
  theme_minimal()

# Save plots
ggsave("inst/gambia/pilot_average_flips.png", p1, width = 8, height = 6)
ggsave("inst/gambia/pilot_max_flips.png", p2, width = 8, height = 6)
ggsave("inst/gambia/pilot_metric.png", p3, width = 8, height = 6)

# Suggest a change_factor
final_flips <- pilot_df %>% filter(iteration == n_pilot_iter)
suggested_cf <- final_flips$change_factor[which.min(abs(final_flips$average_flips - 1))]

cat("\n--- Pilot Run Summary ---\n")
print(final_flips)
cat(paste0("\nSuggested change_factor (closest to 1 flip at end): ", suggested_cf, "\n"))
cat("Check inst/gambia/pilot_average_flips.png and pilot_max_flips.png for details.\n")
