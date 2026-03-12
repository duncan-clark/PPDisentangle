#!/usr/bin/env Rscript
# ============================================================================
# Oklahoma Induced Seismicity Analysis — County-Based Four-Way ETAS
#
# Treatment: OCC directive AOI_20150318 (wastewater injection reduction)
# Tessellation: Oklahoma county boundaries (77 counties)
# State space: Oklahoma state boundary
#
# Four fits:
#   A. Naive independent ETAS   — separate MLE on location-labeled data
#   B. Naive bivariate ETAS     — joint MLE on location-labeled data
#   C. SEM independent ETAS     — adaptive_SEM with model_type="etas"
#   D. SEM bivariate ETAS       — adaptive_SEM with model_type="etas_bivariate"
#
# Each fit produces an ATE: expected per-county savings from treatment.
#
# Usage:
#   Rscript oklahoma_analysis.R
#   Rscript oklahoma_analysis.R --test
# ============================================================================

suppressPackageStartupMessages({
  library(PPDisentangle)
  library(spatstat)
  library(sf)
  library(tigris)
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
TEST_MODE <- "--test" %in% args

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
OUT_DIR  <- file.path(SCRIPT_DIR, "output")
PLOT_DIR <- file.path(OUT_DIR, "plots")
for (d in c(OUT_DIR, PLOT_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ---- Configuration ----
DATA_DIR   <- file.path(SCRIPT_DIR, "oklahoma_induced_seismicity_data_regional20150318")
ETAS_M0    <- 2.5
BETA_GR    <- 2.3
CRS_PROJ   <- 5070

VANILLA_MAXIT <- if (TEST_MODE) 500 else 5000
VANILLA_STARTS <- list(
  list(mu = 1.0, A = 0.2, alpha_m = 0.8, c = 0.05, p = 1.2,
       D = 5.0, gamma = 0.5, q = 1.5),
  list(mu = 2.0, A = 0.1, alpha_m = 0.5, c = 0.1,  p = 1.3,
       D = 2.0, gamma = 0.3, q = 2.0)
)

SEM_N_LABELLINGS  <- if (TEST_MODE) 5  else 10
SEM_N_ITER        <- if (TEST_MODE) 1  else 1
SEM_INNER_ITER    <- if (TEST_MODE) 10 else 100
SEM_INNER_PROPS   <- if (TEST_MODE) 5  else 10
SEM_CHANGE_FACTOR <- 0.01
SEM_PARAM_UPDATE  <- 10
SEM_OUTER_MAXIT       <- if (TEST_MODE) 200 else 1000
SEM_OUTER_MAXIT_BIV   <- if (TEST_MODE) 100 else 200

FIXED_STRUCTURAL <- list(c = 0.05, p = 1.2, D = 5.0, gamma = 0.5, q = 1.5)

ATE_N_SIMS    <- if (TEST_MODE) 20  else 200
ATE_N_TAU_I   <- if (TEST_MODE) 5   else 30
ATE_N_TAU_SIM <- if (TEST_MODE) 5   else 20

etas_names <- c("mu", "A", "alpha_m", "c", "p", "D", "gamma", "q")

cat("=== Oklahoma County-Based ETAS Analysis (4-Way Comparison) ===\n")
cat(sprintf("Mode: %s | SEM iters: %d | Change factor: %.3f\n",
            if (TEST_MODE) "TEST" else "FULL", SEM_N_ITER, SEM_CHANGE_FACTOR))

# ============================================================================
# 1. Data loading
# ============================================================================
cat("\n--- Step 1: Data loading ---\n")

if (!dir.exists(DATA_DIR)) {
  cat("Data not found — running data preparation script...\n")
  source(file.path(SCRIPT_DIR, "Oklahoma_data_and_viz.R"), local = TRUE)
}

ev_all <- fread(file.path(DATA_DIR, "events_all.csv"))
meta   <- jsonlite::fromJSON(readLines(file.path(DATA_DIR, "metadata.json")))

t_star_utc <- as.POSIXct(meta$design$t_star_utc, tz = "UTC")
ev_all[, time_utc := as.POSIXct(time_utc, tz = "UTC")]
ev_all[, t_days := as.numeric(difftime(time_utc, t_star_utc, units = "days"))]
ev_all[, x_km := x_m / 1000]
ev_all[, y_km := y_m / 1000]
if (!"mag" %in% names(ev_all)) stop("No 'mag' column in events data")
ev_all <- ev_all[mag >= ETAS_M0]

post_end_days <- as.numeric(difftime(
  as.POSIXct(meta$design$post_end_utc, tz = "UTC"), t_star_utc, units = "days"))

# ============================================================================
# 2. Oklahoma county tessellation
# ============================================================================
cat("\n--- Step 2: Building county tessellation ---\n")

options(tigris_use_cache = TRUE)
counties_sf <- counties(state = "OK", cb = TRUE, year = 2022)
counties_sf <- st_transform(counties_sf, CRS_PROJ)
counties_sf <- st_make_valid(counties_sf)

ok_boundary <- st_union(counties_sf)
ok_boundary <- st_make_valid(ok_boundary)

bb <- st_bbox(ok_boundary)
win_km <- owin(xrange = c(bb["xmin"], bb["xmax"]) / 1000,
               yrange = c(bb["ymin"], bb["ymax"]) / 1000)

# Convert county polygons to spatstat owin tiles
county_owins <- lapply(seq_len(nrow(counties_sf)), function(i) {
  geom <- st_geometry(counties_sf[i, ])
  coords_list <- st_coordinates(geom)
  x_km <- coords_list[, 1] / 1000
  y_km <- coords_list[, 2] / 1000
  tryCatch(
    owin(poly = list(x = rev(x_km), y = rev(y_km))),
    error = function(e) {
      tryCatch(
        owin(poly = list(x = x_km, y = y_km)),
        error = function(e2) NULL
      )
    }
  )
})

valid_idx <- !sapply(county_owins, is.null)
if (sum(valid_idx) < 50) {
  cat("  Warning: only", sum(valid_idx), "of", nrow(counties_sf), "counties converted.\n")
}

county_owins_valid <- county_owins[valid_idx]
counties_sf_valid  <- counties_sf[valid_idx, ]

names(county_owins_valid) <- counties_sf_valid$NAME
partition <- tess(tiles = county_owins_valid, window = win_km)

cat(sprintf("  Counties in tessellation: %d / %d\n",
            partition$n, nrow(counties_sf)))

# Assign treatment: county centroid inside OCC AOI
aoi_path <- file.path(DATA_DIR, "occ_aoi_layer_2.geojson")
aoi_sf <- st_read(aoi_path, quiet = TRUE)
aoi_sf <- st_transform(aoi_sf, CRS_PROJ)
aoi_sf <- st_make_valid(aoi_sf)
aoi_union <- st_union(aoi_sf)

county_centroids <- st_centroid(counties_sf_valid)
inside_aoi <- lengths(st_within(county_centroids, aoi_union)) > 0

partition_processes <- ifelse(inside_aoi, "treated", "control")
names(partition_processes) <- counties_sf_valid$NAME

treated_idx <- partition_processes == "treated"
treated_names <- names(partition_processes)[treated_idx]
control_ss <- as.owin(partition[!treated_idx])
treated_ss <- as.owin(partition[treated_idx])
state_spaces <- list(control = control_ss, treated = treated_ss)

cat(sprintf("  Treated counties: %d (%s)\n",
            sum(treated_idx),
            paste(treated_names, collapse = ", ")))
cat(sprintf("  Control counties: %d\n", sum(!treated_idx)))

# ---- Assign events to counties ----
assign_county <- function(df) {
  ti <- as.integer(tileindex(df$x, df$y, partition))
  df$location_process <- ifelse(is.na(ti), NA_character_,
                                partition_processes[pmin(pmax(ti, 1), partition$n)])
  df$W <- 1.0
  df$n <- nrow(df)
  df$background <- TRUE
  df
}

pp_pre  <- assign_county(as.data.frame(ev_all[t_days < 0,
  .(x = x_km, y = y_km, t = t_days, mag = mag)]))
pp_post <- assign_county(as.data.frame(ev_all[t_days >= 0 & t_days <= post_end_days,
  .(x = x_km, y = y_km, t = t_days, mag = mag)]))

# Drop events outside Oklahoma
pp_pre  <- pp_pre[!is.na(pp_pre$location_process), ]
pp_post <- pp_post[!is.na(pp_post$location_process), ]

pp_pre$process  <- "control"
pp_post$process <- pp_post$location_process
pp_pre$inferred_process  <- "control"
pp_post$inferred_process <- pp_post$location_process
pp_all <- rbind(pp_pre, pp_post)
pp_all <- pp_all[order(pp_all$t), ]

windowT_post <- c(0, post_end_days)

cat(sprintf("  Events in OK: pre=%d, post=%d, total=%d\n",
            nrow(pp_pre), nrow(pp_post), nrow(pp_all)))
cat(sprintf("  Post-treatment: %d control, %d treated\n",
            sum(pp_post$location_process == "control"),
            sum(pp_post$location_process == "treated")))

# ============================================================================
# 3. Partition and point pattern plots
# ============================================================================
cat("\n--- Step 3: Plots ---\n")

tryCatch({
  counties_plot <- counties_sf_valid
  counties_plot$treatment <- partition_processes

  p_partition <- ggplot() +
    geom_sf(data = counties_plot, aes(fill = treatment),
            color = "grey40", linewidth = 0.4) +
    scale_fill_manual(values = c(control = "#deebf7", treated = "#fc9272"),
                      name = "Assignment") +
    geom_point(data = pp_post,
               aes(x = x * 1000, y = y * 1000),
               size = 0.5, alpha = 0.5) +
    labs(title = "Oklahoma County Partition",
         subtitle = "OCC directive AOI_20150318 (treated = red)",
         x = "Easting (m)", y = "Northing (m)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  ggsave(file.path(PLOT_DIR, "partition_map.png"), p_partition,
         width = 12, height = 7, dpi = 150)
  cat("  Saved partition_map.png\n")
}, error = function(e) cat("  Partition plot error:", e$message, "\n"))

tryCatch({
  p_pre <- ggplot(pp_pre, aes(x = x, y = y, alpha = t)) +
    geom_point(color = "#377eb8", size = 1.5, shape = 21, fill = "#377eb8") +
    scale_alpha_continuous(name = "Time (days)", range = c(0.15, 0.8)) +
    labs(title = "Pre-treatment Events", x = "X (km)", y = "Y (km)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(PLOT_DIR, "pp_pre_treatment.png"), p_pre,
         width = 12, height = 7, dpi = 150)
  cat("  Saved pp_pre_treatment.png\n")

  pp_post_plot <- pp_post
  pp_post_plot$Process <- factor(pp_post_plot$location_process,
                                 levels = c("control", "treated"))
  p_post <- ggplot(pp_post_plot, aes(x = x, y = y, fill = Process, alpha = t)) +
    geom_point(size = 1.8, shape = 21, stroke = 0.3) +
    scale_fill_manual(values = c(control = "#377eb8", treated = "#e41a1c")) +
    scale_alpha_continuous(name = "Time (days)", range = c(0.3, 1)) +
    labs(title = "Post-treatment Events (Location Labels)",
         x = "X (km)", y = "Y (km)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(PLOT_DIR, "pp_post_treatment.png"), p_post,
         width = 12, height = 7, dpi = 150)
  cat("  Saved pp_post_treatment.png\n")
}, error = function(e) cat("  PP plot error:", e$message, "\n"))

# ============================================================================
# 4A. Fit A: Naive independent ETAS
# ============================================================================
cat("\n--- Fit A: Naive independent ETAS ---\n")

fit_best_indep <- function(realiz, zbr, starts, maxit) {
  best_fit <- NULL; best_val <- -Inf
  for (s in starts) {
    fit <- tryCatch(
      fit_etas(params_init = s, realiz = realiz, windowT = windowT_post,
               windowS = win_km, m0 = ETAS_M0, maxit = maxit,
               fixed_params = FIXED_STRUCTURAL, zero_background_region = zbr),
      error = function(e) NULL)
    if (!is.null(fit) && is.finite(fit$value) && fit$value > best_val) {
      best_fit <- fit; best_val <- fit$value
    }
  }
  best_fit
}

naive_control <- pp_post[pp_post$location_process == "control", ]
naive_treated <- pp_post[pp_post$location_process == "treated", ]

cat("  Fitting control...\n")
fitA_ctrl <- fit_best_indep(naive_control, treated_ss, VANILLA_STARTS, VANILLA_MAXIT)
cat("  Fitting treated...\n")
fitA_treat <- fit_best_indep(naive_treated, control_ss, VANILLA_STARTS, VANILLA_MAXIT)

A_ctrl <- if (!is.null(fitA_ctrl)) as.list(fitA_ctrl$par) else VANILLA_STARTS[[1]]
A_treat <- if (!is.null(fitA_treat)) as.list(fitA_treat$par) else VANILLA_STARTS[[1]]
names(A_ctrl) <- etas_names; names(A_treat) <- etas_names

cat("  Control:", paste(etas_names, round(unlist(A_ctrl), 4), sep = "=", collapse = ", "), "\n")
cat("  Treated:", paste(etas_names, round(unlist(A_treat), 4), sep = "=", collapse = ", "), "\n")

# ============================================================================
# 4B. Fit B: Naive bivariate ETAS
# ============================================================================
cat("\n--- Fit B: Naive bivariate ETAS ---\n")

biv_init <- init_bivariate_from_independent(A_ctrl, A_treat)
biv_names <- names(biv_init)

fitB <- tryCatch({
  fit_etas_bivariate(
    params_init = biv_init, realiz = pp_post,
    windowT = windowT_post, windowS = win_km, m0 = ETAS_M0,
    control_state_space = control_ss, treated_state_space = treated_ss,
    maxit = VANILLA_MAXIT, fixed_params = FIXED_STRUCTURAL, trace = 1
  )
}, error = function(e) { cat("  Bivariate fit error:", e$message, "\n"); NULL })

B_params <- if (!is.null(fitB)) fitB$par else biv_init
B_loglik <- if (!is.null(fitB)) fitB$value else NA

cat("  Params:", paste(biv_names, round(B_params, 4), sep = "=", collapse = ", "), "\n")

# ============================================================================
# 4C. Fit C: SEM independent ETAS
# ============================================================================
cat("\n--- Fit C: SEM independent ETAS ---\n")

semC <- tryCatch({
  adaptive_SEM(
    pp_data = pp_all, partition = partition,
    partition_processes = partition_processes,
    statespace = win_km, time_window = windowT_post, treatment_time = 0,
    hawkes_params_control = A_ctrl, hawkes_params_treated = A_treat,
    N_labellings = SEM_N_LABELLINGS, N_iter = SEM_N_ITER, verbose = TRUE,
    model_type = "etas",
    adaptive_control = list(
      param_update_cadence = SEM_PARAM_UPDATE,
      proposal_update_cadence = 1,
      state_spaces = state_spaces,
      iter = SEM_INNER_ITER, n_props = SEM_INNER_PROPS,
      change_factor = SEM_CHANGE_FACTOR, verbose = FALSE,
      update_starting_data = TRUE, include_starting_data = TRUE,
      update_control_params = TRUE, fixed_params = FIXED_STRUCTURAL,
      proposal_method = "simulation",
      outer_maxit = SEM_OUTER_MAXIT
    ),
    m0 = ETAS_M0, beta_gr = BETA_GR
  )
}, error = function(e) { cat("  SEM-indep error:", e$message, "\n"); NULL })

if (!is.null(semC)) {
  C_ctrl  <- semC$hawkes_params_control
  C_treat <- semC$hawkes_params_treated
  if (is.null(names(C_ctrl)))  names(C_ctrl)  <- etas_names
  if (is.null(names(C_treat))) names(C_treat) <- etas_names
  cat("  SEM-indep completed in", round(semC$time, 1), "s\n")
  cat("  Control:", paste(etas_names, round(unlist(C_ctrl), 4), sep = "=", collapse = ", "), "\n")
  cat("  Treated:", paste(etas_names, round(unlist(C_treat), 4), sep = "=", collapse = ", "), "\n")
} else {
  C_ctrl <- A_ctrl; C_treat <- A_treat
  cat("  SEM-indep failed, falling back to naive.\n")
}

# ============================================================================
# 4D. Fit D: SEM bivariate ETAS
# ============================================================================
cat("\n--- Fit D: SEM bivariate ETAS ---\n")

biv_init_D <- init_bivariate_from_independent(A_ctrl, A_treat)
biv_fixed <- list(c = 0.05, p = 1.2, D = 5.0, gamma = 0.5, q = 1.5)

semD <- tryCatch({
  adaptive_SEM(
    pp_data = pp_all, partition = partition,
    partition_processes = partition_processes,
    statespace = win_km, time_window = windowT_post, treatment_time = 0,
    hawkes_params_control = A_ctrl, hawkes_params_treated = A_treat,
    N_labellings = SEM_N_LABELLINGS, N_iter = SEM_N_ITER, verbose = TRUE,
    model_type = "etas_bivariate",
    adaptive_control = list(
      param_update_cadence = SEM_PARAM_UPDATE,
      proposal_update_cadence = 1,
      state_spaces = state_spaces,
      iter = SEM_INNER_ITER, n_props = SEM_INNER_PROPS,
      change_factor = SEM_CHANGE_FACTOR, verbose = FALSE,
      update_starting_data = TRUE, include_starting_data = TRUE,
      update_control_params = TRUE, fixed_params = biv_fixed,
      proposal_method = "simulation",
      outer_maxit = SEM_OUTER_MAXIT, outer_maxit_biv = SEM_OUTER_MAXIT_BIV
    ),
    m0 = ETAS_M0, beta_gr = BETA_GR,
    etas_bivariate_params = biv_init_D
  )
}, error = function(e) { cat("  SEM-biv error:", e$message, "\n"); NULL })

if (!is.null(semD)) {
  D_params <- semD$etas_bivariate_params
  D_ctrl   <- semD$hawkes_params_control
  D_treat  <- semD$hawkes_params_treated
  cat("  SEM-biv completed in", round(semD$time, 1), "s\n")
  cat("  Bivariate params:", paste(biv_names, round(D_params, 4), sep = "=", collapse = ", "), "\n")
} else {
  D_params <- biv_init_D
  D_ctrl <- A_ctrl; D_treat <- A_treat
  cat("  SEM-biv failed, falling back to naive.\n")
}

# ============================================================================
# 5. SEM diagnostic plots
# ============================================================================
cat("\n--- Step 5: SEM plots ---\n")

for (nm in list(list(res = semC, label = "indep", title = "SEM Independent"),
                list(res = semD, label = "biv",   title = "SEM Bivariate"))) {
  tryCatch({
    if (!is.null(nm$res)) {
      sem_lab <- nm$res$adaptive$adaptive_labelling
      sem_post <- if (!is.null(sem_lab)) sem_lab[sem_lab$t >= 0, ] else NULL
      if (!is.null(sem_post) && nrow(sem_post) > 0) {
        sem_post$Process <- factor(sem_post$inferred_process,
                                   levels = c("control", "treated"))
        p_s <- ggplot(sem_post, aes(x = x, y = y, fill = Process, alpha = t)) +
          geom_point(size = 1.8, shape = 21, stroke = 0.3) +
          scale_fill_manual(values = c(control = "#377eb8", treated = "#e41a1c")) +
          scale_alpha_continuous(name = "Time (days)", range = c(0.3, 1)) +
          labs(title = paste0("Post-treatment (", nm$title, " labels)"),
               x = "X (km)", y = "Y (km)") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        ggsave(file.path(PLOT_DIR, sprintf("pp_post_sem_%s.png", nm$label)),
               p_s, width = 12, height = 7, dpi = 150)
        cat(sprintf("  Saved pp_post_sem_%s.png\n", nm$label))
      }
      p_f <- plot_flips(nm$res)
      ggsave(file.path(PLOT_DIR, sprintf("sem_flips_%s.png", nm$label)),
             p_f, width = 8, height = 5, dpi = 150)
      cat(sprintf("  Saved sem_flips_%s.png\n", nm$label))
    }
  }, error = function(e) cat(sprintf("  %s plot error: %s\n", nm$title, e$message)))
}

# ============================================================================
# 6. ATE estimation
# ============================================================================
cat("\n--- Step 6: ATE estimation ---\n")

treated_partitions <- tilenames(partition)[treated_idx]
pp_post_sem_C <- if (!is.null(semC)) semC$adaptive$adaptive_labelling else NULL
pp_post_sem_C <- if (!is.null(pp_post_sem_C)) pp_post_sem_C[pp_post_sem_C$t >= 0, ] else pp_post
pp_post_sem_D <- if (!is.null(semD)) semD$adaptive$adaptive_labelling else NULL
pp_post_sem_D <- if (!is.null(pp_post_sem_D)) pp_post_sem_D[pp_post_sem_D$t >= 0, ] else pp_post

cat("  Computing ATE for Fit A (naive independent)...\n")
ate_A <- tryCatch(
  ATE_estim_etas(
    statespace = win_km, partition = partition,
    observed_data = pp_post, treated_partitions = treated_partitions,
    etas_params = list(control = A_ctrl, treated = A_treat),
    n_sims = ATE_N_SIMS, n_tau_i = ATE_N_TAU_I, n_tau_sims = ATE_N_TAU_SIM,
    windowT = windowT_post, windowS = win_km,
    m0 = ETAS_M0, beta_gr = BETA_GR, fixed_params = FIXED_STRUCTURAL
  ),
  error = function(e) { cat("    ATE-A error:", e$message, "\n"); NULL }
)

cat("  Computing ATE for Fit C (SEM independent)...\n")
ate_C <- tryCatch(
  ATE_estim_etas(
    statespace = win_km, partition = partition,
    observed_data = pp_post_sem_C, treated_partitions = treated_partitions,
    etas_params = list(control = C_ctrl, treated = C_treat),
    n_sims = ATE_N_SIMS, n_tau_i = ATE_N_TAU_I, n_tau_sims = ATE_N_TAU_SIM,
    windowT = windowT_post, windowS = win_km,
    m0 = ETAS_M0, beta_gr = BETA_GR, fixed_params = FIXED_STRUCTURAL
  ),
  error = function(e) { cat("    ATE-C error:", e$message, "\n"); NULL }
)

# For bivariate models, extract marginal params for ATE
B_ctrl_list <- as.list(c(mu = B_params[["mu_0"]], A = B_params[["A_00"]],
  alpha_m = B_params[["alpha_m_00"]], unlist(FIXED_STRUCTURAL)))
B_treat_list <- as.list(c(mu = B_params[["mu_1"]], A = B_params[["A_11"]],
  alpha_m = B_params[["alpha_m_11"]], unlist(FIXED_STRUCTURAL)))

cat("  Computing ATE for Fit B (naive bivariate)...\n")
ate_B <- tryCatch(
  ATE_estim_etas(
    statespace = win_km, partition = partition,
    observed_data = pp_post, treated_partitions = treated_partitions,
    etas_params = list(control = B_ctrl_list, treated = B_treat_list),
    n_sims = ATE_N_SIMS, n_tau_i = ATE_N_TAU_I, n_tau_sims = ATE_N_TAU_SIM,
    windowT = windowT_post, windowS = win_km,
    m0 = ETAS_M0, beta_gr = BETA_GR, fixed_params = FIXED_STRUCTURAL
  ),
  error = function(e) { cat("    ATE-B error:", e$message, "\n"); NULL }
)

D_ctrl_list <- as.list(c(mu = D_params[["mu_0"]], A = D_params[["A_00"]],
  alpha_m = D_params[["alpha_m_00"]], unlist(FIXED_STRUCTURAL)))
D_treat_list <- as.list(c(mu = D_params[["mu_1"]], A = D_params[["A_11"]],
  alpha_m = D_params[["alpha_m_11"]], unlist(FIXED_STRUCTURAL)))

cat("  Computing ATE for Fit D (SEM bivariate)...\n")
ate_D <- tryCatch(
  ATE_estim_etas(
    statespace = win_km, partition = partition,
    observed_data = pp_post_sem_D, treated_partitions = treated_partitions,
    etas_params = list(control = D_ctrl_list, treated = D_treat_list),
    n_sims = ATE_N_SIMS, n_tau_i = ATE_N_TAU_I, n_tau_sims = ATE_N_TAU_SIM,
    windowT = windowT_post, windowS = win_km,
    m0 = ETAS_M0, beta_gr = BETA_GR, fixed_params = FIXED_STRUCTURAL
  ),
  error = function(e) { cat("    ATE-D error:", e$message, "\n"); NULL }
)

# ============================================================================
# 7. Summary tables
# ============================================================================
cat("\n")
cat("===========================================================================\n")
cat("                       PARAMETER COMPARISON TABLE\n")
cat("===========================================================================\n\n")

fmt <- function(x) sprintf("%8.4f", x)

biv_ctrl_val <- function(par, p) {
  if (p == "mu") par[["mu_0"]] else par[[paste0(p, "_00")]]
}
biv_treat_val <- function(par, p) {
  if (p == "mu") par[["mu_1"]] else par[[paste0(p, "_11")]]
}

cat(sprintf("%-12s  %8s  %8s  |  %8s  %8s  |  %8s  %8s  |  %8s  %8s\n",
  "", "A.ctrl", "A.treat", "B.ctrl", "B.treat", "C.ctrl", "C.treat", "D.ctrl", "D.treat"))
cat(paste(rep("-", 106), collapse = ""), "\n")

for (p in c("mu", "A", "alpha_m")) {
  cat(sprintf("%-12s  %8s  %8s  |  %8s  %8s  |  %8s  %8s  |  %8s  %8s\n",
    p,
    fmt(unlist(A_ctrl)[p]), fmt(unlist(A_treat)[p]),
    fmt(biv_ctrl_val(B_params, p)), fmt(biv_treat_val(B_params, p)),
    fmt(unlist(C_ctrl)[p]), fmt(unlist(C_treat)[p]),
    fmt(biv_ctrl_val(D_params, p)), fmt(biv_treat_val(D_params, p))))
}

cat("\nCross-excitation:\n")
cat(sprintf("  B: A_01=%.4f  A_10=%.4f  alpha_m_01=%.4f  alpha_m_10=%.4f\n",
  B_params[["A_01"]], B_params[["A_10"]], B_params[["alpha_m_01"]], B_params[["alpha_m_10"]]))
if (!is.null(semD)) {
  cat(sprintf("  D: A_01=%.4f  A_10=%.4f  alpha_m_01=%.4f  alpha_m_10=%.4f\n",
    D_params[["A_01"]], D_params[["A_10"]], D_params[["alpha_m_01"]], D_params[["alpha_m_10"]]))
}

cat("\n")
cat("===========================================================================\n")
cat("                           ATE COMPARISON\n")
cat("===========================================================================\n\n")
cat(sprintf("%-22s  %12s  %12s  %12s  %12s\n",
            "", "All-Nothing", "Tau-1 (flip)", "ATE-naive", "ATE-spillover"))
cat(paste(rep("-", 80), collapse = ""), "\n")

for (nm in list(list(ate = ate_A, lab = "A: Naive Indep"),
                list(ate = ate_B, lab = "B: Naive Biv"),
                list(ate = ate_C, lab = "C: SEM Indep"),
                list(ate = ate_D, lab = "D: SEM Biv"))) {
  if (!is.null(nm$ate)) {
    an_mean <- mean(nm$ate$all_nothing_sim$ATE, na.rm = TRUE)
    cat(sprintf("%-22s  %12.2f  %12.2f  %12.2f  %12.2f\n",
      nm$lab, an_mean, nm$ate$tau_1_estim,
      nm$ate$ATE_naive, nm$ate$ATE_spillover))
  } else {
    cat(sprintf("%-22s  %12s  %12s  %12s  %12s\n",
      nm$lab, "FAILED", "FAILED", "FAILED", "FAILED"))
  }
}
cat("\nAll-Nothing: mean per-county difference if ALL counties treated vs ALL control.\n")
cat("Tau-1: average one-flip effect (randomly switch one county).\n")
cat("ATE-naive: difference in mean counts between treated and control counties.\n")
cat("ATE-spillover: excess counts in control counties due to cross-triggering.\n")

# ============================================================================
# 8. Save results
# ============================================================================
cat("\n--- Saving results ---\n")

results <- list(
  fitA = list(ctrl = A_ctrl, treat = A_treat, ate = ate_A),
  fitB = list(params = B_params, loglik = B_loglik, fit = fitB, ate = ate_B),
  fitC = list(ctrl = C_ctrl, treat = C_treat, sem = semC, ate = ate_C),
  fitD = list(params = D_params, ctrl = D_ctrl, treat = D_treat, sem = semD, ate = ate_D),
  pp_data = list(pp_pre = pp_pre, pp_post = pp_post),
  counties = list(
    names = counties_sf_valid$NAME,
    treated_names = treated_names,
    n_counties = partition$n,
    n_treated = sum(treated_idx)
  ),
  config = list(
    ETAS_M0 = ETAS_M0, BETA_GR = BETA_GR,
    FIXED_STRUCTURAL = FIXED_STRUCTURAL,
    SEM_N_ITER = SEM_N_ITER, SEM_INNER_ITER = SEM_INNER_ITER,
    SEM_N_LABELLINGS = SEM_N_LABELLINGS,
    SEM_CHANGE_FACTOR = SEM_CHANGE_FACTOR,
    ATE_N_SIMS = ATE_N_SIMS, ATE_N_TAU_I = ATE_N_TAU_I,
    TEST_MODE = TEST_MODE,
    windowT_post = windowT_post,
    n_pre = nrow(pp_pre), n_post = nrow(pp_post),
    n_counties = partition$n, n_treated = sum(treated_idx)
  )
)

out_file <- file.path(OUT_DIR, "oklahoma_results.rds")
saveRDS(results, out_file)
cat(sprintf("Results saved to: %s\n", out_file))
cat(sprintf("Plots saved to:   %s\n", PLOT_DIR))

cat("\n=== Oklahoma Analysis Complete ===\n")
