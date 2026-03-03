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
# Configuration
# =========================================================
#
# Non-IPD placebo analysis: fit Hawkes to non-IPD (pneumococcal
# surveillance) cases using IPD cases for the background rate.
# We expect NO treatment effect on non-IPD cases, so this serves
# as a negative control for the IPD radii study.
#
# To keep computation tractable we randomly sample N_SAMPLE cases
# pre- and post-treatment from the full non-IPD pool.
#
# Two modes controlled by TEMPORAL_MODE:
#   "fixed_beta" - fix beta=0.05, no truncation (default)
#   "truncated"  - free beta, truncate kernel at TRUNC_DAYS
#
# Background: inhomogeneous KDE via bw.diggle on all IPD cases.
#
# Output: RDS with list(results=df, adaptive_by_radius=list, config=list).

QUICK_TEST     <- FALSE  # set TRUE for fast local iteration
MIN_POST_CASES <- 15     # skip radius if either partition has fewer cases
N_SAMPLE       <- 1000   # random sample size pre and post treatment
SAMPLE_SEED    <- 2024

TEMPORAL_MODE   <- "fixed_beta"   # "fixed_beta" or "truncated"
TREATMENT_DATE  <- as_datetime("2011-05-01")
TRUNC_DAYS      <- 30

SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
OUT_DIR <- file.path(SCRIPT_DIR, "output")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

if (TEMPORAL_MODE == "fixed_beta") {
  FIXED_BETA <- 0.05
  T_TRUNC    <- NULL
  OUT_FILE   <- file.path(OUT_DIR, "gambia_radii_results_nonIPD.rds")
} else {
  FIXED_BETA <- NULL
  T_TRUNC    <- TRUNC_DAYS
  OUT_FILE   <- file.path(OUT_DIR, sprintf("gambia_radii_results_nonIPD_trunc%d.rds", TRUNC_DAYS))
}

RADII_KM <- c(seq(0.5, 3.0, by = 0.5), seq(4, 10, by = 1), seq(12, 20, by = 2))

VANILLA_MAXIT  <- if (QUICK_TEST) 100 else 5000
VANILLA_STARTS <- list(
  list(mu = 0.0005, alpha = 1e-7,  beta = 0.05, K = 0.3),
  list(mu = 0.001,  alpha = 1e-6,  beta = 0.03, K = 0.2),
  list(mu = 0.0005, alpha = 1e-8,  beta = 0.1,  K = 0.4),
  list(mu = 0.002,  alpha = 1e-5,  beta = 0.05, K = 0.15),
  list(mu = 0.0003, alpha = 1e-9,  beta = 0.05, K = 0.5)
)

SEM_N_LABELLINGS         <- if (QUICK_TEST) 5 else 20
SEM_N_ITER               <- if (QUICK_TEST) 1 else 2
SEM_INNER_ITER           <- if (QUICK_TEST) 5 else 50
SEM_INNER_N_PROPS        <- if (QUICK_TEST) 5 else 20
SEM_PARAM_UPDATE_CADENCE <- 10
SEM_CHANGE_FACTOR        <- 0.05
SEM_INCLUDE_STARTING     <- TRUE
SEM_UPDATE_STARTING      <- TRUE

# =========================================================
# Data Loading
# =========================================================
gambia_all <- dget(system.file("extdata/gambia/DISCLEAN.txt", package = "PPDisentangle"))
gambia_ipd <- gambia_all %>% filter(IPDEvent == "Y")

# Time origin and study window aligned with the IPD study
TREATMENT_TIME <- as.numeric(as.Date(TREATMENT_DATE) - min(as.Date(gambia_ipd$date_IPD_Pne_Event)))
STUDY_END <- as.numeric(max(as.Date(gambia_ipd$date_IPD_Pne_Event)) - min(as.Date(gambia_ipd$date_IPD_Pne_Event)))
t0 <- min(as.Date(gambia_ipd$date_IPD_Pne_Event))

win_orig <- dget(system.file("extdata/gambia/GAMWIN_PROJ.txt", package = "PPDisentangle"))
S <- diag(c(1000, 1000)); t <- c(0, 0)
win_full <- affine.owin(win_orig, mat = S, vec = t)
hc <- dget(system.file("extdata/gambia/hc.txt", package = "PPDisentangle"))
hc <- affine(hc, mat = S, vec = t)

v_all <- vect(gambia_all, geom = c("long","lat"), crs = "EPSG:4326")
v_all_utm <- project(v_all, "EPSG:32628")
coords_all <- crds(v_all_utm)
gambia_all$easting <- coords_all[,1]
gambia_all$northing <- coords_all[,2]

# =========================================================
# Build non-IPD point pattern with random subsample
# =========================================================
non_ipd_all <- gambia_all %>% filter(IPDEvent == "N")
non_ipd_all <- non_ipd_all[inside.owin(non_ipd_all$easting, non_ipd_all$northing, win_full), ]

pp_nonipd_full <- data.frame(
  x = non_ipd_all$easting,
  y = non_ipd_all$northing,
  t = as.numeric(as.Date(non_ipd_all$date_IPD_Pne_Event)) - as.numeric(t0),
  IPD = FALSE
)
pp_nonipd_full <- pp_nonipd_full[inside.owin(pp_nonipd_full$x, pp_nonipd_full$y, win_full), ]
pp_nonipd_full$t <- pp_nonipd_full$t - min(pp_nonipd_full$t)
pp_nonipd_full <- pp_nonipd_full %>% filter(t <= STUDY_END)

pre_idx  <- which(pp_nonipd_full$t <  TREATMENT_TIME)
post_idx <- which(pp_nonipd_full$t >= TREATMENT_TIME)
cat(sprintf("Non-IPD pool: %d pre-treatment, %d post-treatment\n", length(pre_idx), length(post_idx)))

set.seed(SAMPLE_SEED)
n_pre_sample  <- min(N_SAMPLE, length(pre_idx))
n_post_sample <- min(N_SAMPLE, length(post_idx))
sampled_pre  <- pre_idx[sample.int(length(pre_idx),   n_pre_sample)]
sampled_post <- post_idx[sample.int(length(post_idx), n_post_sample)]
pp_data_full <- pp_nonipd_full[sort(c(sampled_pre, sampled_post)), ]
cat(sprintf("Sampled: %d pre + %d post = %d non-IPD cases\n", n_pre_sample, n_post_sample, nrow(pp_data_full)))

# Jitter exact duplicates
dups <- duplicated(pp_data_full[, c("x", "y", "t")])
if (any(dups)) {
  set.seed(42)
  n_dup <- sum(dups)
  dup_idx <- which(dups)
  pp_data_full$x[dup_idx] <- pp_data_full$x[dup_idx] + runif(n_dup, -0.1, 0.1)
  pp_data_full$y[dup_idx] <- pp_data_full$y[dup_idx] + runif(n_dup, -0.1, 0.1)
  pp_data_full$t[dup_idx] <- pp_data_full$t[dup_idx] + runif(n_dup, -0.01, 0.01)
  cat(sprintf("Jittered %d duplicate non-IPD points\n", n_dup))
}

# =========================================================
# Background KDE (from IPD cases -- reversed role)
# =========================================================
ipd_all <- gambia_all %>% filter(IPDEvent == "Y")
ipd_all <- ipd_all[inside.owin(ipd_all$easting, ipd_all$northing, win_full), ]
X_bg <- ppp(x = ipd_all$easting, y = ipd_all$northing, window = win_full)
bw_sigma <- suppressWarnings(bw.diggle(X_bg))
lambda_im <- density(X_bg, sigma = bw_sigma, edge = TRUE, at = "pixels")
min_nz <- min(lambda_im$v[lambda_im$v > 0], na.rm = TRUE)
lambda_im$v[lambda_im$v <= 0] <- min_nz
cat(sprintf("Background KDE: %d IPD points, bw=%.1f m\n", nrow(ipd_all), as.numeric(bw_sigma)))

# =========================================================
# Helper Functions
# =========================================================
normalize_marks_gambia <- function(df_sub, win_sub, covariate_im, mark_name = "W") {
  if (nrow(df_sub) == 0) return(list(new_df = df_sub, mass = 0, norm = 0))
  cov_in_window <- covariate_im[win_sub, drop = FALSE]
  total_mass_raw <- integral.im(cov_in_window)
  target_area <- spatstat.geom::area(win_sub)
  norm_factor <- target_area / total_mass_raw
  vals_raw <- raster::extract(raster(covariate_im), df_sub[, c("x", "y")])
  vals_raw[is.na(vals_raw)] <- 0
  df_sub[[mark_name]] <- vals_raw * norm_factor
  min_val <- min(df_sub[[mark_name]][df_sub[[mark_name]] > 0], na.rm = TRUE)
  if (is.infinite(min_val) | is.na(min_val)) min_val <- 1e-9
  df_sub[[mark_name]][df_sub[[mark_name]] <= 0] <- min_val
  return(list(new_df = df_sub, mass = total_mass_raw, norm = norm_factor))
}

run_full_fit <- function(df, win, label, fixed_beta = FIXED_BETA,
                         starts = VANILLA_STARTS, maxit = VANILLA_MAXIT,
                         t_trunc = T_TRUNC, verbose_fail = TRUE) {
  if (nrow(df) < 5) {
    if (verbose_fail) cat(sprintf("    [%s] too few points: n=%d\n", label, nrow(df)))
    return(NULL)
  }
  fp <- if (!is.null(fixed_beta)) list(beta = fixed_beta) else NULL
  best_fit <- NULL; best_val <- -Inf
  n_err <- 0L; n_bad_ll <- 0L; n_bad_par <- 0L; last_fit <- NULL
  for (i in seq_along(starts)) {
    fit <- tryCatch(
      suppressWarnings(fit_hawkes(
        params_init = starts[[i]], realiz = df, windowT = range(df$t),
        windowS = win, background_rate_var = "W", maxit = maxit,
        use_fast = TRUE, method = "Nelder-Mead", fixed_params = fp,
        t_trunc = t_trunc)),
      error = function(e) { if (verbose_fail) cat(sprintf("    [%s start %d] error: %s\n", label, i, e$message)); NULL }
    )
    if (is.null(fit)) { n_err <- n_err + 1L; next }
    if (!is.finite(fit$value) || fit$value <= -1e14) {
      n_bad_ll <- n_bad_ll + 1L
      last_fit <- fit
      next
    }
    pars <- as.numeric(unlist(fit$par))
    if (length(pars) != 4L || !all(is.finite(pars)) || pars[4] < 0 || pars[4] >= 0.999) {
      n_bad_par <- n_bad_par + 1L
      last_fit <- fit
      next
    }
    if (fit$value > best_val) { best_val <- fit$value; best_fit <- fit }
  }
  if (is.null(best_fit)) {
    if (verbose_fail)
      cat(sprintf("    [%s] all %d starts failed: %d err, %d bad_ll, %d bad_par%s\n",
          label, length(starts), n_err, n_bad_ll, n_bad_par,
          if (!is.null(last_fit))
            sprintf(" | last: ll=%.2g K=%.4f par_len=%d",
                last_fit$value, as.numeric(unlist(last_fit$par))[4],
                length(as.numeric(unlist(last_fit$par))))
          else ""))
    return(NULL)
  }
  pars <- as.numeric(unlist(best_fit$par))
  names(pars) <- c("mu", "alpha", "beta", "K")
  cat(sprintf("    [%s] best ll=%.1f | %s\n", label,
      best_val, paste(names(pars), signif(pars, 4), sep="=", collapse=" ")))
  return(as.list(pars))
}

as_parlist <- function(p) {
  v <- as.numeric(p); names(v) <- c("mu", "alpha", "beta", "K"); as.list(v)
}

params_valid <- function(p, label = "") {
  K <- p[["K"]]; mu <- p[["mu"]]
  ok <- is.finite(K) && is.finite(mu) && K >= 0 && K < 0.999 && mu > 1e-12
  if (!ok) cat(sprintf("    [%s] degenerate: mu=%.2e K=%.4f\n", label, mu, K))
  ok
}

estimate_savings_simple <- function(p_ctrl, p_trtd, m_ctrl, m_trtd, m_target, duration) {
  denom_c <- 1 - p_ctrl[["K"]]
  denom_t <- 1 - p_trtd[["K"]]
  if (denom_c < 1e-3 || denom_t < 1e-3 || m_ctrl < 1e-10 || m_trtd < 1e-10 ||
      p_ctrl[["mu"]] < 1e-12 || p_trtd[["mu"]] < 1e-12)
    return(list(exp_c = NA, exp_t = NA, savings = NA, pct = NA))
  rate_c <- p_ctrl[["mu"]] / m_ctrl
  rate_t <- p_trtd[["mu"]] / m_trtd
  bg_c <- rate_c * m_target * duration
  bg_t <- rate_t * m_target * duration
  exp_c <- bg_c / denom_c
  exp_t <- bg_t / denom_t
  pct <- if (exp_c > 0) (exp_c - exp_t) / exp_c else NA
  if (!is.na(pct) && (pct < -10 || pct > 10))
    return(list(exp_c = exp_c, exp_t = exp_t, savings = exp_c - exp_t, pct = NA))
  return(list(exp_c = exp_c, exp_t = exp_t, savings = exp_c - exp_t, pct = pct))
}

# =========================================================
# Radii Study
# =========================================================
results_list <- list()
adaptive_by_radius <- list()
louis_by_radius <- list()
skipped <- character()
t_study_start <- proc.time()[3]

cat(sprintf("\nNon-IPD radii study: %s | %d radii | beta=%s | t_trunc=%s | quick=%s | N=%d per period\n",
    TEMPORAL_MODE, length(RADII_KM),
    if (is.null(FIXED_BETA)) "free" else sprintf("%.3f", FIXED_BETA),
    if (is.null(T_TRUNC)) "none" else sprintf("%d days", T_TRUNC),
    QUICK_TEST, N_SAMPLE))

for (r_km in RADII_KM) {
  r_m <- r_km * 1000
  t_radius_start <- proc.time()[3]
  cat(sprintf("\n=== Radius: %.1f km ===\n", r_km))

  discs <- lapply(seq_len(npoints(hc)), function(i)
    disc(radius = r_m, centre = c(hc$x[i], hc$y[i]), npoly = 128))
  treated_ss <- Reduce(union.owin, discs)
  control_ss <- setminus.owin(win_full, treated_ss)

  pp_study <- pp_data_full
  pp_study$location_process <- ifelse(
    inside.owin(pp_study$x, pp_study$y, treated_ss), "treated", "control")
  pp_study$process <- pp_study$location_process
  pp_study$background <- TRUE

  n_pre_t  <- sum(pp_study$t <  TREATMENT_TIME & pp_study$location_process == "treated")
  n_pre_c  <- sum(pp_study$t <  TREATMENT_TIME & pp_study$location_process == "control")
  n_post_t <- sum(pp_study$t >= TREATMENT_TIME & pp_study$location_process == "treated")
  n_post_c <- sum(pp_study$t >= TREATMENT_TIME & pp_study$location_process == "control")
  cat(sprintf("  Pre: trt=%d ctrl=%d | Post: trt=%d ctrl=%d\n", n_pre_t, n_pre_c, n_post_t, n_post_c))

  if (n_post_t < MIN_POST_CASES || n_post_c < MIN_POST_CASES) {
    cat(sprintf("  SKIP: too few post-treatment cases (min=%d)\n", MIN_POST_CASES))
    skipped <- c(skipped, sprintf("%.1f km: too few cases (trt=%d ctrl=%d)", r_km, n_post_t, n_post_c))
    next
  }

  # Normalize marks per partition
  res_pre_trtd  <- normalize_marks_gambia(pp_study %>% filter(t <  TREATMENT_TIME, location_process == "treated"),  treated_ss, lambda_im)
  res_pre_ctrl  <- normalize_marks_gambia(pp_study %>% filter(t <  TREATMENT_TIME, location_process == "control"), control_ss, lambda_im)
  res_post_trtd <- normalize_marks_gambia(pp_study %>% filter(t >= TREATMENT_TIME, location_process == "treated"),  treated_ss, lambda_im)
  res_post_ctrl <- normalize_marks_gambia(pp_study %>% filter(t >= TREATMENT_TIME, location_process == "control"), control_ss, lambda_im)

  # --- Vanilla Hawkes ---
  cat("  Vanilla fits...\n")
  vanilla_trtd <- run_full_fit(res_post_trtd$new_df, treated_ss, "Treated")
  vanilla_ctrl <- run_full_fit(res_post_ctrl$new_df, control_ss, "Control")
  if (is.null(vanilla_trtd) || is.null(vanilla_ctrl)) {
    cat(sprintf("  SKIP: vanilla fit failed (trt=%s ctrl=%s)\n",
                if (is.null(vanilla_trtd)) "FAIL" else "ok",
                if (is.null(vanilla_ctrl)) "FAIL" else "ok"))
    skipped <- c(skipped, sprintf("%.1f km: vanilla fit failed", r_km))
    next
  }

  # --- SEM ---
  cat("  SEM fit...\n")
  pp_sem <- rbind(res_pre_trtd$new_df, res_pre_ctrl$new_df,
                  res_post_trtd$new_df, res_post_ctrl$new_df)

  cov_lookup <- function(x, y) {
    rv <- as.numeric(spatstat.geom::interp.im(lambda_im, list(x = x, y = y)))
    rv[is.na(rv)] <- 0
    is_t <- inside.owin(x, y, treated_ss)
    rv[is_t]  <- rv[is_t]  * res_pre_trtd$norm
    rv[!is_t] <- rv[!is_t] * res_pre_ctrl$norm
    return(rv)
  }

  sem_res <- tryCatch(
    adaptive_SEM(
      pp_data = pp_sem,
      partition = tess(tiles = list("control" = control_ss, "treated" = treated_ss)),
      partition_processes = c("control", "treated"),
      statespace = win_full,
      time_window = c(0, STUDY_END),
      treatment_time = TREATMENT_TIME,
      hawkes_params_control = vanilla_ctrl,
      hawkes_params_treated = vanilla_trtd,
      N_labellings = SEM_N_LABELLINGS,
      N_iter = SEM_N_ITER,
      covariate_lookup = cov_lookup,
      background_rate_var = "W",
      t_trunc = T_TRUNC,
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
        verbose                  = FALSE
      )
    ),
    error = function(e) { cat(sprintf("  SEM ERROR: %s\n", e$message)); NULL }
  )

  if (is.null(sem_res)) {
    cat("  SKIP: SEM failed\n")
    skipped <- c(skipped, sprintf("%.1f km: SEM failed", r_km))
    next
  }

  sem_ctrl <- as_parlist(sem_res$hawkes_params_control)
  sem_trtd <- as_parlist(sem_res$hawkes_params_treated)
  sem_valid_ctrl <- params_valid(sem_ctrl, "SEM ctrl")
  sem_valid_trtd <- params_valid(sem_trtd, "SEM trtd")
  cat(sprintf("  SEM done: K_trt=%.3f K_ctrl=%.3f%s\n", sem_trtd[["K"]], sem_ctrl[["K"]],
      if (!sem_valid_ctrl || !sem_valid_trtd) " [DEGENERATE]" else ""))

  adaptive_by_radius[[as.character(r_km)]] <- list(
    history = sem_res$adaptive_history,
    n_post = n_post_t + n_post_c
  )

  # --- Louis SEs for SEM ---
  if (sem_valid_ctrl && sem_valid_trtd) {
    louis_res <- tryCatch(
      louis_standard_errors(
        sem_res, TREATMENT_TIME, win_full,
        tess(tiles = list("control" = control_ss, "treated" = treated_ss)),
        c("control", "treated"),
        background_rate_var = "W",
        verbose = FALSE,
        t_trunc = T_TRUNC
      ),
      error = function(e) { cat(sprintf("  Louis ERROR: %s\n", e$message)); NULL }
    )
    if (!is.null(louis_res) && !any(is.na(louis_res$se))) {
      louis_by_radius[[as.character(r_km)]] <- list(
        vcov = louis_res$vcov,
        theta_hat = c(unlist(sem_ctrl), unlist(sem_trtd)),
        se = louis_res$se
      )
    }
  }

  # --- Savings ---
  duration_post <- STUDY_END - TREATMENT_TIME
  total_m <- res_pre_trtd$mass + res_pre_ctrl$mass
  sav_vanilla <- estimate_savings_simple(vanilla_ctrl, vanilla_trtd, res_pre_ctrl$mass, res_pre_trtd$mass, total_m, duration_post)
  sav_sem     <- estimate_savings_simple(sem_ctrl, sem_trtd, res_pre_ctrl$mass, res_pre_trtd$mass, total_m, duration_post)

  elapsed_radius <- round((proc.time()[3] - t_radius_start) / 60, 1)
  v_pct_str <- if (is.na(sav_vanilla$pct)) "NA" else sprintf("%.1f%%", sav_vanilla$pct * 100)
  s_pct_str <- if (is.na(sav_sem$pct))     "NA" else sprintf("%.1f%%", sav_sem$pct * 100)
  cat(sprintf("  Savings: vanilla=%s sem=%s (%.1f min)\n", v_pct_str, s_pct_str, elapsed_radius))

  results_list[[as.character(r_km)]] <- data.frame(
    radius_km = r_km, n_post_trt = n_post_t, n_post_ctrl = n_post_c,
    vanilla_K_ctrl = vanilla_ctrl[["K"]], vanilla_K_trtd = vanilla_trtd[["K"]],
    vanilla_mu_ctrl = vanilla_ctrl[["mu"]], vanilla_mu_trtd = vanilla_trtd[["mu"]],
    vanilla_alpha_ctrl = vanilla_ctrl[["alpha"]], vanilla_alpha_trtd = vanilla_trtd[["alpha"]],
    vanilla_beta_ctrl = vanilla_ctrl[["beta"]], vanilla_beta_trtd = vanilla_trtd[["beta"]],
    sem_K_ctrl = sem_ctrl[["K"]], sem_K_trtd = sem_trtd[["K"]],
    sem_mu_ctrl = sem_ctrl[["mu"]], sem_mu_trtd = sem_trtd[["mu"]],
    sem_alpha_ctrl = sem_ctrl[["alpha"]], sem_alpha_trtd = sem_trtd[["alpha"]],
    sem_beta_ctrl = sem_ctrl[["beta"]], sem_beta_trtd = sem_trtd[["beta"]],
    vanilla_savings_pct = sav_vanilla$pct, sem_savings_pct = sav_sem$pct,
    sem_degenerate = !sem_valid_ctrl || !sem_valid_trtd
  )
}

elapsed_total <- round((proc.time()[3] - t_study_start) / 60, 1)

# =========================================================
# Results
# =========================================================
cat(sprintf("\n=== NON-IPD RESULTS: %s | beta=%s | t_trunc=%s | %.1f min | %d/%d radii ===\n",
    TEMPORAL_MODE,
    if (is.null(FIXED_BETA)) "free" else sprintf("%.3f", FIXED_BETA),
    if (is.null(T_TRUNC)) "none" else sprintf("%d", T_TRUNC),
    elapsed_total, length(results_list), length(RADII_KM)))

if (length(skipped) > 0) {
  cat("Skipped:\n")
  for (s in skipped) cat(sprintf("  - %s\n", s))
}

if (length(results_list) == 0) {
  cat("\nNo radii produced results.\n")
  final_results <- data.frame()
} else {
  final_results <- do.call(rbind, results_list)
  final_results$sem_minus_vanilla <- final_results$sem_savings_pct - final_results$vanilla_savings_pct
  print(final_results[, c("radius_km", "n_post_trt", "n_post_ctrl",
                           "vanilla_K_trtd", "vanilla_K_ctrl",
                           "sem_K_trtd", "sem_K_ctrl",
                           "vanilla_savings_pct", "sem_savings_pct", "sem_minus_vanilla")],
        row.names = FALSE)
}

results_full <- list(
  results = final_results,
  adaptive_by_radius = adaptive_by_radius,
  louis_by_radius = louis_by_radius,
  config = list(TEMPORAL_MODE = TEMPORAL_MODE, FIXED_BETA = FIXED_BETA, T_TRUNC = T_TRUNC,
                N_SAMPLE = N_SAMPLE, SAMPLE_SEED = SAMPLE_SEED, data_type = "nonIPD")
)
saveRDS(results_full, OUT_FILE)
cat(sprintf("\nSaved to %s (results + adaptive_history)\n", OUT_FILE))
