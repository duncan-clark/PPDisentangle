##
## Sensitivity analysis: Hawkes K estimates across background specifications
##
## Compares homogeneous, inhomogeneous (bw.diggle, bw*20), and
## non-IPD background sources, all with beta fixed at 0.05 and no jitter.
##
library(PPDisentangle)
library(spatstat)
library(dplyr)
library(terra)
library(raster)
library(lubridate)

TREATMENT_DATE <- as_datetime("2011-05-01")
HEALTH_CENTER_DIST <- 10000
FIXED_BETA <- 0.05
set.seed(42)

# ---- Data prep ----
gambia_data <- dget(system.file("extdata/gambia/DISCLEAN.txt", package = "PPDisentangle"))
ipd_dates   <- as.Date(gambia_data$date_IPD_Pne_Event[gambia_data$IPDEvent == "Y"])
TREATMENT_TIME <- as.numeric(as.Date(TREATMENT_DATE) - min(ipd_dates))
STUDY_END      <- as.numeric(max(ipd_dates) - min(ipd_dates))

win <- dget(system.file("extdata/gambia/GAMWIN_PROJ.txt", package = "PPDisentangle"))
S <- diag(c(1000, 1000)); t_vec <- c(0, 0)
win <- affine.owin(win, mat = S, vec = t_vec)
hc  <- dget(system.file("extdata/gambia/hc.txt", package = "PPDisentangle"))
hc  <- affine(hc, mat = S, vec = t_vec)

v <- vect(gambia_data, geom = c("long", "lat"), crs = "EPSG:4326")
v_utm <- project(v, "EPSG:32628"); coords <- crds(v_utm)
gambia_data$easting_utm  <- coords[, 1]
gambia_data$northing_utm <- coords[, 2]

discs       <- lapply(seq_len(npoints(hc)), function(i) disc(radius = HEALTH_CENTER_DIST,     centre = c(hc$x[i], hc$y[i]), npoly = 200))
discs_outer <- lapply(seq_len(npoints(hc)), function(i) disc(radius = HEALTH_CENTER_DIST * 2, centre = c(hc$x[i], hc$y[i]), npoly = 200))
treated_state_space <- Reduce(union.owin, discs)
control_state_space <- setminus.owin(Reduce(union.owin, discs_outer), treated_state_space)
win_full            <- Reduce(union.owin, discs_outer)

# IPD point pattern (no jitter)
ipd <- gambia_data %>% filter(IPDEvent == "Y")
pp_data <- data.frame(x = ipd$easting_utm, y = ipd$northing_utm,
                      t = as.numeric(as.Date(ipd$date_IPD_Pne_Event)) - min(as.numeric(ipd_dates)))
pp_data <- pp_data[inside.owin(pp_data$x, pp_data$y, win_full), ]
pp_data$t <- pp_data$t - min(pp_data$t)
pp_data <- pp_data %>% filter(t <= STUDY_END)
pp_data$location_process <- ifelse(inside.owin(pp_data$x, pp_data$y, w = treated_state_space), "treated", "control")

# Non-IPD sample for external background
non_ipd <- gambia_data %>% filter(IPDEvent == "N")
non_ipd <- non_ipd[inside.owin(non_ipd$easting_utm, non_ipd$northing_utm, win_full), ]
bg_sample <- non_ipd[sample(nrow(non_ipd), min(1000, nrow(non_ipd))), ]

# ---- Helper functions ----
normalize_marks <- function(df_sub, win_sub, cov_im, mark_name = "W") {
  cov_in <- cov_im[win_sub, drop = FALSE]
  total_mass <- integral.im(cov_in)
  target_area <- spatstat.geom::area(win_sub)
  norm_factor <- target_area / total_mass
  vals_raw <- raster::extract(raster(cov_im), df_sub[, c("x", "y")])
  vals_raw[is.na(vals_raw)] <- 0
  df_sub[[mark_name]] <- vals_raw * norm_factor
  min_val <- min(df_sub[[mark_name]][df_sub[[mark_name]] > 0], na.rm = TRUE)
  if (is.infinite(min_val) | is.na(min_val)) min_val <- 1e-9
  df_sub[[mark_name]][df_sub[[mark_name]] <= 0] <- min_val
  list(new_df = df_sub)
}

fit_segment <- function(df_sub, ws, use_W = TRUE) {
  tw <- range(df_sub$t)
  fit_hawkes(
    params_init = list(mu = 0.0005, alpha = 1e-7, beta = 0.05, K = 0.3),
    realiz = df_sub, windowT = tw, windowS = ws,
    background_rate_var = if (use_W) "W" else NULL,
    maxit = 2000, use_fast = TRUE, method = "Nelder-Mead",
    fixed_params = list(beta = FIXED_BETA)
  )
}

run_all_segments <- function(pp_data, lambda_im = NULL, label = "") {
  seg_specs <- list(
    "Pre-Control"  = list(filt = quote(t < TREATMENT_TIME & location_process == "control"),  ws = control_state_space),
    "Pre-Treated"  = list(filt = quote(t < TREATMENT_TIME & location_process == "treated"),  ws = treated_state_space),
    "Post-Control" = list(filt = quote(t >= TREATMENT_TIME & location_process == "control"), ws = control_state_space),
    "Post-Treated" = list(filt = quote(t >= TREATMENT_TIME & location_process == "treated"), ws = treated_state_space)
  )
  use_W <- !is.null(lambda_im)
  cat(sprintf("\n=== %s ===\n", label))
  results <- list()
  for (nm in names(seg_specs)) {
    df_sub <- pp_data %>% filter(eval(seg_specs[[nm]]$filt))
    if (use_W) df_sub <- normalize_marks(df_sub, seg_specs[[nm]]$ws, lambda_im)$new_df
    fit <- fit_segment(df_sub, seg_specs[[nm]]$ws, use_W)
    p <- fit$par
    cat(sprintf("  %s (N=%d): mu=%.4e  alpha=%.2e  beta=%.4f  K=%.4f  | dist=%4.0fm  ll=%.1f\n",
                nm, nrow(df_sub), p[1], p[2], p[3], p[4], sqrt(1 / p[2]), fit$value))
    results[[nm]] <- list(par = p, ll = fit$value, N = nrow(df_sub))
  }
  invisible(results)
}

make_kde <- function(pts_x, pts_y, sigma) {
  X <- ppp(x = pts_x, y = pts_y, window = win_full)
  im <- density(X, sigma = sigma, edge = TRUE, at = "pixels")
  min_nz <- min(im$v[im$v > 0], na.rm = TRUE)
  im$v[im$v <= 0] <- min_nz
  im
}

# ---- Run fits ----

# 1. Homogeneous
res_homog <- run_all_segments(pp_data, NULL, "Homogeneous, no jitter, beta=0.05")

# 2. IPD-based bw.diggle
X_ipd_pre <- ppp(x = pp_data$x[pp_data$t < TREATMENT_TIME], y = pp_data$y[pp_data$t < TREATMENT_TIME], window = win_full)
bw_ipd <- bw.diggle(X_ipd_pre)
cat(sprintf("\nbw.diggle (IPD pre-treatment): %.1f m\n", bw_ipd))

lambda_ipd_dig <- make_kde(pp_data$x[pp_data$t < TREATMENT_TIME], pp_data$y[pp_data$t < TREATMENT_TIME], bw_ipd)
res_ipd_dig    <- run_all_segments(pp_data, lambda_ipd_dig, "IPD bw.diggle, no jitter, beta=0.05")

lambda_ipd_20  <- make_kde(pp_data$x[pp_data$t < TREATMENT_TIME], pp_data$y[pp_data$t < TREATMENT_TIME], bw_ipd * 20)
res_ipd_20     <- run_all_segments(pp_data, lambda_ipd_20, "IPD bw*20, no jitter, beta=0.05")

# 3. Non-IPD-based bw.diggle
bw_nonipd <- bw.diggle(ppp(x = bg_sample$easting_utm, y = bg_sample$northing_utm, window = win_full))
cat(sprintf("\nbw.diggle (non-IPD 1000): %.1f m\n", bw_nonipd))

lambda_nonipd_dig <- make_kde(bg_sample$easting_utm, bg_sample$northing_utm, bw_nonipd)
res_nonipd_dig    <- run_all_segments(pp_data, lambda_nonipd_dig, "non-IPD bw.diggle, no jitter, beta=0.05")

lambda_nonipd_20  <- make_kde(bg_sample$easting_utm, bg_sample$northing_utm, bw_nonipd * 20)
res_nonipd_20     <- run_all_segments(pp_data, lambda_nonipd_20, "non-IPD bw*20, no jitter, beta=0.05")

# ---- Summary table ----
cat("\n\n=== SUMMARY: K estimates (beta fixed at 0.05, no jitter) ===\n")
cat(sprintf("%-25s  %8s  %8s  %8s  %8s\n", "Setting", "Pre-Ctrl", "Pre-Trt", "Post-Ctrl", "Post-Trt"))
cat(paste0(rep("-", 65), collapse = ""), "\n")
all_res <- list(
  "Homogeneous"      = res_homog,
  "IPD bw.diggle"    = res_ipd_dig,
  "IPD bw*20"        = res_ipd_20,
  "non-IPD bw.diggle"= res_nonipd_dig,
  "non-IPD bw*20"    = res_nonipd_20
)
for (nm in names(all_res)) {
  r <- all_res[[nm]]
  cat(sprintf("%-25s  %8.3f  %8.3f  %8.3f  %8.3f\n", nm,
              r[["Pre-Control"]]$par[4], r[["Pre-Treated"]]$par[4],
              r[["Post-Control"]]$par[4], r[["Post-Treated"]]$par[4]))
}
