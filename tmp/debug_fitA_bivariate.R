#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(spatstat)
  library(sf)
  library(tigris)
  library(data.table)
  library(dplyr)
})

script_dir <- normalizePath(file.path(getwd(), "inst", "oklahoma"), winslash = "/", mustWork = TRUE)
repo_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_dir <- file.path(script_dir, "oklahoma_induced_seismicity_data_regional20150318")
if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required")
pkgload::load_all(repo_dir, quiet = TRUE, export_all = FALSE, helpers = FALSE, attach_testthat = FALSE)

options(error = function() {
  cat("\n--- DEBUG TRACEBACK ---\n")
  traceback(50)
  quit(save = "no", status = 1)
})

ETAS_M0 <- 2.5
BETA_GR <- 2.3
CRS_PROJ <- 5070
VANILLA_MAXIT <- 20
VANILLA_START <- list(mu = 1.0, A = 0.2, alpha_m = 0.8, c = 0.05, p = 1.2, D = 5.0, gamma = 0.5, q = 1.5)

ev_all <- fread(file.path(data_dir, "events_all.csv"))
meta <- jsonlite::fromJSON(readLines(file.path(data_dir, "metadata.json")))
t_star_utc <- as.POSIXct(meta$design$t_star_utc, tz = "UTC")
ev_all[, time_utc := as.POSIXct(time_utc, tz = "UTC")]
ev_all[, t_days := as.numeric(difftime(time_utc, t_star_utc, units = "days"))]
ev_all[, x_km := x_m / 1000]
ev_all[, y_km := y_m / 1000]
ev_all <- ev_all[mag >= ETAS_M0]
post_end_days <- as.numeric(difftime(as.POSIXct(meta$design$post_end_utc, tz = "UTC"), t_star_utc, units = "days"))

counties_sf <- counties(state = "OK", cb = TRUE, year = 2022)
counties_sf <- st_transform(counties_sf, CRS_PROJ)
counties_sf <- st_make_valid(counties_sf)
ok_boundary <- st_make_valid(st_union(counties_sf))
bb <- st_bbox(ok_boundary)
win_km <- owin(xrange = c(bb["xmin"], bb["xmax"]) / 1000, yrange = c(bb["ymin"], bb["ymax"]) / 1000)

county_owins <- lapply(seq_len(nrow(counties_sf)), function(i) {
  geom <- st_geometry(counties_sf[i, ])
  coords <- st_coordinates(geom)
  x_km <- coords[, 1] / 1000
  y_km <- coords[, 2] / 1000
  tryCatch(owin(poly = list(x = rev(x_km), y = rev(y_km))), error = function(e) {
    tryCatch(owin(poly = list(x = x_km, y = y_km)), error = function(e2) NULL)
  })
})
valid_idx <- !sapply(county_owins, is.null)
county_owins_valid <- county_owins[valid_idx]
counties_sf_valid <- counties_sf[valid_idx, ]
names(county_owins_valid) <- counties_sf_valid$NAME
partition <- tess(tiles = county_owins_valid, window = win_km)

aoi_path <- file.path(data_dir, "occ_aoi_layer_2.geojson")
aoi_sf <- st_read(aoi_path, quiet = TRUE)
aoi_sf <- st_make_valid(st_transform(aoi_sf, CRS_PROJ))
aoi_union <- st_union(aoi_sf)
county_centroids <- st_centroid(counties_sf_valid)
inside_aoi <- lengths(st_within(county_centroids, aoi_union)) > 0
partition_processes <- ifelse(inside_aoi, "treated", "control")
treated_idx <- partition_processes == "treated"
control_ss <- as.owin(partition[!treated_idx])
treated_ss <- as.owin(partition[treated_idx])

assign_county <- function(df) {
  ti <- as.integer(tileindex(df$x, df$y, partition))
  df$location_process <- ifelse(is.na(ti), NA_character_, partition_processes[pmin(pmax(ti, 1), partition$n)])
  df$W <- 1.0
  df$n <- nrow(df)
  df$background <- TRUE
  df
}

pp_pre <- assign_county(as.data.frame(ev_all[t_days < 0, .(x = x_km, y = y_km, t = t_days, mag = mag)]))
pp_post <- assign_county(as.data.frame(ev_all[t_days >= 0 & t_days <= post_end_days, .(x = x_km, y = y_km, t = t_days, mag = mag)]))
pp_pre <- pp_pre[!is.na(pp_pre$location_process), ]
pp_post <- pp_post[!is.na(pp_post$location_process), ]
pp_pre$process <- "control"
pp_post$process <- pp_post$location_process
pp_pre$inferred_process <- "control"
pp_post$inferred_process <- pp_post$location_process
pp_all <- rbind(pp_pre, pp_post)
pp_all <- pp_all[order(pp_all$t), ]
windowT_fit <- c(min(pp_pre$t, na.rm = TRUE), post_end_days)

cat("n(pp_all)=", nrow(pp_all), " windowT=", paste(signif(windowT_fit, 6), collapse = ","), "\n", sep = "")
cat("anyNA x/y/t/mag/proc: ",
    anyNA(pp_all$x), anyNA(pp_all$y), anyNA(pp_all$t), anyNA(pp_all$mag), anyNA(pp_all$inferred_process), "\n")
cat("areas control/treated/full: ",
    signif(spatstat.geom::area(control_ss), 6), " / ",
    signif(spatstat.geom::area(treated_ss), 6), " / ",
    signif(spatstat.geom::area(win_km), 6), "\n", sep = "")

A_ctrl <- VANILLA_START
A_treat <- VANILLA_START
biv_init <- init_bivariate_from_independent(A_ctrl, A_treat)
biv_init[c("mu_0", "mu_1")] <- 1.0
biv_init[c("A_00", "A_11")] <- 0.2
biv_init[c("alpha_m_00", "alpha_m_11")] <- 0.8
biv_init[c("c", "p", "D", "gamma", "q")] <- unlist(VANILLA_START[c("c", "p", "D", "gamma", "q")])
fixed_structural <- as.list(biv_init[c("c", "p", "D", "gamma", "q")])

cat("Running one direct loglik evaluation...\n")
ll0 <- loglik_etas_bivariate(
  params = biv_init,
  realiz = pp_all,
  windowT = windowT_fit,
  windowS = win_km,
  m0 = ETAS_M0,
  control_state_space = control_ss,
  treated_state_space = treated_ss,
  treated_background_zero_before = 0,
  t_trunc = 0.556962
)
cat("Initial ll=", ll0, "\n", sep = "")

cat("Running fit_etas_bivariate (debug target: fit-job A)...\n")
fit <- fit_etas_bivariate(
  params_init = biv_init,
  realiz = pp_all,
  windowT = windowT_fit,
  windowS = win_km,
  m0 = ETAS_M0,
  control_state_space = control_ss,
  treated_state_space = treated_ss,
  treated_background_zero_before = 0,
  maxit = VANILLA_MAXIT,
  fixed_params = fixed_structural,
  trace = 1,
  t_trunc = 0.556962
)

cat("fit convergence=", fit$convergence, " value=", fit$value, "\n", sep = "")
print(fit$par)
