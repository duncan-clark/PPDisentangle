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

# 1 year in vaccines were administered
TREATMENT_TIME <- 365
STUDY_END <- 365 * 2 # 2 years of data # to make smaller
# Define treatment as "close to health center" distance in km
HEALTH_CENTER_DIST <- 500

# Get data:
gambia_data <- dget(system.file("extdata/gambia/DISCLEAN.txt", package = "PPDisentangle"))
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
  t = as.numeric(as.Date(gambia_data$date_IPD_Pne_Event)),
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

# jitter the points ~ 1 units
# probably not right but don't want to have to deal with heterogeneity in background rate:
pp_data$x <- pp_data$x + rnorm(nrow(pp_data), mean = 0, sd = 1)
pp_data$y <- pp_data$y + rnorm(nrow(pp_data), mean = 0, sd = 1)
plot_pp(pp_data,partition)
pp_data <- pp_data[inside.owin(pp_data$x, pp_data$y, win),]

# =============================================
# Get inhomogeneous background rate
# =============================================

# function to map northing/easting to lat/long
pop_raster <- raster("gambia_disease/Gambia_100m_Population/GMB_ppp_v2b_2020.tif")
pp_data$W <- raster::extract(pop_raster, data.frame(x=pp_data$long, y=pp_data$lat))
# drop the lat and lon:
# pp_data$lat <- NULL
# pp_data$long <- NULL

# seem odd that most cases are happening in low pop density regions?
# are we sure mapping is happening correctly?
hist(pp_data$W)

compare_to_background <- function(
    pts_df,
    bg_raster,
    coord_cols = c("x","y"),
    plot       = TRUE,
    fact = NULL
) {
  # 1) make SpatialPoints
  if (!all(coord_cols %in% names(pts_df))) {
    warning("coord_cols not all in pts_df")
    return(NULL)
  }
  # aggregate the raster if needed:
  if(!is.null(fact)){
    bg_raster <- raster::aggregate(bg_raster,
                                   fact = fact,            # horizontal and vertical factor
                                   fun  = sum,         # how to combine the background values
                                   expand = TRUE,       # keep full extent
                                   na.rm = TRUE)
  }

  pts_sp <- SpatialPoints(
    coords         = pts_df[, coord_cols],
    proj4string    = crs(bg_raster)
  )

  # 2) rasterize points -> case counts
  case_ras <- rasterize(
    pts_sp,
    bg_raster,
    fun        = "count",
    background = 0
  )

  # 3) extract vectors
  pop_vals    <- getValues(bg_raster)
  case_counts <- getValues(case_ras)

  # 4) build & clean data.frame
  df <- data.frame(
    pop   = pop_vals,
    cases = case_counts
  )
  df <- na.omit(df)

  # 5) optionally plot
  lm_model <- NULL
  if (plot) {
    # par(mfrow = c(1,2))
    # hist(df$pop,   main="Background rate", xlab="rate per cell")
    # hist(df$cases, main="Case counts",     xlab="cases per cell")
    par(mfrow = c(1,1))

    non0 <- df[df$cases > 0, ]
    plot(
      non0$pop, non0$cases,
      log   = "x",
      pch   = 20,
      xlab  = "Background rate (log scale)",
      ylab  = "Case count",
      main  = "Cases vs Background"
    )
    lm_model <- lm(cases ~ pop, data = df)
    abline(lm_model, col = "red", lwd = 2)
  }

  # 6) correlation test
  non0 <- df[df$cases > 0, ]
  cor_test <- cor.test(non0$pop, non0$cases, method = "spearman")

  # 7) return
  out <- list(
    cor_test = cor_test
  )
  if (!is.null(lm_model)) out$lm_model <- lm_model
  return(out)
}

# with raster aggregation - looks a bit better.- maybe working :)
compare_to_background(
  pts_df = pp_data,
  bg_raster = pop_raster,
  coord_cols = c("long", "lat"),
  plot = TRUE,
  fact = 10
)

pop_geo <- rast("gambia_disease/Gambia_100m_Population/GMB_ppp_v2b_2020.tif")
# 2. (Re‑)declare its CRS as WGS84 lon/lat, if not already set
crs(pop_geo) <- "EPSG:4326"
# 3. Project to UTM Zone 28 N (EPSG:32628)
#    – use bilinear interpolation for a continuous surface
#    – set res=c(100,100) to target ~100 m pixels (optional)
pop_utm <- project(
  pop_geo,
  "EPSG:32628",
  method = "bilinear",
)
pop_utm_raster <- raster(pop_utm)
compare_to_background(
  pts_df = pp_data,
  bg_raster = pop_utm_raster,
  coord_cols = c("x", "y"),
  plot = TRUE,
  fact = 10
)

tmp <- raster::extract(pop_utm, data.frame(x=pp_data$x, y=pp_data$y))
tmp1 <- cbind(tmp,pp_data$W)
plot(tmp[,1],tmp1[,2])
# NOT AT ALL CORRELATED OH NO THIS BAD :(

# experiment just replacing - everything on UTM - not as good
# pp_data$W <- raster::extract(pop_utm_raster, data.frame(x=pp_data$x, y=pp_data$y))
# compare_to_background(
#   pts_df = pp_data,
#   bg_raster = pop_utm_raster,
#   coord_cols = c("x", "y"),
#   plot = TRUE,
#   fact = 10
# )


# check we can simulate with the new raster:
tmp <- sim_hawkes(params = list(mu = 1e-10,
                                alpha = 0.1,
                                beta = 0.1,
                                K = 0.1),
                  windowT = c(0, STUDY_END),
                  windowS = win,
                  optimized = TRUE,
                  covariate_lookup =  pop_utm_raster,
                  background_rate_var = 'W'
                  )

# check if the points in tmp are close to "where they should be"?
compare_to_background(
  pts_df = as.data.frame(tmp),
  bg_raster = pop_utm_raster,
  coord_cols = c("x", "y"),
  plot = TRUE,
  fact = 50
)
# looks good!


# function to take in pp data with marks and divide by the total mass so that
# they integrate to one over the mark space.
normalize_marks <- function(pp_data,
                            win,
                            mark_name = "W"
                            ){
  min_non_zero <- min(pp_data[[mark_name]][pp_data[[mark_name]] > 0], na.rm = TRUE)
  pp_data[[mark_name]][is.na(pp_data[[mark_name]]) | pp_data[[mark_name]] ==0] <- min_non_zero
  X = ppp(x = pp_data$x,
          y = pp_data$y,
          marks = pp_data[[mark_name]],
          window = win)
  im_smooth <- Smooth(X,
                      sigma = 2000,
                      edge=TRUE,
                      dimyx=500
                      )
  total_mass <- integral.im(im_smooth)
  pp_data[[mark_name]] <- pp_data[[mark_name]]/total_mass
  min_non_zero <- min(pp_data[[mark_name]][pp_data[[mark_name]] > 0],na.rm = TRUE)
  pp_data[[mark_name]][is.na(pp_data[[mark_name]]) | pp_data[[mark_name]] ==0] <- min_non_zero
  return(list(new_pp_data = pp_data,
              total_mass = total_mass))
}
norm_result <- normalize_marks(pp_data,
                            win = win,
                            mark_name = "W"
                            )
norm_result$total_mass
summary(pp_data$W)
pp_data <- norm_result$new_pp_data
summary(pp_data$W)
# need to divide the raster by the total mass also:
pop_utm_raster <- raster::calc(pop_raster, function(x) x / norm_result$total_mass)

norm_result <- normalize_marks(pp_data,
                               win = win,
                               mark_name = "W"
                               )
norm_result$total_mass


# =============================================

# reset the location process after the jittering:
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

# plot the pre treatment data
plot_pp(pp_data %>% filter(t < TREATMENT_TIME),
        partition,
        "pre_treatment")
pp_data %>%
  group_by(location_process,process) %>%
  summarise(n = n())
# =============================================

# case count over time plot:
# weekly
TIME_PERIOD = 7
pp_data %>%
  mutate(t_bin = floor(t / TIME_PERIOD) * TIME_PERIOD) %>%
  group_by(t_bin, location_process) %>%
  summarise(n = n(), .groups = "drop") %>%
  bind_rows(
    pp_data %>%
      mutate(t_bin = floor(t / TIME_PERIOD) * TIME_PERIOD) %>%
      group_by(t_bin) %>%
      summarise(n = n(), .groups = "drop") %>%
      mutate(location_process = "Total")
  ) %>%
  ggplot(aes(x = t_bin, y = n, color = location_process)) +
  geom_line() +
  geom_vline(xintercept = floor(TREATMENT_TIME / 7) * 7,
             color = "red",
             linetype = "dashed") +
  labs(title = paste0(TIME_PERIOD,"-Day Case Counts Over Time"),
       x = paste0("Time (days, binned into ",TIME_PERIOD,"-day intervals)"),
       y = "Case count",
       color = "Group") +
  theme_minimal()

# IPD cases only - monthly
TIME_PERIOD = 30
keep <- gambia_data$IPDEvent=='Y'
pp_data[keep,] %>%
  mutate(t_bin = floor(t / TIME_PERIOD) * TIME_PERIOD) %>%
  group_by(t_bin, location_process) %>%
  summarise(n = n(), .groups = "drop") %>%
  bind_rows(
    pp_data[keep,] %>%
      mutate(t_bin = floor(t / TIME_PERIOD) * TIME_PERIOD) %>%
      group_by(t_bin) %>%
      summarise(n = n(), .groups = "drop") %>%
      mutate(location_process = "Total")
  ) %>%
  ggplot(aes(x = t_bin, y = n, color = location_process)) +
  geom_line() +
  geom_vline(xintercept = floor(TREATMENT_TIME / TIME_PERIOD) * TIME_PERIOD,
             color = "red",
             linetype = "dashed") +
  labs(title = paste0(TIME_PERIOD,"-Day Case Counts Over Time"),
       x =  paste0("Time (days, binned into ",TIME_PERIOD,"-day intervals)"),
       y = "Case count",
       color = "Group") +
  ylim(c(0,10)) +
  theme_minimal()

# make different data so it is properly normalized:
pre_treated <- normalize_marks(pp_data %>% filter(t < TREATMENT_TIME, location_process == 'treated'),win = treated_state_space,mark_name = "W")$new_pp_data
pre_control <- normalize_marks(pp_data %>% filter(t < TREATMENT_TIME, location_process == 'control'),win = control_state_space,mark_name = "W")$new_pp_data
post_treated <- normalize_marks(pp_data %>% filter(t > TREATMENT_TIME, t < 365*2, location_process == 'treated'),win = treated_state_space,mark_name = "W")$new_pp_data
post_control <- normalize_marks(pp_data %>% filter(t > TREATMENT_TIME, t < 365*2, location_process == 'control'),win = control_state_space,mark_name = "W")$new_pp_data

pre_all <- normalize_marks(pp_data %>% filter(t < TREATMENT_TIME),win = win,mark_name = "W")$new_pp_data
post_all <- normalize_marks(pp_data %>% filter(t > TREATMENT_TIME),win = win,mark_name = "W")$new_pp_data

pre_treatment_fit <- fit_hawkes(params_init = c(10, 0.1, 0.1,0.1),
                                realiz = pre_control,
                                windowT = c(0, TREATMENT_TIME),
                                windowS = win,
                                trace = 1,
                                maxit = 2000,
                                poisson_flag = FALSE,
                                zero_background_region = NULL,
                                background_rate_var = 'W',
                                )
pre_treat_estim <- setNames(as.list(pre_treatment_fit$par), c('mu','alpha', 'beta', 'K'))
pre_treat_estim

# expected squared distance:
dist <- 0.5*sqrt(pi/pre_treat_estim$alpha)
print(dist) # ~300

# only on treated window - so heterogeneity cancels out ?
pre_treatment_fit_treated <- fit_hawkes(params_init = c(10, 0.1, 0.1,0.1),
                                realiz = pre_treated,
                                windowT = c(0, TREATMENT_TIME),
                                windowS = treated_state_space,
                                trace = 1,
                                maxit =1000,
                                poisson_flag = FALSE,
                                zero_background_region = NULL
                                #background_rate_var = 'W'
                                )

pre_treat_estim_treated <- setNames(as.list(pre_treatment_fit_treated$par), c('mu','alpha', 'beta', 'K'))
pre_treat_estim_treated

# only on treated window - so heterogeneity cancels out ?
pre_treatment_fit_control <- fit_hawkes(params_init = c(10, 0.1, 0.1,0.1),
                                        realiz = pre_control,
                                        windowT = c(0, TREATMENT_TIME),
                                        windowS = win,
                                        trace = 1,
                                        maxit =1000,
                                        poisson_flag = FALSE,
                                        zero_background_region = treated_state_space
                                        #background_rate_var = 'W'
)

pre_treat_estim_control <- setNames(as.list(pre_treatment_fit_control$par), c('mu','alpha', 'beta', 'K'))
pre_treat_estim_control

plot_pp(pp_data %>% filter(t > TREATMENT_TIME, t < 365*2),
        partition,
        "post_treatment")

post_treatment_fit_treated <- fit_hawkes(params_init = c(10, 0.1, 0.1,0.1),
                                realiz = post_treated,
                                windowT = c(TREATMENT_TIME,max(pp_data$t)),
                                windowS = win,
                                trace = 1,
                                maxit =1000,
                                poisson_flag = FALSE,
                                zero_background_region = control_state_space
                                #background_rate_var = 'W'
                                )

post_treat_estim_treated <- setNames(as.list(post_treatment_fit_treated$par), c('mu','alpha', 'beta', 'K'))
post_treat_estim_treated

post_treatment_fit_control <- fit_hawkes(params_init = c(10, 0.1, 0.1,0.1),
                                        realiz = post_control,
                                        windowT = c(TREATMENT_TIME,max(pp_data$t)),
                                        windowS = win,
                                        trace = 1,
                                        maxit =1000,
                                        poisson_flag = FALSE,
                                        zero_background_region = treated_state_space
                                        #background_rate_var = 'W'
)

post_treat_estim_control <- setNames(as.list(post_treatment_fit_control$par), c('mu','alpha', 'beta', 'K'))
post_treat_estim_control

post_treat_fit <- fit_hawkes(params_init = c(10, 0.1, 0.1,0.1),
                                realiz = post_treated,
                                windowT = c(TREATMENT_TIME,max(pp_data$t)),
                                windowS = win,
                                trace = 1,
                                maxit =1000,
                                poisson_flag = FALSE,
                                zero_background_region = NULL
                                #background_rate_var = 'W'
)

post_treat_estim <- setNames(as.list(post_treat_fit$par), c('mu','alpha', 'beta', 'K'))
post_treat_estim

results <- data.frame(
  "pre_treat" = unlist(pre_treat_estim),
  "post_treat" = unlist(post_treat_estim),
  "pre_treat_treated" = unlist(pre_treat_estim_treated),
  "post_treat_treated" = unlist(post_treat_estim_treated),
  "pre_treat_control" = unlist(pre_treat_estim_control),
  "post_treat_control" = unlist(post_treat_estim_control)
)
round(results,3)

# GOT TO HERE :
# WHEN WE FIT OVER DIFFERENT REGIONS
# NEED TO MAKE SURE BACKGROUND RATE INTEGRATES TO 1 !

# do the full EM fit:
total_points <- sum(pp_data$location_process == 'treated' & pp_data$t >= TREATMENT_TIME)
mu_start <- (total_points / (STUDY_END))/mean(pp_data$W)
params_init <- list(mu = mu_start,
                   alpha = 0.1,
                   beta = STUDY_END/10,
                   K = 0.1)

# since we don't know the true process can't monitor accuracy in the same way as sim study.
result <- adaptive_SEM(pp_data = pp_data %>% dplyr::select(x,y,t,location_process,process),
                       partition = partition,
                       partition_processes = partition_processes,
                       statespace = win,
                       time_window = c(TREATMENT_TIME,STUDY_END),
                       treatment_time = TREATMENT_TIME,
                       hawkes_params_control = pre_treat_estim_control,
                       hawkes_params_treated = pre_treat_estim_treated,
                       N_labellings = 10,
                       N_iter = 2,
                       verbose = T,
                       # covariate_lookup = pop_utm_raster,
                       # background_rate_var = 'W',
                       adaptive_control = list(
                         update_control_params = TRUE,
                         param_update_cadence = 10,
                         proposal_update_cadence = 1,
                         state_spaces = list(control_state_space,treated_state_space),
                         iter  = 100,
                         n_props = 20,
                         change_factor = 0.02,
                         include_starting_data = F,
                         update_starting_data = T,
                         verbose = T
                         )
                       )

result$adaptive$treated_par[[length(result$adaptive$treated_par)]]
result$adaptive$control_par[[length(result$adaptive$control_par)]]

results$SEM_treated <- unlist(result$adaptive$treated_par[[length(result$adaptive$treated_par)]])
results$SEM_control <- unlist(result$adaptive$control_par[[length(result$adaptive$control_par)]])

result$adaptive$control_par[[length(result$adaptive$control_par)]]


# think about this
# there seems to be little spillover from control to treated - since low background rate
# so I think we are overestimating in control region ?


# round so nicer
results <- round(results, 5)
results

results %>% dplyr::select(pre_treat,
                          pre_treat_treated,
                          pre_treat_control,
                          post_treat_treated,
                          post_treat_control,
                          SEM_treated,
                          SEM_control)

# NEED TO TESSELATE INTO REGIONS - so the thinning works better
# ALSO COMPARE TO ESTIMATED IN POST TREATMENT ON THE TREATED REGIONS

# DO THE LOUIS METHOD:

# =================
# IDEAS:
# ================

# - Fit with no covariates
# - Fit with only near and far from healh center e.g. within 5 miles, 10 - 5 miles.


