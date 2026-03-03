#!/usr/bin/env Rscript
# Deep profiling of EM-style labelling proposal generation.
# Run from package root: Rscript inst/sim_study/profile_em_proposals.R
# Or from inst/sim_study: Rscript profile_em_proposals.R
#
# Measures: per-iteration breakdown (proposals vs precomp vs loglik vs rbind),
#           per-call breakdown inside simulation_labeling_hawkes_hawkes_fast,
#           and Rprof flame data for the full EM loop.
#
# Speedups applied in em_style_labelling path:
#   - statespace converted to owin once at start (avoids repeated as.owin in precomp/GIH).
#   - post_inds = tileindex(post) computed once per proposal-update; passed as
#     points_tile_index to simulation_labeling_hawkes_hawkes_fast (saves n_props-1 tileindex calls).
#   - filt_by_proc = split(pre, pre$location_process) computed once per update; passed
#     to generate_inhomogeneous_hawkes via ... (saves 2*n_props split/tileindex on filtration).
# Standalone calls to simulation_labeling_hawkes_hawkes_fast (e.g. Section 5 or this script)
# do not pass these, so Rprof on "30x simulation_labeling" still shows tileindex/split cost.

suppressPackageStartupMessages({
  library(spatstat)
  library(data.table)
  library(dplyr)
})

if (file.exists("profile_em_proposals.R")) setwd("../..")
devtools::load_all(quiet = TRUE)

# ---------- Setup (match sim_study / test-speedup) ----------
set.seed(42)
OMEGA     <- c(0, 50, 0, 50)
END_TIME  <- 55
TREAT_T   <- 5
NX <- 5; NY <- 5
hp1 <- list(mu = 8, alpha = 0.01, beta = 10, K = 0.5)
hp2 <- list(mu = 8, alpha = 0.01, beta = 10, K = 0.2)

partition <- quadrats(X = OMEGA, nx = NX, ny = NY)
pp <- rep("control", partition$n)
pp[sample(1:25, 12)] <- "treated"
treated_idx <- pp == "treated"
control_ss  <- as.owin(partition[!treated_idx])
treated_ss  <- as.owin(partition[treated_idx])
state_spaces <- if (pp[1] == "treated") list(treated_ss, control_ss) else list(control_ss, treated_ss)

pre_sim <- sim_hawkes(hp1, c(0, TREAT_T), OMEGA)
pre_df  <- data.frame(x = pre_sim$x, y = pre_sim$y, t = pre_sim$t,
                      process = "control", location_process = "control",
                      inferred_process = "control", background = pre_sim$background)

post_sim <- generate_inhomogeneous_hawkes(OMEGA, partition, c(TREAT_T, END_TIME),
              pp, list(control = hp1, treated = hp2),
              state_spaces = state_spaces, filtration = pre_df)
post_df <- as.data.frame(post_sim)
post_df$inferred_process <- post_df$location_process
shared_cols <- c("x","y","t","process","location_process","inferred_process","background")
pre_df  <- pre_df[, shared_cols]
post_df <- post_df[, shared_cols]
obs_data <- rbind(pre_df, post_df)
obs_data$W <- 1

pre  <- obs_data[obs_data$t < TREAT_T, ]
post <- obs_data[obs_data$t >= TREAT_T, ]
post <- post[order(post$t), ]
pre$location_process <- "control"
pre$inferred_process <- NULL

n_post <- nrow(post)
n_pre  <- nrow(pre)
cat("Data: ", n_post, " post-treatment, ", n_pre, " pre-treatment points\n\n")

# ================================================================
# 1. Manual timing: one EM iteration broken into phases
# ================================================================
N_PROPS <- 10
iter <- 1
treatment_time <- TREAT_T
time_window <- c(TREAT_T, END_TIME)
is_pre <- obs_data$t < treatment_time
is_post <- !is_pre
starting_data <- as.data.frame(obs_data)
control_par <- list(hp1)
treated_par <- list(hp2)

# Phase A: proposal generation only (n_props calls to simulation_labeling_hawkes_hawkes_fast)
t_proposals <- system.time({
  post_proposals <- lapply(1:N_PROPS, function(j) {
    set.seed(iter * 10000 + j)
    simulation_labeling_hawkes_hawkes_fast(
      post, partition = partition, partition_process = pp,
      statespace = OMEGA, state_spaces = state_spaces,
      windowT = time_window,
      hawkes_params_control = control_par[[1]],
      hawkes_params_treated = treated_par[[1]],
      change_factor = 0.05, filtration = pre, proximity_weight = 0
    )
  })
})["elapsed"]
cat("=== 1. One EM iteration phase timing (n_props=", N_PROPS, ") ===\n", sep = "")
cat("  A. Proposal generation (n_props x simulation_labeling_fast): ", round(t_proposals, 3), "s\n", sep = "")

# Phase B: precompute_loglik_args (2 calls)
ref_post <- post
t_precomp <- system.time({
  pc_ctrl_all <- precompute_loglik_args(ref_post, OMEGA, treated_ss)
  pc_treat_all <- precompute_loglik_args(ref_post, OMEGA, control_ss)
})["elapsed"]
cat("  B. precompute_loglik_args (2 calls): ", round(t_precomp, 4), "s\n", sep = "")

# Phase C: rbind + build labelling_proposals
t_rbind <- system.time({
  pre$inferred_process <- "control"
  labelling_proposals <- lapply(post_proposals, function(tmp) rbind(pre, tmp))
})["elapsed"]
cat("  C. rbind(pre, tmp) x n_props: ", round(t_rbind, 4), "s\n", sep = "")

# Phase D: metric (vapply loglik over proposals)
ctrl_params_vec <- unlist(hp1)
treat_params_vec <- unlist(hp2)
t_metric <- system.time({
  metric <- vapply(labelling_proposals, function(y) {
    realiz <- y[is_post, ]
    ctrl_rows <- realiz$inferred_process == "control"
    if (!any(ctrl_rows)) return(-Inf)
    control_lik <- loglik_hawk_fast(
      params = ctrl_params_vec, realiz = realiz[ctrl_rows, ],
      windowT = time_window, windowS = OMEGA,
      precomp = list(active_area = pc_ctrl_all$active_area,
                    in_zero_bg = pc_ctrl_all$in_zero_bg_all[ctrl_rows])
    )
    treat_rows <- !ctrl_rows
    if (!any(treat_rows)) return(-Inf)
    treat_lik <- loglik_hawk_fast(
      params = treat_params_vec, realiz = realiz[treat_rows, ],
      windowT = time_window, windowS = OMEGA,
      precomp = list(active_area = pc_treat_all$active_area,
                    in_zero_bg = pc_treat_all$in_zero_bg_all[treat_rows])
    )
    control_lik + treat_lik
  }, numeric(1))
})["elapsed"]
cat("  D. Metric (loglik_hawk_fast x 2 x n_props): ", round(t_metric, 4), "s\n", sep = "")
cat("  Total one-iter (A+B+C+D): ", round(t_proposals + t_precomp + t_rbind + t_metric, 3), "s\n\n", sep = "")

# ================================================================
# 2. Per-call breakdown inside one simulation_labeling_hawkes_hawkes_fast
# ================================================================
cat("=== 2. Single simulation_labeling_hawkes_hawkes_fast breakdown ===\n")
set.seed(1)

# 2a. tileindex only
t_ti <- system.time({
  inds <- as.numeric(tileindex(post$x, post$y, partition))
})["elapsed"]
cat("  tileindex(post): ", round(t_ti * 1000, 2), " ms\n", sep = "")

# 2b. generate_inhomogeneous_hawkes control-only
t_gih_ctrl <- system.time({
  sim_data <- generate_inhomogeneous_hawkes(
    Omega = OMEGA, partition = partition, time_window = time_window,
    partition_processes = pp,
    hawkes_params = list(control = hp1, treated = list(mu = 0, alpha = 0, beta = 0, K = 0)),
    filtration = pre, state_spaces = state_spaces, space_triggering = FALSE
  )
})["elapsed"]
cat("  generate_inhomogeneous_hawkes (control only): ", round(t_gih_ctrl * 1000, 2), " ms\n", sep = "")

# 2c. generate_inhomogeneous_hawkes treated-only
t_gih_treat <- system.time({
  sim_data_t <- generate_inhomogeneous_hawkes(
    Omega = OMEGA, partition = partition, time_window = time_window,
    partition_processes = pp,
    hawkes_params = list(control = list(mu = 0, alpha = 0, beta = 0, K = 0), treated = hp2),
    filtration = pre, state_spaces = state_spaces, space_triggering = FALSE
  )
})["elapsed"]
cat("  generate_inhomogeneous_hawkes (treated only): ", round(t_gih_treat * 1000, 2), " ms\n", sep = "")

# 2d. full call
t_full_label <- system.time({
  out <- simulation_labeling_hawkes_hawkes_fast(
    post, partition = partition, partition_process = pp,
    statespace = OMEGA, state_spaces = state_spaces,
    windowT = time_window,
    hawkes_params_control = hp1, hawkes_params_treated = hp2,
    change_factor = 0.05, filtration = pre, proximity_weight = 0
  )
})["elapsed"]
cat("  Full simulation_labeling_hawkes_hawkes_fast: ", round(t_full_label * 1000, 2), " ms\n", sep = "")
cat("  (Sum of tileindex + 2x GIH approx: ", round((t_ti + t_gih_ctrl + t_gih_treat) * 1000, 2), " ms)\n\n", sep = "")

# ================================================================
# 3. Rprof: full EM run (few iters, many samples)
# ================================================================
cat("=== 3. Rprof on em_style_labelling (iter=3, n_props=10) ===\n")
prof_file <- tempfile()
Rprof(prof_file, interval = 0.002, memory.profiling = FALSE)
set.seed(999)
em_res <- em_style_labelling(
  pp_data = obs_data, partition = partition,
  partition_processes = pp, statespace = OMEGA,
  time_window = time_window, treatment_time = TREAT_T,
  hawkes_params_control = hp1, hawkes_params_treated = hp2,
  update_control_params = FALSE, param_update_cadence = 10,
  iter = 3, n_props = N_PROPS, change_factor = 0.05,
  state_spaces = state_spaces, verbose = FALSE
)
Rprof(NULL)
s <- summaryRprof(prof_file)
cat("Top 25 by self time:\n")
print(head(s$by.self, 25))
cat("\nTop 25 by total time:\n")
print(head(s$by.total, 25))
unlink(prof_file)

# ================================================================
# 4. Rprof: proposal generation only (many calls)
# ================================================================
cat("\n=== 4. Rprof on 30x simulation_labeling_hawkes_hawkes_fast ===\n")
prof_file2 <- tempfile()
Rprof(prof_file2, interval = 0.002, memory.profiling = FALSE)
for (i in 1:30) {
  set.seed(i)
  simulation_labeling_hawkes_hawkes_fast(
    post, partition = partition, partition_process = pp,
    statespace = OMEGA, state_spaces = state_spaces,
    windowT = time_window,
    hawkes_params_control = hp1, hawkes_params_treated = hp2,
    change_factor = 0.05, filtration = pre, proximity_weight = 0
  )
}
Rprof(NULL)
s2 <- summaryRprof(prof_file2)
cat("Top 20 by self time:\n")
print(head(s2$by.self, 20))
cat("\nTop 20 by total time:\n")
print(head(s2$by.total, 20))
unlink(prof_file2)

# ================================================================
# 5. Full EM benchmark (repeat for stability)
# ================================================================
cat("\n=== 5. Full em_style_labelling benchmark (iter=5, n_props=10, 2 runs) ===\n")
em_times <- numeric(2)
for (run in 1:2) {
  set.seed(1000 + run)
  em_times[run] <- system.time({
    em_style_labelling(
      pp_data = obs_data, partition = partition,
      partition_processes = pp, statespace = OMEGA,
      time_window = time_window, treatment_time = TREAT_T,
      hawkes_params_control = hp1, hawkes_params_treated = hp2,
      update_control_params = FALSE, param_update_cadence = 5,
      iter = 5, n_props = N_PROPS, change_factor = 0.05,
      state_spaces = state_spaces, verbose = FALSE
    )
  })["elapsed"]
}
cat("  Run 1: ", round(em_times[1], 2), "s, Run 2: ", round(em_times[2], 2), "s\n", sep = "")
cat("  Mean: ", round(mean(em_times), 2), "s\n", sep = "")
