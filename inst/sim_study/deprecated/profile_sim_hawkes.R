#!/usr/bin/env Rscript
# Profile sim_hawkes_fast and generate_inhomogeneous_hawkes in an EM context.
# Run from package root:  Rscript inst/sim_study/profile_sim_hawkes.R
# Or from inst/sim_study: Rscript profile_sim_hawkes.R
suppressPackageStartupMessages({
  library(spatstat)
  library(data.table)
})

if (file.exists("profile_sim_hawkes.R")) {
  setwd("../..")
}
devtools::load_all(quiet = TRUE)

# ---------- setup (same as test-speedup.R make_test_data) ----------
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

cat("Data size: ", nrow(post), " post-treatment points,",
    nrow(pre), " pre-treatment points\n\n")

# ================================================================
# 1. Micro-benchmark: sim_hawkes_fast in isolation
# ================================================================
cat("=== 1. sim_hawkes_fast micro-benchmark (100 reps) ===\n")
N_REP <- 100

# (a) background-only (K=0) on non-rectangular state space
t_bg <- system.time(for (i in 1:N_REP) {
  set.seed(i)
  sim_hawkes_fast(list(mu = 8, alpha = 0.01, beta = 10, K = 0),
                  c(TREAT_T, END_TIME), control_ss)
})["elapsed"]
cat(sprintf("  bg-only  (K=0, non-rect)   : %.4fs total, %.4fms/call\n",
            t_bg, t_bg / N_REP * 1000))

# (b) full sim on Omega (rectangle)
bg_data <- sim_hawkes_fast(list(mu = 8, alpha = 0.01, beta = 10, K = 0),
                           c(TREAT_T, END_TIME), control_ss)
t_full <- system.time(for (i in 1:N_REP) {
  set.seed(i)
  sim_hawkes_fast(hp1, c(TREAT_T, END_TIME), OMEGA,
                  background_realization = bg_data, filtration = pre)
})["elapsed"]
cat(sprintf("  full sim (rect, bg+filt)    : %.4fs total, %.4fms/call\n",
            t_full, t_full / N_REP * 1000))

# (c) full sim on non-rectangular control_ss
t_full_nr <- system.time(for (i in 1:N_REP) {
  set.seed(i)
  sim_hawkes_fast(hp1, c(TREAT_T, END_TIME), control_ss,
                  background_realization = bg_data, filtration = pre)
})["elapsed"]
cat(sprintf("  full sim (non-rect, bg+filt): %.4fs total, %.4fms/call\n",
            t_full_nr, t_full_nr / N_REP * 1000))

# (d) zero-params (mu=0, K=0)
t_zero <- system.time(for (i in 1:N_REP) {
  sim_hawkes_fast(list(mu = 0, alpha = 0, beta = 0, K = 0),
                  c(TREAT_T, END_TIME), treated_ss)
})["elapsed"]
cat(sprintf("  zero-params (mu=0, K=0)    : %.4fs total, %.4fms/call\n\n",
            t_zero, t_zero / N_REP * 1000))

# ================================================================
# 2. Micro-benchmark: generate_inhomogeneous_hawkes
# ================================================================
cat("=== 2. generate_inhomogeneous_hawkes micro-benchmark (50 reps) ===\n")
N_REP2 <- 50
no_hp <- list(mu = 0, alpha = 0, beta = 0, K = 0)

t_gih_ctrl <- system.time(for (i in 1:N_REP2) {
  set.seed(i)
  generate_inhomogeneous_hawkes(OMEGA, partition, c(TREAT_T, END_TIME), pp,
    list(control = hp1, treated = no_hp),
    state_spaces = state_spaces, filtration = pre)
})["elapsed"]
cat(sprintf("  control-only : %.4fs total, %.4fms/call\n",
            t_gih_ctrl, t_gih_ctrl / N_REP2 * 1000))

t_gih_both <- system.time(for (i in 1:N_REP2) {
  set.seed(i)
  generate_inhomogeneous_hawkes(OMEGA, partition, c(TREAT_T, END_TIME), pp,
    list(control = hp1, treated = hp2),
    state_spaces = state_spaces, filtration = pre)
})["elapsed"]
cat(sprintf("  both procs   : %.4fs total, %.4fms/call\n\n",
            t_gih_both, t_gih_both / N_REP2 * 1000))

# ================================================================
# 3. Micro-benchmark: simulation_labeling_hawkes_hawkes_fast
# ================================================================
cat("=== 3. simulation_labeling_hawkes_hawkes_fast (20 reps) ===\n")
N_REP3 <- 20

t_label_ctrl <- system.time(for (i in 1:N_REP3) {
  set.seed(i)
  simulation_labeling_hawkes_hawkes_fast(
    post, partition = partition, partition_process = pp,
    statespace = OMEGA, state_spaces = state_spaces,
    windowT = c(TREAT_T, END_TIME),
    hawkes_params_control = hp1, hawkes_params_treated = NULL,
    change_factor = 0.05, filtration = pre, proximity_weight = 0)
})["elapsed"]
cat(sprintf("  ctrl-only  : %.4fs total, %.4fms/call\n",
            t_label_ctrl, t_label_ctrl / N_REP3 * 1000))

t_label_both <- system.time(for (i in 1:N_REP3) {
  set.seed(i)
  simulation_labeling_hawkes_hawkes_fast(
    post, partition = partition, partition_process = pp,
    statespace = OMEGA, state_spaces = state_spaces,
    windowT = c(TREAT_T, END_TIME),
    hawkes_params_control = hp1, hawkes_params_treated = hp2,
    change_factor = 0.05, filtration = pre, proximity_weight = 0)
})["elapsed"]
cat(sprintf("  ctrl+treat : %.4fs total, %.4fms/call\n\n",
            t_label_both, t_label_both / N_REP3 * 1000))

# ================================================================
# 4. Mini EM benchmark (3 runs averaged)
# ================================================================
cat("=== 4. Mini em_style_labelling (iter=3, n_props=10, avg of 3 runs) ===\n")
em_times <- numeric(3)
for (run in 1:3) {
  set.seed(999 + run)
  em_times[run] <- system.time({
    em_res <- em_style_labelling(
      pp_data = obs_data, partition = partition,
      partition_processes = pp, statespace = OMEGA,
      time_window = c(TREAT_T, END_TIME), treatment_time = TREAT_T,
      hawkes_params_control = hp1, hawkes_params_treated = hp2,
      update_control_params = FALSE, param_update_cadence = 3,
      iter = 3, n_props = 10, change_factor = 0.05,
      state_spaces = state_spaces, verbose = FALSE)
  })["elapsed"]
}
t_em <- mean(em_times)
cat(sprintf("  runs: %.2fs, %.2fs, %.2fs\n", em_times[1], em_times[2], em_times[3]))
cat(sprintf("  avg total elapsed: %.2fs\n", t_em))
cat(sprintf("  per EM iter  : %.2fs\n", t_em / 3))
cat(sprintf("  per proposal : %.4fs\n\n", t_em / 30))

# ================================================================
# 5. Rprof on labeling call
# ================================================================
cat("=== 5. Rprof breakdown (20 labeling proposals) ===\n")
prof_file <- tempfile()
Rprof(prof_file, interval = 0.005, memory.profiling = FALSE)
for (i in 1:20) {
  set.seed(i + 500)
  simulation_labeling_hawkes_hawkes_fast(
    post, partition = partition, partition_process = pp,
    statespace = OMEGA, state_spaces = state_spaces,
    windowT = c(TREAT_T, END_TIME),
    hawkes_params_control = hp1, hawkes_params_treated = hp2,
    change_factor = 0.05, filtration = pre, proximity_weight = 0)
}
Rprof(NULL)
s <- summaryRprof(prof_file)
cat("\nTop 20 by self time:\n")
print(head(s$by.self, 20))
cat("\nTop 20 by total time:\n")
print(head(s$by.total, 20))

# ================================================================
# 6. Seed-based equivalence verification
# ================================================================
cat("\n=== 6. Seed-based equivalence ===\n")
ref_file <- "/tmp/ppd_reference_outputs.rds"
if (file.exists(ref_file)) {
  ref <- readRDS(ref_file)
  set.seed(100)
  chk1 <- sim_hawkes_fast(hp1, c(TREAT_T, END_TIME), control_ss)
  cat("sim_hawkes_fast x/y/t/bg match: ",
      identical(ref$ref_shf$x, chk1$x) && identical(ref$ref_shf$y, chk1$y) &&
      identical(ref$ref_shf$t, chk1$t) && identical(ref$ref_shf$background, chk1$background), "\n")

  no_hp <- list(mu = 0, alpha = 0, beta = 0, K = 0)
  set.seed(300)
  chk2 <- generate_inhomogeneous_hawkes(OMEGA, partition, c(TREAT_T, END_TIME), pp,
           list(control = hp1, treated = no_hp), state_spaces = state_spaces, filtration = pre)
  cat("generate_inhom(ctrl) x/t/proc match: ",
      identical(ref$ref_gih$x, chk2$x) && identical(ref$ref_gih$t, chk2$t) &&
      identical(ref$ref_gih$process, chk2$process), "\n")

  set.seed(400)
  chk3 <- generate_inhomogeneous_hawkes(OMEGA, partition, c(TREAT_T, END_TIME), pp,
           list(control = hp1, treated = hp2), state_spaces = state_spaces, filtration = pre)
  cat("generate_inhom(both) x/t match: ",
      identical(ref$ref_gih2$x, chk3$x) && identical(ref$ref_gih2$t, chk3$t), "\n")
} else {
  cat("No reference file found. Run once before optimizing to capture reference.\n")
}

# ================================================================
# Baseline comparison (from pre-optimization run)
# ================================================================
cat("\n=== COMPARISON (pre-opt baselines in parens) ===\n")
cat("sim_hawkes_fast:\n")
cat(sprintf("  bg-only(K=0,non-rect) : now=%.2fms  (was 1.03ms)\n", t_bg / N_REP * 1000))
cat(sprintf("  full(rect,bg+filt)    : now=%.2fms  (was 0.80ms)\n", t_full / N_REP * 1000))
cat(sprintf("  full(non-rect,bg+filt): now=%.2fms  (was 1.22ms)\n", t_full_nr / N_REP * 1000))
cat(sprintf("  zero-params           : now=%.2fms  (was 0.27ms)  [%.1fx]\n", t_zero / N_REP * 1000, 0.27 / (t_zero / N_REP * 1000)))
cat("generate_inhomogeneous_hawkes:\n")
cat(sprintf("  control-only : now=%.2fms  (was 8.92ms)\n", t_gih_ctrl / N_REP2 * 1000))
cat(sprintf("  both procs   : now=%.2fms  (was 13.08ms)  [%.1fx]\n", t_gih_both / N_REP2 * 1000, 13.08 / (t_gih_both / N_REP2 * 1000)))
cat("simulation_labeling_hawkes_hawkes_fast:\n")
cat(sprintf("  ctrl-only  : now=%.2fms  (was 18.35ms)\n", t_label_ctrl / N_REP3 * 1000))
cat(sprintf("  ctrl+treat : now=%.2fms  (was 15.95ms)  [%.1fx]\n", t_label_both / N_REP3 * 1000, 15.95 / (t_label_both / N_REP3 * 1000)))
cat(sprintf("EM (avg) : now=%.2fs  (baseline was 1.14s single run)\n", t_em))
