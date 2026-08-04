#!/usr/bin/env Rscript
# Diagnose why bivariate AoN ATE differs so much from marginal ATE.
# Compares several controlled simulators with the same seeds / n_sims.

suppressPackageStartupMessages({
  # run from PPDisentangle/ or inst/oklahoma/
})

args <- commandArgs(trailingOnly = TRUE)
script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- if (length(script_file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[1]), winslash = "/", mustWork = FALSE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
REPO_DIR <- if (basename(SCRIPT_DIR) == "oklahoma" && basename(dirname(SCRIPT_DIR)) == "inst") {
  normalizePath(dirname(dirname(SCRIPT_DIR)), winslash = "/", mustWork = FALSE)
} else if (file.exists(file.path(SCRIPT_DIR, "DESCRIPTION"))) {
  SCRIPT_DIR
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

pkgload::load_all(REPO_DIR, quiet = TRUE, export_all = FALSE, helpers = FALSE, attach_testthat = FALSE)
source(file.path(REPO_DIR, "inst", "oklahoma", "ate_bivariate.R"), local = FALSE)

OUT_DIR <- PPDisentangle::pp_output_path("oklahoma", repo_root = REPO_DIR)
rds_path <- file.path(OUT_DIR, "for_paper.rds")
DATA_DIR <- file.path(REPO_DIR, "inst", "oklahoma", "oklahoma_induced_seismicity_data_regional20150318")

n_sims <- 200L
if ("--quick" %in% args) n_sims <- 50L

res <- readRDS(rds_path)
cfg <- res$config
m0 <- cfg$ETAS_M0
beta_gr <- cfg$BETA_GR
t_trunc <- cfg$SEM_T_TRUNC_DAYS
window_days <- cfg$ATE_WINDOW_DAYS
windowT <- c(0, window_days)
pp_pre <- res$pp_data$pp_pre
filtration <- pp_pre[, c("x", "y", "t", "mag", "inferred_process"), drop = FALSE]

cat("Building geometry...\n")
geom <- oklahoma_build_state_spaces(DATA_DIR)
win <- geom$win_km
ss_obs <- geom$state_spaces
ss_all_ctrl <- list(control = win, treated = NULL)
ss_all_treat <- list(control = NULL, treated = win)

extract_marginals_local <- function(biv_par) {
  structural <- c(c = biv_par[["c"]], p = biv_par[["p"]], D = biv_par[["D"]],
                  gamma = biv_par[["gamma"]], q = biv_par[["q"]])
  list(
    ctrl = as.list(c(mu = biv_par[["mu_0"]], A = biv_par[["A_00"]],
                     alpha_m = biv_par[["alpha_m_00"]], structural)),
    treat = as.list(c(mu = biv_par[["mu_1"]], A = biv_par[["A_11"]],
                      alpha_m = biv_par[["alpha_m_11"]], structural))
  )
}

zero_cross <- function(p) {
  p <- as.list(unlist(p))
  p$A_01 <- 0; p$A_10 <- 0
  p$alpha_m_01 <- 0; p$alpha_m_10 <- 0
  p
}

sim_biv <- function(params, ss, seed) {
  set.seed(seed)
  sim_etas_bivariate(
    params = params, windowT = windowT, windowS = win,
    state_spaces = ss, m0 = m0, beta_gr = beta_gr,
    filtration = filtration, t_trunc = t_trunc
  )
}

sim_univ <- function(params, seed) {
  set.seed(seed)
  sim_etas(params, windowT, windowS = win, m0 = m0, beta_gr = beta_gr,
           filtration = filtration[, c("x", "y", "t", "mag")], t_trunc = t_trunc)
}

mean_n <- function(xs) mean(vapply(xs, function(z) if (is.data.frame(z)) nrow(z) else length(z$t), numeric(1)))

run_scenarios <- function(letter) {
  p <- as.list(unlist(res$fits_named[[letter]]$params))
  marg <- extract_marginals_local(p)
  p0 <- zero_cross(p)
  cat(sprintf("\n========== Fit %s ==========\n", letter))
  cat(sprintf("mu0=%.3f mu1=%.3f A00=%.3f A11=%.3f A01=%.4f A10=%.4f\n",
              p$mu_0, p$mu_1, p$A_00, p$A_11, p$A_01, p$A_10))
  cat(sprintf("Immigrant expectation over %.0f days (Poisson mean, ignoring triggering):\n", window_days))
  cat(sprintf("  only mu0: %.1f | only mu1: %.1f | both mu0+mu1: %.1f\n",
              p$mu_0 * window_days, p$mu_1 * window_days, (p$mu_0 + p$mu_1) * window_days))

  # Collect paired counts
  rows <- list()
  for (s in seq_len(n_sims)) {
    seed <- 100000L + s

    # A) legacy marginal: univ(ctrl) vs univ(treat) on full window
    u_c <- sim_univ(marg$ctrl, seed)
    u_t <- sim_univ(marg$treat, seed)

    # B) bivariate full: all-ctrl vs observed mixed
    b_c <- sim_biv(p, ss_all_ctrl, seed)
    b_o <- sim_biv(p, ss_obs, seed)

    # C) bivariate, cross terms zeroed: all-ctrl vs observed
    z_c <- sim_biv(p0, ss_all_ctrl, seed)
    z_o <- sim_biv(p0, ss_obs, seed)

    # D) bivariate: all-ctrl vs all-treated (single-process worlds, but joint law)
    b_at <- sim_biv(p, ss_all_treat, seed)

    # E) bivariate cross=0: all-ctrl vs all-treated
    z_at <- sim_biv(p0, ss_all_treat, seed)

    # F) independent sum approximation for observed: univ ctrl on ctrl ss + univ treat on treat ss
    #    (approximate by univ on full window is wrong; use bivariate cross=0 components via process_id)
    n_zo_0 <- sum(z_o$process_id == 0)
    n_zo_1 <- sum(z_o$process_id == 1)
    n_bo_0 <- sum(b_o$process_id == 0)
    n_bo_1 <- sum(b_o$process_id == 1)
    n_bc_0 <- sum(b_c$process_id == 0)
    n_bc_1 <- sum(b_c$process_id == 1)

    rows[[s]] <- data.frame(
      s = s,
      marg_c = length(u_c$t),
      marg_t = length(u_t$t),
      marg_saved = length(u_c$t) - length(u_t$t),
      biv_allctrl = nrow(b_c),
      biv_obs = nrow(b_o),
      biv_saved_obs = nrow(b_c) - nrow(b_o),
      biv0_allctrl = nrow(z_c),
      biv0_obs = nrow(z_o),
      biv0_saved_obs = nrow(z_c) - nrow(z_o),
      biv_alltreat = nrow(b_at),
      biv_saved_allt = nrow(b_c) - nrow(b_at),
      biv0_alltreat = nrow(z_at),
      biv0_saved_allt = nrow(z_c) - nrow(z_at),
      biv_obs_n0 = n_bo_0,
      biv_obs_n1 = n_bo_1,
      biv_allctrl_n0 = n_bc_0,
      biv_allctrl_n1 = n_bc_1,
      biv0_obs_n0 = n_zo_0,
      biv0_obs_n1 = n_zo_1
    )
  }
  df <- do.call(rbind, rows)
  sm <- function(x) c(mean = mean(x), sd = sd(x))
  report <- function(nm, x) {
    m <- sm(x)
    cat(sprintf("  %-28s mean=%8.1f  sd=%6.1f\n", nm, m["mean"], m["sd"]))
  }

  cat(sprintf("\n--- Monte Carlo means (n_sims=%d) ---\n", n_sims))
  cat("Legacy marginal (univ ctrl vs univ treat on FULL window):\n")
  report("N(univ ctrl)", df$marg_c)
  report("N(univ treat)", df$marg_t)
  report("saved = ctrl - treat", df$marg_saved)

  cat("\nBivariate FULL (all-ctrl vs OBSERVED mixed allocation):\n")
  report("N(all-ctrl)", df$biv_allctrl)
  report("N(observed)", df$biv_obs)
  report("saved = allctrl - obs", df$biv_saved_obs)
  report("  obs component 0", df$biv_obs_n0)
  report("  obs component 1", df$biv_obs_n1)
  report("  allctrl comp 0", df$biv_allctrl_n0)
  report("  allctrl comp 1 (cross)", df$biv_allctrl_n1)

  cat("\nBivariate CROSS=0 (all-ctrl vs OBSERVED):\n")
  report("N(all-ctrl)", df$biv0_allctrl)
  report("N(observed)", df$biv0_obs)
  report("saved = allctrl - obs", df$biv0_saved_obs)
  report("  obs component 0", df$biv0_obs_n0)
  report("  obs component 1", df$biv0_obs_n1)

  cat("\nBivariate FULL (all-ctrl vs ALL-TREATED):\n")
  report("N(all-treated)", df$biv_alltreat)
  report("saved = allctrl - alltreat", df$biv_saved_allt)

  cat("\nBivariate CROSS=0 (all-ctrl vs ALL-TREATED):\n")
  report("N(all-treated)", df$biv0_alltreat)
  report("saved = allctrl - alltreat", df$biv0_saved_allt)

  cat("\n--- Decomposition of the gap ---\n")
  gap_total <- mean(df$biv_saved_obs) - mean(df$marg_saved)
  gap_estimand_no_cross <- mean(df$biv0_saved_obs) - mean(df$marg_saved)
  gap_cross <- mean(df$biv_saved_obs) - mean(df$biv0_saved_obs)
  gap_vs_allt_full <- mean(df$biv_saved_allt) - mean(df$marg_saved)
  gap_vs_allt_0 <- mean(df$biv0_saved_allt) - mean(df$marg_saved)
  cat(sprintf("Total (biv_obs - marg):              %+.1f\n", gap_total))
  cat(sprintf("  of which estimand/support (cross=0 vs marg): %+.1f\n", gap_estimand_no_cross))
  cat(sprintf("  of which cross terms (biv - biv0) on obs:    %+.1f\n", gap_cross))
  cat(sprintf("If instead all-ctrl vs all-treat (full, with cross): gap vs marg = %+.1f\n", gap_vs_allt_full))
  cat(sprintf("If instead all-ctrl vs all-treat (cross=0):          gap vs marg = %+.1f\n", gap_vs_allt_0))
  cat(sprintf("Univ-ctrl vs biv-allctrl: %.1f vs %.1f (diff %+.1f)\n",
              mean(df$marg_c), mean(df$biv_allctrl), mean(df$biv_allctrl) - mean(df$marg_c)))
  cat(sprintf("Univ-treat vs biv-obs:    %.1f vs %.1f (diff %+.1f)\n",
              mean(df$marg_t), mean(df$biv_obs), mean(df$biv_obs) - mean(df$marg_t)))
  cat(sprintf("Univ-treat vs biv-allt:   %.1f vs %.1f (diff %+.1f)\n",
              mean(df$marg_t), mean(df$biv_alltreat), mean(df$biv_alltreat) - mean(df$marg_t)))

  invisible(df)
}

set.seed(1)
df_E <- run_scenarios("E")
df_F <- run_scenarios("F")

out <- file.path(OUT_DIR, "ate_biv_gap_diagnosis.rds")
saveRDS(list(E = df_E, F = df_F, n_sims = n_sims, note = "controlled ATE gap diagnosis"), out)
cat("\nWrote ", out, "\n", sep = "")
