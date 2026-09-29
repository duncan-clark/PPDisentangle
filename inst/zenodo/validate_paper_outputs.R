#!/usr/bin/env Rscript
# Validate the frozen v0.2.0 paper profile and its regenerated summaries.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1L]] else Sys.getenv("PPDISENTANGLE_OUTPUT_ROOT", "../PPDisentangle-output")
root <- normalizePath(root, mustWork = TRUE)
source("inst/oklahoma/paper/paper_result_helpers.R")
x <- readRDS(file.path(root, "oklahoma/for_paper.rds"))
cfg <- x$config
stopifnot(cfg$n_post == 818L, cfg$n_pre == 1236L, cfg$n_pre_holdout == 1236L,
          cfg$SEM_T_TRUNC_DAYS == 250, cfg$SEM_INNER_ITER == 2000L,
          cfg$BOOT_N_REPS == 512L, cfg$BOOT_SEM_INNER_ITER == 500L,
          cfg$ATE_WINDOW_DAYS == 100, cfg$ATE_N_SIMS == 500L)
expected <- rbind(c(961.576, 1013.084, -51.508), c(1105.422, 930.222, 175.200))
for (i in seq_along(c("C", "D"))) {
  fit <- x$fits_named[[c("C", "D")[[i]]]]
  s <- select_paper_ate(fit$ate, "observed")$all_nothing_sim
  actual <- colMeans(s[, c("c_total", "t_total", "total_saved")])
  stopifnot(max(abs(actual - expected[i, ])) < 1e-8)
}
gen <- file.path(root, "oklahoma/paper/generated")
counts <- read.csv(file.path(gen, "expected_counts.csv"))
stopifnot(nrow(counts) == 10L, all(counts$contrast == "observed"),
          max(abs(as.matrix(counts[1:2, c("control", "intervention", "saved")]) - expected)) < 1e-8,
          identical(readRDS(file.path(gen, "build_manifest.rds"))$contrast, "observed"))
boot <- read.csv(file.path(gen, "ef_bootstrap_summary.csv"))
stopifnot(identical(as.integer(boot$n_boot), c(507L, 512L)),
          all(boot$n_attempted == 512L),
          max(abs(boot$mean_saved - expected[, 3])) < 1e-8)
rob <- file.path(root, "sim_study/paper/robustness_merged_tcal")
raw <- setdiff(list.files(rob, "^robustness_.*\\.rds$"), "robustness_merged_tcal_summary.rds")
stopifnot(length(raw) == 67L)
required <- c("oklahoma/paper/generated/figures/ATE_diff.pdf",
              "oklahoma/paper/generated/figures/ok_partition_county.pdf",
              "oklahoma/paper/generated/figures/ok_point_patterns.pdf",
              "oklahoma/paper/generated/tab_ok_partition_ate.tex",
              "sim_study/generated/figures/simulated_hawkes_hawkes_process.pdf")
stopifnot(all(file.exists(file.path(root, required))))
message("Validated frozen paper profile: observed-vs-none counts; 507/512 and 512/512 bootstrap draws; 67 robustness scenarios.")
