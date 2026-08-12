#!/usr/bin/env Rscript
# Smoke: homogeneous SEM Fit B with t_trunc = 90, allowing truncated Omori p < 1.
# Writes under PPDisentangle-output/smoke_fit_b_p_lt1/ — does not touch
# publication oklahoma_results.rds.

args_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(args_file) > 0) {
  dirname(normalizePath(sub("^--file=", "", args_file[[1]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
pkg_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
smoke_root <- normalizePath(
  file.path(dirname(pkg_root), "PPDisentangle-output", "smoke_fit_b_p_lt1"),
  winslash = "/", mustWork = FALSE
)
dir.create(smoke_root, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(
  PPDISENTANGLE_OUTPUT_ROOT = smoke_root,
  OK_AB_ONLY = "true",
  OK_SKIP_ATE = "true",
  OK_RUN_SENSITIVITY = "false",
  OK_RUN_T_TRUNC_SENSITIVITY = "false",
  OK_RUN_BOOTSTRAP_ATE = "false",
  OK_RUN_FIT_VARIABILITY = "false",
  OK_RUN_SEM_PILOT = "false",
  OK_SEM_T_TRUNC_DAYS = "90",
  OK_SEM_N_ITER = "1",
  OK_SEM_N_LABELLINGS = "5",
  OK_SEM_INNER_ITER = "200",
  OK_SEM_OUTER_MAXIT_BIV = "800",
  OK_VANILLA_MAXIT = "400",
  OK_KDE_VARIANT_MODE = "single",
  OK_REPORT_FORMATS = "",
  OK_SEM_WORKER_LOGS = "false",
  OK_CORES = "1"
)

cat("=== Fit B smoke: t_trunc=90, truncated Omori p>0 allowed ===\n")
cat("Package root:", pkg_root, "\n")
cat("Output root:", smoke_root, "\n")
analysis <- file.path(script_dir, "oklahoma_analysis.R")
if (!file.exists(analysis)) stop("Cannot find ", analysis)
source(analysis, local = FALSE)
