#!/usr/bin/env Rscript
# Recompute observed E/F all-or-nothing ATEs with full bivariate simulation.
# Does NOT rerun SEM fits or bootstrap; preserves bootstrap replicate shape and
# updates observed ATE so HTML Fig (bias-corrected recenter) uses the new center.
#
# Usage (from inst/oklahoma):
#   Rscript recompute_ate_bivariate.R
#   Rscript recompute_ate_bivariate.R --input .../for_paper.rds --contrast observed
#   Rscript recompute_ate_bivariate.R --contrast all_or_nothing \
#       --output .../for_paper_biv_aon.rds
#
# contrast:
#   observed       — control-everywhere vs observed mixed allocation (default)
#   all_or_nothing — control-everywhere vs treated-everywhere (original AoN estimand)
#
# Default output: for_paper_biv_ate.rds (observed) or for_paper_biv_aon.rds (AoN)

args <- commandArgs(trailingOnly = TRUE)

script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- if (length(script_file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[1]), winslash = "/", mustWork = FALSE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
REPO_DIR <- if (basename(SCRIPT_DIR) == "oklahoma" &&
                basename(dirname(SCRIPT_DIR)) == "inst") {
  normalizePath(dirname(dirname(SCRIPT_DIR)), winslash = "/", mustWork = FALSE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("pkgload is required. Install with install.packages('pkgload').")
}
pkgload::load_all(REPO_DIR, quiet = TRUE, export_all = FALSE, helpers = FALSE, attach_testthat = FALSE)

source(file.path(SCRIPT_DIR, "ate_bivariate.R"), local = FALSE)

get_arg_val <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) < 1L || idx[[1]] >= length(args)) return(default)
  args[[idx[[1]] + 1L]]
}

OUT_DIR <- PPDisentangle::pp_output_path("oklahoma", repo_root = REPO_DIR)
input_arg <- get_arg_val("--input", NULL)
input_rds <- if (!is.null(input_arg) && nzchar(input_arg)) {
  normalizePath(input_arg, winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(OUT_DIR, "for_paper.rds"), winslash = "/", mustWork = TRUE)
}
contrast <- tolower(trimws(get_arg_val("--contrast", "observed")))
if (!contrast %in% c("observed", "all_or_nothing")) {
  stop("--contrast must be 'observed' or 'all_or_nothing', got: ", contrast)
}
default_out <- if (identical(contrast, "all_or_nothing")) {
  file.path(OUT_DIR, "for_paper_biv_aon.rds")
} else {
  file.path(OUT_DIR, "for_paper_biv_ate.rds")
}
output_rds <- normalizePath(
  get_arg_val("--output", default_out),
  winslash = "/",
  mustWork = FALSE
)
cat("Contrast: ", contrast, "\n", sep = "")

DATA_DIR <- file.path(SCRIPT_DIR, "oklahoma_induced_seismicity_data_regional20150318")
if (!dir.exists(DATA_DIR)) {
  stop("Data directory not found: ", DATA_DIR)
}

cat("Loading: ", input_rds, "\n", sep = "")
res <- readRDS(input_rds)
cfg <- res$config
if (is.null(cfg)) stop("RDS missing config")

n_sims <- suppressWarnings(as.integer(cfg$ATE_N_SIMS))
if (!is.finite(n_sims) || is.na(n_sims) || n_sims < 1L) n_sims <- 500L
window_days <- suppressWarnings(as.numeric(cfg$ATE_WINDOW_DAYS))
if (!is.finite(window_days) || is.na(window_days) || window_days <= 0) window_days <- 100
m0 <- suppressWarnings(as.numeric(cfg$ETAS_M0))
if (!is.finite(m0)) m0 <- 2.5
beta_gr <- suppressWarnings(as.numeric(cfg$BETA_GR))
if (!is.finite(beta_gr)) beta_gr <- 2.3
t_trunc <- suppressWarnings(as.numeric(cfg$SEM_T_TRUNC_DAYS))
if (!is.finite(t_trunc)) t_trunc <- NULL
use_crn <- isTRUE(cfg$ATE_USE_CRN)
crn_pair <- isTRUE(cfg$ATE_CRN_PAIR)
crn_base <- suppressWarnings(as.integer(cfg$ATE_CRN_BASE))
if (length(crn_base) != 1L || !is.finite(crn_base) || is.na(crn_base)) {
  gs <- suppressWarnings(as.integer(cfg$OK_GLOBAL_SEED))
  crn_base <- if (length(gs) == 1L && is.finite(gs) && !is.na(gs)) {
    as.integer(100000L + 1000L * gs)
  } else {
    100000L
  }
}

cat("Building county supports from data dir...\n")
geom <- oklahoma_build_state_spaces(DATA_DIR, crs_proj = 5070L)
pp_pre <- res$pp_data$pp_pre
if (is.null(pp_pre)) stop("RDS missing pp_data$pp_pre")

# Match analysis convention: condition ATE on pre-treatment history.
filtration <- pp_pre[, c("x", "y", "t", "mag", "inferred_process"), drop = FALSE]

recompute_one <- function(letter) {
  fit <- res$fits_named[[letter]]
  if (is.null(fit) || is.null(fit$params)) {
    stop("Missing fits_named$", letter, "$params")
  }
  old_ate <- fit$ate
  old_mean <- if (!is.null(old_ate$all_nothing_sim$total_saved)) {
    mean(old_ate$all_nothing_sim$total_saved, na.rm = TRUE)
  } else {
    NA_real_
  }
  cat(sprintf("\n=== Fit %s (old marginal mean saved = %.3f) ===\n", letter, old_mean))
  new_ate <- ate_estim_bivariate(
    biv_params = fit$params,
    windowT = c(0, window_days),
    windowS = geom$win_km,
    state_spaces_obs = geom$state_spaces,
    label = paste0(letter, ": bivariate ", contrast),
    n_sims = n_sims,
    m0 = m0,
    beta_gr = beta_gr,
    filtration_history = filtration,
    t_trunc = t_trunc,
    n_tiles = geom$n_tiles,
    crn_base_seed = crn_base,
    use_crn = use_crn,
    crn_pair = crn_pair,
    quiet = FALSE,
    contrast = contrast
  )
  new_mean <- mean(new_ate$all_nothing_sim$total_saved, na.rm = TRUE)
  list(old_ate = old_ate, new_ate = new_ate, old_mean = old_mean, new_mean = new_mean)
}

out_E <- recompute_one("E")
out_F <- recompute_one("F")

# Preserve marginal ATEs for comparison; swap in bivariate as primary.
res$fits_named$E$ate_marginal <- out_E$old_ate
res$fits_named$F$ate_marginal <- out_F$old_ate
res$fits_named$E$ate <- out_E$new_ate
res$fits_named$F$ate <- out_F$new_ate

if (!is.null(res$fitE)) {
  res$fitE$ate_marginal <- out_E$old_ate
  res$fitE$ate <- out_E$new_ate
}
if (!is.null(res$fitF)) {
  res$fitF$ate_marginal <- out_F$old_ate
  res$fitF$ate <- out_F$new_ate
}

contrast_note <- if (identical(contrast, "all_or_nothing")) {
  "control-everywhere vs treated-everywhere (original AoN estimand, bivariate law with cross terms)"
} else {
  "control-everywhere vs observed mixed allocation"
}
res$ate_bivariate_patch <- list(
  created = as.character(Sys.time()),
  source_rds = input_rds,
  method = paste0("bivariate_", contrast),
  contrast = contrast,
  note = paste0(
    "E/F ATEs recomputed with full bivariate forward simulation (",
    contrast_note, "). Bootstrap replicates unchanged; ",
    "HTML bias-corrected figure recenters old bootstrap shape to new observed means."
  ),
  n_sims = n_sims,
  window_days = window_days,
  t_trunc = t_trunc,
  E = list(old_marginal_mean = out_E$old_mean, new_bivariate_mean = out_E$new_mean),
  F = list(old_marginal_mean = out_F$old_mean, new_bivariate_mean = out_F$new_mean)
)

if (!is.null(res$config)) {
  res$config$ATE_METHOD <- paste0("bivariate_", contrast)
  res$config$ATE_BIVARIATE_PATCHED <- TRUE
  res$config$ATE_CONTRAST <- contrast
}

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(res, output_rds)
cat("\nWrote: ", output_rds, "\n", sep = "")
cat(sprintf(
  "Summary:\n  E: marginal %.1f -> bivariate %.1f\n  F: marginal %.1f -> bivariate %.1f\n",
  out_E$old_mean, out_E$new_mean, out_F$old_mean, out_F$new_mean
))
cat("Bootstrap replicates left unchanged.\n")
cat("Render HTML with:\n")
cat("  Rscript render_oklahoma_report.R --results ", output_rds, " --to html\n", sep = "")
