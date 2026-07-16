#!/usr/bin/env Rscript
# Robustness-suite launcher for PPDisentangle simulation studies.
#
# This script runs the existing Hawkes simulation study across a grid of
# K-separation and signal-to-noise scenarios. Each scenario writes the normal
# sim_study.R result object with extra scenario metadata in config.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(PPDisentangle)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg_val <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) < 1L || idx[[1]] >= length(args)) return(default)
  args[[idx[[1]] + 1L]]
}

parse_num_vec <- function(raw, default) {
  if (is.null(raw) || !nzchar(raw)) return(default)
  toks <- unlist(strsplit(raw, "[,;[:space:]]+", perl = TRUE), use.names = FALSE)
  vals <- suppressWarnings(as.numeric(toks[nzchar(toks)]))
  vals <- vals[is.finite(vals)]
  if (length(vals) < 1L) return(default)
  unique(vals)
}

format_num_tag <- function(x) {
  gsub("[^0-9A-Za-z]+", "p", sprintf("%.3f", as.numeric(x)))
}

# Control branching ratio (high); treatment lowers K relative to control.
ROBUSTNESS_CONTROL_K <- 0.8
ROBUSTNESS_TREATED_K_DEFAULT <- 0.2

# SNR grid: background-rate multiplier with mu recalibrated toward target_points.
DEFAULT_MU_SCALES <- c(0.25, 0.5, 1, 1.5, 2)
DEFAULT_MU_SCALES_STR <- paste(DEFAULT_MU_SCALES, collapse = ",")
format_mu_scale_grid <- function(scales = DEFAULT_MU_SCALES, for_tex = FALSE) {
  vals <- vapply(scales, function(x) {
    s <- format(x, trim = TRUE, scientific = FALSE)
    sub("\\.?0+$", "", s)
  }, character(1L))
  body <- paste(vals, collapse = ", ")
  if (for_tex) paste0("\\{", body, "\\}") else paste0("{", body, "}")
}

# Treated-K grid for k_separation with control fixed at k_anchor (default 0.8).
default_k_separation_values <- function(k_anchor = ROBUSTNESS_CONTROL_K) {
  vals <- seq(0.1, 0.7, by = 0.1)
  vals <- vals[vals >= 0 & vals < 1 & abs(vals - k_anchor) > 1e-10]
  sort(unique(vals))
}

script_dir <- {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)
  if (length(file_arg) > 0L) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }
}
repo_dir <- if (basename(script_dir) == "sim_study" && basename(dirname(script_dir)) == "inst") {
  normalizePath(dirname(dirname(script_dir)), winslash = "/", mustWork = FALSE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
sim_root <- PPDisentangle::pp_output_path("sim_study", repo_root = repo_dir)
out_dir <- sim_root
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pp_sims <- suppressWarnings(as.integer(get_arg_val("--sims", Sys.getenv("PP_SIMS", "32"))))
if (!is.finite(pp_sims) || is.na(pp_sims) || pp_sims < 1L) pp_sims <- 32L
allocated_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", Sys.getenv("PP_CPUS", as.character(pp_sims)))))
if (!is.finite(allocated_cpus) || is.na(allocated_cpus) || allocated_cpus < 1L) allocated_cpus <- pp_sims
scenario_workers <- suppressWarnings(as.integer(get_arg_val("--scenario-workers", Sys.getenv("PP_ROBUSTNESS_SCENARIO_WORKERS", ""))))
if (!is.finite(scenario_workers) || is.na(scenario_workers) || scenario_workers < 1L) {
  scenario_workers <- max(1L, floor(allocated_cpus / pp_sims))
}
target_points <- suppressWarnings(as.numeric(get_arg_val("--target-points", Sys.getenv("PP_TARGET_POINTS", "2500"))))
if (!is.finite(target_points) || is.na(target_points) || target_points <= 0) target_points <- 2500
test_mode <- "--test" %in% args
decay_reps <- suppressWarnings(as.integer(get_arg_val("--decay-reps", Sys.getenv("PP_DECAY_REPS", "2000"))))
if (!is.finite(decay_reps) || is.na(decay_reps) || decay_reps < 1L) decay_reps <- 2000L
refresh_decay <- "--refresh-decay" %in% args ||
  tolower(Sys.getenv("PP_REFRESH_DECAY", "false")) %in% c("1", "true", "yes", "y")

k_anchor <- suppressWarnings(as.numeric(get_arg_val("--k-anchor", Sys.getenv("PP_K_ANCHOR", as.character(ROBUSTNESS_CONTROL_K)))))
if (!is.finite(k_anchor) || is.na(k_anchor) || k_anchor < 0 || k_anchor >= 1) k_anchor <- ROBUSTNESS_CONTROL_K
k_values_default <- default_k_separation_values(k_anchor)
k_values_raw <- get_arg_val("--k-values", Sys.getenv("PP_K_VALUES", ""))
k_values <- if (is.null(k_values_raw) || !nzchar(k_values_raw)) {
  k_values_default
} else {
  parse_num_vec(k_values_raw, default = k_values_default)
}
k_values <- k_values[k_values >= 0 & k_values < 1 & abs(k_values - k_anchor) > 1e-10]
if (length(k_values) < 1L) {
  stop(sprintf(
    "No valid treated K values after filtering anchor=%.3f. Provide --k-values or PP_K_VALUES.",
    k_anchor
  ))
}
message(sprintf(
  "[robustness] k_separation anchor=%.3f treated K grid (%d): %s",
  k_anchor, length(k_values), paste(signif(k_values, 4), collapse = ", ")
))
mu_scales <- parse_num_vec(get_arg_val("--mu-scales", Sys.getenv("PP_MU_SCALES", DEFAULT_MU_SCALES_STR)),
                           default = DEFAULT_MU_SCALES)
mu_scales <- mu_scales[mu_scales > 0]
message(sprintf(
  "[robustness] snr mu_scale grid (%d): %s",
  length(mu_scales), paste(signif(mu_scales, 4), collapse = ", ")
))

scenario_filter <- get_arg_val("--scenario-set", Sys.getenv("PP_ROBUSTNESS_SCENARIO_SET", "all"))
if (!nzchar(scenario_filter)) scenario_filter <- "all"
replot_basename <- get_arg_val("--replot", Sys.getenv("PP_ROBUSTNESS_REPLOT", ""))
if (!nzchar(replot_basename)) replot_basename <- NULL
resume_from <- get_arg_val("--resume-from", Sys.getenv("PP_ROBUSTNESS_RESUME_FROM", ""))
if (!nzchar(resume_from)) resume_from <- NULL
run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_basename <- if (!is.null(replot_basename)) {
  replot_basename
} else {
  paste0("robustness_", run_stamp)
}
if (is.null(replot_basename) && nzchar(Sys.getenv("SLURM_JOB_ID", ""))) {
  output_basename <- paste0("robustness_", Sys.getenv("SLURM_JOB_ID"))
}

# Archived paper runs live under sim_study/paper/<basename>/; new cluster runs
# still write flat under sim_study/. Prefer the paper folder when replotting.
# Figure/tex outputs always go under sim_study/generated/robustness/ (not under paper/).
resolve_replot_out_dir <- function(base_out_dir, basename) {
  candidates <- c(
    file.path(base_out_dir, "paper", basename),
    file.path(base_out_dir, basename),
    base_out_dir
  )
  for (cand in candidates) {
    manifest_hit <- file.path(cand, paste0(basename, "_manifest.csv"))
    if (file.exists(manifest_hit)) return(normalizePath(cand, winslash = "/", mustWork = FALSE))
  }
  base_out_dir
}
if (!is.null(replot_basename)) {
  out_dir <- resolve_replot_out_dir(sim_root, replot_basename)
  message("[robustness] replot data dir: ", out_dir)
}

k_separation <- data.frame(
  scenario_family = "k_separation",
  control_k = k_anchor,
  treated_k = k_values,
  mu_scale = 1,
  sim_kernel = "exponential",
  fit_kernel = "exponential",
  sim_spatial_kernel = "exponential",
  fit_spatial_kernel = "exponential",
  treatment_assignment = "random",
  stringsAsFactors = FALSE
)

pretreatment_assignment <- expand.grid(
  scenario_family = "pretreatment_assignment",
  control_k = ROBUSTNESS_CONTROL_K,
  treated_k = ROBUSTNESS_TREATED_K_DEFAULT,
  mu_scale = 1,
  sim_kernel = "exponential",
  fit_kernel = "exponential",
  sim_spatial_kernel = "exponential",
  fit_spatial_kernel = "exponential",
  treatment_assignment = c(
    "highest_count_50pct",
    "lowest_count_50pct",
    "count_propensity_50pct",
    "contiguous_aoi_50pct",
    "voronoi_random_50pct"
  ),
  stringsAsFactors = FALSE
)

snr_scale <- data.frame(
  scenario_family = "snr_scale",
  control_k = ROBUSTNESS_CONTROL_K,
  treated_k = ROBUSTNESS_TREATED_K_DEFAULT,
  mu_scale = mu_scales,
  sim_kernel = "exponential",
  fit_kernel = "exponential",
  sim_spatial_kernel = "exponential",
  fit_spatial_kernel = "exponential",
  treatment_assignment = "random",
  stringsAsFactors = FALSE
)

spatiotemporal_kernel_mismatch <- expand.grid(
  scenario_family = "spatiotemporal_kernel_mismatch",
  control_k = ROBUSTNESS_CONTROL_K,
  treated_k = ROBUSTNESS_TREATED_K_DEFAULT,
  mu_scale = 1,
  sim_kernel = c("exponential", "power_law"),
  fit_kernel = c("exponential", "power_law"),
  sim_spatial_kernel = c("exponential", "power_law"),
  fit_spatial_kernel = c("exponential", "power_law"),
  treatment_assignment = "random",
  stringsAsFactors = FALSE
)

scenarios <- bind_rows(
  k_separation,
  pretreatment_assignment,
  snr_scale,
  spatiotemporal_kernel_mismatch
) %>%
  mutate(
    k_delta = abs(.data$treated_k - .data$control_k),
    scenario_id = paste0(
      .data$scenario_family,
      "_kc", format_num_tag(.data$control_k),
      "_kt", format_num_tag(.data$treated_k),
      "_mu", format_num_tag(.data$mu_scale),
      "_sim", .data$sim_kernel,
      "_fit", .data$fit_kernel,
      "_ssim", .data$sim_spatial_kernel,
      "_sfit", .data$fit_spatial_kernel,
      "_assign", .data$treatment_assignment
    )
  ) %>%
  distinct(.data$scenario_id, .keep_all = TRUE)

if (!identical(tolower(scenario_filter), "all")) {
  keep_families <- unlist(strsplit(scenario_filter, "[,;[:space:]]+", perl = TRUE), use.names = FALSE)
  keep_families <- keep_families[nzchar(keep_families)]
  # Alias legacy family names used in older jobs / docs.
  if ("high_count_assignment" %in% keep_families) {
    keep_families <- unique(c(keep_families, "pretreatment_assignment"))
  }
  if (any(c("kernel_mismatch", "spatial_kernel_mismatch") %in% keep_families)) {
    keep_families <- unique(c(keep_families, "spatiotemporal_kernel_mismatch"))
  }
  scenarios <- scenarios %>% filter(.data$scenario_family %in% keep_families)
}
if (nrow(scenarios) < 1L) stop("No robustness scenarios selected.")
scenario_workers <- max(1L, min(as.integer(scenario_workers), nrow(scenarios)))
scenario_cpus <- max(1L, floor(allocated_cpus / scenario_workers))
message(sprintf(
  "[robustness] scenario workers=%d | allocated_cpus=%d | cpus_per_scenario=%d | sims_per_scenario=%d",
  scenario_workers, allocated_cpus, scenario_cpus, pp_sims
))

robustness_dir <- file.path(sim_root, "generated", "robustness")
fig_dir <- file.path(robustness_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(robustness_dir, recursive = TRUE, showWarnings = FALSE)
message("[robustness] generated dir: ", robustness_dir)
source(file.path(script_dir, "decay_validation_utils.R"), local = FALSE)

run_one <- function(row_id) {
  sc <- scenarios[row_id, , drop = FALSE]
  run_basename <- paste0(output_basename, "_", sc$scenario_id)
  rds_path <- file.path(out_dir, paste0(run_basename, ".rds"))

  # Resume: reuse a completed RDS from a prior job with the same scenario_id.
  if (!is.null(resume_from)) {
    resume_rds <- file.path(out_dir, paste0(resume_from, "_", sc$scenario_id, ".rds"))
    if (!file.exists(resume_rds) && identical(out_dir, sim_root)) {
      # Also look in paper/ archive layout.
      resume_rds_alt <- file.path(sim_root, "paper", resume_from, paste0(resume_from, "_", sc$scenario_id, ".rds"))
      if (file.exists(resume_rds_alt)) resume_rds <- resume_rds_alt
    }
    if (file.exists(resume_rds)) {
      message(sprintf(
        "[robustness] %d/%d %s | reusing %s",
        row_id, nrow(scenarios), sc$scenario_id, basename(resume_rds)
      ))
      rds_path <- resume_rds
      run_basename <- sub("\\.rds$", "", basename(resume_rds))
      return(data.frame(
        scenario_id = sc$scenario_id,
        scenario_family = sc$scenario_family,
        control_k = sc$control_k,
        treated_k = sc$treated_k,
        k_delta = sc$k_delta,
        mu_scale = sc$mu_scale,
        sim_kernel = sc$sim_kernel,
        fit_kernel = sc$fit_kernel,
        sim_spatial_kernel = sc$sim_spatial_kernel,
        fit_spatial_kernel = sc$fit_spatial_kernel,
        treatment_assignment = sc$treatment_assignment,
        target_points = target_points,
        run_basename = run_basename,
        rds_path = rds_path,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Skip recompute if this job already wrote the RDS (safe restart).
  if (file.exists(rds_path)) {
    message(sprintf(
      "[robustness] %d/%d %s | existing RDS found, skipping recompute",
      row_id, nrow(scenarios), sc$scenario_id
    ))
    return(data.frame(
      scenario_id = sc$scenario_id,
      scenario_family = sc$scenario_family,
      control_k = sc$control_k,
      treated_k = sc$treated_k,
      k_delta = sc$k_delta,
      mu_scale = sc$mu_scale,
      sim_kernel = sc$sim_kernel,
      fit_kernel = sc$fit_kernel,
      sim_spatial_kernel = sc$sim_spatial_kernel,
      fit_spatial_kernel = sc$fit_spatial_kernel,
      treatment_assignment = sc$treatment_assignment,
      target_points = target_points,
      run_basename = run_basename,
      rds_path = rds_path,
      stringsAsFactors = FALSE
    ))
  }

  message(sprintf(
    "[robustness] %d/%d %s | K=(%.3f, %.3f) mu_scale=%.3f target_points=%.0f",
    row_id, nrow(scenarios), sc$scenario_id, sc$control_k, sc$treated_k,
    sc$mu_scale, target_points
  ))
  env <- c(
    PP_OUTPUT_BASENAME = run_basename,
    PP_OUTPUT_TAG = "",
    PP_SCENARIO_ID = sc$scenario_id,
    PP_CPUS = as.character(scenario_cpus),
    SLURM_CPUS_PER_TASK = as.character(scenario_cpus),
    PP_CONTROL_K = as.character(sc$control_k),
    PP_TREATED_K = as.character(sc$treated_k),
    PP_MU_SCALE = as.character(sc$mu_scale),
    PP_TARGET_POINTS = as.character(target_points),
    PP_SIM_KERNEL = sc$sim_kernel,
    PP_FIT_KERNEL = sc$fit_kernel,
    PP_SIM_SPATIAL_KERNEL = sc$sim_spatial_kernel,
    PP_FIT_SPATIAL_KERNEL = sc$fit_spatial_kernel,
    PP_TREATMENT_ASSIGNMENT = sc$treatment_assignment,
    PP_DECAY_REPS = as.character(decay_reps)
  )
  cmd_args <- c(file.path("inst", "sim_study", "sim_study.R"), "--sims", as.character(pp_sims))
  if (test_mode) cmd_args <- c(cmd_args, "--test")
  old_wd <- getwd()
  setwd(repo_dir)
  on.exit(setwd(old_wd), add = TRUE)
  status <- system2("Rscript", args = cmd_args, env = sprintf("%s=%s", names(env), unname(env)))
  if (!identical(status, 0L)) stop(sprintf("sim_study.R failed for %s", sc$scenario_id))
  if (!file.exists(rds_path)) stop(sprintf("Expected output missing: %s", rds_path))
  data.frame(
    scenario_id = sc$scenario_id,
    scenario_family = sc$scenario_family,
    control_k = sc$control_k,
    treated_k = sc$treated_k,
    k_delta = sc$k_delta,
    mu_scale = sc$mu_scale,
    sim_kernel = sc$sim_kernel,
    fit_kernel = sc$fit_kernel,
    sim_spatial_kernel = sc$sim_spatial_kernel,
    fit_spatial_kernel = sc$fit_spatial_kernel,
    treatment_assignment = sc$treatment_assignment,
    target_points = target_points,
    run_basename = run_basename,
    rds_path = rds_path,
    stringsAsFactors = FALSE
  )
}

run_one_safe <- function(row_id) {
  tryCatch(
    run_one(row_id),
    error = function(e) {
      sc <- scenarios[row_id, , drop = FALSE]
      warning(sprintf(
        "[robustness] scenario %s failed: %s",
        sc$scenario_id, conditionMessage(e)
      ), call. = FALSE)
      NULL
    }
  )
}

is_manifest_row <- function(x) {
  is.data.frame(x) && nrow(x) >= 1L && "scenario_id" %in% names(x)
}

manifest <- if (!is.null(replot_basename)) {
  manifest_path <- file.path(out_dir, paste0(output_basename, "_manifest.csv"))
  if (!file.exists(manifest_path)) {
    stop(sprintf("Replot manifest not found: %s", manifest_path))
  }
  message("[robustness] replot from existing summaries: ", output_basename)
  read.csv(manifest_path, stringsAsFactors = FALSE)
} else {
  if (!is.null(resume_from)) {
    message("[robustness] resume-from prior basename: ", resume_from)
  }
  scenario_ids <- seq_len(nrow(scenarios))
  raw_results <- if (scenario_workers > 1L) {
    parallel::mclapply(
      scenario_ids,
      run_one_safe,
      mc.cores = scenario_workers,
      mc.preschedule = FALSE
    )
  } else {
    lapply(scenario_ids, run_one_safe)
  }
  # Drop NULLs and any try-error leftovers from parallel failures.
  kept <- Filter(is_manifest_row, raw_results)
  n_fail <- length(scenario_ids) - length(kept)
  if (n_fail > 0L) {
    warning(sprintf(
      "[robustness] %d/%d scenarios failed; continuing with successful rows",
      n_fail, length(scenario_ids)
    ), call. = FALSE)
  }
  if (!length(kept)) {
    stop("[robustness] all scenarios failed; no manifest rows to aggregate")
  }
  bind_rows(kept)
}

resolve_manifest_rds_path <- function(row, out_dir) {
  rds_path <- row$rds_path[[1]]
  if (file.exists(rds_path)) return(rds_path)
  if ("run_basename" %in% names(row) && nzchar(row$run_basename[[1]])) {
    local_path <- file.path(out_dir, paste0(row$run_basename[[1]], ".rds"))
    if (file.exists(local_path)) return(local_path)
  }
  if ("scenario_id" %in% names(row) && nzchar(row$scenario_id[[1]])) {
    alt <- list.files(
      out_dir,
      pattern = paste0(row$scenario_id[[1]], "\\.rds$"),
      full.names = TRUE,
      recursive = TRUE
    )
    if (length(alt) > 0L) return(alt[[1L]])
  }
  rds_path
}

load_label_summary_from_manifest <- function(manifest_df, out_dir) {
  rows <- lapply(seq_len(nrow(manifest_df)), function(i) {
    row <- manifest_df[i, , drop = FALSE]
    rds_path <- resolve_manifest_rds_path(row, out_dir)
    if (!file.exists(rds_path)) {
      warning(sprintf("Missing RDS for label recovery replot: %s", basename(rds_path)))
      return(NULL)
    }
    res <- readRDS(rds_path)
    cls <- rebuild_class_metrics_from_rds(res)
    if (is.null(cls) || nrow(cls) < 1L) return(NULL)
    transform(cls, scenario_id = row$scenario_id[[1]])
  })
  bind_rows(Filter(Negate(is.null), rows)) %>%
    left_join(manifest_df, by = "scenario_id")
}

# Label-recovery plots use a looser high-mu guard than ATE reporting.
# Spatial power-law fits often inflate control mu for oracle/naive above the
# default 2x threshold even when labels themselves are well-defined.
LABEL_MU_FAIL_MULTIPLIER <- {
  raw <- suppressWarnings(as.numeric(Sys.getenv("PP_LABEL_MU_FAIL_MULTIPLIER", "5")))
  if (!is.finite(raw) || raw < 1) 5 else raw
}

rebuild_class_metrics_from_rds <- function(res, mu_fail_multiplier = LABEL_MU_FAIL_MULTIPLIER) {
  lr <- res$label_recovery_prefilter
  if (is.null(lr) || nrow(lr) < 1L) {
    # Fall back to stored post-filter metrics only if prefilter is unavailable.
    if (!is.null(res$class_metrics) && nrow(res$class_metrics) > 0L) {
      return(res$class_metrics)
    }
    return(NULL)
  }

  fs <- res$fit_status
  true_mu <- NA_real_
  cfg <- res$config
  if (!is.null(cfg)) {
    true_mu <- suppressWarnings(as.numeric(cfg$TRUE_MU %||% cfg$true_mu %||% NA_real_))
    if (!is.finite(true_mu) && !is.null(cfg$hawkes_par_1$mu)) {
      true_mu <- suppressWarnings(as.numeric(cfg$hawkes_par_1$mu))
    }
  }
  mu_mult <- suppressWarnings(as.numeric(mu_fail_multiplier))
  if (!is.finite(mu_mult) || mu_mult < 1) mu_mult <- 5

  lr_kept <- lapply(unique(lr$method), function(method_nm) {
    fit_label <- if (identical(method_nm, "SEM_adaptive")) "SEM_full" else method_nm
    method_fs <- if (!is.null(fs) && nrow(fs) > 0L) {
      fs %>% filter(.data$labelling == fit_label)
    } else {
      NULL
    }
    if (is.null(method_fs) || nrow(method_fs) < 1L) {
      return(lr %>% filter(.data$method == method_nm))
    }

    keep_mask <- rep(TRUE, nrow(method_fs))
    if ("crazy_fit" %in% names(method_fs)) {
      crazy <- method_fs$crazy_fit
      crazy[is.na(crazy)] <- FALSE
      keep_mask <- keep_mask & !crazy
    }
    if (is.finite(true_mu) && all(c("control_mu", "treated_mu") %in% names(method_fs))) {
      thr <- mu_mult * true_mu
      mu_ok <- is.finite(method_fs$control_mu) & is.finite(method_fs$treated_mu) &
        method_fs$control_mu <= thr & method_fs$treated_mu <= thr
      keep_mask <- keep_mask & mu_ok
    } else if ("include_result" %in% names(method_fs)) {
      incl <- method_fs$include_result
      incl[is.na(incl)] <- FALSE
      keep_mask <- keep_mask & incl
    }
    kept_sims <- unique(method_fs$sim_id[keep_mask])
    if (length(kept_sims) < 1L) return(NULL)
    lr %>% filter(.data$method == method_nm, .data$sim_id %in% kept_sims)
  })
  lr_all <- bind_rows(Filter(Negate(is.null), lr_kept))
  if (nrow(lr_all) < 1L) return(NULL)
  lr_all %>%
    mutate(method = normalize_labelling_method(.data$method)) %>%
    group_by(.data$method) %>%
    summarize(
      n_kept = n(),
      mean_accuracy = mean(.data$accuracy, na.rm = TRUE),
      mean_balanced_accuracy = mean(.data$balanced_accuracy, na.rm = TRUE),
      mean_precision_treated = mean(.data$precision_treated, na.rm = TRUE),
      mean_recall_treated = mean(.data$recall_treated, na.rm = TRUE),
      mean_specificity_control = mean(.data$specificity_control, na.rm = TRUE),
      mean_f1_treated = mean(.data$f1_treated, na.rm = TRUE),
      .groups = "drop"
    )
}

rebuild_support_contrast_summary_from_rds <- function(res, mu_fail_multiplier = LABEL_MU_FAIL_MULTIPLIER) {
  cfg <- res$config
  if (is.null(cfg) || is.null(cfg$hawkes_par_1) || is.null(cfg$hawkes_par_2)) return(NULL)
  ctrl_true <- cfg$hawkes_par_1
  treat_true <- cfg$hawkes_par_2
  r_c_true <- as.numeric(ctrl_true$mu) / max(1e-8, 1 - as.numeric(ctrl_true$K))
  r_t_true <- as.numeric(treat_true$mu) / max(1e-8, 1 - as.numeric(treat_true$K))
  denom <- r_c_true - r_t_true
  if (!is.finite(denom) || abs(denom) < 1e-8) return(NULL)

  global_ref <- if (!is.null(res$all_nothing_ATE)) res$all_nothing_ATE else NA_real_
  specs <- res$support_contrast_specs
  rf <- res$results_flat
  tasks <- res$tasks
  fs <- res$fit_status
  if (is.null(rf) || length(rf) < 1L || is.null(tasks) || length(tasks) < 1L) return(NULL)
  if (is.null(specs) || nrow(specs) < 1L) return(NULL)

  # Recompute regime contrasts as (k treated cells) - (all control) so replots
  # pick up the new definitions even when RDS specs still store older
  # flip-vs-observed contrasts.
  n_cells <- sum(specs$contrast_family == "single_cell_flip")
  if (!is.finite(n_cells) || n_cells < 1L) {
    nx <- suppressWarnings(as.integer(cfg$NX %||% NA_integer_))
    ny <- suppressWarnings(as.integer(cfg$NY %||% NA_integer_))
    n_cells <- if (is.finite(nx) && is.finite(ny) && nx > 0L && ny > 0L) nx * ny else 100L
  }
  treat_prop <- suppressWarnings(as.numeric(cfg$TREAT_PROP %||% 0.5))
  if (!is.finite(treat_prop) || treat_prop <= 0 || treat_prop >= 1) treat_prop <- 0.5
  n_treated_50 <- as.integer(round(n_cells * treat_prop))
  n_treated_50 <- max(1L, min(n_cells - 1L, n_treated_50))
  psi_global <- specs$psi_truth[specs$contrast_family == "global_1_0"]
  psi_global <- if (length(psi_global) < 1L) NA_real_ else psi_global[[1L]]
  if (is.finite(psi_global)) {
    psi_single <- (1 / n_cells) * psi_global
    psi_random_50 <- (n_treated_50 / n_cells) * psi_global
    hit_single <- which(specs$contrast_family == "single_cell_flip")
    if (length(hit_single) >= 1L) {
      specs$psi_truth[hit_single] <- psi_single
      if ("contrast_id" %in% names(specs)) {
        specs$contrast_id[hit_single] <- paste0(
          "single_cell_vs_all_control_",
          if ("cell_id" %in% names(specs)) specs$cell_id[hit_single] else seq_along(hit_single)
        )
      }
      if ("hamming_distance_from_zobs" %in% names(specs)) {
        specs$hamming_distance_from_zobs[hit_single] <- NA_integer_
      }
    }
    hit_50 <- which(specs$contrast_family == "random_50pct_flip")
    if (length(hit_50) >= 1L) {
      specs$psi_truth[hit_50] <- psi_random_50
      if ("contrast_id" %in% names(specs)) {
        specs$contrast_id[hit_50] <- "random_50pct_vs_all_control"
      }
      if ("hamming_distance_from_zobs" %in% names(specs)) {
        specs$hamming_distance_from_zobs[hit_50] <- NA_integer_
      }
    }
  }

  true_mu <- suppressWarnings(as.numeric(cfg$TRUE_MU %||% cfg$true_mu %||% ctrl_true$mu))
  mu_mult <- suppressWarnings(as.numeric(mu_fail_multiplier))
  if (!is.finite(mu_mult) || mu_mult < 1) mu_mult <- 5
  thr <- if (is.finite(true_mu)) mu_mult * true_mu else Inf

  task_ok <- vapply(seq_along(tasks), function(k) {
    if (!is.null(fs) && nrow(fs) > 0L && "task_idx" %in% names(fs)) {
      row <- fs[fs$task_idx == k, , drop = FALSE]
      if (nrow(row) >= 1L) {
        if ("crazy_fit" %in% names(row) && isTRUE(row$crazy_fit[[1]])) return(FALSE)
        ctrl_mu <- suppressWarnings(as.numeric(row$control_mu[[1]]))
        treat_mu <- suppressWarnings(as.numeric(row$treated_mu[[1]]))
        if (!is.finite(ctrl_mu) || !is.finite(treat_mu)) return(FALSE)
        return(ctrl_mu <= thr && treat_mu <= thr)
      }
    }
    r <- rf[[k]]
    if (is.null(r) || is.null(r$treated_pp) || is.null(r$control_pp)) return(FALSE)
    ctrl_mu <- suppressWarnings(as.numeric(r$control_pp$mu))
    treat_mu <- suppressWarnings(as.numeric(r$treated_pp$mu))
    if (!is.finite(ctrl_mu) || !is.finite(treat_mu)) return(FALSE)
    ctrl_mu <= thr && treat_mu <= thr
  }, logical(1))

  kept_idx <- which(task_ok)
  if (length(kept_idx) < 1L) return(NULL)

  # Productivity-rescaled plug-in contrasts under true control + fitted treated.
  # Skip SEM_adaptive (SEM_full is the reported SEM labelling).
  row_chunks <- lapply(kept_idx, function(k) {
    lab <- tasks[[k]]$labelling_name
    if (identical(lab, "SEM_adaptive")) return(NULL)
    r <- rf[[k]]
    tp <- r$treated_pp
    if (is.null(tp)) return(NULL)
    r_t_fit <- as.numeric(tp$mu) / max(1e-8, 1 - as.numeric(tp$K))
    if (!is.finite(r_t_fit)) return(NULL)
    scale <- (r_c_true - r_t_fit) / denom
    data.frame(
      labelling = lab,
      contrast_family = specs$contrast_family,
      psi_truth = specs$psi_truth,
      psi_estimate = specs$psi_truth * scale,
      stringsAsFactors = FALSE
    )
  })

  df <- bind_rows(Filter(Negate(is.null), row_chunks))
  if (is.null(df) || nrow(df) < 1L) return(NULL)

  df %>%
    mutate(
      labelling = normalize_labelling_method(.data$labelling),
      bias = .data$psi_estimate - .data$psi_truth,
      abs_error = abs(.data$bias),
      pct_error = ifelse(
        is.finite(.data$psi_truth) & abs(.data$psi_truth) > 1e-8,
        100 * .data$abs_error / abs(.data$psi_truth),
        ifelse(
          .data$contrast_family == "single_cell_flip" &
            is.finite(global_ref) & abs(global_ref) > 1e-8,
          100 * .data$abs_error / abs(global_ref),
          NA_real_
        )
      ),
      pct_bias = ifelse(
        is.finite(.data$psi_truth) & abs(.data$psi_truth) > 1e-8,
        100 * .data$bias / abs(.data$psi_truth),
        ifelse(
          .data$contrast_family == "single_cell_flip" &
            is.finite(global_ref) & abs(global_ref) > 1e-8,
          100 * .data$bias / abs(global_ref),
          NA_real_
        )
      )
    ) %>%
    group_by(.data$labelling, .data$contrast_family) %>%
    summarize(
      n = n(),
      mean_psi_estimate = mean(.data$psi_estimate, na.rm = TRUE),
      mean_psi_truth = mean(.data$psi_truth, na.rm = TRUE),
      mean_bias = mean(.data$bias, na.rm = TRUE),
      se_bias = ifelse(
        sum(is.finite(.data$bias)) > 1L,
        stats::sd(.data$bias, na.rm = TRUE) / sqrt(sum(is.finite(.data$bias))),
        NA_real_
      ),
      mean_pct_bias = mean(.data$pct_bias, na.rm = TRUE),
      se_pct_bias = ifelse(
        sum(is.finite(.data$pct_bias)) > 1L,
        stats::sd(.data$pct_bias, na.rm = TRUE) / sqrt(sum(is.finite(.data$pct_bias))),
        NA_real_
      ),
      mean_abs_error = mean(.data$abs_error, na.rm = TRUE),
      mean_pct_error = mean(.data$pct_error, na.rm = TRUE),
      .groups = "drop"
    )
}

load_support_summary_from_manifest <- function(manifest_df, out_dir) {
  rows <- lapply(seq_len(nrow(manifest_df)), function(i) {
    row <- manifest_df[i, , drop = FALSE]
    rds_path <- resolve_manifest_rds_path(row, out_dir)
    if (!file.exists(rds_path)) {
      warning(sprintf("Missing RDS for support contrast replot: %s", basename(rds_path)))
      return(NULL)
    }
    res <- readRDS(rds_path)
    rebuilt <- rebuild_support_contrast_summary_from_rds(res)
    if (is.null(rebuilt) || nrow(rebuilt) < 1L) return(NULL)
    rebuilt$scenario_id <- row$scenario_id[[1]]
    rebuilt
  })
  bind_rows(Filter(Negate(is.null), rows)) %>%
    left_join(manifest_df, by = "scenario_id")
}

refresh_decay_summaries <- function(manifest_df, reps, out_dir) {
  message(sprintf("[robustness] refreshing decay validation with %d reps per scenario", reps))
  rows <- lapply(seq_len(nrow(manifest_df)), function(i) {
    row <- manifest_df[i, , drop = FALSE]
    rds_path <- resolve_manifest_rds_path(row, out_dir)
    message(sprintf("[robustness] decay refresh %d/%d %s", i, nrow(manifest_df), row$scenario_id))
    out <- refresh_decay_from_rds(rds_path, reps)
    if (is.null(out)) return(NULL)
    summary_df <- as.data.frame(out$summary)
    summary_df$scenario_id <- row$scenario_id[[1]]
    rate_df <- if (!is.null(out$rate_lines) && nrow(out$rate_lines) > 0) {
      as.data.frame(out$rate_lines)
    } else {
      NULL
    }
    if (!is.null(rate_df)) {
      rate_df$scenario_id <- row$scenario_id[[1]]
    }
    list(summary = summary_df, rate_lines = rate_df)
  })
  summary_bind <- bind_rows(lapply(rows, `[[`, "summary"))
  rate_bind <- bind_rows(lapply(rows, function(x) x$rate_lines))
  list(summary = summary_bind, rate_lines = rate_bind)
}

default_spatial_decay_showcase_specs <- function() {
  # Three spatial kernels at fixed (K_0,K_1)=(0.8,0.2):
  # 1) exponential (Gaussian radial factor)
  # 2) mean-matched power-law (same expected distance, heavier tails)
  # 3) power-law with deliberately larger mean distance (smaller q, near q=1.5)
  data.frame(
    showcase_id = c("spatial_exp", "spatial_pl_matched", "spatial_pl_fat"),
    decay_label = c(
      "Exponential spatial",
      "Power-law spatial (mean-matched)",
      "Power-law spatial (larger mean distance)"
    ),
    control_k = ROBUSTNESS_CONTROL_K,
    treated_k = ROBUSTNESS_TREATED_K_DEFAULT,
    sim_kernel = c("exponential", "exponential", "exponential"),
    sim_spatial_kernel = c("exponential", "power_law", "power_law"),
    # q -> 1.5+ for heavier tails; fat case also uses inflated mean distance below.
    spatial_q = c(NA_real_, 1.8, 1.55),
    match_spatial_mean = c(FALSE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
}

default_temporal_decay_showcase_specs <- function() {
  # Three temporal kernels at fixed (K_0,K_1)=(0.8,0.2):
  # 1) exponential
  # 2) mean-matched power-law (same expected lag, heavier tails; needs p>2)
  # 3) power-law with deliberately larger mean lag (p closer to 2)
  data.frame(
    showcase_id = c("temporal_exp", "temporal_pl_matched", "temporal_pl_fat"),
    decay_label = c(
      "Exponential temporal",
      "Power-law temporal (mean-matched)",
      "Power-law temporal (larger mean lag)"
    ),
    control_k = ROBUSTNESS_CONTROL_K,
    treated_k = ROBUSTNESS_TREATED_K_DEFAULT,
    sim_kernel = c("exponential", "power_law", "power_law"),
    sim_spatial_kernel = c("exponential", "exponential", "exponential"),
    # p -> 2+ for heavier Omori tails; fat case also uses inflated mean lag below.
    power_p = c(NA_real_, 2.4, 2.05),
    match_temporal_mean = c(FALSE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
}

pick_decay_template_rds <- function(manifest_df, out_dir) {
  preferred <- manifest_df %>%
    dplyr::filter(
      .data$scenario_family == "k_separation",
      abs(.data$treated_k - ROBUSTNESS_TREATED_K_DEFAULT) < 1e-8,
      .data$sim_kernel == "exponential"
    )
  candidates <- if (nrow(preferred) > 0L) preferred else manifest_df
  for (i in seq_len(nrow(candidates))) {
    path <- resolve_manifest_rds_path(candidates[i, , drop = FALSE], out_dir)
    if (file.exists(path)) return(path)
  }
  stop("No local RDS template found for decay showcase refresh.")
}

refresh_decay_showcase <- function(manifest_df, specs, reps, out_dir) {
  template_rds <- pick_decay_template_rds(manifest_df, out_dir)
  message(sprintf("[robustness] decay showcase refresh (%d specs, %d reps) from %s",
                  nrow(specs), reps, basename(template_rds)))
  res0 <- readRDS(template_rds)
  cfg0 <- res0$config
  target_spatial_mean <- hawkes_exponential_spatial_mean(cfg0$HAWKES_ALPHA)
  bind_rows(lapply(seq_len(nrow(specs)), function(i) {
    spec <- specs[i, , drop = FALSE]
    message(sprintf("[robustness] decay showcase %d/%d %s",
                    i, nrow(specs), spec$decay_label))
    spatial_kernel <- if ("sim_spatial_kernel" %in% names(spec)) {
      as.character(spec$sim_spatial_kernel[[1]])
    } else {
      "exponential"
    }
    spatial_q <- if ("spatial_q" %in% names(spec)) as.numeric(spec$spatial_q[[1]]) else NA_real_
    spatial_d <- NULL
    if (identical(spatial_kernel, "power_law") && is.finite(target_spatial_mean) && is.finite(spatial_q)) {
      if (isTRUE(as.logical(spec$match_spatial_mean[[1]]))) {
        spatial_d <- hawkes_power_law_spatial_d_for_mean(target_spatial_mean, spatial_q)
      } else {
        # Fat-tailed: keep heavier tails and deliberately inflate mean distance.
        spatial_d <- hawkes_power_law_spatial_d_for_mean(5 * target_spatial_mean, spatial_q)
      }
    }
    refresh_decay_for_spec(
      rds_path = template_rds,
      control_k = spec$control_k,
      treated_k = spec$treated_k,
      sim_kernel = spec$sim_kernel,
      decay_reps = reps,
      showcase_id = spec$showcase_id,
      decay_label = spec$decay_label,
      sim_spatial_kernel = spatial_kernel,
      spatial_q = if (is.finite(spatial_q)) spatial_q else NULL,
      spatial_d = spatial_d
    )
  }))
}

refresh_temporal_decay_showcase <- function(manifest_df, specs, reps, out_dir) {
  template_rds <- pick_decay_template_rds(manifest_df, out_dir)
  message(sprintf("[robustness] temporal decay showcase refresh (%d specs, %d reps) from %s",
                  nrow(specs), reps, basename(template_rds)))
  res0 <- readRDS(template_rds)
  cfg0 <- res0$config
  # Match the main-study temporal scale (HAWKES_BETA, typically 10 => mean lag 0.1).
  target_temporal_mean <- hawkes_exponential_temporal_mean(cfg0$HAWKES_BETA)
  bind_rows(lapply(seq_len(nrow(specs)), function(i) {
    spec <- specs[i, , drop = FALSE]
    message(sprintf("[robustness] temporal decay showcase %d/%d %s",
                    i, nrow(specs), spec$decay_label))
    power_p <- if ("power_p" %in% names(spec)) as.numeric(spec$power_p[[1]]) else NA_real_
    power_c <- NULL
    sim_kernel <- as.character(spec$sim_kernel[[1]])
    if (identical(sim_kernel, "power_law") &&
        is.finite(target_temporal_mean) && is.finite(power_p) && power_p > 2) {
      if (isTRUE(as.logical(spec$match_temporal_mean[[1]]))) {
        power_c <- hawkes_power_law_c_for_mean(target_temporal_mean, power_p)
      } else {
        # Fat-tailed: keep heavier tails and deliberately inflate mean lag.
        power_c <- hawkes_power_law_c_for_mean(5 * target_temporal_mean, power_p)
      }
    }
    refresh_temporal_decay_for_spec(
      rds_path = template_rds,
      control_k = spec$control_k,
      treated_k = spec$treated_k,
      sim_kernel = sim_kernel,
      decay_reps = reps,
      showcase_id = spec$showcase_id,
      decay_label = spec$decay_label,
      power_c = power_c,
      power_p = if (is.finite(power_p)) power_p else NULL,
      time_bin_width = 0.2
    )
  }))
}

refresh_label_flip_decay_showcase <- function(manifest_df, specs, reps, out_dir) {
  template_rds <- pick_decay_template_rds(manifest_df, out_dir)
  message(sprintf("[robustness] label-flip spatial decay showcase refresh (%d specs, %d reps) from %s",
                  nrow(specs), reps, basename(template_rds)))
  res0 <- readRDS(template_rds)
  cfg0 <- res0$config
  target_spatial_mean <- hawkes_exponential_spatial_mean(cfg0$HAWKES_ALPHA)
  bind_rows(lapply(seq_len(nrow(specs)), function(i) {
    spec <- specs[i, , drop = FALSE]
    message(sprintf("[robustness] label-flip spatial decay showcase %d/%d %s",
                    i, nrow(specs), spec$decay_label))
    spatial_kernel <- if ("sim_spatial_kernel" %in% names(spec)) {
      as.character(spec$sim_spatial_kernel[[1]])
    } else {
      "exponential"
    }
    spatial_q <- if ("spatial_q" %in% names(spec)) as.numeric(spec$spatial_q[[1]]) else NA_real_
    spatial_d <- NULL
    if (identical(spatial_kernel, "power_law") && is.finite(target_spatial_mean) && is.finite(spatial_q)) {
      if (isTRUE(as.logical(spec$match_spatial_mean[[1]]))) {
        spatial_d <- hawkes_power_law_spatial_d_for_mean(target_spatial_mean, spatial_q)
      } else {
        spatial_d <- hawkes_power_law_spatial_d_for_mean(5 * target_spatial_mean, spatial_q)
      }
    }
    refresh_label_flip_decay_for_spec(
      rds_path = template_rds,
      control_k = spec$control_k,
      treated_k = spec$treated_k,
      sim_kernel = spec$sim_kernel,
      decay_reps = reps,
      showcase_id = paste0("label_", spec$showcase_id),
      decay_label = spec$decay_label,
      sim_spatial_kernel = spatial_kernel,
      spatial_q = if (is.finite(spatial_q)) spatial_q else NULL,
      spatial_d = spatial_d
    )
  }))
}

refresh_label_flip_temporal_decay_showcase <- function(manifest_df, specs, reps, out_dir) {
  template_rds <- pick_decay_template_rds(manifest_df, out_dir)
  message(sprintf("[robustness] label-flip temporal decay showcase refresh (%d specs, %d reps) from %s",
                  nrow(specs), reps, basename(template_rds)))
  res0 <- readRDS(template_rds)
  cfg0 <- res0$config
  target_temporal_mean <- hawkes_exponential_temporal_mean(cfg0$HAWKES_BETA)
  bind_rows(lapply(seq_len(nrow(specs)), function(i) {
    spec <- specs[i, , drop = FALSE]
    message(sprintf("[robustness] label-flip temporal decay showcase %d/%d %s",
                    i, nrow(specs), spec$decay_label))
    power_p <- if ("power_p" %in% names(spec)) as.numeric(spec$power_p[[1]]) else NA_real_
    power_c <- NULL
    sim_kernel <- as.character(spec$sim_kernel[[1]])
    if (identical(sim_kernel, "power_law") &&
        is.finite(target_temporal_mean) && is.finite(power_p) && power_p > 2) {
      if (isTRUE(as.logical(spec$match_temporal_mean[[1]]))) {
        power_c <- hawkes_power_law_c_for_mean(target_temporal_mean, power_p)
      } else {
        power_c <- hawkes_power_law_c_for_mean(5 * target_temporal_mean, power_p)
      }
    }
    refresh_label_flip_temporal_decay_for_spec(
      rds_path = template_rds,
      control_k = spec$control_k,
      treated_k = spec$treated_k,
      sim_kernel = sim_kernel,
      decay_reps = reps,
      showcase_id = paste0("label_", spec$showcase_id),
      decay_label = spec$decay_label,
      power_c = power_c,
      power_p = if (is.finite(power_p)) power_p else NULL,
      time_bin_width = 0.2
    )
  }))
}

manifest_path <- file.path(out_dir, paste0(output_basename, "_manifest.csv"))
if (is.null(replot_basename)) {
  write.csv(manifest, manifest_path, row.names = FALSE)
}

read_summary_csv <- function(suffix) {
  path <- file.path(out_dir, paste0(output_basename, suffix))
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE)
}

safe_metric_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) < 1L) return(NULL)
  hit[[1L]]
}

safe_label_col <- function(df) {
  safe_metric_col(df, c("method", "labelling"))
}

core_labelling_ids <- c("oracle", "naive", "SEM_adaptive", "SEM_full", "SEM")
core_labelling_levels <- c("oracle", "naive", "SEM")

normalize_labelling_method <- function(x) {
  x <- as.character(x)
  x[x %in% c("SEM_full", "SEM_adaptive")] <- "SEM"
  x
}

collapse_sem_labelling_rows <- function(df, label_col = NULL) {
  if (is.null(df) || nrow(df) < 1L) return(df)
  if (is.null(label_col)) label_col <- safe_label_col(df)
  if (is.null(label_col) || !(label_col %in% names(df))) return(df)

  key_cols <- intersect(
    c(
      label_col, "scenario_id", "contrast_family", "scenario_family", "control_k",
      "treated_k", "k_delta", "mu_scale", "sim_kernel", "fit_kernel",
      "sim_spatial_kernel", "fit_spatial_kernel", "treatment_assignment", "target_points"
    ),
    names(df)
  )
  if (length(key_cols) < 1L) return(df)

  out <- df %>%
    mutate(!!label_col := normalize_labelling_method(.data[[label_col]]))
  if (nrow(out) == nrow(dplyr::distinct(out, dplyr::across(dplyr::all_of(key_cols))))) {
    out[[label_col]] <- factor(out[[label_col]], levels = core_labelling_levels)
    return(out)
  }

  count_col <- if ("n" %in% names(out)) {
    "n"
  } else if ("n_kept" %in% names(out)) {
    "n_kept"
  } else {
    NULL
  }
  mean_cols <- intersect(
    c(
      "mean_accuracy", "mean_balanced_accuracy", "mean_precision_treated",
      "mean_recall_treated", "mean_specificity_control", "mean_f1_treated",
      "mean_psi_estimate", "mean_psi_truth", "mean_bias", "mean_pct_bias", "mean_abs_error", "mean_pct_error",
      "all_nothing_theory", "tau_1_estim", "ATE_total", "ATE_treatment",
      "ATE_spillover", "ATE_naive", "mean_all_nothing_true_control",
      "tau_1_true_control", "true_all_nothing_ATE", "true_tau_1"
    ),
    names(out)
  )
  passthrough_cols <- setdiff(
    names(out),
    c(key_cols, count_col, mean_cols)
  )

  collapsed <- out %>%
    group_by(dplyr::across(dplyr::all_of(key_cols))) %>%
    summarize(
      dplyr::across(dplyr::any_of(count_col), ~ sum(.x, na.rm = TRUE)),
      dplyr::across(dplyr::any_of(mean_cols), ~ {
        if (!is.null(count_col) && count_col %in% names(cur_data_all())) {
          w <- cur_data_all()[[count_col]]
          if (sum(w, na.rm = TRUE) > 0) {
            return(stats::weighted.mean(.x, w, na.rm = TRUE))
          }
        }
        mean(.x, na.rm = TRUE)
      }),
      dplyr::across(dplyr::any_of(passthrough_cols), ~ dplyr::first(.x)),
      .groups = "drop"
    )
  collapsed[[label_col]] <- factor(collapsed[[label_col]], levels = core_labelling_levels)
  collapsed
}

subset_core_labelling <- function(df, label_col = NULL) {
  if (is.null(df) || nrow(df) < 1L) return(df[0, , drop = FALSE])
  if (is.null(label_col)) label_col <- safe_label_col(df)
  if (is.null(label_col) || !(label_col %in% names(df))) return(df[0, , drop = FALSE])
  out <- df %>% filter(.data[[label_col]] %in% core_labelling_ids)
  if (nrow(out) < 1L) return(out)
  collapse_sem_labelling_rows(out, label_col)
}

apply_core_labelling_subset <- function(df) {
  lbl_col <- safe_label_col(df)
  if (is.null(lbl_col)) return(df)
  subset_core_labelling(df, lbl_col)
}

summary_rows <- if (!is.null(replot_basename)) {
  NULL
} else {
  lapply(seq_len(nrow(manifest)), function(i) {
    res <- readRDS(manifest$rds_path[[i]])
    cls <- res$class_metrics
    sup <- res$support_contrast_summary
    ate <- res$summary_df
    decay <- if (!is.null(res$decay_validation)) res$decay_validation$summary else NULL
    true_all_nothing <- if (!is.null(res$all_nothing_ATE)) res$all_nothing_ATE else NA_real_
    true_tau_1 <- if (!is.null(res$true_tau_1)) res$true_tau_1 else NA_real_
    list(
      label = if (!is.null(cls)) transform(cls, scenario_id = manifest$scenario_id[[i]]) else NULL,
      support = if (!is.null(sup)) transform(sup, scenario_id = manifest$scenario_id[[i]]) else NULL,
      ate = if (!is.null(ate)) transform(
        ate,
        scenario_id = manifest$scenario_id[[i]],
        true_all_nothing_ATE = true_all_nothing,
        true_tau_1 = true_tau_1
      ) else NULL,
      decay = if (!is.null(decay)) transform(decay, scenario_id = manifest$scenario_id[[i]]) else NULL
    )
  })
}
spatial_decay_showcase_specs <- default_spatial_decay_showcase_specs()
temporal_decay_showcase_specs <- default_temporal_decay_showcase_specs()
manifest_has_k_separation <- any(manifest$scenario_family == "k_separation", na.rm = TRUE)
auto_refresh_decay_showcase <- is.null(replot_basename) && isTRUE(manifest_has_k_separation)
decay_showcase_summary <- if (isTRUE(refresh_decay) || isTRUE(auto_refresh_decay_showcase)) {
  if (isTRUE(auto_refresh_decay_showcase) && !isTRUE(refresh_decay)) {
    message("[robustness] auto-refreshing spatial decay showcase for appendix figures")
  }
  refresh_decay_showcase(manifest, spatial_decay_showcase_specs, decay_reps, out_dir)
} else if (!is.null(replot_basename)) {
  showcase_csv <- read_summary_csv("_decay_showcase_summary.csv")
  if (!is.null(showcase_csv) && nrow(showcase_csv) > 0L) {
    showcase_csv
  } else {
    message("[robustness] decay showcase CSV missing; refreshing showcase specs")
    refresh_decay_showcase(manifest, spatial_decay_showcase_specs, decay_reps, out_dir)
  }
} else {
  NULL
}
if ((is.null(decay_showcase_summary) || nrow(decay_showcase_summary) < 1L) &&
    (isTRUE(refresh_decay) || isTRUE(auto_refresh_decay_showcase))) {
  stop("Decay showcase refresh produced no rows.")
}
decay_temporal_showcase_summary <- if (isTRUE(refresh_decay) || isTRUE(auto_refresh_decay_showcase)) {
  if (isTRUE(auto_refresh_decay_showcase) && !isTRUE(refresh_decay)) {
    message("[robustness] auto-refreshing temporal decay showcase for appendix figures")
  }
  refresh_temporal_decay_showcase(manifest, temporal_decay_showcase_specs, decay_reps, out_dir)
} else if (!is.null(replot_basename)) {
  temporal_csv <- read_summary_csv("_decay_temporal_showcase_summary.csv")
  if (!is.null(temporal_csv) && nrow(temporal_csv) > 0L) {
    temporal_csv
  } else {
    message("[robustness] temporal decay showcase CSV missing; refreshing showcase specs")
    refresh_temporal_decay_showcase(manifest, temporal_decay_showcase_specs, decay_reps, out_dir)
  }
} else {
  NULL
}
if ((is.null(decay_temporal_showcase_summary) || nrow(decay_temporal_showcase_summary) < 1L) &&
    (isTRUE(refresh_decay) || isTRUE(auto_refresh_decay_showcase))) {
  stop("Temporal decay showcase refresh produced no rows.")
}
decay_label_flip_showcase_summary <- if (isTRUE(refresh_decay) || isTRUE(auto_refresh_decay_showcase)) {
  if (isTRUE(auto_refresh_decay_showcase) && !isTRUE(refresh_decay)) {
    message("[robustness] auto-refreshing label-flip spatial decay showcase for appendix figures")
  }
  refresh_label_flip_decay_showcase(manifest, spatial_decay_showcase_specs, decay_reps, out_dir)
} else if (!is.null(replot_basename)) {
  label_flip_csv <- read_summary_csv("_decay_label_flip_showcase_summary.csv")
  if (!is.null(label_flip_csv) && nrow(label_flip_csv) > 0L) {
    label_flip_csv
  } else {
    message("[robustness] label-flip spatial decay showcase CSV missing; refreshing")
    refresh_label_flip_decay_showcase(manifest, spatial_decay_showcase_specs, decay_reps, out_dir)
  }
} else {
  NULL
}
if ((is.null(decay_label_flip_showcase_summary) || nrow(decay_label_flip_showcase_summary) < 1L) &&
    (isTRUE(refresh_decay) || isTRUE(auto_refresh_decay_showcase))) {
  stop("Label-flip spatial decay showcase refresh produced no rows.")
}
decay_label_flip_temporal_showcase_summary <- if (isTRUE(refresh_decay) || isTRUE(auto_refresh_decay_showcase)) {
  if (isTRUE(auto_refresh_decay_showcase) && !isTRUE(refresh_decay)) {
    message("[robustness] auto-refreshing label-flip temporal decay showcase for appendix figures")
  }
  refresh_label_flip_temporal_decay_showcase(manifest, temporal_decay_showcase_specs, decay_reps, out_dir)
} else if (!is.null(replot_basename)) {
  label_flip_temporal_csv <- read_summary_csv("_decay_label_flip_temporal_showcase_summary.csv")
  if (!is.null(label_flip_temporal_csv) && nrow(label_flip_temporal_csv) > 0L) {
    label_flip_temporal_csv
  } else {
    message("[robustness] label-flip temporal decay showcase CSV missing; refreshing")
    refresh_label_flip_temporal_decay_showcase(manifest, temporal_decay_showcase_specs, decay_reps, out_dir)
  }
} else {
  NULL
}
if ((is.null(decay_label_flip_temporal_showcase_summary) || nrow(decay_label_flip_temporal_showcase_summary) < 1L) &&
    (isTRUE(refresh_decay) || isTRUE(auto_refresh_decay_showcase))) {
  stop("Label-flip temporal decay showcase refresh produced no rows.")
}
label_summary <- if (!is.null(replot_basename)) {
  rebuilt <- load_label_summary_from_manifest(manifest, out_dir)
  if (!is.null(rebuilt) && nrow(rebuilt) > 0L) {
    message(sprintf(
      "[robustness] label recovery rebuilt from per-scenario RDS (label mu filter=%.1fx)",
      LABEL_MU_FAIL_MULTIPLIER
    ))
    rebuilt
  } else {
    read_summary_csv("_label_recovery_summary.csv")
  }
} else {
  bind_rows(lapply(summary_rows, `[[`, "label")) %>%
    left_join(manifest, by = "scenario_id")
}
support_summary <- if (!is.null(replot_basename)) {
  rebuilt <- load_support_summary_from_manifest(manifest, out_dir)
  if (!is.null(rebuilt) && nrow(rebuilt) > 0L) {
    message(sprintf(
      "[robustness] support contrast summary rebuilt (true control; label mu filter=%.1fx)",
      LABEL_MU_FAIL_MULTIPLIER
    ))
    rebuilt
  } else {
    read_summary_csv("_support_contrast_summary.csv")
  }
} else {
  bind_rows(lapply(summary_rows, `[[`, "support")) %>%
    left_join(manifest, by = "scenario_id")
}
ate_summary <- if (!is.null(replot_basename)) {
  read_summary_csv("_ate_summary.csv")
} else {
  bind_rows(lapply(summary_rows, `[[`, "ate")) %>%
    left_join(manifest, by = "scenario_id")
}
decay_summary <- if (!is.null(replot_basename)) {
  read_summary_csv("_decay_validation_summary.csv")
} else {
  decay_rows <- Filter(Negate(is.null), lapply(summary_rows, `[[`, "decay"))
  if (length(decay_rows) > 0L) {
    bind_rows(decay_rows) %>%
      left_join(manifest, by = "scenario_id")
  } else {
    NULL
  }
}

if (!is.null(label_summary) && nrow(label_summary) > 0) {
  label_summary <- apply_core_labelling_subset(label_summary)
}
if (!is.null(support_summary) && nrow(support_summary) > 0) {
  support_summary <- apply_core_labelling_subset(support_summary)
}
if (!is.null(ate_summary) && nrow(ate_summary) > 0) {
  ate_summary <- apply_core_labelling_subset(ate_summary)
}

if (is.null(replot_basename) || (!is.null(label_summary) && nrow(label_summary) > 0L)) {
  write.csv(label_summary, file.path(out_dir, paste0(output_basename, "_label_recovery_summary.csv")), row.names = FALSE)
}
if (is.null(replot_basename) || (!is.null(support_summary) && nrow(support_summary) > 0L)) {
  write.csv(support_summary, file.path(out_dir, paste0(output_basename, "_support_contrast_summary.csv")), row.names = FALSE)
}
if (is.null(replot_basename)) {
  write.csv(ate_summary, file.path(out_dir, paste0(output_basename, "_ate_summary.csv")), row.names = FALSE)
  if (!is.null(decay_summary)) {
    write.csv(decay_summary, file.path(out_dir, paste0(output_basename, "_decay_validation_summary.csv")), row.names = FALSE)
  }
  if (!is.null(decay_showcase_summary) && nrow(decay_showcase_summary) > 0) {
    write.csv(decay_showcase_summary, file.path(out_dir, paste0(output_basename, "_decay_showcase_summary.csv")), row.names = FALSE)
  }
  if (!is.null(decay_temporal_showcase_summary) && nrow(decay_temporal_showcase_summary) > 0) {
    write.csv(decay_temporal_showcase_summary, file.path(out_dir, paste0(output_basename, "_decay_temporal_showcase_summary.csv")), row.names = FALSE)
  }
  if (!is.null(decay_label_flip_showcase_summary) && nrow(decay_label_flip_showcase_summary) > 0) {
    write.csv(decay_label_flip_showcase_summary, file.path(out_dir, paste0(output_basename, "_decay_label_flip_showcase_summary.csv")), row.names = FALSE)
  }
  if (!is.null(decay_label_flip_temporal_showcase_summary) && nrow(decay_label_flip_temporal_showcase_summary) > 0) {
    write.csv(decay_label_flip_temporal_showcase_summary, file.path(out_dir, paste0(output_basename, "_decay_label_flip_temporal_showcase_summary.csv")), row.names = FALSE)
  }
} else if (isTRUE(refresh_decay) || isTRUE(auto_refresh_decay_showcase)) {
  if (!is.null(decay_showcase_summary) && nrow(decay_showcase_summary) > 0) {
    write.csv(decay_showcase_summary, file.path(out_dir, paste0(output_basename, "_decay_showcase_summary.csv")), row.names = FALSE)
  }
  if (!is.null(decay_temporal_showcase_summary) && nrow(decay_temporal_showcase_summary) > 0) {
    write.csv(decay_temporal_showcase_summary, file.path(out_dir, paste0(output_basename, "_decay_temporal_showcase_summary.csv")), row.names = FALSE)
  }
  if (!is.null(decay_label_flip_showcase_summary) && nrow(decay_label_flip_showcase_summary) > 0) {
    write.csv(decay_label_flip_showcase_summary, file.path(out_dir, paste0(output_basename, "_decay_label_flip_showcase_summary.csv")), row.names = FALSE)
  }
  if (!is.null(decay_label_flip_temporal_showcase_summary) && nrow(decay_label_flip_temporal_showcase_summary) > 0) {
    write.csv(decay_label_flip_temporal_showcase_summary, file.path(out_dir, paste0(output_basename, "_decay_label_flip_temporal_showcase_summary.csv")), row.names = FALSE)
  }
} else {
  # Persist showcase CSVs refreshed because they were missing during replot.
  if (!is.null(decay_showcase_summary) && nrow(decay_showcase_summary) > 0) {
    write.csv(decay_showcase_summary, file.path(out_dir, paste0(output_basename, "_decay_showcase_summary.csv")), row.names = FALSE)
  }
  if (!is.null(decay_temporal_showcase_summary) && nrow(decay_temporal_showcase_summary) > 0) {
    write.csv(decay_temporal_showcase_summary, file.path(out_dir, paste0(output_basename, "_decay_temporal_showcase_summary.csv")), row.names = FALSE)
  }
  if (!is.null(decay_label_flip_showcase_summary) && nrow(decay_label_flip_showcase_summary) > 0) {
    write.csv(decay_label_flip_showcase_summary, file.path(out_dir, paste0(output_basename, "_decay_label_flip_showcase_summary.csv")), row.names = FALSE)
  }
  if (!is.null(decay_label_flip_temporal_showcase_summary) && nrow(decay_label_flip_temporal_showcase_summary) > 0) {
    write.csv(decay_label_flip_temporal_showcase_summary, file.path(out_dir, paste0(output_basename, "_decay_label_flip_temporal_showcase_summary.csv")), row.names = FALSE)
  }
}

save_plot_pair <- function(plot, stem, width = 7.2, height = 4.6) {
  if (is.null(plot) || !inherits(plot, "ggplot")) return(NULL)
  # Simulation details belong in LaTeX captions, not plot subtitles.
  plot <- plot + ggplot2::labs(subtitle = NULL)
  pdf_path <- file.path(fig_dir, paste0(stem, ".pdf"))
  png_path <- file.path(fig_dir, paste0(stem, ".png"))
  ggplot2::ggsave(pdf_path, plot, width = width, height = height, units = "in")
  ggplot2::ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 300)
  list(stem = stem, pdf = pdf_path, png = png_path)
}

make_k_separation_line_plot <- function(df, x_col, y_col, lbl_col, title, subtitle, ylab,
                                        family = "k_separation") {
  plot_df <- df %>%
    filter(.data$scenario_family == family)
  if (nrow(plot_df) < 1L) return(NULL)
  ggplot(plot_df, aes(x = .data[[x_col]], y = .data[[y_col]],
                      color = .data[[lbl_col]], group = .data[[lbl_col]])) +
    geom_point(size = 1.8, alpha = 0.85) +
    {if (dplyr::n_distinct(plot_df[[x_col]]) > 1L) geom_line(linewidth = 0.7, alpha = 0.8)} +
    labs(
      x = expression(abs(K[treated] - K[control])),
      y = ylab,
      color = "Labelling",
      title = title,
      subtitle = subtitle
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

kernel_pair_label <- function(sim_kernel, fit_kernel) {
  short_kernel <- function(k) {
    if (identical(k, "power_law")) "power-law" else as.character(k)
  }
  paste0(short_kernel(sim_kernel), " sim / ", short_kernel(fit_kernel), " fit")
}

# Compact truth/fit labels for the joint spatiotemporal kernel matrix.
st_kernel_short <- function(kernel) {
  k <- tolower(as.character(kernel))
  ifelse(k %in% c("power_law", "power-law", "powerlaw"), "PL", "Exp")
}

st_kernel_combo_label <- function(temporal_kernel, spatial_kernel) {
  paste0("T:", st_kernel_short(temporal_kernel), ", S:", st_kernel_short(spatial_kernel))
}

st_kernel_combo_levels <- c(
  "T:Exp, S:Exp",
  "T:PL, S:Exp",
  "T:Exp, S:PL",
  "T:PL, S:PL"
)

canonical_kernel_pair_levels <- function() {
  expand.grid(
    sim_kernel = c("exponential", "power_law"),
    fit_kernel = c("exponential", "power_law"),
    stringsAsFactors = FALSE
  ) %>%
    mutate(kernel_pair = kernel_pair_label(.data$sim_kernel, .data$fit_kernel)) %>%
    pull(.data$kernel_pair)
}

spatial_kernel_pair_label <- function(sim_spatial_kernel, fit_spatial_kernel) {
  kernel_pair_label(sim_spatial_kernel, fit_spatial_kernel)
}

canonical_spatial_kernel_pair_levels <- function() {
  canonical_kernel_pair_levels()
}

spatial_decay_showcase_label_levels <- function() {
  c(
    "Exponential spatial",
    "Power-law spatial (mean-matched)",
    "Power-law spatial (larger mean distance)"
  )
}

temporal_decay_showcase_label_levels <- function() {
  c(
    "Exponential temporal",
    "Power-law temporal (mean-matched)",
    "Power-law temporal (larger mean lag)"
  )
}

relabel_decay_showcase_labels <- function(labels) {
  dplyr::recode(
    as.character(labels),
    "Fat-tailed spatial (larger mean distance)" = "Power-law spatial (larger mean distance)",
    "Fat-tailed temporal (larger mean lag)" = "Power-law temporal (larger mean lag)",
    .default = as.character(labels)
  )
}

prepare_decay_showcase_plot_df <- function(showcase_df, max_distance = NULL) {
  if (is.null(showcase_df) || nrow(showcase_df) < 1L) return(NULL)
  label_levels <- spatial_decay_showcase_label_levels()
  y_col <- if ("mean_abs_delta_plot" %in% names(showcase_df)) {
    "mean_abs_delta_plot"
  } else if ("mean_abs_delta" %in% names(showcase_df)) {
    "mean_abs_delta"
  } else {
    "mean_delta"
  }
  if (is.null(max_distance) || !is.finite(max_distance)) {
    if ("max_unclipped_distance" %in% names(showcase_df)) {
      max_distance <- min(showcase_df$max_unclipped_distance, na.rm = TRUE)
    } else {
      # Central cell on the default 100x100 / 10x10 grid: exterior buffer
      # hits Omega at distance 40; older showcase CSVs lack this column.
      max_distance <- 40
    }
  }
  out <- showcase_df %>%
    dplyr::filter(is.finite(.data$d_mid), is.finite(.data[[y_col]])) %>%
    dplyr::mutate(
      decay_label = factor(
        relabel_decay_showcase_labels(.data$decay_label),
        levels = label_levels
      ),
      mean_abs_delta_plot = .data[[y_col]]
    ) %>%
    dplyr::filter(.data$d_mid <= max_distance, .data$mean_abs_delta_plot > 0)
  if ("area_keep_frac" %in% names(out)) {
    out <- out %>%
      dplyr::filter(is.na(.data$area_keep_frac) | .data$area_keep_frac >= 0.95)
  }
  if (nrow(out) < 1L) return(NULL)
  pos <- out$mean_abs_delta_plot[is.finite(out$mean_abs_delta_plot) & out$mean_abs_delta_plot > 0]
  eps <- if (length(pos) > 0L) min(pos, na.rm = TRUE) / 2 else 1e-6
  if (!is.finite(eps) || eps <= 0) eps <- 1e-6
  out <- out %>%
    dplyr::mutate(mean_abs_delta_plot = pmax(.data$mean_abs_delta_plot, eps)) %>%
    dplyr::arrange(.data$decay_label, .data$d_mid)
  attr(out, "max_distance") <- max_distance
  attr(out, "y_floor") <- eps
  out
}

make_decay_validation_plot <- function(
    showcase_df,
    decay_reps_used = NA_integer_,
    title = "Single-cell allocation flip: spatial decay",
    subtitle = "One treated cell flipped to control under CRN; area-normalized |Delta N|; clipped to untruncated annuli",
    xlab = "Distance from flipped cell") {
  decay_plot_df <- prepare_decay_showcase_plot_df(showcase_df)
  if (is.null(decay_plot_df) || nrow(decay_plot_df) < 1L) return(NULL)
  xmax <- attr(decay_plot_df, "max_distance")
  if (!is.finite(xmax)) xmax <- max(decay_plot_df$d_mid, na.rm = TRUE)
  y_floor <- attr(decay_plot_df, "y_floor")
  if (!is.finite(y_floor) || y_floor <= 0) {
    y_floor <- min(decay_plot_df$mean_abs_delta_plot, na.rm = TRUE) / 2
  }
  y_hi <- max(decay_plot_df$mean_abs_delta_plot, na.rm = TRUE)
  y_breaks <- c(0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3)
  y_breaks <- y_breaks[y_breaks >= y_floor * 0.9 & y_breaks <= y_hi * 1.1]
  if (length(y_breaks) < 2L) {
    y_breaks <- scales::log_breaks(4)(c(y_floor, y_hi))
  }
  ggplot(decay_plot_df, aes(x = .data$d_mid, y = .data$mean_abs_delta_plot,
                            color = .data$decay_label,
                            group = .data$decay_label)) +
    geom_line(linewidth = 0.8, alpha = 0.95) +
    geom_point(size = 1.2, alpha = 0.9) +
    scale_y_log10(
      breaks = y_breaks,
      limits = c(y_floor, y_hi * 1.12),
      labels = scales::label_number(accuracy = 0.01)
    ) +
    coord_cartesian(xlim = c(0, xmax)) +
    labs(
      x = xlab,
      y = expression(Mean~"|"*Delta*N*"|"/area~per~annulus~(log[10])),
      color = "Spatial kernel",
      title = title,
      subtitle = subtitle
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

prepare_decay_showcase_cumulative_plot_df <- function(showcase_df, max_distance = NULL) {
  if (is.null(showcase_df) || nrow(showcase_df) < 1L) return(NULL)
  label_levels <- spatial_decay_showcase_label_levels()
  y_col <- if ("mean_abs_delta_plot" %in% names(showcase_df)) {
    "mean_abs_delta_plot"
  } else if ("mean_abs_delta" %in% names(showcase_df)) {
    "mean_abs_delta"
  } else {
    "mean_delta"
  }
  if (is.null(max_distance) || !is.finite(max_distance)) {
    if ("max_unclipped_distance" %in% names(showcase_df)) {
      max_distance <- min(showcase_df$max_unclipped_distance, na.rm = TRUE)
    } else {
      max_distance <- 40
    }
  }
  out <- showcase_df %>%
    dplyr::filter(is.finite(.data$d_mid), is.finite(.data[[y_col]])) %>%
    dplyr::mutate(
      decay_label = factor(
        relabel_decay_showcase_labels(.data$decay_label),
        levels = label_levels
      ),
      mean_abs_delta_plot = .data[[y_col]],
      annulus_area_plot = if ("annulus_area" %in% names(.)) .data$annulus_area else 1
    ) %>%
    dplyr::filter(.data$d_mid <= max_distance, .data$mean_abs_delta_plot > 0)
  if ("area_keep_frac" %in% names(out)) {
    out <- out %>%
      dplyr::filter(is.na(.data$area_keep_frac) | .data$area_keep_frac >= 0.95)
  }
  if (nrow(out) < 1L) return(NULL)
  out <- out %>%
    dplyr::group_by(.data$decay_label) %>%
    dplyr::arrange(.data$d_mid, .by_group = TRUE) %>%
    dplyr::mutate(
      bin_mass = .data$mean_abs_delta_plot * .data$annulus_area_plot,
      cumulative_abs_delta = cumsum(.data$bin_mass)
    ) %>%
    dplyr::ungroup()
  attr(out, "max_distance") <- max_distance
  out
}

make_decay_validation_cumulative_plot <- function(
    showcase_df,
    decay_reps_used = NA_integer_,
    title = "Single-cell allocation flip: spatial decay (cumulative)",
    subtitle = "One treated cell flipped to control; cumulative |Delta N| within radius of that cell",
    xlab = "Distance from flipped cell") {
  decay_plot_df <- prepare_decay_showcase_cumulative_plot_df(showcase_df)
  if (is.null(decay_plot_df) || nrow(decay_plot_df) < 1L) return(NULL)
  xmax <- attr(decay_plot_df, "max_distance")
  if (!is.finite(xmax)) xmax <- max(decay_plot_df$d_mid, na.rm = TRUE)
  ggplot(decay_plot_df, aes(x = .data$d_mid, y = .data$cumulative_abs_delta,
                            color = .data$decay_label,
                            group = .data$decay_label)) +
    geom_line(linewidth = 0.8, alpha = 0.95) +
    geom_point(size = 1.2, alpha = 0.9) +
    coord_cartesian(xlim = c(0, xmax)) +
    labs(
      x = xlab,
      y = "Cumulative mean |Delta N| within radius",
      color = "Spatial kernel",
      title = title,
      subtitle = subtitle
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

prepare_temporal_decay_showcase_plot_df <- function(showcase_df, min_lag = 0, max_lag = 10) {
  if (is.null(showcase_df) || nrow(showcase_df) < 1L) return(NULL)
  x_col <- if ("t_mid" %in% names(showcase_df)) "t_mid" else "d_mid"
  y_col <- if ("mean_abs_delta" %in% names(showcase_df)) "mean_abs_delta" else "mean_delta"
  label_levels <- temporal_decay_showcase_label_levels()
  out <- showcase_df %>%
    dplyr::filter(is.finite(.data[[x_col]]), is.finite(.data[[y_col]])) %>%
    dplyr::mutate(
      decay_label = factor(
        relabel_decay_showcase_labels(.data$decay_label),
        levels = label_levels
      ),
      mean_abs_delta_plot = .data[[y_col]],
      x_val = .data[[x_col]]
    ) %>%
    dplyr::filter(
      .data$x_val >= min_lag,
      .data$x_val <= max_lag,
      .data$mean_abs_delta_plot > 0
    )
  if (nrow(out) < 1L) return(NULL)
  out
}

make_temporal_decay_validation_plot <- function(
    showcase_df,
    decay_reps_used = NA_integer_,
    title = "Single-cell allocation flip: temporal decay",
    subtitle = "One treated cell flipped to control under CRN; per-lag-bin |Delta N| after t*",
    xlab = "Lag since treatment time") {
  decay_plot_df <- prepare_temporal_decay_showcase_plot_df(showcase_df)
  if (is.null(decay_plot_df) || nrow(decay_plot_df) < 1L) return(NULL)
  x_col <- if ("t_mid" %in% names(decay_plot_df)) "t_mid" else "d_mid"
  ggplot(decay_plot_df, aes(x = .data[[x_col]], y = .data$mean_abs_delta_plot,
                            color = .data$decay_label,
                            group = .data$decay_label)) +
    geom_line(linewidth = 0.8, alpha = 0.95) +
    geom_point(size = 1.2, alpha = 0.9) +
    scale_y_log10() +
    coord_cartesian(xlim = c(0, 10)) +
    labs(
      x = xlab,
      y = expression(Mean~"|"*Delta*N*"|"~per~lag~bin~(log[10])),
      color = "Temporal kernel",
      title = title,
      subtitle = subtitle
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

prepare_temporal_decay_showcase_cumulative_plot_df <- function(showcase_df, min_lag = 0, max_lag = 10) {
  base_df <- prepare_temporal_decay_showcase_plot_df(
    showcase_df,
    min_lag = min_lag,
    max_lag = max_lag
  )
  if (is.null(base_df) || nrow(base_df) < 1L) return(NULL)
  x_col <- if ("t_mid" %in% names(base_df)) "t_mid" else "d_mid"
  base_df %>%
    dplyr::group_by(.data$decay_label) %>%
    dplyr::arrange(.data[[x_col]], .by_group = TRUE) %>%
    dplyr::mutate(cumulative_abs_delta = cumsum(.data$mean_abs_delta_plot)) %>%
    dplyr::ungroup()
}

make_temporal_decay_validation_cumulative_plot <- function(
    showcase_df,
    decay_reps_used = NA_integer_,
    title = "Single-cell allocation flip: temporal decay (cumulative)",
    subtitle = "One treated cell flipped to control; cumulative |Delta N| up to each lag after t*",
    xlab = "Lag since treatment time") {
  decay_plot_df <- prepare_temporal_decay_showcase_cumulative_plot_df(showcase_df)
  if (is.null(decay_plot_df) || nrow(decay_plot_df) < 1L) return(NULL)
  x_col <- if ("t_mid" %in% names(decay_plot_df)) "t_mid" else "d_mid"
  ggplot(decay_plot_df, aes(x = .data[[x_col]], y = .data$cumulative_abs_delta,
                            color = .data$decay_label,
                            group = .data$decay_label)) +
    geom_line(linewidth = 0.8, alpha = 0.95) +
    geom_point(size = 1.2, alpha = 0.9) +
    coord_cartesian(xlim = c(0, 10)) +
    labs(
      x = xlab,
      y = "Cumulative mean |Delta N| up to lag",
      color = "Temporal kernel",
      title = title,
      subtitle = subtitle
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

assignment_display_label <- function(assignment) {
  vapply(as.character(assignment), function(a) {
    if (is.na(a) || !nzchar(a)) return("Highest count")
    switch(
      a,
      highest_count_50pct = "Highest count",
      highest_count = "Highest count",
      lowest_count_50pct = "Lowest count",
      lowest_count = "Lowest count",
      count_propensity_50pct = "Count propensity",
      count_propensity = "Count propensity",
      propensity_50pct = "Count propensity",
      contiguous_aoi_50pct = "Contiguous AOI",
      contiguous_aoi = "Contiguous AOI",
      contiguous_high_count_50pct = "Contiguous AOI",
      voronoi_random_50pct = "Voronoi random",
      voronoi_random = "Voronoi random",
      voronoi = "Voronoi random",
      pretreatment_voronoi_random = "Voronoi random",
      random = "Random",
      random_50pct = "Random",
      a
    )
  }, character(1L), USE.NAMES = FALSE)
}

assignment_display_levels <- c(
  "Highest count",
  "Lowest count",
  "Count propensity",
  "Contiguous AOI",
  "Voronoi random",
  "Random"
)

make_high_count_assignment_bar_plot <- function(df, y_col, lbl_col, title, subtitle, ylab) {
  plot_df <- df %>%
    filter(.data$scenario_family %in% c("pretreatment_assignment", "high_count_assignment"))
  if (nrow(plot_df) < 1L) return(NULL)
  plot_df <- plot_df %>%
    mutate(scenario_label = assignment_display_label(.data$treatment_assignment))
  plot_df$scenario_label <- factor(
    plot_df$scenario_label,
    levels = intersect(assignment_display_levels, unique(plot_df$scenario_label))
  )
  ggplot(plot_df, aes(x = .data$scenario_label, y = .data[[y_col]], fill = .data[[lbl_col]])) +
    geom_col(position = position_dodge(width = 0.8), width = 0.55) +
    labs(
      x = NULL,
      y = ylab,
      fill = "Labelling",
      title = title,
      subtitle = subtitle
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "bottom")
}

make_snr_scale_line_plot <- function(df, y_col, lbl_col, title, subtitle, ylab,
                                     family = "snr_scale") {
  plot_df <- df %>%
    filter(.data$scenario_family == family)
  if (nrow(plot_df) < 1L) return(NULL)
  ggplot(plot_df, aes(x = .data$mu_scale, y = .data[[y_col]],
                      color = .data[[lbl_col]], group = .data[[lbl_col]])) +
    geom_point(size = 1.8, alpha = 0.85) +
    {if (dplyr::n_distinct(plot_df$mu_scale) > 1L) geom_line(linewidth = 0.7, alpha = 0.8)} +
    labs(
      x = expression(paste(mu, " scale")),
      y = ylab,
      color = "Labelling",
      title = title,
      subtitle = subtitle
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

make_kernel_mismatch_bar_plot <- function(df, y_col, lbl_col, title, subtitle, ylab) {
  plot_df <- df %>%
    filter(.data$scenario_family == "kernel_mismatch") %>%
    mutate(kernel_pair = kernel_pair_label(.data$sim_kernel, .data$fit_kernel))
  if (nrow(plot_df) < 1L) return(NULL)
  plot_df$kernel_pair <- factor(plot_df$kernel_pair, levels = canonical_kernel_pair_levels())
  ggplot(plot_df, aes(x = .data$kernel_pair, y = .data[[y_col]], fill = .data[[lbl_col]])) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    labs(
      x = "Kernel scenario",
      y = ylab,
      fill = "Labelling",
      title = title,
      subtitle = subtitle
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
}

make_spatial_kernel_mismatch_bar_plot <- function(df, y_col, lbl_col, title, subtitle, ylab) {
  if (!all(c("sim_spatial_kernel", "fit_spatial_kernel") %in% names(df))) return(NULL)
  plot_df <- df %>%
    filter(.data$scenario_family == "spatial_kernel_mismatch") %>%
    mutate(kernel_pair = spatial_kernel_pair_label(.data$sim_spatial_kernel, .data$fit_spatial_kernel))
  if (nrow(plot_df) < 1L) return(NULL)
  plot_df$kernel_pair <- factor(plot_df$kernel_pair, levels = canonical_spatial_kernel_pair_levels())
  ggplot(plot_df, aes(x = .data$kernel_pair, y = .data[[y_col]], fill = .data[[lbl_col]])) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    labs(
      x = "Spatial kernel scenario",
      y = ylab,
      fill = "Labelling",
      title = title,
      subtitle = subtitle
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
}

# Short display names for the all-scenario labelling summary plot.
scenario_display_name <- function(scenario_family, treated_k, k_delta, mu_scale,
                                  sim_kernel, fit_kernel,
                                  sim_spatial_kernel, fit_spatial_kernel,
                                  scenario_id, treatment_assignment = NULL) {
  fam <- as.character(scenario_family)
  n <- length(fam)
  out <- character(n)
  assign <- if (is.null(treatment_assignment)) {
    rep(NA_character_, n)
  } else {
    as.character(treatment_assignment)
  }
  for (i in seq_len(n)) {
    if (fam[[i]] %in% c("pretreatment_assignment", "high_count_assignment")) {
      out[[i]] <- paste0("Assign: ", assignment_display_label(assign[[i]]))
    } else if (identical(fam[[i]], "snr_scale")) {
      out[[i]] <- sprintf("SNR: mu x %s", format(mu_scale[[i]], trim = TRUE, scientific = FALSE))
    } else if (identical(fam[[i]], "k_separation")) {
      out[[i]] <- sprintf(
        "K-sep: K_1=%s (dK=%s)",
        format(treated_k[[i]], trim = TRUE, scientific = FALSE),
        format(k_delta[[i]], trim = TRUE, scientific = FALSE)
      )
    } else if (fam[[i]] %in% c(
      "spatiotemporal_kernel_mismatch", "kernel_mismatch", "spatial_kernel_mismatch"
    )) {
      out[[i]] <- sprintf(
        "Kernels: truth %s | fit %s",
        st_kernel_combo_label(sim_kernel[[i]], sim_spatial_kernel[[i]]),
        st_kernel_combo_label(fit_kernel[[i]], fit_spatial_kernel[[i]])
      )
    } else {
      out[[i]] <- as.character(scenario_id[[i]])
    }
  }
  out
}

scenario_family_order <- c(
  "pretreatment_assignment",
  "high_count_assignment",
  "snr_scale",
  "k_separation",
  "spatiotemporal_kernel_mismatch",
  "kernel_mismatch",
  "spatial_kernel_mismatch",
  "effect_modification",
  "geometry_transport"
)

# Structured DTAITE studies write label summaries under generated/robustness/tables/.
# Fold them into the all-scenario overview scatter (Figure 1).
load_structured_dtaite_label_rows <- function(robustness_dir) {
  table_dir <- file.path(robustness_dir, "tables")
  rows <- list()

  effect_csv <- file.path(table_dir, "effect_modification_label_summary.csv")
  if (file.exists(effect_csv)) {
    ef <- utils::read.csv(effect_csv, stringsAsFactors = FALSE)
    keep <- ef$stratum == "all" & ef$method %in% c("naive", "SEM")
    ef <- ef[keep, , drop = FALSE]
    if (nrow(ef) > 0L) {
      acc <- if ("accuracy" %in% names(ef)) ef$accuracy else ef$balanced_accuracy
      rows[[length(rows) + 1L]] <- data.frame(
        method = as.character(ef$method),
        mean_accuracy = as.numeric(acc),
        mean_balanced_accuracy = as.numeric(ef$balanced_accuracy),
        scenario_id = sprintf(
          "effect_modification_h%s",
          format(ef$h_true, trim = TRUE, scientific = FALSE)
        ),
        scenario_family = "effect_modification",
        stringsAsFactors = FALSE
      )
    }
  }

  geo_csv <- file.path(table_dir, "geometry_transport_label_summary.csv")
  if (file.exists(geo_csv)) {
    geo <- utils::read.csv(geo_csv, stringsAsFactors = FALSE)
    if ("stratum" %in% names(geo)) {
      geo <- geo[geo$stratum == "all", , drop = FALSE]
    }
    geo <- geo[geo$method %in% c("naive", "SEM"), , drop = FALSE]
    if (nrow(geo) > 0L) {
      acc <- if ("accuracy" %in% names(geo)) geo$accuracy else geo$balanced_accuracy
      rows[[length(rows) + 1L]] <- data.frame(
        method = as.character(geo$method),
        mean_accuracy = as.numeric(acc),
        mean_balanced_accuracy = as.numeric(geo$balanced_accuracy),
        scenario_id = "geometry_transport",
        scenario_family = "geometry_transport",
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(rows)) return(NULL)
  bind_rows(rows)
}

make_all_scenarios_label_scatter <- function(df, y_col, lbl_col,
                                             title = "Labelling accuracy across robustness scenarios",
                                             subtitle = "Each point is one fitted scenario.") {
  if (is.null(df) || nrow(df) < 1L) return(NULL)
  plot_df <- df %>%
    mutate(!!lbl_col := normalize_labelling_method(.data[[lbl_col]])) %>%
    filter(.data[[lbl_col]] %in% c("naive", "SEM"))
  if (nrow(plot_df) < 1L) return(NULL)
  if (!"treatment_assignment" %in% names(plot_df)) {
    plot_df$treatment_assignment <- NA_character_
  }

  # Alias legacy family names for a compact legend.
  family_labels <- c(
    pretreatment_assignment = "Pretreatment assignment",
    high_count_assignment = "Pretreatment assignment",
    snr_scale = "Signal-to-noise",
    k_separation = "K-separation",
    spatiotemporal_kernel_mismatch = "Kernel misspecification",
    kernel_mismatch = "Kernel misspecification",
    spatial_kernel_mismatch = "Kernel misspecification",
    effect_modification = "Covariate effect",
    geometry_transport = "Coarseness"
  )
  family_levels <- c(
    "Pretreatment assignment",
    "Signal-to-noise",
    "K-separation",
    "Kernel misspecification",
    "Covariate effect",
    "Coarseness"
  )
  family_colors <- c(
    "Pretreatment assignment" = "#3D405B",
    "Signal-to-noise" = "#E07A5F",
    "K-separation" = "#81B29A",
    "Kernel misspecification" = "#F2CC8F",
    "Covariate effect" = "#9B5DE5",
    "Coarseness" = "#00BBF9"
  )

  wide <- plot_df %>%
    mutate(
      family_label = {
        fam <- as.character(.data$scenario_family)
        mapped <- unname(family_labels[fam])
        ifelse(is.na(mapped), fam, mapped)
      }
    ) %>%
    select(
      "scenario_id",
      "family_label",
      dplyr::all_of(c(lbl_col, y_col))
    ) %>%
    rename(accuracy = !!y_col) %>%
    tidyr::pivot_wider(names_from = !!lbl_col, values_from = "accuracy") %>%
    filter(is.finite(.data$naive), is.finite(.data$SEM))
  if (nrow(wide) < 1L) return(NULL)

  # Drop unused family levels so the legend only shows families present.
  present <- intersect(family_levels, unique(as.character(wide$family_label)))
  wide$family_label <- factor(wide$family_label, levels = present)
  lims <- range(c(wide$naive, wide$SEM), finite = TRUE)
  pad <- max(0.02, 0.05 * diff(lims))
  lims <- pmin(1, pmax(0, c(lims[1] - pad, lims[2] + pad)))

  ggplot(wide, aes(x = .data$naive, y = .data$SEM, color = .data$family_label)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.45) +
    geom_point(size = 2.6, alpha = 0.9) +
    coord_fixed(xlim = lims, ylim = lims) +
    scale_color_manual(values = family_colors[present], drop = FALSE) +
    labs(
      x = "Naive mean raw accuracy",
      y = "SEM mean raw accuracy",
      color = "Scenario family",
      title = title,
      subtitle = subtitle
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}

plot_files <- list()
if (!is.null(label_summary) && nrow(label_summary) > 0) {
  structured_label_rows <- load_structured_dtaite_label_rows(robustness_dir)
  if (!is.null(structured_label_rows) && nrow(structured_label_rows) > 0L) {
    message(sprintf(
      "[robustness] adding %d structured DTAITE label rows to all-scenario overview",
      nrow(structured_label_rows)
    ))
    label_summary <- bind_rows(label_summary, structured_label_rows)
  }
  lbl_col <- safe_label_col(label_summary)
  acc_col <- safe_metric_col(label_summary, c("mean_accuracy", "mean_balanced_accuracy", "accuracy"))
  if (!is.null(lbl_col) && !is.null(acc_col)) {
    plot_files$all_scenarios_label_recovery <- save_plot_pair(
      make_all_scenarios_label_scatter(
        label_summary, acc_col, lbl_col,
        title = "Labelling accuracy across all robustness scenarios",
        subtitle = "SEM vs naive; dashed line is equal accuracy. Points above the line favour SEM."
      ),
      "robustness_all_scenarios_label_recovery",
      width = 6.8, height = 6.4
    )
  }
}

support_contrast_display_labels <- c(
  single_cell_flip = "Single-cell vs all-control",
  random_50pct_flip = "Random 50% vs all-control",
  global_1_0 = "All-or-nothing DTAITE"
)

# Bias support figures focus on the all-or-nothing DTAITE; other contrasts
# are retained in the summary CSVs but omitted from plots.
SUPPORT_PLOT_CONTRAST <- "global_1_0"

filter_support_plot_contrast <- function(df) {
  if (is.null(df) || nrow(df) < 1L) return(df)
  if (!"contrast_family" %in% names(df)) return(df)
  df %>% dplyr::filter(.data$contrast_family == SUPPORT_PLOT_CONTRAST)
}

prepare_support_plot_df <- function(df, family, keep_kernel_pair = FALSE, keep_spatial_pair = FALSE) {
  if (is.null(df) || nrow(df) < 1L) return(NULL)
  lbl_col <- safe_label_col(df)
  if (is.null(lbl_col)) return(NULL)
  out <- df %>%
    filter(.data$scenario_family == family)
  if (nrow(out) < 1L) return(NULL)
  if (!"mean_pct_error" %in% names(out)) {
    out$mean_pct_error <- NA_real_
  }
  global_ref <- out %>%
    filter(.data$contrast_family == "global_1_0") %>%
    distinct(.data$scenario_id, .data$mean_psi_truth) %>%
    rename(global_ate_truth = .data$mean_psi_truth)
  out <- out %>%
    left_join(global_ref, by = "scenario_id") %>%
    mutate(
      mean_pct_error = ifelse(
        is.finite(.data$mean_pct_error),
        .data$mean_pct_error,
        ifelse(
          is.finite(.data$mean_psi_truth) & abs(.data$mean_psi_truth) > 1e-8,
          100 * .data$mean_abs_error / abs(.data$mean_psi_truth),
          ifelse(
            .data$contrast_family == "single_cell_flip" &
              is.finite(.data$global_ate_truth) & abs(.data$global_ate_truth) > 1e-8,
            100 * .data$mean_abs_error / abs(.data$global_ate_truth),
            NA_real_
          )
        )
      )
    ) %>%
    select(-.data$global_ate_truth) %>%
    filter(is.finite(.data$mean_pct_error))
  if (nrow(out) < 1L) return(NULL)
  if (isTRUE(keep_spatial_pair) && all(c("sim_spatial_kernel", "fit_spatial_kernel") %in% names(out))) {
    out <- out %>%
      mutate(kernel_pair = spatial_kernel_pair_label(.data$sim_spatial_kernel, .data$fit_spatial_kernel)) %>%
      group_by(.data$contrast_family, .data[[lbl_col]], .data$kernel_pair) %>%
      summarize(mean_pct_error = mean(.data$mean_pct_error, na.rm = TRUE), .groups = "drop")
    if (nrow(out) < 1L) return(NULL)
    out$contrast_display <- support_contrast_display_labels[out$contrast_family]
    out$contrast_display <- factor(
      out$contrast_display,
      levels = unname(support_contrast_display_labels)
    )
    pair_levels <- canonical_spatial_kernel_pair_levels()
    out$kernel_pair <- factor(out$kernel_pair, levels = pair_levels)
    return(out)
  }
  if (isTRUE(keep_kernel_pair) && all(c("sim_kernel", "fit_kernel") %in% names(out))) {
    out <- out %>%
      mutate(kernel_pair = kernel_pair_label(.data$sim_kernel, .data$fit_kernel)) %>%
      group_by(.data$contrast_family, .data[[lbl_col]], .data$kernel_pair) %>%
      summarize(mean_pct_error = mean(.data$mean_pct_error, na.rm = TRUE), .groups = "drop")
    if (nrow(out) < 1L) return(NULL)
    out$contrast_display <- support_contrast_display_labels[out$contrast_family]
    out$contrast_display <- factor(
      out$contrast_display,
      levels = unname(support_contrast_display_labels)
    )
    pair_levels <- canonical_kernel_pair_levels()
    out$kernel_pair <- factor(out$kernel_pair, levels = pair_levels)
    return(out)
  }
  out <- out %>%
    group_by(.data$contrast_family, .data[[lbl_col]]) %>%
    summarize(mean_pct_error = mean(.data$mean_pct_error, na.rm = TRUE), .groups = "drop")
  if (nrow(out) < 1L) return(NULL)
  out$contrast_family <- factor(
    out$contrast_family,
    levels = names(support_contrast_display_labels),
    labels = unname(support_contrast_display_labels)
  )
  out
}

make_high_count_support_plot <- function(df, title, stem, subtitle = NULL) {
  lbl_col <- safe_label_col(df)
  if (is.null(lbl_col)) return(NULL)
  plot_labelling_levels <- c("naive", "SEM")
  out <- df %>%
    filter(.data$scenario_family %in% c("pretreatment_assignment", "high_count_assignment")) %>%
    filter_support_plot_contrast()
  if (nrow(out) < 1L) return(NULL)
  if (!"mean_bias" %in% names(out)) out$mean_bias <- NA_real_
  if (!"se_bias" %in% names(out)) out$se_bias <- NA_real_
  if (!"treatment_assignment" %in% names(out)) {
    out$treatment_assignment <- "highest_count_50pct"
  }
  support_plot_df <- out %>%
    mutate(
      mean_bias = ifelse(
        is.finite(.data$mean_bias),
        .data$mean_bias,
        ifelse(
          is.finite(.data$mean_psi_estimate) & is.finite(.data$mean_psi_truth),
          .data$mean_psi_estimate - .data$mean_psi_truth,
          NA_real_
        )
      ),
      se_bias = ifelse(is.finite(.data$se_bias), .data$se_bias, NA_real_),
      !!lbl_col := normalize_labelling_method(.data[[lbl_col]]),
      assignment_label = assignment_display_label(.data$treatment_assignment)
    ) %>%
    filter(.data[[lbl_col]] %in% plot_labelling_levels, is.finite(.data$mean_bias)) %>%
    group_by(.data$assignment_label, .data[[lbl_col]]) %>%
    summarize(
      mean_bias = mean(.data$mean_bias, na.rm = TRUE),
      se_bias = dplyr::coalesce(
        if (dplyr::n() == 1L) .data$se_bias[[1]] else NA_real_,
        {
          s <- stats::sd(.data$mean_bias, na.rm = TRUE)
          ifelse(is.finite(s) && dplyr::n() > 1L, s / sqrt(dplyr::n()), NA_real_)
        }
      ),
      .groups = "drop"
    ) %>%
    mutate(
      bias_lo = .data$mean_bias - .data$se_bias,
      bias_hi = .data$mean_bias + .data$se_bias
    )
  if (nrow(support_plot_df) < 1L) return(NULL)

  # Prefer the five pretreatment rules in fixed order when present.
  assign_levels <- intersect(
    c("Highest count", "Lowest count", "Count propensity", "Contiguous AOI", "Voronoi random"),
    unique(as.character(support_plot_df$assignment_label))
  )
  if (length(assign_levels) < 1L) {
    assign_levels <- intersect(assignment_display_levels, unique(as.character(support_plot_df$assignment_label)))
  }
  support_plot_df$assignment_label <- factor(support_plot_df$assignment_label, levels = assign_levels)
  support_plot_df[[lbl_col]] <- factor(support_plot_df[[lbl_col]], levels = plot_labelling_levels)

  has_se <- any(is.finite(support_plot_df$se_bias))
  dodge <- position_dodge(width = 0.35)
  p_support <- ggplot(
    support_plot_df,
    aes(
      x = .data$assignment_label,
      y = .data$mean_bias,
      color = .data[[lbl_col]],
      shape = .data[[lbl_col]],
      group = .data$assignment_label
    )
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_line(color = "grey70", linewidth = 0.45, alpha = 0.9) +
    {if (has_se) {
      geom_errorbar(
        aes(ymin = .data$bias_lo, ymax = .data$bias_hi),
        width = 0.12,
        linewidth = 0.35,
        alpha = 0.8,
        position = dodge
      )
    }} +
    geom_point(size = 2.8, position = dodge) +
    scale_color_manual(values = c(naive = "#E07A5F", SEM = "#3D405B")) +
    scale_shape_manual(values = c(naive = 16, SEM = 17)) +
    labs(
      x = "Pretreatment assignment rule",
      y = "DTAITE bias (estimate - truth)",
      color = "Labelling",
      shape = "Labelling",
      title = title,
      subtitle = subtitle
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "bottom")
  save_plot_pair(p_support, stem, width = 8.0, height = 4.6)
}

make_support_plot <- function(df, family, title, stem, subtitle, keep_kernel_pair = FALSE) {
  lbl_col <- safe_label_col(df)
  if (is.null(lbl_col)) return(NULL)
  support_plot_df <- prepare_support_plot_df(df, family, keep_kernel_pair = keep_kernel_pair)
  if (is.null(support_plot_df) || nrow(support_plot_df) < 1L) return(NULL)
  p_support <- ggplot(support_plot_df, aes(x = .data$contrast_family, y = .data$mean_pct_error,
                                             fill = .data[[lbl_col]])) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    labs(
      x = "Off-support contrast",
      y = "Mean percentage error (%)",
      fill = "Labelling",
      title = title,
      subtitle = subtitle
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
  save_plot_pair(p_support, stem, width = 8.5, height = 5.0)
}

prepare_parameter_sweep_support_df <- function(df) {
  if (is.null(df) || nrow(df) < 1L) return(NULL)
  lbl_col <- safe_label_col(df)
  if (is.null(lbl_col)) return(NULL)
  plot_labelling_levels <- c("naive", "SEM")
  panel_levels <- c(
    "K-separation: |K_1 - K_0|",
    "Signal-to-noise: mu scale"
  )

  finish_panel <- function(out, panel_label, param_col) {
    if (is.null(out) || nrow(out) < 1L) return(NULL)
    if (!param_col %in% names(out)) return(NULL)
    if (!"mean_bias" %in% names(out)) out$mean_bias <- NA_real_
    if (!"se_bias" %in% names(out)) out$se_bias <- NA_real_
    out %>%
      mutate(
        mean_bias = ifelse(
          is.finite(.data$mean_bias),
          .data$mean_bias,
          ifelse(
            is.finite(.data$mean_psi_estimate) & is.finite(.data$mean_psi_truth),
            .data$mean_psi_estimate - .data$mean_psi_truth,
            NA_real_
          )
        ),
        se_bias = ifelse(is.finite(.data$se_bias), .data$se_bias, NA_real_),
        !!lbl_col := normalize_labelling_method(.data[[lbl_col]]),
        param_value = as.numeric(.data[[param_col]]),
        panel = panel_label
      ) %>%
      filter(.data[[lbl_col]] %in% plot_labelling_levels, is.finite(.data$mean_bias), is.finite(.data$param_value)) %>%
      group_by(.data$panel, .data$param_value, .data[[lbl_col]]) %>%
      summarize(
        mean_bias = mean(.data$mean_bias, na.rm = TRUE),
        se_bias = dplyr::coalesce(
          if (dplyr::n() == 1L) .data$se_bias[[1]] else NA_real_,
          {
            s <- stats::sd(.data$mean_bias, na.rm = TRUE)
            ifelse(is.finite(s) && dplyr::n() > 1L, s / sqrt(dplyr::n()), NA_real_)
          }
        ),
        .groups = "drop"
      ) %>%
      mutate(
        bias_lo = .data$mean_bias - .data$se_bias,
        bias_hi = .data$mean_bias + .data$se_bias
      )
  }

  k_out <- NULL
  if ("k_delta" %in% names(df)) {
    k_out <- finish_panel(
      df %>%
        filter(.data$scenario_family == "k_separation") %>%
        filter_support_plot_contrast(),
      panel_levels[[1]],
      "k_delta"
    )
  }
  snr_out <- NULL
  if ("mu_scale" %in% names(df)) {
    snr_out <- finish_panel(
      df %>%
        filter(.data$scenario_family == "snr_scale") %>%
        filter_support_plot_contrast(),
      panel_levels[[2]],
      "mu_scale"
    )
  }
  out <- dplyr::bind_rows(k_out, snr_out)
  if (is.null(out) || nrow(out) < 1L) return(NULL)
  out$panel <- factor(out$panel, levels = panel_levels)
  out[[lbl_col]] <- factor(out[[lbl_col]], levels = plot_labelling_levels)
  out
}

make_parameter_sweep_support_plot <- function(df, title, stem, subtitle) {
  lbl_col <- safe_label_col(df)
  if (is.null(lbl_col)) return(NULL)
  support_plot_df <- prepare_parameter_sweep_support_df(df)
  if (is.null(support_plot_df) || nrow(support_plot_df) < 1L) return(NULL)
  has_se <- any(is.finite(support_plot_df$se_bias))
  dodge <- ggplot2::position_dodge(width = 0.04)
  p_support <- ggplot(
    support_plot_df,
    aes(
      x = .data$param_value,
      y = .data$mean_bias,
      color = .data[[lbl_col]],
      shape = .data[[lbl_col]],
      group = .data[[lbl_col]]
    )
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 0.7, alpha = 0.8) +
    {if (has_se) {
      geom_errorbar(
        aes(ymin = .data$bias_lo, ymax = .data$bias_hi),
        width = 0.03,
        linewidth = 0.35,
        alpha = 0.7,
        position = dodge
      )
    }} +
    geom_point(size = 1.8, alpha = 0.85, position = dodge) +
    facet_wrap(~ panel, scales = "free_x", nrow = 1) +
    scale_color_manual(values = c(naive = "#E07A5F", SEM = "#3D405B")) +
    scale_shape_manual(values = c(naive = 16, SEM = 17)) +
    labs(
      x = NULL,
      y = "DTAITE bias (estimate - truth)",
      color = "Labelling",
      shape = "Labelling",
      title = title,
      subtitle = subtitle
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  save_plot_pair(p_support, stem, width = 9.5, height = 4.6)
}

make_spatiotemporal_kernel_heatmap <- function(df, title, stem, subtitle = NULL) {
  lbl_col <- safe_label_col(df)
  if (is.null(lbl_col)) return(NULL)
  needed <- c("sim_kernel", "fit_kernel", "sim_spatial_kernel", "fit_spatial_kernel")
  if (!all(needed %in% names(df))) return(NULL)
  heatmap_labelling_levels <- c("naive", "SEM")
  out <- df %>%
    filter(
      .data$scenario_family %in% c(
        "spatiotemporal_kernel_mismatch",
        "kernel_mismatch",
        "spatial_kernel_mismatch"
      )
    ) %>%
    filter_support_plot_contrast()
  if (nrow(out) < 1L) return(NULL)
  if (!"mean_bias" %in% names(out)) out$mean_bias <- NA_real_
  support_plot_df <- out %>%
    mutate(
      mean_bias = ifelse(
        is.finite(.data$mean_bias),
        .data$mean_bias,
        ifelse(
          is.finite(.data$mean_psi_estimate) & is.finite(.data$mean_psi_truth),
          .data$mean_psi_estimate - .data$mean_psi_truth,
          NA_real_
        )
      ),
      !!lbl_col := normalize_labelling_method(.data[[lbl_col]]),
      truth_label = st_kernel_combo_label(.data$sim_kernel, .data$sim_spatial_kernel),
      fit_label = st_kernel_combo_label(.data$fit_kernel, .data$fit_spatial_kernel)
    ) %>%
    filter(.data[[lbl_col]] %in% heatmap_labelling_levels, is.finite(.data$mean_bias)) %>%
    group_by(.data[[lbl_col]], .data$truth_label, .data$fit_label) %>%
    summarize(mean_bias = mean(.data$mean_bias, na.rm = TRUE), .groups = "drop")
  if (nrow(support_plot_df) < 1L) return(NULL)

  # Complete the 4x4 grid so missing cells still appear (e.g. partial replot).
  methods <- unique(as.character(support_plot_df[[lbl_col]]))
  methods <- methods[methods %in% heatmap_labelling_levels]
  if (length(methods) < 1L) methods <- heatmap_labelling_levels
  grid <- expand.grid(
    truth_label = st_kernel_combo_levels,
    fit_label = st_kernel_combo_levels,
    method_tmp = methods,
    stringsAsFactors = FALSE
  )
  names(grid)[names(grid) == "method_tmp"] <- lbl_col
  support_plot_df <- grid %>%
    left_join(support_plot_df, by = c("truth_label", "fit_label", lbl_col))

  support_plot_df$truth_label <- factor(support_plot_df$truth_label, levels = rev(st_kernel_combo_levels))
  support_plot_df$fit_label <- factor(support_plot_df$fit_label, levels = st_kernel_combo_levels)
  support_plot_df[[lbl_col]] <- factor(
    as.character(support_plot_df[[lbl_col]]),
    levels = heatmap_labelling_levels
  )

  lim <- max(abs(support_plot_df$mean_bias), na.rm = TRUE)
  if (!is.finite(lim) || lim <= 0) lim <- 1
  p_heat <- ggplot(
    support_plot_df,
    aes(x = .data$fit_label, y = .data$truth_label, fill = .data$mean_bias)
  ) +
    geom_tile(color = "grey85", linewidth = 0.3) +
    facet_wrap(stats::as.formula(paste("~", lbl_col)), nrow = 1) +
    scale_fill_gradient2(
      low = "#3D405B",
      mid = "white",
      high = "#E07A5F",
      midpoint = 0,
      limits = c(-lim, lim),
      na.value = "grey95",
      name = "DTAITE bias"
    ) +
    labs(
      x = "Fitted kernels (temporal, spatial)",
      y = "True kernels (temporal, spatial)",
      title = title,
      subtitle = subtitle
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank(),
      legend.position = "bottom"
    )
  save_plot_pair(p_heat, stem, width = 8.5, height = 4.8)
}

if (!is.null(support_summary) && nrow(support_summary) > 0 &&
    all(c("contrast_family", "mean_abs_error") %in% names(support_summary))) {
  lbl_col <- safe_label_col(support_summary)
  if (!is.null(lbl_col)) {
    plot_files$parameter_sweep_support <- make_parameter_sweep_support_plot(
      support_summary,
      title = "K-separation and signal-to-noise: all-or-nothing DTAITE bias",
      stem = "robustness_parameter_sweep_support_contrasts",
      subtitle = "Naive vs SEM only. Left: bias vs |K_1 - K_0|; right: bias vs mu scale at (K_0,K_1)=(0.8,0.2). Error bars are +/- 1 SE across replications."
    )
    plot_files$spatiotemporal_kernel_heatmap <- make_spatiotemporal_kernel_heatmap(
      support_summary,
      title = "Spatiotemporal kernel misspecification: all-or-nothing DTAITE bias",
      stem = "robustness_spatiotemporal_kernel_mismatch_heatmap",
      subtitle = "Rows = true (temporal, spatial) kernels; columns = fitted kernels. Cells are mean DTAITE bias at (K_0,K_1)=(0.8,0.2)."
    )
    plot_files$high_count_assignment_support <- make_high_count_support_plot(
      support_summary,
      title = "Pretreatment-informed assignment: all-or-nothing DTAITE bias",
      stem = "robustness_pretreatment_assignment_support_contrasts",
      subtitle = "Naive vs SEM DTAITE bias by assignment rule at (K_0,K_1)=(0.8,0.2); grey lines connect the two methods. Error bars are +/- 1 SE across replications."
    )
  }
}

if (!is.null(decay_showcase_summary) && nrow(decay_showcase_summary) > 0 &&
    all(c("d_mid", "mean_abs_delta", "decay_label") %in% names(decay_showcase_summary))) {
  decay_reps_used <- if (isTRUE(refresh_decay) || isTRUE(auto_refresh_decay_showcase)) decay_reps else 2000L
  p_decay <- make_decay_validation_plot(decay_showcase_summary, decay_reps_used)
  plot_files$decay <- save_plot_pair(p_decay, "robustness_decay_validation", width = 8.5, height = 5.2)
  p_decay_cumulative <- make_decay_validation_cumulative_plot(decay_showcase_summary, decay_reps_used)
  plot_files$decay_cumulative <- save_plot_pair(
    p_decay_cumulative,
    "robustness_decay_validation_cumulative",
    width = 8.5,
    height = 5.2
  )
}

if (!is.null(decay_temporal_showcase_summary) && nrow(decay_temporal_showcase_summary) > 0 &&
    all(c("mean_abs_delta", "decay_label") %in% names(decay_temporal_showcase_summary)) &&
    any(c("t_mid", "d_mid") %in% names(decay_temporal_showcase_summary))) {
  decay_reps_used <- if (isTRUE(refresh_decay) || isTRUE(auto_refresh_decay_showcase)) decay_reps else 2000L
  p_decay_temporal <- make_temporal_decay_validation_plot(decay_temporal_showcase_summary, decay_reps_used)
  plot_files$decay_temporal <- save_plot_pair(
    p_decay_temporal,
    "robustness_temporal_decay_validation",
    width = 8.5,
    height = 5.2
  )
  p_decay_temporal_cumulative <- make_temporal_decay_validation_cumulative_plot(
    decay_temporal_showcase_summary,
    decay_reps_used
  )
  plot_files$decay_temporal_cumulative <- save_plot_pair(
    p_decay_temporal_cumulative,
    "robustness_temporal_decay_validation_cumulative",
    width = 8.5,
    height = 5.2
  )
}

if (!is.null(decay_label_flip_showcase_summary) && nrow(decay_label_flip_showcase_summary) > 0 &&
    all(c("d_mid", "mean_abs_delta", "decay_label") %in% names(decay_label_flip_showcase_summary))) {
  decay_reps_used <- if (isTRUE(refresh_decay) || isTRUE(auto_refresh_decay_showcase)) decay_reps else 2000L
  p_label_flip <- make_decay_validation_plot(
    decay_label_flip_showcase_summary,
    decay_reps_used,
    title = "Single-event label flip: spatial decay",
    subtitle = "One treated event's latent label flipped to control under CRN; area-normalized |Delta N|; clipped to untruncated annuli",
    xlab = "Distance from flipped event"
  )
  plot_files$decay_label_flip <- save_plot_pair(
    p_label_flip,
    "robustness_label_flip_decay_validation",
    width = 8.5,
    height = 5.2
  )
  p_label_flip_cum <- make_decay_validation_cumulative_plot(
    decay_label_flip_showcase_summary,
    decay_reps_used,
    title = "Single-event label flip: spatial decay (cumulative)",
    subtitle = "One treated event's latent label flipped to control; cumulative |Delta N| within radius of that event",
    xlab = "Distance from flipped event"
  )
  plot_files$decay_label_flip_cumulative <- save_plot_pair(
    p_label_flip_cum,
    "robustness_label_flip_decay_validation_cumulative",
    width = 8.5,
    height = 5.2
  )
}

if (!is.null(decay_label_flip_temporal_showcase_summary) &&
    nrow(decay_label_flip_temporal_showcase_summary) > 0 &&
    all(c("mean_abs_delta", "decay_label") %in% names(decay_label_flip_temporal_showcase_summary)) &&
    any(c("t_mid", "d_mid") %in% names(decay_label_flip_temporal_showcase_summary))) {
  decay_reps_used <- if (isTRUE(refresh_decay) || isTRUE(auto_refresh_decay_showcase)) decay_reps else 2000L
  p_label_flip_temporal <- make_temporal_decay_validation_plot(
    decay_label_flip_temporal_showcase_summary,
    decay_reps_used,
    title = "Single-event label flip: temporal decay",
    subtitle = "One treated event's latent label flipped to control under CRN; per-lag-bin |Delta N| after that event",
    xlab = "Lag since flipped event"
  )
  plot_files$decay_label_flip_temporal <- save_plot_pair(
    p_label_flip_temporal,
    "robustness_label_flip_temporal_decay_validation",
    width = 8.5,
    height = 5.2
  )
  p_label_flip_temporal_cum <- make_temporal_decay_validation_cumulative_plot(
    decay_label_flip_temporal_showcase_summary,
    decay_reps_used,
    title = "Single-event label flip: temporal decay (cumulative)",
    subtitle = "One treated event's latent label flipped to control; cumulative |Delta N| up to each lag after that event",
    xlab = "Lag since flipped event"
  )
  plot_files$decay_label_flip_temporal_cumulative <- save_plot_pair(
    p_label_flip_temporal_cum,
    "robustness_label_flip_temporal_decay_validation_cumulative",
    width = 8.5,
    height = 5.2
  )
}

latex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x
}

ROBUSTNESS_FIG_STEMS <- list(
  all_scenarios_label = "robustness_all_scenarios_label_recovery",
  parameter_sweep_support = "robustness_parameter_sweep_support_contrasts",
  spatiotemporal_kernel_heatmap = "robustness_spatiotemporal_kernel_mismatch_heatmap",
  high_count_support = "robustness_pretreatment_assignment_support_contrasts",
  effect_modification = "robustness_effect_modification",
  geometry_transport = "robustness_geometry_transport",
  decay_spatial = "robustness_decay_validation",
  decay_spatial_cumulative = "robustness_decay_validation_cumulative",
  decay_temporal = "robustness_temporal_decay_validation",
  decay_temporal_cumulative = "robustness_temporal_decay_validation_cumulative",
  decay_label_flip_spatial = "robustness_label_flip_decay_validation",
  decay_label_flip_spatial_cumulative = "robustness_label_flip_decay_validation_cumulative",
  decay_label_flip_temporal = "robustness_label_flip_temporal_decay_validation",
  decay_label_flip_temporal_cumulative = "robustness_label_flip_temporal_decay_validation_cumulative"
)

prune_stale_robustness_figures <- function(fig_dir, keep_stems) {
  keep_files <- unique(c(
    unlist(keep_stems, use.names = FALSE),
    "simulation_robustness_appendix"
  ))
  keep_paths <- as.character(unlist(lapply(
    keep_files,
    function(stem) file.path(fig_dir, paste0(stem, c(".pdf", ".png")))
  )))
  stale <- list.files(fig_dir, pattern = "\\.(pdf|png)$", full.names = TRUE)
  stale <- stale[!stale %in% keep_paths]
  if (length(stale) > 0L) {
    unlink(stale)
    message("[robustness] pruned stale figures: ", length(stale))
  }
  invisible(stale)
}

tex_fig <- function(stem, label, caption) {
  paste0(
    "\\begin{figure}[H]\n",
    "\\centering\n",
    "\\RobustnessIncludeFig{", stem, ".pdf}\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\end{figure}\n"
  )
}

tex_subsubsection <- function(title, label, intro_lines, fig_lines) {
  c(
    paste0("\\subsubsection{", title, "}\\label{", label, "}"),
    "",
    intro_lines,
    "",
    fig_lines,
    "\\FloatBarrier",
    ""
  )
}

scenario_counts <- manifest %>%
  count(.data$scenario_family, name = "n_scenarios")
scenario_table <- paste(
  apply(scenario_counts, 1, function(row) {
    paste0(latex_escape(row[["scenario_family"]]), " (", row[["n_scenarios"]], ")")
  }),
  collapse = ", "
)

tex_lines <- c(
  "% Auto-generated by inst/sim_study/sim_study_robustness.R",
  "%",
  "% Single deliverable for Overleaf and local PDF preview.",
  "%",
  "% Overleaf:",
  "%   1. Copy this file and figures/*.pdf to plots/sim_study/robustness/",
  "%   2. In revision.tex, inside \\section{Additional simulation results}:",
  "%        \\input{plots/sim_study/robustness/simulation_robustness_appendix.tex}",
  "%",
  "% Local preview (from generated/robustness/):",
  "%   pdflatex -jobname=robustness \"\\\\def\\\\robustnessstandalone{}\\\\input{simulation_robustness_appendix.tex}\"",
  "%",
  paste0(
    "% Run: ", latex_escape(basename(paste0(output_basename, "_summary.rds"))),
    " | scenario families: ", scenario_table
  ),
  "% Missing figure PDFs render as placeholders until the production run completes.",
  "",
  "\\ifdefined\\robustnessstandalone",
  "\\documentclass[11pt]{article}",
  "\\usepackage[margin=1in]{geometry}",
  "\\usepackage{graphicx}",
  "\\usepackage{amsmath,amssymb}",
  "\\usepackage{booktabs}",
  "\\usepackage{longtable}",
  "\\usepackage{float}",
  "\\usepackage{placeins}",
  "\\usepackage{hyperref}",
  "\\usepackage{cleveref}",
  "\\begin{document}",
  "\\appendix",
  "\\section{Additional simulation results}\\label{app:additional_simulation}",
  "\\fi",
  "",
  "% Overleaf parent document needs \\usepackage{float,placeins} in the preamble.",
  "",
  "\\subsection{Alternative simulation specifications}\\label{app:simulation_robustness}",
  "",
  "The main simulation in \\cref{sec:simulation_study} uses one baseline parameterisation---control $K_0=0.8$, treated $K_1=0.2$---and random treatment assignment. This subsection reports a structured robustness suite designed around the finite-sample assumptions in \\cref{sec:finite-sample-theory} and the identification discussion in \\cref{sec:identification}. Unless stated otherwise, all scenarios use the same window $[0,100]\\times[0,100]$, treatment time $t^\\star=10$, end time $T=110$, a $10\\times10$ cell tessellation, exponential temporal triggering with $(\\alpha,\\beta)=(0.01,10)$, and exponential models at fit time. Treatment lowers the branching ratio relative to control and is assigned to $50\\%$ of cells.",
  "",
  "For each scenario we simulate an ensemble of independent replications and report two classes of outcome. First, \\emph{label recovery}: raw accuracy of the naive and SEM labellings of post-$t^\\star$ events (as in the main simulation study). Second, \\emph{off-support allocation contrasts}: mean percentage error in plug-in expected-count contrasts for (i) a single treated cell, (ii) a random $50\\%$ treated regime, and (iii) all-treated, each versus all-control. All contrast estimates fix control parameters at truth and use fitted treated parameters from each labelling, matching \\cref{fig:all_nothing}. These contrasts probe extrapolation across treatment intensity (SEM fits with explosive parameter estimates are excluded throughout).",
  "",
  "The branching-ratio grid, spatiotemporal kernel misspecification matrix, assignment designs, signal-to-noise scalings, and decay diagnostics below are generated automatically by \\texttt{sim\\_study\\_robustness.R} in the \\texttt{PPDisentangle} package. Figures are reported in \\cref{app:robustness-parameter-sweep,app:robustness-effect-modification,app:robustness-geometry-transport,app:robustness-kernel,app:robustness-high-count,app:robustness-decay-cell,app:robustness-decay-label}.",
  "",
  "\\makeatletter",
  "\\ifdefined\\robustnessstandalone",
  "  \\def\\robustness@figdir{figures/}",
  "  \\def\\robustness@tabledir{tables/}",
  "\\else",
  "  \\providecommand{\\robustness@figdir}{plots/sim_study/robustness/}",
  "  \\providecommand{\\robustness@tabledir}{plots/sim_study/robustness/tables/}",
  "\\fi",
  "\\providecommand{\\RobustnessIncludeFig}[1]{%",
  "  \\IfFileExists{\\robustness@figdir #1}{%",
  "    \\includegraphics[width=\\textwidth]{\\robustness@figdir #1}%",
  "  }{%",
  "    \\fbox{\\begin{minipage}{0.95\\textwidth}%",
  "      \\vspace{0.12\\textheight}%",
  "      \\centering",
  "      {\\large\\itshape Placeholder figure}\\\\[0.75em]",
  "      {\\begingroup\\catcode`\\_=12\\relax\\ttfamily #1\\endgroup}\\\\[0.75em]",
  "      {\\footnotesize Awaiting robustness simulation output}%",
  "      \\vspace{0.12\\textheight}%",
  "    \\end{minipage}}%",
  "  }%",
  "}",
  "\\providecommand{\\RobustnessInputTable}[1]{%",
  "  \\IfFileExists{\\robustness@tabledir #1}{%",
  "    \\input{\\robustness@tabledir #1}%",
  "  }{}%",
  "}",
  "\\makeatother",
  "",
  "\\Cref{fig:robustness-all-scenarios-label} summarises mean raw labelling accuracy for naive and SEM across every fitted robustness scenario with a DTAITE contrast (K-separation, signal-to-noise, kernel misspecification, pretreatment assignment, binary-covariate effect modification, and coarseness / geometry transport).",
  "",
  tex_fig(
    ROBUSTNESS_FIG_STEMS$all_scenarios_label,
    "fig:robustness-all-scenarios-label",
    "Mean raw labelling accuracy for SEM versus naive across all robustness scenarios with DTAITE contrasts. Each point is one fitted scenario; colour marks the scenario family (including covariate-effect and coarseness studies). The dashed line is equal accuracy (points above favour SEM)."
  ),
  "\\FloatBarrier",
  "",
  tex_subsubsection(
    "K-separation and signal-to-noise",
    "app:robustness-parameter-sweep",
    c(
      "The main study in \\cref{sec:simulation_study} fixes a large separation between control and treated branching ratios, with treatment lowering the branching ratio from $K_0=0.8$ to $K_1=0.2$. Here we stress-test performance as that separation shrinks. Throughout, control and treated processes share $(\\alpha,\\beta)=(0.01,10)$ and exponential temporal triggering; treatment is assigned independently to $50\\%$ of cells. The control branching ratio is fixed at $K_0=0.8$ and the treated ratio varies over $K_1\\in\\{0.1,0.2,0.3,0.4,0.5,0.6,0.7\\}$, always below control. For each pair $(K_0,K_1)$ the background rate $\\mu$ is recalibrated so the expected post-treatment point count remains fixed; small $K_0-K_1$ therefore isolates weak component separation at comparable event abundance.",
      sprintf(
        "Finite-sample behaviour also depends on how many post-treatment events are available relative to the ambiguity created by spillover. Holding $(K_0,K_1)=(0.8,0.2)$ and the random assignment design fixed, we vary a background-rate multiplier $\\mu_{\\mathrm{scale}}\\in%s$ while calibrating $\\mu$ from a common target abundance; realized point counts therefore vary across the grid. Lower $\\mu_{\\mathrm{scale}}$ increases the relative contribution of triggered offspring to the observed pattern; higher values push the process toward a Poisson-like regime with weaker history dependence.",
        format_mu_scale_grid(mu_scales, for_tex = TRUE)
      ),
      "\\Cref{fig:robustness-parameter-sweep-support} reports all-or-nothing DTAITE bias for both grids (naive and SEM). All contrast estimates use true control parameters with fitted treated parameters (as in \\cref{fig:all_nothing}). Labelling accuracy for these scenarios is summarised with the other robustness checks in \\cref{fig:robustness-all-scenarios-label}."
    ),
    c(
      tex_fig(
        ROBUSTNESS_FIG_STEMS$parameter_sweep_support,
        "fig:robustness-parameter-sweep-support",
        "All-or-nothing DTAITE bias ($\\widehat\\psi-\\psi$) for naive and SEM across two parameter sweeps. Left: K-separation vs $|K_1-K_0|$. Right: signal-to-noise vs $\\mu$ scale at $(K_0,K_1)=(0.8,0.2)$. Error bars are $\\pm 1$ SE across replications."
      )
    )
  ),
  tex_subsubsection(
    "Binary-covariate treatment-effect modification",
    "app:robustness-effect-modification",
    c(
      "We next consider a minimal structured extension of the homogeneous treatment process. On the fixed $10\\times10$ grid let $X_{r,c}=1$ when $r+c$ is even and $X_{r,c}=-1$ otherwise. A single recorded allocation $z_{\\mathrm{obs}}$ treats exactly $25$ cells in each covariate stratum, and both $X$ and $z_{\\mathrm{obs}}$ are held fixed across replications, so the experiment isolates point-process and estimation variability rather than redrawing the design.",
      "Only the treatment-component source rate changes:",
      "\\[",
      "m_h(X)=\\frac{\\exp(hX)}{\\cosh(h)},\\qquad",
      "\\lambda_1^z(t,x)=\\mu_1 m_h(X_{c(x)})\\mathbf{1}\\{z_{c(x)}=1\\}",
      "+K_1\\sum_{\\substack{i:t_i<t\\\\r_i=1}}g_t(t-t_i)g_s(x-x_i).",
      "\\]",
      "The balanced strata imply $\\{m_h(1)+m_h(-1)\\}/2=1$, while the stratum source-rate ratio is $\\exp(2h)$; exact $X$-balance of $z_{\\mathrm{obs}}$ likewise keeps the treated-area mean of $m_h$ equal to one. We use $h\\in\\{0,0.3,0.6\\}$ and estimate unrestricted $h\\in\\mathbb R$ jointly with the treated Hawkes parameters $(\\mu_1,\\alpha_1,\\beta_1,K_1)$, holding control parameters at their true values. The correctly specified heterogeneous model is fitted under naive and SEM labels; a nested homogeneous model with $h=0$ is fitted under a separate homogeneous SEM.",
      "Before the full experiment, a true-label pilot gates the study. For each $h$ it requires convergence rate at least $0.8$, absolute bias of $\\widehat h$ at most $0.25$, and RMSE at most $0.4$, and it checks profile likelihoods, Hessian Wald coverage for $h$, correlation of $\\widehat h$ with $\\widehat\\mu_1$, stratum-specific event counts and optimizer diagnostics.",
      "We report three DTAITE contrasts versus all-control: the all-or-nothing contrast $\\Psi_{\\mathrm{global}}=\\Lambda^{1}(D)-\\Lambda^{0}(D)$, and two fixed single-cell interventions from the all-control baseline---one $X=+1$ cell and one $X=-1$ cell. Within every outer replication, truths and fitted counterfactuals condition on that replication's pre-treatment history and use common random numbers across regimes. Causal summaries report bias, RMSE, MAE, empirical variation, forward-simulation Monte Carlo error, Monte Carlo uncertainty of the summaries, and fitting failures. Label summaries use balanced accuracy and component-specific recall, including within each $X$ stratum."
    ),
    c(
      tex_fig(
        ROBUSTNESS_FIG_STEMS$effect_modification,
        "fig:robustness-effect-modification",
        "Binary-covariate effect-modification study. Panel A shows the fixed covariate and observed allocation, with the two single-cell flip locations outlined. Panel B reports estimated DTAITE bias versus $h$ for all-or-nothing and the single $X=+1$ and $X=-1$ cell flips (naive and SEM)."
      )
    )
  ),
  tex_subsubsection(
    "Transport across balanced allocation geometries",
    "app:robustness-geometry-transport",
    c(
      "Starting at the checkerboard $z^{\\mathrm{cb}}_{r,c}=\\mathbf{1}\\{r+c\\text{ is even}\\}$, a fixed sequence reverses $m\\in\\{0,5,10,15,20,25\\}$ discrepant reflected pairs $(r,c)$ and $(r,11-c)$, ending at the left--right block $z^{\\mathrm{block}}_{r,c}=\\mathbf{1}\\{c\\leq5\\}$ while retaining exactly $50\\%$ treatment.",
      "\\paragraph{Rook adjacency and coarseness.}",
      "Index the $10\\times10$ cells by $(r,c)$ with $r,c\\in\\{1,\\ldots,10\\}$. Two distinct cells $j=(r,c)$ and $\\ell=(r',c')$ are rook-adjacent if they share a side, $|r-r'|+|c-c'|=1$, and are not diagonal neighbours. Let $E$ be the undirected edge set of all such pairs. On a $10\\times10$ grid there are $90$ horizontal and $90$ vertical edges, so $|E|=180$.",
      "For an allocation $z\\in\\{0,1\\}^{100}$, the cut count is $D(z)=\\sum_{(j,\\ell)\\in E}\\mathbf{1}\\{z_j\\neq z_\\ell\\}$. The checkerboard maximises the cut ($D(z^{\\mathrm{cb}})=180$), while the left--right block has only the ten centreline edges discordant ($D(z^{\\mathrm{block}})=10$). Normalised coarseness is $C(z)=\\{180-D(z)\\}/170$, so $C(z^{\\mathrm{cb}})=0$ and $C(z^{\\mathrm{block}})=1$. Larger $C(z)$ means greater spatial aggregation. Figures use realised $C(z^{(m)})$ rather than $m/25$, and mark $C(z_{\\mathrm{obs}})$.",
      "The five-cell focal band in column 5 remains treated throughout. Each naive or SEM model is fitted once under $z_{\\mathrm{obs}}$ and evaluated under all six targets. We report estimated bias in the DTAITE for each target allocation."
    ),
    c(
      tex_fig(
        ROBUSTNESS_FIG_STEMS$geometry_transport,
        "fig:robustness-geometry-transport",
        "Fixed balanced-allocation path and estimated DTAITE bias against realized coarseness for naive and SEM. The vertical line marks the observed allocation."
      )
    )
  ),
  tex_subsubsection(
    "Mis-specified spatiotemporal kernels",
    "app:robustness-kernel",
    c(
      "Assumption~\\ref{ass:G3} and the identification discussion in \\cref{sec:identification} treat the triggering kernel as part of the causal specification: misspecification changes both fit and the transmitted spillover geometry. We therefore cross a $2\\times 2$ truth grid of temporal and spatial kernels (exponential or mean-matched power-law in each margin) with the same $2\\times 2$ fitted grid, yielding a $4\\times 4$ misspecification matrix at $(K_0,K_1)=(0.8,0.2)$ with random $50\\%$ treatment assignment.",
      "Power-law factors are mean-matched to their exponential counterparts so that misspecification isolates heavier tails rather than a change in expected triggering distance or lag. The resulting matrix includes correctly specified fits, temporal-only and spatial-only misspecification, and joint misspecification of both margins.",
      "\\Cref{fig:robustness-kernel-heatmap} reports all-or-nothing DTAITE bias across this matrix for naive and SEM; labelling accuracy appears in \\cref{fig:robustness-all-scenarios-label}."
    ),
    c(
      tex_fig(ROBUSTNESS_FIG_STEMS$spatiotemporal_kernel_heatmap, "fig:robustness-kernel-heatmap",
              "Spatiotemporal kernel misspecification matrix for all-or-nothing DTAITE bias ($\\widehat\\psi-\\psi$). Rows are true (temporal, spatial) kernels; columns are fitted kernels. Facets are naive and SEM. Empty cells indicate scenarios not yet available in the plotted run.")
    )
  ),
  tex_subsubsection(
    "Pre-treatment data informed treatment allocation",
    "app:robustness-high-count",
    c(
      "The baseline design randomises treated cells. Many applications instead assign treatment using pre-treatment activity. We therefore simulate a reference pre-treatment catalogue from the control Hawkes law on $[0,t^\\star]$, count points per cell, and compare five assignment rules that each treat $50\\%$ of cells at $K_0=0.8$, $K_1=0.2$:",
      "\\begin{enumerate}\\item \\emph{Highest count}: treat the cells with the largest reference counts (anti-aligned intensity targeting).",
      "\\item \\emph{Lowest count}: treat the quietest cells (positively aligned intensity targeting).",
      "\\item \\emph{Count propensity}: draw a soft logistic propensity on $\\log(1+N_j^{\\mathrm{pre}})$, calibrated so the expected treated count is $50\\%$, then sample exactly half the cells without replacement.",
      "\\item \\emph{Contiguous AOI}: grow a rook-contiguous treated region from the hottest seed cell by repeatedly adjoining the highest-count frontier neighbour (Oklahoma-style geometry).",
      "\\item \\emph{Voronoi random}: take the first $100$ pretreatment events as Dirichlet seeds (yielding up to $100$ Voronoi areas), randomly treat $50\\%$ of those tiles, and assign each analysis-grid cell the treatment of the tile containing its centroid.",
      "\\end{enumerate}",
      "Highest-count assignment is especially demanding for naive location-based labelling, because cell location and latent process type are deliberately misaligned. Contiguous AOI assignment additionally concentrates spillover along a spatial boundary. Voronoi-random assignment uses pretreatment geometry to define irregular contiguous units while keeping treatment random at the tile level. \\Cref{fig:robustness-high-count-support} summarises all-or-nothing DTAITE bias across the five rules; labelling accuracy appears in \\cref{fig:robustness-all-scenarios-label}."
    ),
    c(
      tex_fig(ROBUSTNESS_FIG_STEMS$high_count_support, "fig:robustness-high-count-support",
              "Pretreatment-informed assignment all-or-nothing DTAITE bias ($\\widehat\\psi-\\psi$) by assignment rule for naive and SEM. Grey lines connect the two methods within each rule. Error bars are $\\pm 1$ SE across replications.")
    )
  ),
  tex_subsubsection(
    "Single-cell allocation flip under alternative kernels",
    "app:robustness-decay-cell",
    c(
      "Supplement Proposition~\\ref{prop:decay-snippet} (and Corollary~\\ref{cor:ani-snippet}) shows that the influence of a single allocation coordinate flip $z_j$ on $\\Lambda^z_\\theta(B)$ decays geometrically for finite-range kernels, exponentially for light-tailed kernels, and polynomially for regularly varying spatial kernels. We illustrate these allocation-influence decay laws by forward simulation under true Hawkes parameters.",
      "In each replicate we flip the same fixed central \\emph{cell} from treated to control (one allocation coordinate), simulate the two allocations with common random numbers from a shared pre-treatment history, and report mean $|\\Delta N|$ per unit area in distance annuli around that cell (spatial) or mean $|\\Delta N|$ in lag bins after $t^\\star$ (temporal). Companion cumulative panels sum the same annulus or lag-bin masses out to each distance or lag. Spatial curves are restricted to untruncated annuli. Both axes fix $(K_0,K_1)=(0.8,0.2)$ and compare exponential, mean-matched power-law, and larger-mean power-law kernels. Temporal kernels use the main-study mean lag ($\\beta=10$). \\Cref{fig:robustness-decay,fig:robustness-decay-cumulative,fig:robustness-decay-temporal,fig:robustness-decay-temporal-cumulative}."
    ),
    c(
      tex_fig(ROBUSTNESS_FIG_STEMS$decay_spatial, "fig:robustness-decay",
              "Single-cell allocation flip: spatial decay at $(K_0,K_1)=(0.8,0.2)$. One fixed central treated \\emph{cell} is flipped to control under CRN; curves show mean $|\\Delta N|$ per unit area in each distance annulus around that cell on a $\\log_{10}$ $y$-scale, restricted to untruncated annuli."),
      tex_fig(ROBUSTNESS_FIG_STEMS$decay_spatial_cumulative, "fig:robustness-decay-cumulative",
              "Single-cell allocation flip: cumulative spatial decay. Curves show cumulative mean $|\\Delta N|$ within each distance radius of the flipped cell."),
      tex_fig(ROBUSTNESS_FIG_STEMS$decay_temporal, "fig:robustness-decay-temporal",
              "Single-cell allocation flip: temporal decay at $(K_0,K_1)=(0.8,0.2)$. One fixed central treated cell is flipped to control under CRN; curves show mean $|\\Delta N|$ per lag bin after $t^\\star$ on a $\\log_{10}$ $y$-scale ($\\beta=10$)."),
      tex_fig(ROBUSTNESS_FIG_STEMS$decay_temporal_cumulative, "fig:robustness-decay-temporal-cumulative",
              "Single-cell allocation flip: cumulative temporal decay. Curves show cumulative mean $|\\Delta N|$ up to each lag after $t^\\star$.")
    )
  ),
  tex_subsubsection(
    "Single-event label flip under alternative kernels",
    "app:robustness-decay-label",
    c(
      "A complementary diagnostic flips the latent process label of a single post-treatment \\emph{event}, not an allocation cell. In each replicate we place one treated-labelled seed event in the fixed central treated cell at $t^\\star$, then generate offspring under treated versus control Hawkes parameters with common random numbers (full cascade). Mean $|\\Delta N|$ is reported in distance annuli around that event or in lag bins after the event time, with companion cumulative panels. Spatial curves are again restricted to untruncated annuli. The same three-kernel comparisons and $(K_0,K_1)=(0.8,0.2)$ design are used. \\Cref{fig:robustness-label-flip-decay,fig:robustness-label-flip-decay-cumulative,fig:robustness-label-flip-decay-temporal,fig:robustness-label-flip-decay-temporal-cumulative}."
    ),
    c(
      tex_fig(ROBUSTNESS_FIG_STEMS$decay_label_flip_spatial, "fig:robustness-label-flip-decay",
              "Single-event label flip: spatial decay at $(K_0,K_1)=(0.8,0.2)$. One treated seed \\emph{event} in the central cell has its latent label flipped to control under CRN; curves show mean $|\\Delta N|$ per unit area in each distance annulus around that event on a $\\log_{10}$ $y$-scale, restricted to untruncated annuli."),
      tex_fig(ROBUSTNESS_FIG_STEMS$decay_label_flip_spatial_cumulative, "fig:robustness-label-flip-decay-cumulative",
              "Single-event label flip: cumulative spatial decay. Curves show cumulative mean $|\\Delta N|$ within each distance radius of the flipped event."),
      tex_fig(ROBUSTNESS_FIG_STEMS$decay_label_flip_temporal, "fig:robustness-label-flip-decay-temporal",
              "Single-event label flip: temporal decay at $(K_0,K_1)=(0.8,0.2)$. One treated seed event has its latent label flipped to control under CRN; curves show mean $|\\Delta N|$ per lag bin after that event on a $\\log_{10}$ $y$-scale ($\\beta=10$)."),
      tex_fig(ROBUSTNESS_FIG_STEMS$decay_label_flip_temporal_cumulative, "fig:robustness-label-flip-decay-temporal-cumulative",
              "Single-event label flip: cumulative temporal decay. Curves show cumulative mean $|\\Delta N|$ up to each lag after the flipped event.")
    )
  ),
  "",
  "\\ifdefined\\robustnessstandalone",
  "\\end{document}",
  "\\fi"
)

tex_path <- file.path(robustness_dir, "simulation_robustness_appendix.tex")
writeLines(tex_lines[nzchar(tex_lines)], tex_path)
prune_stale_robustness_figures(fig_dir, ROBUSTNESS_FIG_STEMS)

saveRDS(
  list(
    output_basename = output_basename,
    manifest = manifest,
    label_summary = label_summary,
    support_summary = support_summary,
    ate_summary = ate_summary,
    decay_summary = decay_summary,
    plot_files = plot_files,
    tex_path = tex_path
  ),
  file.path(out_dir, paste0(output_basename, "_summary.rds"))
)

message("[robustness] wrote manifest: ", manifest_path)
message("[robustness] wrote figures: ", fig_dir)
message("[robustness] wrote tex: ", tex_path)
