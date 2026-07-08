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
DEFAULT_MU_SCALES <- c(0.5, 1, 1.5)
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
out_dir <- PPDisentangle::pp_output_path("sim_study", repo_root = repo_dir)
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
run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_basename <- if (!is.null(replot_basename)) {
  replot_basename
} else {
  paste0("robustness_", run_stamp)
}
if (is.null(replot_basename) && nzchar(Sys.getenv("SLURM_JOB_ID", ""))) {
  output_basename <- paste0("robustness_", Sys.getenv("SLURM_JOB_ID"))
}

k_separation <- data.frame(
  scenario_family = "k_separation",
  control_k = k_anchor,
  treated_k = k_values,
  mu_scale = 1,
  sim_kernel = "exponential",
  fit_kernel = "exponential",
  treatment_assignment = "random",
  stringsAsFactors = FALSE
)

high_count_assignment <- data.frame(
  scenario_family = "high_count_assignment",
  control_k = ROBUSTNESS_CONTROL_K,
  treated_k = ROBUSTNESS_TREATED_K_DEFAULT,
  mu_scale = 1,
  sim_kernel = "exponential",
  fit_kernel = "exponential",
  treatment_assignment = "highest_count_50pct",
  stringsAsFactors = FALSE
)

snr_scale <- data.frame(
  scenario_family = "snr_scale",
  control_k = ROBUSTNESS_CONTROL_K,
  treated_k = ROBUSTNESS_TREATED_K_DEFAULT,
  mu_scale = mu_scales,
  sim_kernel = "exponential",
  fit_kernel = "exponential",
  treatment_assignment = "random",
  stringsAsFactors = FALSE
)

kernel_mismatch <- expand.grid(
  scenario_family = "kernel_mismatch",
  control_k = ROBUSTNESS_CONTROL_K,
  treated_k = ROBUSTNESS_TREATED_K_DEFAULT,
  mu_scale = 1,
  sim_kernel = c("exponential", "power_law"),
  fit_kernel = c("exponential", "power_law"),
  treatment_assignment = "random",
  stringsAsFactors = FALSE
)

scenarios <- bind_rows(k_separation, high_count_assignment, snr_scale, kernel_mismatch) %>%
  mutate(
    k_delta = abs(.data$treated_k - .data$control_k),
    scenario_id = paste0(
      .data$scenario_family,
      "_kc", format_num_tag(.data$control_k),
      "_kt", format_num_tag(.data$treated_k),
      "_mu", format_num_tag(.data$mu_scale),
      "_sim", .data$sim_kernel,
      "_fit", .data$fit_kernel,
      "_assign", .data$treatment_assignment
    )
  ) %>%
  distinct(.data$scenario_id, .keep_all = TRUE)

if (!identical(tolower(scenario_filter), "all")) {
  keep_families <- unlist(strsplit(scenario_filter, "[,;[:space:]]+", perl = TRUE), use.names = FALSE)
  keep_families <- keep_families[nzchar(keep_families)]
  scenarios <- scenarios %>% filter(.data$scenario_family %in% keep_families)
}
if (nrow(scenarios) < 1L) stop("No robustness scenarios selected.")
scenario_workers <- max(1L, min(as.integer(scenario_workers), nrow(scenarios)))
scenario_cpus <- max(1L, floor(allocated_cpus / scenario_workers))
message(sprintf(
  "[robustness] scenario workers=%d | allocated_cpus=%d | cpus_per_scenario=%d | sims_per_scenario=%d",
  scenario_workers, allocated_cpus, scenario_cpus, pp_sims
))

robustness_dir <- file.path(out_dir, "generated", "robustness")
fig_dir <- file.path(robustness_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(robustness_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(script_dir, "decay_validation_utils.R"), local = FALSE)

run_one <- function(row_id) {
  sc <- scenarios[row_id, , drop = FALSE]
  run_basename <- paste0(output_basename, "_", sc$scenario_id)
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
  rds_path <- file.path(out_dir, paste0(run_basename, ".rds"))
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
    treatment_assignment = sc$treatment_assignment,
    target_points = target_points,
    run_basename = run_basename,
    rds_path = rds_path,
    stringsAsFactors = FALSE
  )
}

manifest <- if (!is.null(replot_basename)) {
  manifest_path <- file.path(out_dir, paste0(output_basename, "_manifest.csv"))
  if (!file.exists(manifest_path)) {
    stop(sprintf("Replot manifest not found: %s", manifest_path))
  }
  message("[robustness] replot from existing summaries: ", output_basename)
  read.csv(manifest_path, stringsAsFactors = FALSE)
} else {
  scenario_ids <- seq_len(nrow(scenarios))
  if (scenario_workers > 1L) {
    bind_rows(parallel::mclapply(
      scenario_ids,
      run_one,
      mc.cores = scenario_workers,
      mc.preschedule = FALSE
    ))
  } else {
    bind_rows(lapply(scenario_ids, run_one))
  }
}

resolve_manifest_rds_path <- function(row, out_dir) {
  rds_path <- row$rds_path[[1]]
  if (file.exists(rds_path)) return(rds_path)
  if ("run_basename" %in% names(row) && nzchar(row$run_basename[[1]])) {
    local_path <- file.path(out_dir, paste0(row$run_basename[[1]], ".rds"))
    if (file.exists(local_path)) return(local_path)
  }
  if ("scenario_id" %in% names(row) && nzchar(row$scenario_id[[1]])) {
    alt <- list.files(out_dir, pattern = paste0(row$scenario_id[[1]], "\\.rds$"), full.names = TRUE)
    if (length(alt) > 0L) return(alt[[1L]])
  }
  rds_path
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
  data.frame(
    showcase_id = c("k02_exp", "k05_exp", "k07_exp"),
    decay_label = c(
      "K_c=0.8, K_t=0.2 (exp)",
      "K_c=0.8, K_t=0.5 (exp)",
      "K_c=0.8, K_t=0.7 (exp)"
    ),
    control_k = ROBUSTNESS_CONTROL_K,
    treated_k = c(0.2, 0.5, 0.7),
    sim_kernel = c("exponential", "exponential", "exponential"),
    stringsAsFactors = FALSE
  )
}

default_temporal_decay_showcase_specs <- function() {
  data.frame(
    showcase_id = c("k02_exp", "k05_exp", "k02_power"),
    decay_label = c(
      "K_c=0.8, K_t=0.2 (exp)",
      "K_c=0.8, K_t=0.5 (exp)",
      "K_c=0.8, K_t=0.2 (power-law)"
    ),
    control_k = ROBUSTNESS_CONTROL_K,
    treated_k = c(0.2, 0.5, 0.2),
    sim_kernel = c("exponential", "exponential", "power_law"),
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
  bind_rows(lapply(seq_len(nrow(specs)), function(i) {
    spec <- specs[i, , drop = FALSE]
    message(sprintf("[robustness] decay showcase %d/%d %s",
                    i, nrow(specs), spec$decay_label))
    refresh_decay_for_spec(
      rds_path = template_rds,
      control_k = spec$control_k,
      treated_k = spec$treated_k,
      sim_kernel = spec$sim_kernel,
      decay_reps = reps,
      showcase_id = spec$showcase_id,
      decay_label = spec$decay_label
    )
  }))
}

refresh_temporal_decay_showcase <- function(manifest_df, specs, reps, out_dir) {
  template_rds <- pick_decay_template_rds(manifest_df, out_dir)
  message(sprintf("[robustness] temporal decay showcase refresh (%d specs, %d reps) from %s",
                  nrow(specs), reps, basename(template_rds)))
  bind_rows(lapply(seq_len(nrow(specs)), function(i) {
    spec <- specs[i, , drop = FALSE]
    message(sprintf("[robustness] temporal decay showcase %d/%d %s",
                    i, nrow(specs), spec$decay_label))
    refresh_temporal_decay_for_spec(
      rds_path = template_rds,
      control_k = spec$control_k,
      treated_k = spec$treated_k,
      sim_kernel = spec$sim_kernel,
      decay_reps = reps,
      showcase_id = spec$showcase_id,
      decay_label = spec$decay_label
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
label_summary <- if (!is.null(replot_basename)) {
  read_summary_csv("_label_recovery_summary.csv")
} else {
  bind_rows(lapply(summary_rows, `[[`, "label")) %>%
    left_join(manifest, by = "scenario_id")
}
support_summary <- if (!is.null(replot_basename)) {
  read_summary_csv("_support_contrast_summary.csv")
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

if (is.null(replot_basename)) {
  write.csv(label_summary, file.path(out_dir, paste0(output_basename, "_label_recovery_summary.csv")), row.names = FALSE)
  write.csv(support_summary, file.path(out_dir, paste0(output_basename, "_support_contrast_summary.csv")), row.names = FALSE)
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
} else if (isTRUE(refresh_decay)) {
  if (!is.null(decay_showcase_summary) && nrow(decay_showcase_summary) > 0) {
    write.csv(decay_showcase_summary, file.path(out_dir, paste0(output_basename, "_decay_showcase_summary.csv")), row.names = FALSE)
  }
  if (!is.null(decay_temporal_showcase_summary) && nrow(decay_temporal_showcase_summary) > 0) {
    write.csv(decay_temporal_showcase_summary, file.path(out_dir, paste0(output_basename, "_decay_temporal_showcase_summary.csv")), row.names = FALSE)
  }
}

safe_metric_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) < 1L) return(NULL)
  hit[[1L]]
}

safe_label_col <- function(df) {
  safe_metric_col(df, c("method", "labelling"))
}

save_plot_pair <- function(plot, stem, width = 7.2, height = 4.6) {
  if (is.null(plot) || !inherits(plot, "ggplot")) return(NULL)
  pdf_path <- file.path(fig_dir, paste0(stem, ".pdf"))
  png_path <- file.path(fig_dir, paste0(stem, ".png"))
  ggplot2::ggsave(pdf_path, plot, width = width, height = height, units = "in")
  ggplot2::ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 300)
  list(stem = stem, pdf = pdf_path, png = png_path)
}

standard_labelling_levels <- c("SEM_adaptive", "SEM_full", "best", "naive", "oracle")

make_k_separation_line_plot <- function(df, x_col, y_col, lbl_col, title, subtitle, ylab,
                                        family = "k_separation") {
  plot_df <- df %>%
    filter(.data$scenario_family == family) %>%
    filter(.data[[lbl_col]] %in% standard_labelling_levels)
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

spatial_decay_showcase_label_levels <- function() {
  c(
    "K_c=0.8, K_t=0.2 (exp)",
    "K_c=0.8, K_t=0.5 (exp)",
    "K_c=0.8, K_t=0.7 (exp)"
  )
}

temporal_decay_showcase_label_levels <- function() {
  c(
    "K_c=0.8, K_t=0.2 (exp)",
    "K_c=0.8, K_t=0.5 (exp)",
    "K_c=0.8, K_t=0.2 (power-law)"
  )
}

prepare_decay_showcase_plot_df <- function(showcase_df, max_distance = 35) {
  if (is.null(showcase_df) || nrow(showcase_df) < 1L) return(NULL)
  label_levels <- spatial_decay_showcase_label_levels()
  out <- showcase_df %>%
    dplyr::filter(is.finite(.data$d_mid), is.finite(.data$mean_abs_delta)) %>%
    dplyr::mutate(
      decay_label = factor(.data$decay_label, levels = label_levels)
    )
  if (nrow(out) < 1L) return(NULL)
  positive <- out$mean_abs_delta[out$mean_abs_delta > 0]
  eps <- if (length(positive) > 0L) min(positive, na.rm = TRUE) / 2 else 1e-6
  out$mean_abs_delta_plot <- pmax(out$mean_abs_delta, eps)
  out %>% dplyr::filter(.data$d_mid <= max_distance)
}

make_decay_validation_plot <- function(showcase_df, decay_reps_used = NA_integer_) {
  decay_plot_df <- prepare_decay_showcase_plot_df(showcase_df)
  if (is.null(decay_plot_df) || nrow(decay_plot_df) < 1L) return(NULL)
  reps_txt <- if (is.finite(decay_reps_used)) as.character(decay_reps_used) else "stored"
  ggplot(decay_plot_df, aes(x = .data$d_mid, y = .data$mean_abs_delta_plot,
                            color = .data$decay_label,
                            group = .data$decay_label)) +
    geom_line(linewidth = 0.8, alpha = 0.95) +
    geom_point(size = 1.2, alpha = 0.9) +
    scale_y_log10() +
    labs(
      x = "Distance from flipped event",
      y = "Mean |Delta N| per annulus",
      color = "Scenario",
      title = "Forward-simulation spatial decay validation (true parameters)",
      subtitle = paste0(
        "Single latent label flip; three exponential-kernel K-separation scenarios (spatial kernel shared). ",
        "Decay reps per scenario: ", reps_txt, "."
      )
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

prepare_temporal_decay_showcase_plot_df <- function(showcase_df) {
  if (is.null(showcase_df) || nrow(showcase_df) < 1L) return(NULL)
  x_col <- if ("t_mid" %in% names(showcase_df)) "t_mid" else "d_mid"
  label_levels <- temporal_decay_showcase_label_levels()
  out <- showcase_df %>%
    dplyr::filter(is.finite(.data[[x_col]]), is.finite(.data$mean_abs_delta)) %>%
    dplyr::mutate(
      decay_label = factor(.data$decay_label, levels = label_levels)
    )
  if (nrow(out) < 1L) return(NULL)
  positive <- out$mean_abs_delta[out$mean_abs_delta > 0]
  eps <- if (length(positive) > 0L) min(positive, na.rm = TRUE) / 2 else 1e-6
  out$mean_abs_delta_plot <- pmax(out$mean_abs_delta, eps)
  out
}

make_temporal_decay_validation_plot <- function(showcase_df, decay_reps_used = NA_integer_) {
  decay_plot_df <- prepare_temporal_decay_showcase_plot_df(showcase_df)
  if (is.null(decay_plot_df) || nrow(decay_plot_df) < 1L) return(NULL)
  x_col <- if ("t_mid" %in% names(decay_plot_df)) "t_mid" else "d_mid"
  reps_txt <- if (is.finite(decay_reps_used)) as.character(decay_reps_used) else "stored"
  ggplot(decay_plot_df, aes(x = .data[[x_col]], y = .data$mean_abs_delta_plot,
                            color = .data$decay_label,
                            group = .data$decay_label)) +
    geom_line(linewidth = 0.8, alpha = 0.95) +
    geom_point(size = 1.2, alpha = 0.9) +
    scale_y_log10() +
    labs(
      x = "Lag since flipped event",
      y = "Mean |Delta N| per time bin",
      color = "Scenario",
      title = "Forward-simulation temporal decay validation (true parameters)",
      subtitle = paste0(
        "Single latent label flip on one event in a treated cell (control/treated K differ by scenario). ",
        "Offspring counted by lag since flip time. Decay reps per scenario: ",
        reps_txt, "."
      )
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

make_high_count_assignment_bar_plot <- function(df, y_col, lbl_col, title, subtitle, ylab) {
  plot_df <- df %>%
    filter(.data$scenario_family == "high_count_assignment") %>%
    filter(.data[[lbl_col]] %in% standard_labelling_levels)
  if (nrow(plot_df) < 1L) return(NULL)
  plot_df <- plot_df %>%
    mutate(
      scenario_label = sprintf(
        "K_c=%.1f, K_t=%.1f",
        .data$control_k,
        .data$treated_k
      )
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
    theme(legend.position = "bottom")
}

make_snr_scale_line_plot <- function(df, y_col, lbl_col, title, subtitle, ylab,
                                     family = "snr_scale") {
  plot_df <- df %>%
    filter(.data$scenario_family == family) %>%
    filter(.data[[lbl_col]] %in% standard_labelling_levels)
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
    filter(.data[[lbl_col]] %in% standard_labelling_levels) %>%
    mutate(kernel_pair = kernel_pair_label(.data$sim_kernel, .data$fit_kernel))
  if (nrow(plot_df) < 1L) return(NULL)
  pair_levels <- unique(plot_df$kernel_pair)
  plot_df$kernel_pair <- factor(plot_df$kernel_pair, levels = pair_levels)
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

plot_files <- list()
if (!is.null(label_summary) && nrow(label_summary) > 0) {
  lbl_col <- safe_label_col(label_summary)
  acc_col <- safe_metric_col(label_summary, c("mean_balanced_accuracy", "mean_accuracy", "accuracy"))
  if (!is.null(lbl_col) && !is.null(acc_col)) {
    ylab_acc <- gsub("_", " ", acc_col)
    plot_files$k_separation_label_recovery <- save_plot_pair(
      make_k_separation_line_plot(
        label_summary, "k_delta", acc_col, lbl_col,
        title = "K-separation sensitivity: label recovery",
        subtitle = "Control K fixed at 0.8; treated K in {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7}. Treatment lowers K; higher is better; small separation is the hard regime.",
        ylab = ylab_acc
      ),
      "robustness_k_separation_label_recovery"
    )
    plot_files$snr_scale_label_recovery <- save_plot_pair(
      make_snr_scale_line_plot(
        label_summary, acc_col, lbl_col,
        title = "Signal-to-noise scaling: label recovery",
        subtitle = sprintf(
          "K fixed at (0.8, 0.2); mu scale in %s.",
          format_mu_scale_grid(mu_scales)
        ),
        ylab = ylab_acc
      ),
      "robustness_snr_scale_label_recovery"
    )
    plot_files$high_count_assignment_label_recovery <- save_plot_pair(
      make_high_count_assignment_bar_plot(
        label_summary, acc_col, lbl_col,
        title = "High-count treatment assignment: label recovery",
        subtitle = "Treatment assigned to the 50% of cells with most points in a reference pre-treatment catalogue; K_c=0.8, K_t=0.2.",
        ylab = ylab_acc
      ),
      "robustness_high_count_assignment_label_recovery",
      width = 7.2, height = 4.8
    )
    plot_files$kernel_mismatch_label_recovery <- save_plot_pair(
      make_kernel_mismatch_bar_plot(
        label_summary, acc_col, lbl_col,
        title = "Kernel mismatch: label recovery",
        subtitle = "All four sim/fit kernel pairs at K=(0.8, 0.2).",
        ylab = ylab_acc
      ),
      "robustness_kernel_mismatch_label_recovery",
      width = 8.0, height = 4.8
    )
  }
}

if (!is.null(ate_summary) && nrow(ate_summary) > 0) {
  lbl_col <- safe_label_col(ate_summary)
  ate_col <- safe_metric_col(ate_summary, c("mean_all_nothing", "mean_all_nothing_true_control"))
  if (!is.null(lbl_col) && !is.null(ate_col) && "true_all_nothing_ATE" %in% names(ate_summary)) {
    ate_plot_base <- ate_summary %>%
      mutate(abs_all_nothing_error = abs(.data[[ate_col]] - .data$true_all_nothing_ATE)) %>%
      filter(is.finite(.data$abs_all_nothing_error))
    ate_k_separation_plot_base <- ate_plot_base %>%
      filter(.data[[lbl_col]] != "SEM_full")
    plot_files$k_separation_ate_error <- save_plot_pair(
      make_k_separation_line_plot(
        ate_k_separation_plot_base, "k_delta", "abs_all_nothing_error", lbl_col,
        title = "K-separation sensitivity: ATE error",
        subtitle = "Absolute all-or-nothing ATE error relative to the known simulation truth.",
        ylab = "Absolute all-or-nothing ATE error"
      ),
      "robustness_k_separation_ate_error"
    )
    plot_files$snr_scale_ate_error <- save_plot_pair(
      make_snr_scale_line_plot(
        ate_plot_base %>% filter(.data[[lbl_col]] != "SEM_full"),
        "abs_all_nothing_error", lbl_col,
        title = "Signal-to-noise scaling: ATE error",
        subtitle = "Absolute all-or-nothing ATE error at K=(0.8, 0.2) across mu scales.",
        ylab = "Absolute all-or-nothing ATE error"
      ),
      "robustness_snr_scale_ate_error"
    )
    plot_files$high_count_assignment_ate_error <- save_plot_pair(
      make_high_count_assignment_bar_plot(
        ate_plot_base %>% filter(.data[[lbl_col]] != "SEM_full"),
        "abs_all_nothing_error", lbl_col,
        title = "High-count treatment assignment: ATE error",
        subtitle = "Treatment assigned to the 50% of cells with most points in a reference pre-treatment catalogue; K_c=0.8, K_t=0.2.",
        ylab = "Absolute all-or-nothing ATE error"
      ),
      "robustness_high_count_assignment_ate_error",
      width = 7.2, height = 4.8
    )
    plot_files$kernel_mismatch_ate_error <- save_plot_pair(
      make_kernel_mismatch_bar_plot(
        ate_plot_base %>% filter(.data[[lbl_col]] != "SEM_full"),
        "abs_all_nothing_error", lbl_col,
        title = "Kernel mismatch: ATE error",
        subtitle = "Absolute all-or-nothing ATE error at K=(0.8, 0.2) for each sim/fit kernel pair.",
        ylab = "Absolute all-or-nothing ATE error"
      ),
      "robustness_kernel_mismatch_ate_error",
      width = 8.0, height = 4.8
    )
  }
}

support_contrast_display_labels <- c(
  single_cell_flip = "Single-cell flip",
  random_50pct_flip = "Random 50% flip",
  global_1_0 = "All-or-nothing ATE"
)

prepare_support_plot_df <- function(df, family, keep_kernel_pair = FALSE) {
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
    pair_levels <- unique(out$kernel_pair)
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

make_snr_scale_support_plot <- function(df, title, stem, subtitle) {
  lbl_col <- safe_label_col(df)
  if (is.null(lbl_col)) return(NULL)
  out <- df %>%
    filter(.data$scenario_family == "snr_scale")
  if (nrow(out) < 1L) return(NULL)
  if (!"mean_pct_error" %in% names(out)) out$mean_pct_error <- NA_real_
  global_ref <- out %>%
    filter(.data$contrast_family == "global_1_0") %>%
    distinct(.data$scenario_id, .data$mean_psi_truth) %>%
    rename(global_ate_truth = .data$mean_psi_truth)
  support_plot_df <- out %>%
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
    filter(is.finite(.data$mean_pct_error)) %>%
    group_by(.data$contrast_family, .data[[lbl_col]], .data$mu_scale) %>%
    summarize(mean_pct_error = mean(.data$mean_pct_error, na.rm = TRUE), .groups = "drop")
  if (nrow(support_plot_df) < 1L) return(NULL)
  support_plot_df$contrast_display <- support_contrast_display_labels[support_plot_df$contrast_family]
  support_plot_df$contrast_display <- factor(
    support_plot_df$contrast_display,
    levels = unname(support_contrast_display_labels)
  )
  p_support <- ggplot(support_plot_df, aes(x = .data$mu_scale, y = .data$mean_pct_error,
                                           color = .data[[lbl_col]], group = .data[[lbl_col]])) +
    geom_point(size = 1.6, alpha = 0.85) +
    {if (dplyr::n_distinct(support_plot_df$mu_scale) > 1L) geom_line(linewidth = 0.7, alpha = 0.8)} +
    facet_wrap(~ .data$contrast_display, scales = "free_y") +
    labs(
      x = expression(paste(mu, " scale")),
      y = "Mean percentage error (%)",
      color = "Labelling",
      title = title,
      subtitle = subtitle
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  save_plot_pair(p_support, stem, width = 9.5, height = 5.5)
}

make_kernel_support_plot <- function(df, title, stem, subtitle) {
  lbl_col <- safe_label_col(df)
  if (is.null(lbl_col)) return(NULL)
  support_plot_df <- prepare_support_plot_df(df, "kernel_mismatch", keep_kernel_pair = TRUE)
  if (is.null(support_plot_df) || nrow(support_plot_df) < 1L) return(NULL)
  p_support <- ggplot(support_plot_df, aes(x = .data$kernel_pair, y = .data$mean_pct_error,
                                             fill = .data[[lbl_col]])) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~ .data$contrast_display, scales = "free_y") +
    labs(
      x = "Kernel scenario",
      y = "Mean percentage error (%)",
      fill = "Labelling",
      title = title,
      subtitle = subtitle
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
  save_plot_pair(p_support, stem, width = 10.0, height = 5.5)
}

if (!is.null(support_summary) && nrow(support_summary) > 0 &&
    all(c("contrast_family", "mean_abs_error") %in% names(support_summary))) {
  lbl_col <- safe_label_col(support_summary)
  if (!is.null(lbl_col)) {
    plot_files$k_separation_support <- make_support_plot(
      support_summary, "k_separation",
      title = "K-separation: off-support allocation contrasts",
      stem = "robustness_k_separation_support_contrasts",
      subtitle = "Mean percentage error in psi estimates. Single-cell flips use |all-or-nothing ATE truth| as the reference scale when the local contrast is ~0.",
      keep_kernel_pair = FALSE
    )
    plot_files$snr_scale_support <- make_snr_scale_support_plot(
      support_summary,
      title = "Signal-to-noise scaling: off-support allocation contrasts",
      stem = "robustness_snr_scale_support_contrasts",
      subtitle = "Mean percentage error in psi estimates at K=(0.8, 0.2) by mu scale. Single-cell flips use |all-or-nothing ATE truth| as the reference scale when the local contrast is ~0."
    )
    plot_files$kernel_mismatch_support <- make_kernel_support_plot(
      support_summary,
      title = "Kernel mismatch: off-support allocation contrasts",
      stem = "robustness_kernel_mismatch_support_contrasts",
      subtitle = "Mean percentage error in psi estimates by sim/fit kernel pair at K=(0.8, 0.2). Single-cell flips use |all-or-nothing ATE truth| as the reference scale when the local contrast is ~0."
    )
    plot_files$high_count_assignment_support <- make_support_plot(
      support_summary, "high_count_assignment",
      title = "High-count treatment assignment: off-support allocation contrasts",
      stem = "robustness_high_count_assignment_support_contrasts",
      subtitle = "Mean percentage error in psi estimates at K_c=0.8, K_t=0.2 under high-count assignment. Single-cell flips use |all-or-nothing ATE truth| as the reference scale when the local contrast is ~0.",
      keep_kernel_pair = FALSE
    )
  }
}

if (!is.null(decay_showcase_summary) && nrow(decay_showcase_summary) > 0 &&
    all(c("d_mid", "mean_abs_delta", "decay_label") %in% names(decay_showcase_summary))) {
  decay_reps_used <- if (isTRUE(refresh_decay) || isTRUE(auto_refresh_decay_showcase)) decay_reps else 2000L
  p_decay <- make_decay_validation_plot(decay_showcase_summary, decay_reps_used)
  plot_files$decay <- save_plot_pair(p_decay, "robustness_decay_validation", width = 8.5, height = 5.2)
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
}

latex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x
}

ROBUSTNESS_FIG_STEMS <- list(
  k_separation_label = "robustness_k_separation_label_recovery",
  k_separation_ate = "robustness_k_separation_ate_error",
  k_separation_support = "robustness_k_separation_support_contrasts",
  kernel_label = "robustness_kernel_mismatch_label_recovery",
  kernel_ate = "robustness_kernel_mismatch_ate_error",
  kernel_support = "robustness_kernel_mismatch_support_contrasts",
  high_count_label = "robustness_high_count_assignment_label_recovery",
  high_count_ate = "robustness_high_count_assignment_ate_error",
  high_count_support = "robustness_high_count_assignment_support_contrasts",
  snr_label = "robustness_snr_scale_label_recovery",
  snr_ate = "robustness_snr_scale_ate_error",
  snr_support = "robustness_snr_scale_support_contrasts",
  decay_spatial = "robustness_decay_validation",
  decay_temporal = "robustness_temporal_decay_validation"
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
  "%   pdflatex -jobname=robustness_standalone \"\\\\def\\\\robustnessstandalone{}\\\\input{simulation_robustness_appendix.tex}\"",
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
  "The main simulation in \\cref{sec:simulation_study} uses one baseline parameterisation---control $K_c=0.8$, treated $K_t=0.2$---and random treatment assignment. This subsection reports a structured robustness suite designed around the finite-sample assumptions in \\cref{sec:finite-sample-theory} and the identification discussion in \\cref{sec:identification}. Unless stated otherwise, all scenarios use the same window $[0,100]\\times[0,100]$, treatment time $t^\\star=10$, end time $T=110$, a $10\\times10$ cell tessellation, exponential temporal triggering with $(\\alpha,\\beta)=(0.01,10)$, and exponential models at fit time. Treatment lowers the branching ratio relative to control and is assigned to $50\\%$ of cells.",
  "",
  "For each scenario we simulate an ensemble of independent replications and report three classes of outcome. First, \\emph{label recovery}: balanced accuracy of the oracle, naive, best-proposal, adaptive SEM, and full SEM labellings of post-$t^\\star$ events. Second, \\emph{all-or-nothing ATE error}: absolute bias in the plug-in all-treated versus all-control contrast relative to the simulation truth (SEM fits with explosive parameter estimates are excluded). Third, \\emph{off-support allocation contrasts}: mean percentage error in the expected count contrast under (i) flipping a single cell relative to the observed assignment, (ii) a random $50\\%$ relabelling of cells, and (iii) the global all-treated versus all-control allocation. The first two contrasts are local perturbations of the observed design; the third coincides with the all-or-nothing estimand and probes extrapolation to an unsupported allocation.",
  "",
  "The branching-ratio grid, kernel misspecification, assignment designs, signal-to-noise scalings, and single-flip decay diagnostics below are generated automatically by \\texttt{sim\\_study\\_robustness.R} in the \\texttt{PPDisentangle} package. Figures are reported in \\cref{app:robustness-k-separation,app:robustness-kernel,app:robustness-high-count,app:robustness-snr,app:robustness-decay}.",
  "",
  "\\makeatletter",
  "\\ifdefined\\robustnessstandalone",
  "  \\def\\robustness@figdir{figures/}",
  "\\else",
  "  \\providecommand{\\robustness@figdir}{plots/sim_study/robustness/}",
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
  "\\makeatother",
  "",
  tex_subsubsection(
    "Control and treatment process K-separation",
    "app:robustness-k-separation",
    c(
      "The main study in \\cref{sec:simulation_study} fixes a large separation between control and treated branching ratios, with treatment lowering the branching ratio from $K_c=0.8$ to $K_t=0.2$. Here we stress-test performance as that separation shrinks. Throughout, control and treated processes share $(\\alpha,\\beta)=(0.01,10)$ and exponential temporal triggering; treatment is assigned independently to $50\\%$ of cells. The control branching ratio is fixed at $K_c=0.8$ and the treated ratio varies over $K_t\\in\\{0.1,0.2,0.3,0.4,0.5,0.6,0.7\\}$, always below control. For each pair $(K_c,K_t)$ the background rate $\\mu$ is recalibrated so the expected post-treatment point count remains fixed; small $K_c-K_t$ therefore isolates weak component separation at comparable event abundance.",
      "\\Cref{fig:robustness-k-separation-label,fig:robustness-k-separation-ate,fig:robustness-k-separation-support} report balanced labelling accuracy, absolute all-or-nothing ATE error (excluding explosive SEM fits), and mean percentage error in three off-support allocation contrasts: a single-cell flip from the observed assignment, a random $50\\%$ relabelling of cells, and the global all-treated versus all-control contrast. These diagnostics are directly motivated by the finite-sample theory: label recovery feeds the E-step surrogates in Assumptions~\\ref{ass:G4}--\\ref{ass:G8}, while off-support contrasts probe the regime-stable extrapolations in \\cref{sec:identification}."
    ),
    c(
      tex_fig(ROBUSTNESS_FIG_STEMS$k_separation_label, "fig:robustness-k-separation-label",
              "K-separation sensitivity for label recovery. Control $K=0.8$; treated $K \\in \\{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7\\}$ with treatment lowering $K$; $\\mu$ is calibrated to hold expected point count fixed."),
      tex_fig(ROBUSTNESS_FIG_STEMS$k_separation_ate, "fig:robustness-k-separation-ate",
              "K-separation sensitivity for all-or-nothing ATE error."),
      tex_fig(ROBUSTNESS_FIG_STEMS$k_separation_support, "fig:robustness-k-separation-support",
              "K-separation off-support allocation contrast accuracy (mean percentage error). The rightmost contrast is the all-or-nothing ATE.")
    )
  ),
  tex_subsubsection(
    "Mis-specified temporal kernels",
    "app:robustness-kernel",
    c(
      "Assumption~\\ref{ass:G3} and the identification discussion in \\cref{sec:identification} treat the triggering kernel as part of the causal specification: misspecification changes both fit and the transmitted spillover geometry. We therefore simulate from exponential and power-law temporal kernels and fit exponential or power-law models, covering all four sim/fit pairings at $K=(0.8,0.2)$ with random $50\\%$ treatment assignment and the same $\\mu$-calibration scheme as above.",
      "Power-law temporal decay generates heavier tails than the fitted exponential, so this exercise is a structured check on whether SEM labelling and plug-in ATE estimation remain stable when the true memory of the process is longer than the fitted kernel assumes.",
      "\\Cref{fig:robustness-kernel-label,fig:robustness-kernel-ate,fig:robustness-kernel-support} report label recovery, all-or-nothing ATE error, and off-support contrast accuracy across the four sim/fit kernel pairings."
    ),
    c(
      tex_fig(ROBUSTNESS_FIG_STEMS$kernel_label, "fig:robustness-kernel-label",
              "Kernel specification check for label recovery across all four sim/fit kernel pairs at $K=(0.8, 0.2)$."),
      tex_fig(ROBUSTNESS_FIG_STEMS$kernel_ate, "fig:robustness-kernel-ate",
              "Kernel specification check for all-or-nothing ATE error across all four sim/fit kernel pairs."),
      tex_fig(ROBUSTNESS_FIG_STEMS$kernel_support, "fig:robustness-kernel-support",
              "Kernel mismatch off-support allocation contrast accuracy (mean percentage error) by sim/fit kernel pair, faceted by contrast type including all-or-nothing ATE.")
    )
  ),
  tex_subsubsection(
    "Pre-treatment data informed treatment allocation",
    "app:robustness-high-count",
    c(
      "The baseline design randomises treated cells. Many applications instead concentrate treatment on historically active locations---for example, when a policy targets high-activity areas while the intervention reduces subsequent productivity. To mimic this structure, we simulate a reference pre-treatment catalogue from the control Hawkes law on $[0,t^\\star]$, count points per cell, and assign treatment to the $50\\%$ of cells with the largest reference counts. The post-treatment law uses $K_c=0.8$ and $K_t=0.2$, so treated cells are less self-exciting than controls even though they were selected for prior activity.",
      "This scenario is especially demanding for naive location-based labelling, because cell location and latent process type are deliberately misaligned. \\Cref{fig:robustness-high-count-label,fig:robustness-high-count-ate,fig:robustness-high-count-support} summarise labelling accuracy, all-or-nothing ATE error, and off-support allocation contrast accuracy under this assignment mechanism."
    ),
    c(
      tex_fig(ROBUSTNESS_FIG_STEMS$high_count_label, "fig:robustness-high-count-label",
              "Label recovery when treatment is assigned to the 50\\% of cells with the most points in a reference pre-treatment catalogue at $K_c=0.8$, $K_t=0.2$."),
      tex_fig(ROBUSTNESS_FIG_STEMS$high_count_ate, "fig:robustness-high-count-ate",
              "All-or-nothing ATE error under high-count treatment assignment at $K_c=0.8$, $K_t=0.2$. Treatment lowers the branching ratio relative to control, as in the main simulation study."),
      tex_fig(ROBUSTNESS_FIG_STEMS$high_count_support, "fig:robustness-high-count-support",
              "High-count assignment off-support allocation contrast accuracy (mean percentage error). The rightmost contrast is the all-or-nothing ATE.")
    )
  ),
  tex_subsubsection(
    "Signal to noise ratio",
    "app:robustness-snr",
    c(
      sprintf(
        "Finite-sample behaviour also depends on how many post-treatment events are available relative to the ambiguity created by spillover. Holding $K=(0.8,0.2)$ and the random assignment design fixed, we vary a background-rate multiplier $\\mu_{\\mathrm{scale}}\\in%s$ while calibrating $\\mu$ from a common target abundance; realized point counts therefore vary across the grid. Lower $\\mu_{\\mathrm{scale}}$ increases the relative contribution of triggered offspring to the observed pattern; higher values push the process toward a Poisson-like regime with weaker history dependence.",
        format_mu_scale_grid(mu_scales, for_tex = TRUE)
      ),
      "\\Cref{fig:robustness-snr-label,fig:robustness-snr-ate,fig:robustness-snr-support} trace label recovery, ATE error, and off-support contrast accuracy across this signal-to-noise grid."
    ),
    c(
      tex_fig(ROBUSTNESS_FIG_STEMS$snr_label, "fig:robustness-snr-label",
              sprintf(
                "Signal-to-noise sensitivity for label recovery at $K=(0.8, 0.2)$ with $\\mu$ scale in %s.",
                format_mu_scale_grid(mu_scales, for_tex = TRUE)
              )),
      tex_fig(ROBUSTNESS_FIG_STEMS$snr_ate, "fig:robustness-snr-ate",
              "Signal-to-noise sensitivity for all-or-nothing ATE error at $K=(0.8, 0.2)$ across $\\mu$ scales."),
      tex_fig(ROBUSTNESS_FIG_STEMS$snr_support, "fig:robustness-snr-support",
              "Signal-to-noise off-support allocation contrast accuracy (mean percentage error) across $\\mu$ scales at $K=(0.8, 0.2)$.")
    )
  ),
  tex_subsubsection(
    "Validating temporal and spatial decay",
    "app:robustness-decay",
    c(
      "Assumption~\\ref{ass:G3} requires that a single latent label flip has localized influence on the intensity path. We validate this forward-simulation property directly under true Hawkes parameters, separate from the estimation experiments above. In each replicate we sample one post-treatment event in a treated cell, flip only its latent label (control $\\leftrightarrow$ treated), and simulate the two resulting catalogues under the same parameter values. Spatial decay is measured by the mean absolute count difference $|\\Delta N|$ in distance annuli around the flipped event; temporal decay uses the same flip but bins offspring by lag since the flip time.",
      "Three exponential-kernel $K$-separation scenarios are shown in the spatial panel: $K_c=0.8$ with $K_t\\in\\{0.2,0.5,0.7\\}$. The temporal panel overlays $K_t\\in\\{0.2,0.5\\}$ under exponential triggering with an additional power-law temporal kernel at $K=(0.8,0.2)$. \\Cref{fig:robustness-decay,fig:robustness-decay-temporal} confirm geometric decay in space and algebraic decay in time for the exponential specification, with slower temporal persistence under power-law triggering."
    ),
    c(
      tex_fig(ROBUSTNESS_FIG_STEMS$decay_spatial, "fig:robustness-decay",
              "Forward-simulation spatial decay validation for three exponential-kernel $K$-separation scenarios: $K_c=0.8$ with $K_t\\in\\{0.2,0.5,0.7\\}$. One event in a cell receives a latent label flip; curves show mean $|\\Delta N|$ per spatial annulus around the flip point."),
      tex_fig(ROBUSTNESS_FIG_STEMS$decay_temporal, "fig:robustness-decay-temporal",
              "Forward-simulation temporal decay validation for $K_c=0.8$ with $K_t\\in\\{0.2,0.5\\}$ (exponential) and $K_t=0.2$ (power-law temporal kernel). One event in a cell receives a latent label flip and curves show mean $|\\Delta N|$ per lag bin since the flip time.")
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
