#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg_val <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) < 1L || idx[[1]] >= length(args)) return(default)
  args[[idx[[1]] + 1L]]
}

parse_multiplier_vec <- function(raw, default = c(0.1, 0.5, 1, 2)) {
  if (is.null(raw) || !nzchar(raw)) return(default)
  vals <- suppressWarnings(as.numeric(strsplit(raw, ",", fixed = TRUE)[[1]]))
  vals <- vals[is.finite(vals) & vals > 0]
  if (length(vals) < 1L) return(default)
  unique(vals)
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
out_dir <- file.path(repo_dir, "output", "sim_study")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
old_wd <- getwd()
setwd(repo_dir)
on.exit(setwd(old_wd), add = TRUE)

pp_sims <- suppressWarnings(as.integer(get_arg_val("--sims", Sys.getenv("PP_SIMS", "32"))))
if (!is.finite(pp_sims) || is.na(pp_sims) || pp_sims < 1L) pp_sims <- 32L
test_mode <- "--test" %in% args
multipliers <- parse_multiplier_vec(
  raw = get_arg_val("--multipliers", Sys.getenv("PP_POST_TIME_MULTIPLIERS", "0.1,0.5,1,2"))
)
run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_basename <- paste0("time_sweep_", run_stamp)
if (nzchar(Sys.getenv("SLURM_JOB_ID", ""))) {
  output_basename <- paste0("time_sweep_", Sys.getenv("SLURM_JOB_ID"))
}

format_tag <- function(mult) {
  gsub("[^0-9A-Za-z]+", "p", sprintf("tm_%0.3f", mult))
}

build_core_method <- function(df) {
  if (is.null(df) || nrow(df) < 1L || !("labelling" %in% names(df))) return(NULL)
  out <- df %>%
    filter(.data$labelling %in% c("oracle", "naive", "SEM_full")) %>%
    mutate(method = ifelse(.data$labelling == "SEM_full", "SEM", as.character(.data$labelling)))
  out$method <- factor(out$method, levels = c("oracle", "naive", "SEM"))
  out
}

slim_sim_result <- function(res) {
  if (is.null(res) || !is.list(res)) return(NULL)
  out <- res
  out$results_flat <- NULL
  out$tasks <- NULL
  out$EM_results <- NULL
  out$sem_diagnostics_all <- NULL
  out
}

runs_all <- list()
runs_core <- list()
runs_tc_all <- list()
runs_tc_core <- list()
per_run_results <- list()
run_manifest <- list()
processed_mults <- numeric(0)
for (mult in multipliers) {
  tag <- format_tag(mult)
  run_basename <- paste0(output_basename, "_", tag)
  message(sprintf("[time-sweep] running multiplier=%s (tag=%s)", signif(mult, 4), tag))
  env <- c(
    PP_POST_TIME_MULTIPLIER = as.character(mult),
    PP_OUTPUT_BASENAME = run_basename,
    PP_OUTPUT_TAG = ""
  )
  cmd_args <- c(file.path("inst", "sim_study", "sim_study.R"), "--sims", as.character(pp_sims))
  if (test_mode) cmd_args <- c(cmd_args, "--test")
  status <- system2(command = "Rscript", args = cmd_args, env = sprintf("%s=%s", names(env), unname(env)))
  if (!identical(status, 0L)) stop(sprintf("sim_study.R failed for multiplier=%s", mult))

  rds_path <- file.path(out_dir, paste0(run_basename, ".rds"))
  if (!file.exists(rds_path)) stop(sprintf("Expected output missing: %s", rds_path))
  res <- readRDS(rds_path)
  df <- if (!is.null(res$results_df)) res$results_df else NULL
  if (is.null(df) || nrow(df) < 1L) {
    stop(sprintf("results_df missing/empty for multiplier=%s (%s)", mult, rds_path))
  }
  if (!("all_nothing_theory" %in% names(df))) {
    stop(sprintf("results_df lacks all_nothing_theory for multiplier=%s (%s)", mult, rds_path))
  }
  got_mult <- if (!is.null(res$config) && !is.null(res$config$POST_TIME_MULTIPLIER)) {
    suppressWarnings(as.numeric(res$config$POST_TIME_MULTIPLIER))
  } else {
    NA_real_
  }
  if (is.finite(got_mult) && !isTRUE(all.equal(got_mult, as.numeric(mult), tolerance = 1e-8))) {
    stop(sprintf("Multiplier mismatch: requested=%s but file has POST_TIME_MULTIPLIER=%s (%s)", mult, got_mult, rds_path))
  }
  all_df <- df
  all_df$method <- as.character(all_df$labelling)
  core <- build_core_method(df)
  df_tc <- if (!is.null(res$results_df_true_control)) res$results_df_true_control else NULL
  tc_core <- if (!is.null(df_tc)) build_core_method(df_tc) else NULL
  end_time <- if (!is.null(res$config) && !is.null(res$config$END_TIME)) as.numeric(res$config$END_TIME) else NA_real_
  treat_time <- if (!is.null(res$config) && !is.null(res$config$TREATMENT_TIME)) as.numeric(res$config$TREATMENT_TIME) else NA_real_
  window_len <- if (is.finite(end_time) && is.finite(treat_time)) end_time - treat_time else NA_real_
  all_df$time_multiplier <- as.numeric(mult)
  all_df$time_window <- window_len
  runs_all[[length(runs_all) + 1L]] <- all_df
  per_run_results[[tag]] <- list(
    multiplier = as.numeric(mult),
    tag = tag,
    run_basename = run_basename,
    rds_path = rds_path,
    result = slim_sim_result(res)
  )
  run_manifest[[length(run_manifest) + 1L]] <- data.frame(
    multiplier = as.numeric(mult),
    tag = tag,
    run_basename = run_basename,
    rds_path = rds_path,
    rows_results_df = nrow(df),
    stringsAsFactors = FALSE
  )
  processed_mults <- c(processed_mults, as.numeric(mult))
  if (!is.null(core) && nrow(core) > 0) {
    core$time_multiplier <- as.numeric(mult)
    core$time_window <- window_len
    runs_core[[length(runs_core) + 1L]] <- core
  }
  if (!is.null(df_tc) && nrow(df_tc) > 0 && "all_nothing_true_control" %in% names(df_tc)) {
    tc_df <- df_tc
    tc_df$method <- as.character(tc_df$labelling)
    tc_df$time_multiplier <- as.numeric(mult)
    tc_df$time_window <- window_len
    runs_tc_all[[length(runs_tc_all) + 1L]] <- tc_df
    if (!is.null(tc_core) && nrow(tc_core) > 0) {
      tc_core$time_multiplier <- as.numeric(mult)
      tc_core$time_window <- window_len
      runs_tc_core[[length(runs_tc_core) + 1L]] <- tc_core
    }
  }
}

missing_mults <- setdiff(as.numeric(multipliers), unique(processed_mults))
if (length(missing_mults) > 0) {
  stop(sprintf("Time sweep missing multipliers: %s", paste(signif(missing_mults, 5), collapse = ", ")))
}

if (length(runs_all) < 1L) stop("No valid runs found for time sweep.")
sweep_df_all <- bind_rows(runs_all)
sweep_df_all$time_mult_f <- factor(
  paste0(signif(sweep_df_all$time_multiplier, 3), "x"),
  levels = paste0(signif(multipliers, 3), "x")
)
method_levels_all <- unique(sweep_df_all$method)
sweep_df_all$method <- factor(sweep_df_all$method, levels = method_levels_all)
p_all <- ggplot(sweep_df_all, aes(x = .data$method, y = .data$all_nothing_theory, fill = .data$time_mult_f)) +
  geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single")) +
  labs(x = "Method", y = "All-Nothing ATE Estimate", fill = "Post-treatment window") +
  theme_minimal() +
  theme(plot.title = element_blank(), plot.subtitle = element_blank())

png_all <- file.path(out_dir, paste0(output_basename, "_all_nothing_grouped_boxplot_all_methods.png"))
pdf_all <- file.path(out_dir, paste0(output_basename, "_all_nothing_grouped_boxplot_all_methods.pdf"))
ggsave(filename = png_all, plot = p_all, width = 10.5, height = 5.5, dpi = 300)
ggsave(filename = pdf_all, plot = p_all, width = 10.5, height = 5.5)

sweep_df_core <- if (length(runs_core) > 0) bind_rows(runs_core) else data.frame()
if (nrow(sweep_df_core) > 0) {
  sweep_df_core$time_mult_f <- factor(
    paste0(signif(sweep_df_core$time_multiplier, 3), "x"),
    levels = paste0(signif(multipliers, 3), "x")
  )
  sweep_df_core$method <- factor(sweep_df_core$method, levels = c("oracle", "naive", "SEM"))
  p_core <- ggplot(sweep_df_core, aes(x = .data$method, y = .data$all_nothing_theory, fill = .data$time_mult_f)) +
    geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single")) +
    labs(x = "Method", y = "All-Nothing ATE Estimate", fill = "Post-treatment window") +
    theme_minimal() +
    theme(plot.title = element_blank(), plot.subtitle = element_blank())
  ggsave(
    filename = file.path(out_dir, paste0(output_basename, "_all_nothing_grouped_boxplot_core_methods.png")),
    plot = p_core, width = 8.5, height = 5.5, dpi = 300
  )
  ggsave(
    filename = file.path(out_dir, paste0(output_basename, "_all_nothing_grouped_boxplot_core_methods.pdf")),
    plot = p_core, width = 8.5, height = 5.5
  )
}

sweep_df_true_control_all <- if (length(runs_tc_all) > 0) bind_rows(runs_tc_all) else data.frame()
sweep_df_true_control_core <- if (length(runs_tc_core) > 0) bind_rows(runs_tc_core) else data.frame()
p_all_true_control <- NULL
p_core_true_control <- NULL
png_all_true_control <- NA_character_
pdf_all_true_control <- NA_character_
png_core_true_control <- NA_character_
pdf_core_true_control <- NA_character_
if (nrow(sweep_df_true_control_all) > 0) {
  sweep_df_true_control_all$time_mult_f <- factor(
    paste0(signif(sweep_df_true_control_all$time_multiplier, 3), "x"),
    levels = paste0(signif(multipliers, 3), "x")
  )
  method_levels_tc_all <- unique(sweep_df_true_control_all$method)
  sweep_df_true_control_all$method <- factor(sweep_df_true_control_all$method, levels = method_levels_tc_all)
  p_all_true_control <- ggplot(
    sweep_df_true_control_all,
    aes(x = .data$method, y = .data$all_nothing_true_control, fill = .data$time_mult_f)
  ) +
    geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single")) +
    labs(x = "Method", y = "All-Nothing ATE Estimate (True/Control-Fixed)", fill = "Post-treatment window") +
    theme_minimal() +
    theme(plot.title = element_blank(), plot.subtitle = element_blank())
  png_all_true_control <- file.path(out_dir, paste0(output_basename, "_all_nothing_grouped_boxplot_all_methods_true_control.png"))
  pdf_all_true_control <- file.path(out_dir, paste0(output_basename, "_all_nothing_grouped_boxplot_all_methods_true_control.pdf"))
  ggsave(filename = png_all_true_control, plot = p_all_true_control, width = 10.5, height = 5.5, dpi = 300)
  ggsave(filename = pdf_all_true_control, plot = p_all_true_control, width = 10.5, height = 5.5)
}
if (nrow(sweep_df_true_control_core) > 0) {
  sweep_df_true_control_core$time_mult_f <- factor(
    paste0(signif(sweep_df_true_control_core$time_multiplier, 3), "x"),
    levels = paste0(signif(multipliers, 3), "x")
  )
  sweep_df_true_control_core$method <- factor(sweep_df_true_control_core$method, levels = c("oracle", "naive", "SEM"))
  p_core_true_control <- ggplot(
    sweep_df_true_control_core,
    aes(x = .data$method, y = .data$all_nothing_true_control, fill = .data$time_mult_f)
  ) +
    geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single")) +
    labs(x = "Method", y = "All-Nothing ATE Estimate (True/Control-Fixed)", fill = "Post-treatment window") +
    theme_minimal() +
    theme(plot.title = element_blank(), plot.subtitle = element_blank())
  png_core_true_control <- file.path(out_dir, paste0(output_basename, "_all_nothing_grouped_boxplot_core_methods_true_control.png"))
  pdf_core_true_control <- file.path(out_dir, paste0(output_basename, "_all_nothing_grouped_boxplot_core_methods_true_control.pdf"))
  ggsave(filename = png_core_true_control, plot = p_core_true_control, width = 8.5, height = 5.5, dpi = 300)
  ggsave(filename = pdf_core_true_control, plot = p_core_true_control, width = 8.5, height = 5.5)
}

csv_out <- file.path(out_dir, paste0(output_basename, "_all_nothing_grouped_data.csv"))
write.csv(sweep_df_all, csv_out, row.names = FALSE)
if (nrow(sweep_df_core) > 0) {
  write.csv(sweep_df_core, file.path(out_dir, paste0(output_basename, "_all_nothing_grouped_data_core_methods.csv")), row.names = FALSE)
}
if (nrow(sweep_df_true_control_all) > 0) {
  write.csv(sweep_df_true_control_all, file.path(out_dir, paste0(output_basename, "_all_nothing_grouped_data_all_methods_true_control.csv")), row.names = FALSE)
}
if (nrow(sweep_df_true_control_core) > 0) {
  write.csv(sweep_df_true_control_core, file.path(out_dir, paste0(output_basename, "_all_nothing_grouped_data_core_methods_true_control.csv")), row.names = FALSE)
}

rds_out <- file.path(out_dir, paste0(output_basename, "_summary.rds"))
saveRDS(list(
  output_basename = output_basename,
  multipliers = multipliers,
  run_manifest = if (length(run_manifest) > 0) bind_rows(run_manifest) else data.frame(),
  # Master list of the four sim_study runs (one per multiplier)
  per_run_results = per_run_results,
  sweep_df_all = sweep_df_all,
  sweep_df_core = sweep_df_core,
  sweep_df_true_control_all = sweep_df_true_control_all,
  sweep_df_true_control_core = sweep_df_true_control_core,
  grouped_plots = list(
    all_methods = p_all,
    core_methods = if (exists("p_core")) p_core else NULL,
    all_methods_true_control = p_all_true_control,
    core_methods_true_control = p_core_true_control
  ),
  grouped_plot_files = list(
    all_methods_png = png_all,
    all_methods_pdf = pdf_all,
    core_methods_png = if (nrow(sweep_df_core) > 0) file.path(out_dir, paste0(output_basename, "_all_nothing_grouped_boxplot_core_methods.png")) else NA_character_,
    core_methods_pdf = if (nrow(sweep_df_core) > 0) file.path(out_dir, paste0(output_basename, "_all_nothing_grouped_boxplot_core_methods.pdf")) else NA_character_,
    all_methods_true_control_png = png_all_true_control,
    all_methods_true_control_pdf = pdf_all_true_control,
    core_methods_true_control_png = png_core_true_control,
    core_methods_true_control_pdf = pdf_core_true_control
  )
), rds_out)

message("[time-sweep] wrote grouped boxplot:")
message("  ", png_all)
message("  ", pdf_all)
message("[time-sweep] wrote data: ", csv_out)
message("[time-sweep] wrote summary: ", rds_out)
