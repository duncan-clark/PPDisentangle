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
target_points <- suppressWarnings(as.numeric(get_arg_val("--target-points", Sys.getenv("PP_TARGET_POINTS", "2500"))))
if (!is.finite(target_points) || is.na(target_points) || target_points <= 0) target_points <- 2500
test_mode <- "--test" %in% args

k_anchor <- suppressWarnings(as.numeric(get_arg_val("--k-anchor", Sys.getenv("PP_K_ANCHOR", "0.2"))))
if (!is.finite(k_anchor) || is.na(k_anchor) || k_anchor < 0 || k_anchor >= 1) k_anchor <- 0.2
k_values <- parse_num_vec(get_arg_val("--k-values", Sys.getenv("PP_K_VALUES", "0.3,0.4,0.5,0.6,0.7,0.8")),
                          default = seq(0.3, 0.8, by = 0.1))
k_values <- k_values[k_values >= 0 & k_values < 1 & abs(k_values - k_anchor) > 1e-10]
mu_scales <- parse_num_vec(get_arg_val("--mu-scales", Sys.getenv("PP_MU_SCALES", "0.5,2")),
                           default = c(0.5, 2))
mu_scales <- mu_scales[mu_scales > 0]

scenario_filter <- get_arg_val("--scenario-set", Sys.getenv("PP_ROBUSTNESS_SCENARIO_SET", "all"))
run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_basename <- paste0("robustness_", run_stamp)
if (nzchar(Sys.getenv("SLURM_JOB_ID", ""))) {
  output_basename <- paste0("robustness_", Sys.getenv("SLURM_JOB_ID"))
}

k_separation <- data.frame(
  scenario_family = "k_separation",
  control_k = k_anchor,
  treated_k = k_values,
  mu_scale = 1,
  sim_kernel = "exponential",
  fit_kernel = "exponential",
  stringsAsFactors = FALSE
)

snr_scale <- data.frame(
  scenario_family = "snr_scale",
  control_k = 0.2,
  treated_k = 0.8,
  mu_scale = mu_scales,
  sim_kernel = "exponential",
  fit_kernel = "exponential",
  stringsAsFactors = FALSE
)

joint_shrink <- expand.grid(
  scenario_family = "joint_k_mu_shrink",
  control_k = k_anchor,
  treated_k = k_values[k_values <= min(max(k_values), k_anchor + 0.3)],
  mu_scale = mu_scales,
  sim_kernel = "exponential",
  fit_kernel = "exponential",
  stringsAsFactors = FALSE
)

kernel_mismatch <- expand.grid(
  scenario_family = "kernel_mismatch",
  control_k = 0.2,
  treated_k = 0.8,
  mu_scale = 1,
  sim_kernel = c("exponential", "power_law"),
  fit_kernel = c("exponential", "power_law"),
  stringsAsFactors = FALSE
)

scenarios <- bind_rows(k_separation, snr_scale, joint_shrink, kernel_mismatch) %>%
  mutate(
    k_delta = abs(.data$treated_k - .data$control_k),
    scenario_id = paste0(
      .data$scenario_family,
      "_kc", format_num_tag(.data$control_k),
      "_kt", format_num_tag(.data$treated_k),
      "_mu", format_num_tag(.data$mu_scale),
      "_sim", .data$sim_kernel,
      "_fit", .data$fit_kernel
    )
  ) %>%
  distinct(.data$scenario_id, .keep_all = TRUE)

if (!identical(tolower(scenario_filter), "all")) {
  keep_families <- unlist(strsplit(scenario_filter, "[,;[:space:]]+", perl = TRUE), use.names = FALSE)
  keep_families <- keep_families[nzchar(keep_families)]
  scenarios <- scenarios %>% filter(.data$scenario_family %in% keep_families)
}
if (nrow(scenarios) < 1L) stop("No robustness scenarios selected.")

fig_dir <- file.path(out_dir, "generated", "robustness", "figures")
tex_dir <- file.path(out_dir, "generated", "robustness", "tex")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tex_dir, recursive = TRUE, showWarnings = FALSE)

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
    PP_CONTROL_K = as.character(sc$control_k),
    PP_TREATED_K = as.character(sc$treated_k),
    PP_MU_SCALE = as.character(sc$mu_scale),
    PP_TARGET_POINTS = as.character(target_points),
    PP_SIM_KERNEL = sc$sim_kernel,
    PP_FIT_KERNEL = sc$fit_kernel
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
    target_points = target_points,
    run_basename = run_basename,
    rds_path = rds_path,
    stringsAsFactors = FALSE
  )
}

manifest <- bind_rows(lapply(seq_len(nrow(scenarios)), run_one))
manifest_path <- file.path(out_dir, paste0(output_basename, "_manifest.csv"))
write.csv(manifest, manifest_path, row.names = FALSE)

summary_rows <- lapply(seq_len(nrow(manifest)), function(i) {
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
label_summary <- bind_rows(lapply(summary_rows, `[[`, "label")) %>%
  left_join(manifest, by = "scenario_id")
support_summary <- bind_rows(lapply(summary_rows, `[[`, "support")) %>%
  left_join(manifest, by = "scenario_id")
ate_summary <- bind_rows(lapply(summary_rows, `[[`, "ate")) %>%
  left_join(manifest, by = "scenario_id")
decay_summary <- bind_rows(lapply(summary_rows, `[[`, "decay")) %>%
  left_join(manifest, by = "scenario_id")

write.csv(label_summary, file.path(out_dir, paste0(output_basename, "_label_recovery_summary.csv")), row.names = FALSE)
write.csv(support_summary, file.path(out_dir, paste0(output_basename, "_support_contrast_summary.csv")), row.names = FALSE)
write.csv(ate_summary, file.path(out_dir, paste0(output_basename, "_ate_summary.csv")), row.names = FALSE)
write.csv(decay_summary, file.path(out_dir, paste0(output_basename, "_decay_validation_summary.csv")), row.names = FALSE)

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

plot_files <- list()
if (!is.null(label_summary) && nrow(label_summary) > 0) {
  lbl_col <- safe_label_col(label_summary)
  acc_col <- safe_metric_col(label_summary, c("mean_balanced_accuracy", "mean_accuracy", "accuracy"))
  if (!is.null(lbl_col) && !is.null(acc_col)) {
    label_plot_df <- label_summary %>%
      filter(.data[[lbl_col]] %in% c("SEM_adaptive", "SEM_full", "best", "naive", "oracle"))
    p_label <- ggplot(label_plot_df, aes(x = .data$k_delta, y = .data[[acc_col]],
                                         color = .data[[lbl_col]], group = .data[[lbl_col]])) +
      geom_point(size = 1.8, alpha = 0.85) +
      {if (dplyr::n_distinct(label_plot_df$k_delta) > 1L) geom_line(linewidth = 0.7, alpha = 0.8)} +
      facet_wrap(~ scenario_family, scales = "free_x") +
      labs(
        x = "|K treated - K control|",
        y = gsub("_", " ", acc_col),
        color = "Labelling",
        title = "Label recovery across robustness scenarios",
        subtitle = "Higher curves indicate better recovery; small K separation is the hard regime."
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
    plot_files$label_recovery <- save_plot_pair(p_label, "robustness_label_recovery")
  }
}

if (!is.null(ate_summary) && nrow(ate_summary) > 0) {
  lbl_col <- safe_label_col(ate_summary)
  ate_col <- safe_metric_col(ate_summary, c("mean_all_nothing", "mean_all_nothing_true_control"))
  if (!is.null(lbl_col) && !is.null(ate_col) && "true_all_nothing_ATE" %in% names(ate_summary)) {
    ate_plot_df <- ate_summary %>%
      mutate(abs_all_nothing_error = abs(.data[[ate_col]] - .data$true_all_nothing_ATE)) %>%
      filter(is.finite(.data$abs_all_nothing_error))
    p_ate <- ggplot(ate_plot_df, aes(x = .data$k_delta, y = .data$abs_all_nothing_error,
                                     color = .data[[lbl_col]], group = .data[[lbl_col]])) +
      geom_point(size = 1.8, alpha = 0.85) +
      {if (dplyr::n_distinct(ate_plot_df$k_delta) > 1L) geom_line(linewidth = 0.7, alpha = 0.8)} +
      facet_wrap(~ scenario_family, scales = "free_x") +
      labs(
        x = "|K treated - K control|",
        y = "Absolute all-or-nothing ATE error",
        color = "Labelling",
        title = "ATE robustness across scenario families",
        subtitle = "Error is relative to the known all-or-nothing truth from the scenario parameters."
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
    plot_files$ate_error <- save_plot_pair(p_ate, "robustness_ate_error")
  }
}

if (!is.null(support_summary) && nrow(support_summary) > 0 &&
    all(c("contrast_family", "mean_abs_error") %in% names(support_summary))) {
  lbl_col <- safe_label_col(support_summary)
  if (!is.null(lbl_col)) {
    support_plot_df <- support_summary %>%
      filter(is.finite(.data$mean_abs_error))
    p_support <- ggplot(support_plot_df, aes(x = .data$contrast_family, y = .data$mean_abs_error,
                                             fill = .data[[lbl_col]])) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      facet_wrap(~ scenario_family, scales = "free_y") +
      labs(
        x = "Off-support contrast",
        y = "Mean absolute error",
        fill = "Labelling",
        title = "Off-support allocation contrast robustness",
        subtitle = "Includes single-cell flips, random 50% flips, and the global all-treated/all-control contrast."
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
    plot_files$support <- save_plot_pair(p_support, "robustness_support_contrasts", width = 8.5, height = 5.0)
  }
}

if (!is.null(label_summary) && nrow(label_summary) > 0 &&
    any(label_summary$scenario_family == "kernel_mismatch")) {
  lbl_col <- safe_label_col(label_summary)
  acc_col <- safe_metric_col(label_summary, c("mean_balanced_accuracy", "mean_accuracy", "accuracy"))
  if (!is.null(lbl_col) && !is.null(acc_col)) {
    kernel_df <- label_summary %>%
      filter(.data$scenario_family == "kernel_mismatch") %>%
      mutate(kernel_pair = paste0(.data$sim_kernel, " sim / ", .data$fit_kernel, " fit"))
    p_kernel <- ggplot(kernel_df, aes(x = .data$kernel_pair, y = .data[[acc_col]], fill = .data[[lbl_col]])) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      labs(
        x = "Kernel scenario",
        y = gsub("_", " ", acc_col),
        fill = "Labelling",
        title = "Kernel misspecification robustness",
        subtitle = "Power-law and exponential Hawkes simulations fitted under both kernels."
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
    plot_files$kernel <- save_plot_pair(p_kernel, "robustness_kernel_mismatch", width = 8.0, height = 4.8)
  }
}

if (!is.null(decay_summary) && nrow(decay_summary) > 0 &&
    all(c("d_mid", "mean_abs_delta") %in% names(decay_summary))) {
  decay_plot_df <- decay_summary %>%
    filter(is.finite(.data$d_mid), is.finite(.data$mean_abs_delta)) %>%
    group_by(.data$scenario_family, .data$sim_kernel, .data$fit_kernel, .data$d_mid) %>%
    summarize(mean_abs_delta = mean(.data$mean_abs_delta, na.rm = TRUE), .groups = "drop")
  positive <- decay_plot_df$mean_abs_delta[decay_plot_df$mean_abs_delta > 0]
  eps <- if (length(positive) > 0L) min(positive, na.rm = TRUE) / 2 else 1e-6
  decay_plot_df$mean_abs_delta_plot <- pmax(decay_plot_df$mean_abs_delta, eps)
  p_decay <- ggplot(decay_plot_df, aes(x = .data$d_mid, y = .data$mean_abs_delta_plot,
                                       color = .data$scenario_family)) +
    geom_line(linewidth = 0.7, alpha = 0.85) +
    geom_point(size = 1.2, alpha = 0.8) +
    scale_y_log10() +
    labs(
      x = "Distance from flipped cell",
      y = "Mean |Delta N| per annulus",
      color = "Scenario family",
      title = "Forward-simulation decay validation",
      subtitle = "Coupled single-cell flips; outside-cell clusters cancel by construction."
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  plot_files$decay <- save_plot_pair(p_decay, "robustness_decay_validation")
}

latex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x
}

tex_fig <- function(file_info, label, caption) {
  if (is.null(file_info)) return(character(0))
  paste0(
    "\\begin{figure}[t]\n",
    "\\centering\n",
    "\\includegraphics[width=0.95\\linewidth]{figures/", file_info$stem, ".pdf}\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\end{figure}\n"
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
  "% Copy this file and the generated figures/ directory to Overleaf, or \\input{} it from a wrapper.",
  "\\section{Robustness simulation checks}\\label{sec:robustness-simulations}",
  "",
  paste0(
    "This robustness suite re-runs the simulation study over scenario families ",
    scenario_table,
    ". Each scenario stores the full simulation object, label-recovery summaries, ",
    "off-support allocation contrasts, all-or-nothing ATE summaries, and the forward-only ",
    "decay-validation diagnostic. The figures below are generated from ",
    "\\texttt{", latex_escape(basename(paste0(output_basename, "_summary.rds"))), "}."
  ),
  "",
  tex_fig(plot_files$label_recovery, "fig:robustness-label-recovery",
          "Label recovery across robustness scenarios. Curves are plotted against the separation in branching ratios, with facets for the scenario family."),
  tex_fig(plot_files$ate_error, "fig:robustness-ate-error",
          "Absolute error of the all-or-nothing ATE estimator under the robustness scenarios. The reference truth is known from the simulation parameters."),
  tex_fig(plot_files$support, "fig:robustness-support",
          "Accuracy of off-support allocation contrasts, including single-cell flips, random 50\\% flips, and the global all-treated/all-control contrast."),
  tex_fig(plot_files$kernel, "fig:robustness-kernel",
          "Kernel misspecification checks: exponential and power-law Hawkes simulations fitted with both exponential and power-law kernels."),
  tex_fig(plot_files$decay, "fig:robustness-decay",
          "Forward-simulation decay validation for coupled single-cell flips. The diagnostic records annular $|\\Delta N|$ with fitting disabled.")
)

tex_path <- file.path(tex_dir, "robustness.tex")
writeLines(tex_lines[nzchar(tex_lines)], tex_path)

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
