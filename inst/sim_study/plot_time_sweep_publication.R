#!/usr/bin/env Rscript
#
# Builds publication boxplots + optional parameter tables from the time-sweep summary.
# Default input: <repo>/output/sim_study/time_sweep_5228509_summary.rds
# Default outputs (same layout as inst/oklahoma/paper/generated/):
#   <repo>/inst/sim_study/generated/figures/*.pdf|.png
#   <repo>/inst/sim_study/generated/tab_sim_time_sweep_param_tables.tex
# Optional: --output-prefix STEM (STEM.pdf, STEM_estimated_control.pdf, STEM_param_tables.tex)
#           --figures-dir --tex-dir --rds-search-dir

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# ---------------------------------------------------------------------
# User-facing defaults (interactive-friendly)
# ---------------------------------------------------------------------
DEFAULT_TARGET_MULTIPLIERS <- c(0.1, 0.5, 1.0)
DEFAULT_Y_QUANTILES <- c(0.02, 0.98)  # for robust y-limits
# Default summary from time sweep (written under output/sim_study/ by the sweep pipeline)
DEFAULT_SUMMARY_RDS <- "time_sweep_5228509_summary.rds"
# Default figure/table basenames under inst/sim_study/generated/ (matches paper.tex)
DEFAULT_STEM_TRUE_CONTROL <- "results_true_control"
DEFAULT_STEM_ESTIMATED_CONTROL <- "results_estimated_control"
DEFAULT_TEX_PARAM_TABLES <- "tab_sim_time_sweep_param_tables.tex"

default_summary_rds <- function(repo_root = NULL) {
  root <- if (is.null(repo_root)) find_repo_root() else repo_root
  normalizePath(
    file.path(root, "output", "sim_study", DEFAULT_SUMMARY_RDS),
    winslash = "/",
    mustWork = FALSE
  )
}

default_sim_figures_dir <- function(repo_root = NULL) {
  root <- if (is.null(repo_root)) find_repo_root() else repo_root
  normalizePath(
    file.path(root, "inst", "sim_study", "generated", "figures"),
    winslash = "/",
    mustWork = FALSE
  )
}

default_sim_tex_dir <- function(repo_root = NULL) {
  root <- if (is.null(repo_root)) find_repo_root() else repo_root
  normalizePath(
    file.path(root, "inst", "sim_study", "generated"),
    winslash = "/",
    mustWork = FALSE
  )
}

default_rds_search_dir <- function(repo_root = NULL) {
  root <- if (is.null(repo_root)) find_repo_root() else repo_root
  normalizePath(
    file.path(root, "output", "sim_study"),
    winslash = "/",
    mustWork = FALSE
  )
}

# ---------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------
safe_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) < 1L) return(NA_real_)
  x[[1]]
}

get_arg_val <- function(args, flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) < 1L || idx[[1]] >= length(args)) return(default)
  args[[idx[[1]] + 1L]]
}

find_repo_root <- function(start_dir = getwd()) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  repeat {
    if (dir.exists(file.path(cur, ".git"))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) {
      return(normalizePath(start_dir, winslash = "/", mustWork = FALSE))
    }
    cur <- parent
  }
}

validate_summary_structure <- function(summary_obj) {
  if (is.null(summary_obj$sweep_df_true_control_core) ||
      nrow(summary_obj$sweep_df_true_control_core) < 1L) {
    stop("sweep_df_true_control_core missing or empty in input summary.")
  }
  if (is.null(summary_obj$sweep_df_core) || nrow(summary_obj$sweep_df_core) < 1L) {
    stop("sweep_df_core missing or empty in input summary.")
  }
}

prepare_core_true_control_df <- function(summary_obj, target_multipliers) {
  df <- summary_obj$sweep_df_true_control_core
  if (!"method" %in% names(df) && "labelling" %in% names(df)) {
    df$method <- ifelse(df$labelling == "SEM_full", "SEM", as.character(df$labelling))
  }

  required_cols <- c(
    "time_multiplier", "method",
    "all_nothing_true_control_per_post_time", "true_ate_per_post_time"
  )
  miss <- setdiff(required_cols, names(df))
  if (length(miss) > 0L) {
    stop("Missing columns in sweep_df_true_control_core: ", paste(miss, collapse = ", "))
  }

  out <- df %>%
    filter(vapply(.data$time_multiplier, function(z) any(abs(z - target_multipliers) < 1e-9), logical(1))) %>%
    mutate(
      method = factor(.data$method, levels = c("oracle", "naive", "SEM")),
      mult_chr = formatC(.data$time_multiplier, format = "f", digits = 1)
    )
  if (nrow(out) < 1L) stop("No rows left after multiplier filtering.")
  out
}

prepare_core_estimated_control_df <- function(summary_obj, target_multipliers) {
  df <- summary_obj$sweep_df_core
  if (!"method" %in% names(df) && "labelling" %in% names(df)) {
    df$method <- ifelse(df$labelling == "SEM_full", "SEM", as.character(df$labelling))
  }
  required_cols <- c(
    "time_multiplier", "method",
    "all_nothing_theory_per_post_time", "true_ate_per_post_time"
  )
  miss <- setdiff(required_cols, names(df))
  if (length(miss) > 0L) {
    stop("Missing columns in sweep_df_core: ", paste(miss, collapse = ", "))
  }
  out <- df %>%
    filter(vapply(.data$time_multiplier, function(z) any(abs(z - target_multipliers) < 1e-9), logical(1))) %>%
    mutate(
      method = factor(.data$method, levels = c("oracle", "naive", "SEM")),
      mult_chr = formatC(.data$time_multiplier, format = "f", digits = 1)
    )
  if (nrow(out) < 1L) stop("No rows left after multiplier filtering in sweep_df_core.")
  out
}

avg_post_points_for_multiplier <- function(mult, summary_obj, out_dir) {
  runs <- summary_obj$per_run_results
  if (is.null(runs) || length(runs) < 1L) return(NA_real_)

  keys <- names(runs)
  dist_to_target <- vapply(keys, function(k) abs(safe_num(runs[[k]]$multiplier) - mult), numeric(1))
  key <- keys[[which.min(dist_to_target)]]
  run <- runs[[key]]
  res <- run$result

  # 1) Best path: use empirical post-treatment counts if obs_data exists.
  if (!is.null(res$obs_data) && length(res$obs_data) > 0L) {
    t_star <- safe_num(res$config$TREATMENT_TIME)
    post_counts <- vapply(
      res$obs_data,
      function(z) sum(as.data.frame(z)$t >= t_star),
      numeric(1)
    )
    return(mean(post_counts, na.rm = TRUE))
  }

  # 2) If summary is slim, try loading the full per-run RDS by basename.
  if (!is.null(run$rds_path) && nzchar(run$rds_path)) {
    local_rds <- file.path(out_dir, basename(run$rds_path))
    if (file.exists(local_rds)) {
      full_res <- tryCatch(readRDS(local_rds), error = function(e) NULL)
      if (!is.null(full_res) && !is.null(full_res$obs_data) && length(full_res$obs_data) > 0L) {
        t_star <- safe_num(full_res$config$TREATMENT_TIME)
        post_counts <- vapply(
          full_res$obs_data,
          function(z) sum(as.data.frame(z)$t >= t_star),
          numeric(1)
        )
        return(mean(post_counts, na.rm = TRUE))
      }
    }
  }

  # 3) Last-resort fallback (no tile multiplier; returns total expected post points).
  cfg <- res$config
  if (is.null(cfg) || is.null(cfg$hawkes_par_1) || is.null(cfg$hawkes_par_2)) return(NA_real_)
  dt <- safe_num(cfg$END_TIME) - safe_num(cfg$TREATMENT_TIME)
  if (!is.finite(dt) || dt <= 0) return(NA_real_)
  e_ctrl <- safe_num(cfg$hawkes_par_1$mu) * dt / max(1e-8, 1 - safe_num(cfg$hawkes_par_1$K))
  e_treat <- safe_num(cfg$hawkes_par_2$mu) * dt / max(1e-8, 1 - safe_num(cfg$hawkes_par_2$K))
  e_ctrl + e_treat
}

build_multiplier_legend_info <- function(summary_obj, target_multipliers, out_dir) {
  info <- data.frame(
    time_multiplier = target_multipliers,
    avg_post_points = vapply(
      target_multipliers,
      avg_post_points_for_multiplier,
      numeric(1),
      summary_obj = summary_obj,
      out_dir = out_dir
    ),
    stringsAsFactors = FALSE
  )
  info$mult_chr <- formatC(info$time_multiplier, format = "f", digits = 1)
  info$legend_label <- ifelse(
    is.finite(info$avg_post_points),
    sprintf("%sx (~%.0f points)", info$mult_chr, info$avg_post_points),
    sprintf("%sx", info$mult_chr)
  )
  # Shorter than legend_label but keeps ~N (post-treatment points).
  info$legend_label_compact <- ifelse(
    is.finite(info$avg_post_points),
    sprintf("%s×, ~%.0f", info$mult_chr, info$avg_post_points),
    sprintf("%s×", info$mult_chr)
  )
  info
}

build_publication_plot <- function(
    df,
    legend_info,
    y_col,
    y_quantiles = DEFAULT_Y_QUANTILES,
    coord_ylim = NULL,
    font_scale = 1,
    legend_fill_col = "legend_label",
    fill_title = "Post treatment window",
    true_dtaite_label_y = NULL,
    true_dtaite_label_x = 3.8,
    true_dtaite_label_hjust = 1,
    true_dtaite_label_vjust = -0.5) {
  join_cols <- c("time_multiplier", legend_fill_col)
  df <- df %>% left_join(legend_info[, join_cols, drop = FALSE], by = "time_multiplier")
  ord <- order(legend_info$time_multiplier)
  lvl <- legend_info[[legend_fill_col]][ord]
  df$legend_fill <- factor(df[[legend_fill_col]], levels = lvl)

  true_ate <- mean(df$true_ate_per_post_time, na.rm = TRUE)
  y_vals <- df[[y_col]]
  if (is.null(coord_ylim)) {
    y_q <- as.numeric(stats::quantile(y_vals, probs = y_quantiles, na.rm = TRUE))
    y_pad <- max(0.02, 0.08 * diff(y_q))
    y_lims <- c(y_q[[1]] - y_pad, y_q[[2]] + y_pad)
  } else {
    y_lims <- coord_ylim
  }

  bs <- 18.2 * font_scale
  ax_txt <- 15.4 * font_scale
  ann_size <- 5.88 * font_scale
  out_sz <- 0.98 * font_scale
  ann_y <- if (!is.null(true_dtaite_label_y)) {
    true_dtaite_label_y
  } else if (is.null(coord_ylim)) {
    true_ate + 0.01
  } else {
    y_lims[[2]] - 0.12 * diff(y_lims)
  }

  ggplot(
    df,
    aes(x = .data$method, y = .data[[y_col]], fill = .data$legend_fill)
  ) +
    geom_boxplot(
      position = position_dodge2(width = 0.8, preserve = "single"),
      outlier.alpha = 0.5, outlier.size = out_sz
    ) +
    geom_hline(yintercept = true_ate, linetype = "dashed", linewidth = 1.25, colour = "black") +
    annotate(
      "text",
      x = true_dtaite_label_x,
      y = ann_y,
      label = "True DTAITE",
      hjust = true_dtaite_label_hjust,
      vjust = true_dtaite_label_vjust,
      size = ann_size,
      fontface = "bold"
    ) +
    coord_cartesian(ylim = y_lims) +
    labs(x = NULL, y = "All-Nothing DTAITE", fill = fill_title) +
    theme_minimal(base_size = bs) +
    theme(
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.text = element_text(size = ax_txt),
      axis.title.y = element_text(size = bs),
      legend.title = element_text(size = bs),
      legend.text = element_text(size = ax_txt),
      legend.position = "right"
    )
}

collect_param_summary <- function(summary_obj, target_multipliers, process = c("treated", "control")) {
  process <- match.arg(process)
  runs <- summary_obj$per_run_results
  if (is.null(runs) || length(runs) < 1L) return(data.frame())
  keys <- names(runs)
  out <- lapply(target_multipliers, function(mult) {
    d <- vapply(keys, function(k) abs(safe_num(runs[[k]]$multiplier) - mult), numeric(1))
    key <- keys[[which.min(d)]]
    res <- runs[[key]]$result
    tab <- if (identical(process, "treated")) res$treated_param_summary else res$control_param_summary
    if (is.null(tab) || nrow(tab) < 1L) return(NULL)
    tab <- as.data.frame(tab)
    tab$time_multiplier <- mult
    tab
  })
  out <- Filter(Negate(is.null), out)
  if (length(out) < 1L) return(data.frame())
  dplyr::bind_rows(out)
}

latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x
}

build_param_tables_tex_lines <- function(summary_obj, target_multipliers) {
  treated_df <- collect_param_summary(summary_obj, target_multipliers, process = "treated")
  control_df <- collect_param_summary(summary_obj, target_multipliers, process = "control")
  if (nrow(treated_df) < 1L && nrow(control_df) < 1L) return(character(0))

  fmt_num <- function(x, d = 3) ifelse(is.finite(x), sprintf(paste0("%.", d, "f"), x), "NA")
  fmt_mean_se <- function(mean_x, sd_x, n_x, d = 3) {
    mean_v <- as.numeric(mean_x)
    sd_v <- as.numeric(sd_x)
    n_v <- as.numeric(n_x)
    se_v <- ifelse(is.finite(sd_v) & is.finite(n_v) & n_v > 0, sd_v / sqrt(n_v), NA_real_)
    ifelse(
      is.finite(mean_v) & is.finite(se_v),
      paste0(fmt_num(mean_v, d), " (", fmt_num(se_v, d), ")"),
      "NA"
    )
  }

  to_display_rows <- function(df) {
    if (nrow(df) < 1L) return(data.frame())
    max_mult <- max(as.numeric(df$time_multiplier), na.rm = TRUE)
    x <- df %>%
      dplyr::filter(abs(as.numeric(.data$time_multiplier) - max_mult) < 1e-9) %>%
      dplyr::mutate(
        method = dplyr::case_when(
          .data$labelling == "oracle" ~ "Oracle",
          .data$labelling == "naive" ~ "Naive",
          .data$labelling %in% c("SEM_full", "SEM_adaptive") ~ "SEM",
          TRUE ~ NA_character_
        ),
        sem_priority = dplyr::case_when(
          .data$labelling == "SEM_full" ~ 1L,
          .data$labelling == "SEM_adaptive" ~ 2L,
          TRUE ~ 3L
        )
      ) %>%
      dplyr::filter(!is.na(.data$method)) %>%
      dplyr::arrange(.data$method, .data$sem_priority) %>%
      dplyr::group_by(.data$method) %>%
      dplyr::slice(1L) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        mu_col = fmt_mean_se(.data$mu_mean, .data$mu_sd, .data$n),
        alpha_col = fmt_mean_se(.data$alpha_mean, .data$alpha_sd, .data$n),
        beta_col = fmt_mean_se(.data$beta_mean, .data$beta_sd, .data$n),
        K_col = fmt_mean_se(.data$K_mean, .data$K_sd, .data$n),
        method_order = match(.data$method, c("Oracle", "Naive", "SEM"))
      ) %>%
      dplyr::arrange(.data$method_order)
    x
  }

  build_rows <- function(df) {
    x <- to_display_rows(df)
    if (nrow(x) < 1L) return(character(0))
    apply(x, 1, function(r) {
      paste(
        latex_escape(r[["method"]]),
        r[["mu_col"]],
        r[["alpha_col"]],
        r[["beta_col"]],
        r[["K_col"]],
        sep = " & "
      )
    })
  }

  treated_rows <- build_rows(treated_df)
  control_rows <- build_rows(control_df)

  c(
    "% Auto-generated by inst/sim_study/plot_time_sweep_publication.R",
    "Figures \\ref{fig:treated_params} and \\ref{fig:control_params} summarise the treated and control parameter estimates across simulations respectively.",
    "",
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{Treated parameter mean estimates across simulations for largest post treatment time window, with standard errors in parentheses.}",
    "\\label{tab:treated_params_summary}",
    "\\begin{tabular}{lrrrr}",
    "\\toprule",
    "Method & $\\mu$ mean (SE) & $\\alpha$ mean (SE) & $\\beta$ mean (SE) & $K$ mean (SE) \\\\",
    "\\midrule",
    if (length(treated_rows) > 0) paste0(treated_rows, " \\\\") else "No data \\\\",
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    "",
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{Control parameter mean estimates across simulations for largest post treatment time window, with standard errors in parentheses.}",
    "\\label{tab:control_params_summary}",
    "\\begin{tabular}{lrrrr}",
    "\\toprule",
    "Method & $\\mu$ mean (SE) & $\\alpha$ mean (SE) & $\\beta$ mean (SE) & $K$ mean (SE) \\\\",
    "\\midrule",
    if (length(control_rows) > 0) paste0(control_rows, " \\\\") else "No data \\\\",
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}",
    ""
  )
}

print_param_tables_tex <- function(summary_obj, target_multipliers) {
  lines <- build_param_tables_tex_lines(summary_obj, target_multipliers)
  if (length(lines) < 1L) {
    message("No treated/control parameter summary rows found.")
    return(invisible(NULL))
  }
  cat(paste(lines, collapse = "\n"), "\n")
  invisible(lines)
}

write_param_tables_tex <- function(summary_obj, target_multipliers, tex_file) {
  lines <- build_param_tables_tex_lines(summary_obj, target_multipliers)
  if (length(lines) < 1L) return(invisible(NULL))
  writeLines(lines, con = tex_file, useBytes = TRUE)
}

# ---------------------------------------------------------------------
# Main entry point (call this interactively)
# ---------------------------------------------------------------------
run_publication_plot <- function(
    input_rds = NULL,
    output_prefix = NULL,
    target_multipliers = DEFAULT_TARGET_MULTIPLIERS,
    y_quantiles = DEFAULT_Y_QUANTILES,
    figures_dir = NULL,
    tex_dir = NULL,
    rds_search_dir = NULL) {
  repo_root <- find_repo_root()
  rds_path <- if (is.null(input_rds) || !nzchar(input_rds)) {
    default_summary_rds(repo_root)
  } else {
    normalizePath(input_rds, winslash = "/", mustWork = FALSE)
  }
  fig_dir <- if (is.null(figures_dir)) default_sim_figures_dir(repo_root) else normalizePath(figures_dir, winslash = "/", mustWork = FALSE)
  tx_dir <- if (is.null(tex_dir)) default_sim_tex_dir(repo_root) else normalizePath(tex_dir, winslash = "/", mustWork = FALSE)
  search_dir <- if (is.null(rds_search_dir)) default_rds_search_dir(repo_root) else normalizePath(rds_search_dir, winslash = "/", mustWork = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tx_dir, recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(rds_path)) {
    stop("Input RDS not found: ", rds_path, call. = FALSE)
  }
  summary_obj <- readRDS(rds_path)
  validate_summary_structure(summary_obj)

  df <- prepare_core_true_control_df(summary_obj, target_multipliers = target_multipliers)
  df_est <- prepare_core_estimated_control_df(summary_obj, target_multipliers = target_multipliers)
  legend_info <- build_multiplier_legend_info(
    summary_obj = summary_obj,
    target_multipliers = target_multipliers,
    out_dir = search_dir
  )
  message("Read summary: ", rds_path)
  p_true <- build_publication_plot(
    df,
    legend_info,
    y_col = "all_nothing_true_control_per_post_time",
    y_quantiles = y_quantiles,
    coord_ylim = c(-0.35, -0.05),
    font_scale = 2,
    legend_fill_col = "legend_label_compact",
    fill_title = "Post window",
    true_dtaite_label_y = -0.35,
    true_dtaite_label_x = 2,
    true_dtaite_label_hjust = 0.5,
    true_dtaite_label_vjust = -0.5
  )
  p_est <- build_publication_plot(df_est, legend_info, y_col = "all_nothing_theory_per_post_time", y_quantiles = y_quantiles)

  use_custom <- !is.null(output_prefix) && nzchar(output_prefix)
  stem_true <- if (use_custom) output_prefix else DEFAULT_STEM_TRUE_CONTROL
  stem_est <- if (use_custom) paste0(output_prefix, "_estimated_control") else DEFAULT_STEM_ESTIMATED_CONTROL
  tex_base <- if (use_custom) paste0(output_prefix, "_param_tables.tex") else DEFAULT_TEX_PARAM_TABLES

  png_file_true <- file.path(fig_dir, paste0(stem_true, ".png"))
  pdf_file_true <- file.path(fig_dir, paste0(stem_true, ".pdf"))
  png_file_est <- file.path(fig_dir, paste0(stem_est, ".png"))
  pdf_file_est <- file.path(fig_dir, paste0(stem_est, ".pdf"))
  tex_file <- file.path(tx_dir, tex_base)
  ggsave(png_file_true, p_true, width = 12.6, height = 7.56, dpi = 320)
  ggsave(pdf_file_true, p_true, width = 12.6, height = 7.56)
  ggsave(png_file_est, p_est, width = 12.6, height = 7.56, dpi = 320)
  ggsave(pdf_file_est, p_est, width = 12.6, height = 7.56)
  write_param_tables_tex(summary_obj, target_multipliers = target_multipliers, tex_file = tex_file)

  message("Wrote: ", png_file_true)
  message("Wrote: ", pdf_file_true)
  message("Wrote: ", png_file_est)
  message("Wrote: ", pdf_file_est)
  message("Wrote: ", tex_file)
  message("Legend info:")
  print(legend_info[, c("time_multiplier", "avg_post_points", "legend_label")], row.names = FALSE)
  invisible(list(
    plot_true_control = p_true,
    plot_estimated_control = p_est,
    legend_info = legend_info,
    png_true_control = png_file_true,
    pdf_true_control = pdf_file_true,
    png_estimated_control = png_file_est,
    pdf_estimated_control = pdf_file_est,
    tex = tex_file
  ))
}

run_param_tables_only <- function(
    input_rds = NULL,
    target_multipliers = DEFAULT_TARGET_MULTIPLIERS,
    output_prefix = NULL,
    tex_dir = NULL,
    print_console = TRUE) {
  repo_root <- find_repo_root()
  rds_path <- if (is.null(input_rds) || !nzchar(input_rds)) {
    default_summary_rds(repo_root)
  } else {
    normalizePath(input_rds, winslash = "/", mustWork = FALSE)
  }
  tx_dir <- if (is.null(tex_dir)) default_sim_tex_dir(repo_root) else normalizePath(tex_dir, winslash = "/", mustWork = FALSE)
  dir.create(tx_dir, recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(rds_path)) stop("Input RDS not found: ", rds_path, call. = FALSE)
  message("Read summary: ", rds_path)
  summary_obj <- readRDS(rds_path)
  validate_summary_structure(summary_obj)

  use_custom <- !is.null(output_prefix) && nzchar(output_prefix)
  tex_base <- if (use_custom) paste0(output_prefix, "_param_tables.tex") else DEFAULT_TEX_PARAM_TABLES
  tex_file <- file.path(tx_dir, tex_base)
  write_param_tables_tex(summary_obj, target_multipliers = target_multipliers, tex_file = tex_file)
  if (isTRUE(print_console)) {
    print_param_tables_tex(summary_obj, target_multipliers = target_multipliers)
  }
  message("Wrote: ", tex_file)
  invisible(tex_file)
}

# ---------------------------------------------------------------------
# CLI mode (Rscript ...)
# ---------------------------------------------------------------------
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  repo_root <- find_repo_root()
  default_in <- default_summary_rds(repo_root)
  input_rds <- get_arg_val(args, "--input", default_in)
  op_raw <- get_arg_val(args, "--output-prefix", "")
  output_prefix <- if (nzchar(op_raw)) op_raw else NULL
  fig_dir_arg <- get_arg_val(args, "--figures-dir", "")
  tex_dir_arg <- get_arg_val(args, "--tex-dir", "")
  rds_search_arg <- get_arg_val(args, "--rds-search-dir", "")
  figures_dir <- if (nzchar(fig_dir_arg)) fig_dir_arg else NULL
  tex_dir <- if (nzchar(tex_dir_arg)) tex_dir_arg else NULL
  rds_search_dir <- if (nzchar(rds_search_arg)) rds_search_arg else NULL
  mode <- tolower(get_arg_val(args, "--mode", "both"))
  if (mode == "plot") {
    run_publication_plot(
      input_rds = input_rds,
      output_prefix = output_prefix,
      figures_dir = figures_dir,
      tex_dir = tex_dir,
      rds_search_dir = rds_search_dir
    )
  } else if (mode == "tables") {
    run_param_tables_only(
      input_rds = input_rds,
      output_prefix = output_prefix,
      tex_dir = tex_dir,
      print_console = TRUE
    )
  } else {
    run_publication_plot(
      input_rds = input_rds,
      output_prefix = output_prefix,
      figures_dir = figures_dir,
      tex_dir = tex_dir,
      rds_search_dir = rds_search_dir
    )
  }
}
