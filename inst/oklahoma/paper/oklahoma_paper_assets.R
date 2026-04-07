# Oklahoma application figures and LaTeX table fragments for the paper.
#
# Usage (from repo root):
#   Rscript inst/oklahoma/paper/oklahoma_paper_assets.R
#   (defaults to output/oklahoma/for_paper.rds, then inst/oklahoma/paper/for_paper.rds)
#   Rscript inst/oklahoma/paper/oklahoma_paper_assets.R --input /path/to/other.rds
#
# Or in R: setwd("<repo>"); source("inst/oklahoma/paper/oklahoma_paper_assets.R")
#
# Writes (default: figures next to table fragments under generated/):
#   inst/oklahoma/paper/generated/figures/*.pdf:
#     ATE_diff, ok_bootstrap_*, ok_sim_* — require bootstrap_ate in RDS
#     cumulative_count — events CSV + jsonlite
#     ok_partition_county, ok_point_patterns — sf + tigris + AOI geojson under --data-dir
#     ok_sem_f_confusion_trace, ok_sem_f_metric_trace, ok_sem_f_flips — fits_named$F$fit$adaptive
#   inst/oklahoma/paper/generated/*.tex:
#     tab_ok_counts, tab_ok_sem_config, tab_ok_confusion_F, tab_ok_ef_params
#     tab_ok_bootstrap_summary — only if bootstrap present
#   Override PDF location: --plots-dir <path>
#
# Optional packages: jsonlite (cumulative plot), tigris (county maps; same as oklahoma_report.qmd).
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(data.table)
})

find_repo_root <- function(start_dir = getwd()) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(cur, ".git"))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) return(normalizePath(start_dir, winslash = "/", mustWork = TRUE))
    cur <- parent
  }
}

get_arg_val <- function(args, flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) < 1L || idx[[1]] >= length(args)) return(default)
  args[[idx[[1]] + 1L]]
}

args <- commandArgs(trailingOnly = TRUE)
repo_root <- find_repo_root()
input_arg <- get_arg_val(args, "--input", NULL)
if (!is.null(input_arg) && nzchar(input_arg)) {
  input_rds <- normalizePath(input_arg, winslash = "/", mustWork = TRUE)
} else {
  input_candidates <- c(
    file.path(repo_root, "output", "oklahoma", "for_paper.rds"),
    file.path(repo_root, "inst", "oklahoma", "paper", "for_paper.rds")
  )
  input_hit <- input_candidates[file.exists(input_candidates)][1L]
  if (is.na(input_hit)) {
    stop(
      "No results RDS found. Place for_paper.rds in output/oklahoma/ or inst/oklahoma/paper/, ",
      "or pass --input /path/to/results.rds"
    )
  }
  input_rds <- normalizePath(input_hit, winslash = "/", mustWork = TRUE)
}
plots_dir <- normalizePath(
  get_arg_val(
    args,
    "--plots-dir",
    file.path(repo_root, "inst", "oklahoma", "paper", "generated", "figures")
  ),
  winslash = "/",
  mustWork = FALSE
)
tex_dir <- normalizePath(
  get_arg_val(args, "--tex-dir", file.path(repo_root, "inst", "oklahoma", "paper", "generated")),
  winslash = "/",
  mustWork = FALSE
)
data_dir <- normalizePath(
  get_arg_val(
    args,
    "--data-dir",
    file.path(repo_root, "inst", "oklahoma", "oklahoma_induced_seismicity_data_regional20150318")
  ),
  winslash = "/",
  mustWork = FALSE
)

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tex_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_rds)) stop("Input RDS not found: ", input_rds)

res <- readRDS(input_rds)
if (is.null(res$fits_named$E) || is.null(res$fits_named$F)) {
  stop("RDS must contain fits_named$E and fits_named$F.")
}

boot_obj <- res$bootstrap_ate
have_boot <- !is.null(boot_obj) &&
  !is.null(boot_obj$fit_E) && !is.null(boot_obj$fit_F) &&
  !is.null(boot_obj$fit_E$replicate_summary) && !is.null(boot_obj$fit_F$replicate_summary)
if (!have_boot) {
  message("No bootstrap_ate replicate summaries in RDS; skipping bootstrap figures and tab_ok_bootstrap_summary.")
}

ate_days <- suppressWarnings(as.numeric(res$config$ATE_WINDOW_DAYS))
if (!is.finite(ate_days)) ate_days <- 100

fit_E <- res$fits_named$E
fit_F <- res$fits_named$F

tex_escape_cfg <- function(s) {
  s <- as.character(s)
  s <- gsub("&", intToUtf8(c(92L, 38L)), s, fixed = TRUE)
  s <- gsub("%", intToUtf8(c(92L, 37L)), s, fixed = TRUE)
  s <- gsub("_", intToUtf8(c(92L, 95L)), s, fixed = TRUE)
  s <- gsub("#", intToUtf8(c(92L, 35L)), s, fixed = TRUE)
  s
}

fmt_cfg_val <- function(x, default = "---") {
  if (is.null(x) || length(x) < 1L) return(default)
  v <- x[[1L]]
  if (is.logical(v)) return(if (isTRUE(v)) "yes" else if (identical(v, FALSE)) "no" else default)
  if (is.numeric(v)) {
    if (!is.finite(v)) return(default)
    if (abs(v - round(v)) < 1e-9) return(as.character(as.integer(round(v))))
    return(as.character(signif(v, 6)))
  }
  if (is.character(v) && !nzchar(v)) return(default)
  tex_escape_cfg(as.character(v)[1])
}

# ---- Expected-count rows (same quantities as HTML report / former Rmd appendix) ----
extract_totals <- function(ate_obj) {
  out <- list(c_total = numeric(0), t_total = numeric(0), saved_total = numeric(0))
  if (is.null(ate_obj) || is.null(ate_obj$all_nothing_sim)) return(out)
  sim_df <- as.data.frame(ate_obj$all_nothing_sim)
  if ("c_total" %in% names(sim_df)) out$c_total <- suppressWarnings(as.numeric(sim_df$c_total))
  if ("t_total" %in% names(sim_df)) out$t_total <- suppressWarnings(as.numeric(sim_df$t_total))
  if ("total_saved" %in% names(sim_df)) out$saved_total <- suppressWarnings(as.numeric(sim_df$total_saved))
  if (length(out$saved_total) < 1L && length(out$c_total) > 0 && length(out$t_total) > 0) {
    out$saved_total <- out$c_total - out$t_total
  }
  out
}

summ_one <- function(fit_obj) {
  tt <- extract_totals(if (!is.null(fit_obj)) fit_obj$ate else NULL)
  c_mean <- if (length(tt$c_total) > 0) mean(tt$c_total, na.rm = TRUE) else NA_real_
  t_mean <- if (length(tt$t_total) > 0) mean(tt$t_total, na.rm = TRUE) else NA_real_
  saved_mean <- if (length(tt$saved_total) > 0) mean(tt$saved_total, na.rm = TRUE) else NA_real_
  c(round(c_mean, 0), round(t_mean, 0), round(saved_mean, 0))
}

vals_E <- summ_one(fit_E)
vals_F <- summ_one(fit_F)

write_tex_ok_counts <- function(path) {
  fmt <- function(x) if (is.finite(x)) as.character(as.integer(x)) else "---"
  lines <- c(
    "% Auto-generated by oklahoma_paper_assets.R — do not edit by hand",
    "\\begin{table}",
    "\\caption{\\label{tab:ok_counts}Observed and model-implied all-or-nothing expected counts naive and SEM approach.}",
    "\\centering",
    "\\begin{tabular}[t]{lrrr}",
    "\\toprule",
    "Method & $E[N_C]$ & $E[N_A]$ & $\\Delta=E[N_C]-E[N_A]$\\\\",
    "\\midrule",
    sprintf("Naive & %s & %s & %s\\\\", fmt(vals_E[1]), fmt(vals_E[2]), fmt(vals_E[3])),
    sprintf("SEM & %s & %s & %s\\\\", fmt(vals_F[1]), fmt(vals_F[2]), fmt(vals_F[3])),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
  writeLines(lines, path)
  message("Wrote ", path)
}

# Figures included at column width in paper.tex are scaled down for legibility (13 pt ggplot default as reference).
theme_ok_boot_hist <- theme_minimal(base_size = 13 * 1.5) +
  theme(
    legend.position = "top",
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    strip.text = element_text(face = "bold")
  )
theme_paper_shrunk <- theme_minimal(base_size = 13 * 1.8) +
  theme(
    legend.position = "top",
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    strip.text = element_text(face = "bold")
  )

# ---- Bootstrap (bias-corrected replicates, same as oklahoma_publication_ef.R) ----
if (have_boot) {
boot_E <- as.data.frame(boot_obj$fit_E$replicate_summary)
boot_F <- as.data.frame(boot_obj$fit_F$replicate_summary)
boot_E$model <- "Naive Bivaraite ETAS with KDE"
boot_F$model <- "SEM Bivariate ETAS with KDE"
boot_df <- bind_rows(boot_E, boot_F)

if (!"ate_total_mean" %in% names(boot_df)) {
  stop("Expected bootstrap column ate_total_mean is missing.")
}
boot_df <- boot_df %>% filter(is.finite(.data$ate_total_mean))

get_total_effect_sim <- function(ate_obj) {
  if (is.null(ate_obj) || is.null(ate_obj$all_nothing_sim)) return(numeric(0))
  sim_df <- as.data.frame(ate_obj$all_nothing_sim)
  if ("total_saved" %in% names(sim_df)) return(suppressWarnings(as.numeric(sim_df$total_saved)))
  if ("total_effect" %in% names(sim_df)) return(suppressWarnings(as.numeric(sim_df$total_effect)))
  numeric(0)
}
obs_ate_mean <- function(fit_obj) {
  if (is.null(fit_obj) || is.null(fit_obj$ate)) return(NA_real_)
  vals <- get_total_effect_sim(fit_obj$ate)
  vals <- vals[is.finite(vals)]
  if (length(vals) < 1L) return(NA_real_)
  mean(vals, na.rm = TRUE)
}
recenter_boot <- function(boot_vals, target_mean) {
  bv <- suppressWarnings(as.numeric(boot_vals))
  if (!is.finite(target_mean) || is.na(target_mean)) return(bv)
  bmean <- mean(bv, na.rm = TRUE)
  if (!is.finite(bmean) || is.na(bmean)) return(bv)
  bv + (target_mean - bmean)
}

obs_E <- obs_ate_mean(fit_E)
obs_F <- obs_ate_mean(fit_F)
target_mean_by_model <- c(
  "Naive Bivaraite ETAS with KDE" = obs_E,
  "SEM Bivariate ETAS with KDE" = obs_F
)
boot_df <- boot_df %>%
  group_by(.data$model) %>%
  mutate(ate_total_mean = recenter_boot(.data$ate_total_mean, target_mean_by_model[.data$model[[1]]])) %>%
  ungroup() %>%
  filter(is.finite(.data$ate_total_mean))

boot_df <- boot_df %>%
  group_by(.data$model) %>%
  arrange(.data$rep, .by_group = TRUE) %>%
  mutate(
    rep_idx = row_number(),
    running_mean = cumsum(.data$ate_total_mean) / row_number()
  ) %>%
  ungroup()

sim_E <- as.data.frame(fit_E$ate$all_nothing_sim)
sim_F <- as.data.frame(fit_F$ate$all_nothing_sim)
sim_E$model <- "Naive Bivaraite ETAS with KDE"
sim_F$model <- "SEM Bivariate ETAS with KDE"
sim_df <- bind_rows(sim_E, sim_F) %>% filter(is.finite(.data$total_saved))

theme_pub <- theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    strip.text = element_text(face = "bold")
  )

p_boot_hist <- ggplot(boot_df, aes(x = .data$ate_total_mean, fill = .data$model)) +
  geom_histogram(bins = 30, alpha = 0.55, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.9) +
  scale_fill_manual(values = c(
    "Naive Bivaraite ETAS with KDE" = "#e6550d",
    "SEM Bivariate ETAS with KDE" = "#3182bd"
  )) +
  labs(
    x = "bootstrap replicate total saved / 100 day",
    y = "Count",
    fill = "Model"
  ) +
  xlim(-400, 700) +
  theme_ok_boot_hist

p_boot_ecdf <- ggplot(boot_df, aes(x = .data$ate_total_mean, colour = .data$model)) +
  stat_ecdf(linewidth = 1.0) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.9) +
  scale_colour_manual(values = c(
    "Naive Bivaraite ETAS with KDE" = "#e6550d",
    "SEM Bivariate ETAS with KDE" = "#3182bd"
  )) +
  labs(
    x = "bootstrap replicate total saved / 100 day",
    y = "Cumulative proportion",
    colour = "Model"
  ) +
  theme_pub

p_boot_running <- ggplot(boot_df, aes(x = .data$rep_idx, y = .data$running_mean, colour = .data$model)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.9) +
  facet_wrap(~model, scales = "free_y") +
  scale_colour_manual(values = c(
    "Naive Bivaraite ETAS with KDE" = "#e6550d",
    "SEM Bivariate ETAS with KDE" = "#3182bd"
  )) +
  labs(
    x = "Bootstrap replicate index",
    y = bquote("Running mean of " ~ hat(Delta)[AoN]),
    colour = "Model"
  ) +
  theme_pub +
  theme(legend.position = "none")

p_sim_hist <- ggplot(sim_df, aes(x = .data$total_saved, fill = .data$model)) +
  geom_histogram(bins = 30, alpha = 0.55, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.9) +
  scale_fill_manual(values = c(
    "Naive Bivaraite ETAS with KDE" = "#e6550d",
    "SEM Bivariate ETAS with KDE" = "#3182bd"
  )) +
  labs(
    x = bquote(hat(Delta)[AoN] ~ "(Monte Carlo total saved /" ~ .(ate_days) ~ " days)"),
    y = "Count",
    fill = "Model"
  ) +
  theme_pub

# Paper-facing model labels (used in figure legends)
ggsave(file.path(plots_dir, "ATE_diff.pdf"), p_boot_hist, width = 8.8, height = 5.2)
ggsave(file.path(plots_dir, "ok_bootstrap_ecdf.pdf"), p_boot_ecdf, width = 8.8, height = 5.2)
ggsave(file.path(plots_dir, "ok_bootstrap_running_mean.pdf"), p_boot_running, width = 8.8, height = 5.2)
ggsave(file.path(plots_dir, "ok_sim_total_saved_hist.pdf"), p_sim_hist, width = 8.8, height = 5.2)
message("Wrote bootstrap-related figures under ", plots_dir)
} else {
  boot_df <- NULL
}

# Cumulative catalog plot
ev_csv <- file.path(data_dir, "events_all.csv")
if (file.exists(ev_csv) && requireNamespace("jsonlite", quietly = TRUE)) {
  meta_ok <- jsonlite::fromJSON(file.path(data_dir, "metadata.json"), simplifyVector = TRUE)
  m0_c <- res$config$ETAS_M0
  if (is.null(m0_c) || !length(m0_c) || all(is.na(m0_c))) m0_c <- meta_ok$design$min_mag
  ev_all <- fread(ev_csv)
  ev_all[, time_utc := as.POSIXct(time_utc, tz = "UTC")]
  t_star_utc <- as.POSIXct(meta_ok$design$t_star_utc, tz = "UTC")
  ev_all[, t_days := as.numeric(difftime(time_utc, t_star_utc, units = "days"))]
  setorder(ev_all, time_utc)
  ev_all <- ev_all[mag >= m0_c]
  ev_all[, cum_n := seq_len(.N)]
  t_star_lbl <- format(t_star_utc, "%d %b %Y", tz = "UTC")
  p_cum <- ggplot(ev_all, aes(x = .data$t_days, y = .data$cum_n)) +
    geom_step(color = "#2171b5", linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    annotate(
      "text",
      x = 5,
      y = max(ev_all$cum_n, na.rm = TRUE) * 0.1,
      label = paste0("OCC directive\n(", t_star_lbl, ")"),
      color = "red",
      hjust = 0,
      size = 3.5 * 1.8
    ) +
    labs(x = "Days since treatment (t*)", y = "Cumulative events") +
    theme_paper_shrunk +
    theme(legend.position = "none")
  ggsave(file.path(plots_dir, "cumulative_count.pdf"), p_cum, width = 8.8, height = 5.2)
  message("Wrote cumulative_count.pdf")
} else {
  message(
    "Skipped cumulative_count.pdf (need ", ev_csv,
    " and package jsonlite)."
  )
}

# ---- Design maps (HTML report § County partition / Point patterns) ----
coerce_pp_df <- function(obj) {
  if (is.null(obj)) return(data.frame())
  if (inherits(obj, "ppp")) {
    out <- data.frame(x = obj$x, y = obj$y)
    if (!is.null(obj$marks)) out <- cbind(out, as.data.frame(obj$marks), stringsAsFactors = FALSE)
    return(out)
  }
  if (is.data.frame(obj)) return(as.data.frame(obj))
  if (is.matrix(obj)) return(as.data.frame(obj))
  if (is.list(obj) && all(c("x", "y") %in% names(obj))) {
    n <- min(length(obj$x), length(obj$y))
    out <- data.frame(x = as.numeric(obj$x)[seq_len(n)], y = as.numeric(obj$y)[seq_len(n)])
    if ("t" %in% names(obj)) out$t <- as.numeric(obj$t)[seq_len(min(n, length(obj$t)))]
    for (nm in setdiff(names(obj), c("x", "y", "t"))) {
      v <- obj[[nm]]
      if (length(v) >= n && is.atomic(v)) out[[nm]] <- v[seq_len(n)]
    }
    return(out)
  }
  out <- tryCatch(as.data.frame(obj), error = function(e) NULL)
  if (is.null(out)) stop("Could not coerce point-pattern object to data frame.")
  out
}

cfg_paper <- if (!is.null(res$config)) res$config else list()
if (is.null(cfg_paper$ETAS_M0) || !length(cfg_paper$ETAS_M0)) {
  if (requireNamespace("jsonlite", quietly = TRUE) && file.exists(file.path(data_dir, "metadata.json"))) {
    meta_cfg <- jsonlite::fromJSON(file.path(data_dir, "metadata.json"), simplifyVector = TRUE)
    cfg_paper$ETAS_M0 <- meta_cfg$design$min_mag
  }
}

if (is.null(res$pp_data) || is.null(res$counties)) {
  support_path <- normalizePath(file.path(repo_root, "inst", "oklahoma", "oklahoma_report_support.R"), winslash = "/")
  if (!file.exists(support_path)) {
    stop("RDS lacks pp_data/counties; install support file or use a full results RDS: ", support_path)
  }
  source(support_path, local = TRUE)
  rebuilt <- oklahoma_report_rebuild_pp_and_counties(data_dir, etas_m0 = cfg_paper$ETAS_M0)
  res$pp_data <- rebuilt$pp_data
  res$counties <- rebuilt$counties
}

cty <- res$counties
pp_pre <- coerce_pp_df(res$pp_data$pp_pre)
pp_post <- coerce_pp_df(res$pp_data$pp_post)
if (!("t" %in% names(pp_pre))) pp_pre$t <- numeric(nrow(pp_pre))
if (!("t" %in% names(pp_post))) pp_post$t <- numeric(nrow(pp_post))

geojson <- file.path(data_dir, "occ_aoi_layer_2.geojson")
if (requireNamespace("sf", quietly = TRUE) && requireNamespace("tigris", quietly = TRUE) && file.exists(geojson)) {
  options(tigris_use_cache = TRUE)
  counties_sf <- tigris::counties(state = "OK", cb = TRUE, year = 2022)
  counties_sf <- sf::st_transform(counties_sf, 5070)
  counties_sf <- sf::st_make_valid(counties_sf)
  aoi_sf <- sf::st_read(geojson, quiet = TRUE)
  aoi_sf <- sf::st_transform(aoi_sf, 5070)
  aoi_sf <- sf::st_make_valid(aoi_sf)
  aoi_union <- sf::st_union(aoi_sf)
  county_centroids <- sf::st_centroid(counties_sf)
  inside_aoi <- lengths(sf::st_within(county_centroids, aoi_union)) > 0
  counties_sf$treatment <- ifelse(inside_aoi, "treated", "control")
  p_part <- ggplot() +
    geom_sf(data = counties_sf, aes(fill = .data$treatment), color = "grey40", linewidth = 0.4) +
    scale_fill_manual(values = c(control = "#deebf7", treated = "#fc9272"), name = "Assignment") +
    geom_sf(data = aoi_sf, fill = NA, color = "red", linewidth = 1, linetype = "dashed") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_blank(), plot.subtitle = element_blank())
  ggsave(file.path(plots_dir, "ok_partition_county.pdf"), p_part, width = 8.8, height = 6.2)
  message("Wrote ok_partition_county.pdf")

  counties_km <- counties_sf
  sf::st_geometry(counties_km) <- sf::st_geometry(counties_sf) / 1000
  sf::st_crs(counties_km) <- NA

  xs <- c(pp_pre$x, pp_post$x)
  ys <- c(pp_pre$y, pp_post$y)
  xs <- xs[is.finite(xs)]
  ys <- ys[is.finite(ys)]
  if (length(xs) == 0L || length(ys) == 0L) {
    bb_map <- sf::st_bbox(counties_km)
    xlims <- as.numeric(bb_map[c("xmin", "xmax")])
    ylims <- as.numeric(bb_map[c("ymin", "ymax")])
  } else {
    xlims <- range(xs, na.rm = TRUE)
    ylims <- range(ys, na.rm = TRUE)
  }

  pp_post_plot <- pp_post
  process_col <- dplyr::case_when(
    "location_process" %in% names(pp_post_plot) ~ "location_process",
    "inferred_process" %in% names(pp_post_plot) ~ "inferred_process",
    "process" %in% names(pp_post_plot) ~ "process",
    TRUE ~ NA_character_
  )
  nr_post <- nrow(pp_post_plot)
  if (nr_post == 0L) {
    pp_post_plot$Process <- factor(character(0), levels = c("control", "treated"))
  } else if (is.na(process_col)) {
    pp_post_plot$Process <- factor(rep("control", nr_post), levels = c("control", "treated"))
  } else {
    pp_post_plot$Process <- factor(as.character(pp_post_plot[[process_col]]), levels = c("control", "treated"))
  }

  period_levels <- c("Pre-treatment", "Post-treatment")
  pre_plot_df <- data.frame(
    x = pp_pre$x,
    y = pp_pre$y,
    t = pp_pre$t,
    period = factor(period_levels[[1]], levels = period_levels),
    stringsAsFactors = FALSE
  )
  post_plot_df <- data.frame(
    x = pp_post_plot$x,
    y = pp_post_plot$y,
    t = pp_post_plot$t,
    Process = pp_post_plot$Process,
    period = factor(period_levels[[2]], levels = period_levels),
    stringsAsFactors = FALSE
  )

  p_point_patterns <- ggplot() +
    geom_sf(
      data = counties_km,
      mapping = aes(fill = .data$treatment),
      colour = "grey60",
      linewidth = 0.25,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = pre_plot_df,
      mapping = aes(x = .data$x, y = .data$y, alpha = .data$t),
      colour = "#2166ac",
      size = 0.8,
      inherit.aes = FALSE,
      na.rm = TRUE
    ) +
    geom_point(
      data = post_plot_df,
      mapping = aes(x = .data$x, y = .data$y, colour = .data$Process, alpha = .data$t),
      size = 0.8,
      inherit.aes = FALSE,
      na.rm = TRUE
    ) +
    scale_fill_manual(
      name = "Assignment",
      values = c(control = "#deebf7", treated = "#fc9272")
    ) +
    scale_colour_manual(
      name = "Process",
      values = c(control = "#2166ac", treated = "#b2182b")
    ) +
    scale_alpha_continuous(name = "Time (days)", range = c(0.2, 1), na.value = NA) +
    facet_wrap(~period, nrow = 1L, drop = FALSE) +
    coord_sf(xlim = xlims, ylim = ylims, datum = NA, expand = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      panel.spacing.x = grid::unit(3, "mm"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.spacing.y = grid::unit(3, "pt"),
      plot.margin = margin(4, 4, 4, 4)
    )

  # coord_sf fixes geographic aspect; if the PDF page is too wide for the chosen
  # height, ggplot adds large left/right gutters. Size from bbox and two columns.
  x_span <- abs(diff(xlims))
  y_span <- abs(diff(ylims))
  fig_w_in <- 14
  if (!is.finite(x_span) || x_span <= 0 || !is.finite(y_span) || y_span <= 0) {
    fig_h_in <- 5.8
  } else {
    geo_asp <- x_span / y_span
    inner_w <- fig_w_in - 0.65
    facet_gap_in <- 0.18
    h_panel_in <- (inner_w - facet_gap_in) / (2 * geo_asp)
    fig_h_in <- h_panel_in + 1.05
    fig_h_in <- min(max(fig_h_in, 4.5), 10)
  }

  ggsave(
    filename = file.path(plots_dir, "ok_point_patterns.pdf"),
    plot = p_point_patterns,
    width = fig_w_in,
    height = fig_h_in
  )
  message("Wrote ok_point_patterns.pdf")
} else {
  message(
    "Skipped ok_partition_county.pdf / ok_point_patterns.pdf (need packages sf and tigris, and ",
    geojson, ")."
  )
}

# ---- SEM diagnostics: model F only (confusion trace, score trace, flips) ----
fit_F_blk <- res$fits_named$F
sem_fit_F <- if (!is.null(fit_F_blk) && !is.null(fit_F_blk$fit)) fit_F_blk$fit else NULL

get_sem_metric_trace <- function(fit_obj) {
  if (is.null(fit_obj)) return(NULL)
  sem_obj <- if (!is.null(fit_obj$fit) && !is.null(fit_obj$fit$adaptive)) {
    fit_obj$fit
  } else if (!is.null(fit_obj$sem) && !is.null(fit_obj$sem$adaptive)) {
    fit_obj$sem
  } else {
    NULL
  }
  if (is.null(sem_obj) || is.null(sem_obj$adaptive)) return(NULL)
  ad <- sem_obj$adaptive
  if (!is.null(ad$metrics) && length(ad$metrics) > 0) return(as.numeric(ad$metrics))
  if (!is.null(ad$all_metrics) && length(ad$all_metrics) > 0) {
    vals <- vapply(ad$all_metrics, function(v) {
      if (is.null(v) || length(v) == 0) return(NA_real_)
      max(as.numeric(v), na.rm = TRUE)
    }, numeric(1))
    vals[is.infinite(vals)] <- NA_real_
    return(vals)
  }
  NULL
}

extract_selected_confusion_trace <- function(sem_obj, model_label) {
  ad <- sem_obj$adaptive
  if (is.null(ad) || is.null(ad$class_results) || is.null(ad$all_metrics)) return(NULL)
  if (length(ad$class_results) < 1 || length(ad$all_metrics) < 1) return(NULL)
  offset <- 0L
  rows <- list()
  for (i in seq_along(ad$all_metrics)) {
    metric_i <- ad$all_metrics[[i]]
    n_i <- length(metric_i)
    if (n_i < 1L) next
    best_j <- which.max(as.numeric(metric_i))
    class_idx <- offset + best_j
    if (class_idx > length(ad$class_results)) break
    cm <- ad$class_results[[class_idx]]
    if (is.null(cm)) {
      offset <- offset + n_i
      next
    }
    rows[[length(rows) + 1L]] <- data.frame(
      iteration = i,
      Model = model_label,
      control_as_control = as.numeric(if (!is.null(cm$control_as_control)) cm$control_as_control else 0),
      control_as_treated = as.numeric(if (!is.null(cm$control_as_treated)) cm$control_as_treated else 0),
      treated_as_control = as.numeric(if (!is.null(cm$treated_as_control)) cm$treated_as_control else 0),
      treated_as_treated = as.numeric(if (!is.null(cm$treated_as_treated)) cm$treated_as_treated else 0)
    )
    offset <- offset + n_i
  }
  if (length(rows) < 1L) return(NULL)
  do.call(rbind, rows)
}

if (!is.null(sem_fit_F) && !is.null(sem_fit_F$adaptive)) {
  trace_df <- extract_selected_confusion_trace(sem_fit_F, "F: SEM Biv+KDE (all-free)")
  if (!is.null(trace_df) && nrow(trace_df) > 0) {
    trace_long <- pivot_longer(
      trace_df,
      cols = c(control_as_control, control_as_treated, treated_as_control, treated_as_treated),
      names_to = "entry",
      values_to = "count"
    )
    trace_long$entry <- factor(
      trace_long$entry,
      levels = c("control_as_control", "control_as_treated", "treated_as_control", "treated_as_treated"),
      labels = c("Control->Control", "Control->Treated", "Treated->Control", "Treated->Treated")
    )
    p_conf <- ggplot(trace_long, aes(x = .data$iteration, y = .data$count, colour = .data$entry)) +
      geom_line(linewidth = 0.9) +
      scale_color_manual(
        name = "Confusion entry",
        values = c(
          "Control->Control" = "#1f78b4",
          "Control->Treated" = "#e31a1c",
          "Treated->Control" = "#ff7f00",
          "Treated->Treated" = "#33a02c"
        )
      ) +
      labs(x = "Inner iteration", y = "Count") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "top", plot.title = element_blank(), plot.subtitle = element_blank())
    ggsave(file.path(plots_dir, "ok_sem_f_confusion_trace.pdf"), p_conf, width = 8.8, height = 5.2)
    message("Wrote ok_sem_f_confusion_trace.pdf")
  } else {
    message("Skipped ok_sem_f_confusion_trace.pdf (no class_results trace).")
  }

  tr_m <- get_sem_metric_trace(fit_F_blk)
  if (!is.null(tr_m) && length(tr_m) > 0) {
    tr_m <- as.numeric(tr_m)
    ok <- is.finite(tr_m)
    dfm <- data.frame(Iteration = seq_along(tr_m)[ok], LogLik = tr_m[ok])
    if (nrow(dfm) > 0) {
      dfm <- dfm %>% arrange(.data$Iteration) %>% mutate(RunningBest = cummax(.data$LogLik))
      last_lbl <- sprintf("final = %.2f", dfm$LogLik[nrow(dfm)])
      p_met <- ggplot(dfm, aes(x = .data$Iteration, y = .data$LogLik)) +
        geom_line(linewidth = 0.9, colour = "#1f78b4") +
        geom_line(aes(y = .data$RunningBest), linewidth = 0.7, linetype = "dashed", alpha = 0.7, colour = "#1f78b4") +
        annotate("text", x = Inf, y = -Inf, label = last_lbl, hjust = 1.05, vjust = -0.5, size = 3.2) +
        labs(
          x = "Adaptive SEM iteration",
          y = "Selected score (solid) / running best (dashed)"
        ) +
        theme_minimal(base_size = 12) +
        theme(plot.title = element_blank(), plot.subtitle = element_blank())
      ggsave(file.path(plots_dir, "ok_sem_f_metric_trace.pdf"), p_met, width = 8.8, height = 4.8)
      message("Wrote ok_sem_f_metric_trace.pdf")
    }
  } else {
    message("Skipped ok_sem_f_metric_trace.pdf (no metrics trace).")
  }

  adF <- sem_fit_F$adaptive
  if (!is.null(adF$average_flips) && !is.null(adF$max_metric_flips)) {
    n <- min(length(adF$average_flips), length(adF$max_metric_flips))
    if (n >= 1L) {
      dff <- data.frame(
        Iteration = seq_len(n),
        Average = as.numeric(adF$average_flips[seq_len(n)]),
        Accepted = as.numeric(adF$max_metric_flips[seq_len(n)])
      )
      dffl <- pivot_longer(dff, cols = c(Average, Accepted), names_to = "series", values_to = "flips")
      p_fl <- ggplot(dffl, aes(x = .data$Iteration, y = .data$flips, colour = .data$series)) +
        geom_line(linewidth = 0.9) +
        scale_colour_manual(values = c(Average = "#6a3d9a", Accepted = "#e31a1c"), name = NULL) +
        labs(x = "Iteration", y = "Flips") +
        theme_minimal(base_size = 12) +
        theme(legend.position = "top", plot.title = element_blank(), plot.subtitle = element_blank())
      ggsave(file.path(plots_dir, "ok_sem_f_flips.pdf"), p_fl, width = 8.8, height = 4.8)
      message("Wrote ok_sem_f_flips.pdf")
    }
  } else {
    message("Skipped ok_sem_f_flips.pdf (no flip vectors).")
  }
} else {
  message("Skipped SEM F PDFs (fits_named$F$fit or adaptive block missing).")
}

fmt_scaled_num <- function(x, digits = 5) {
  if (is.null(x)) return("---")
  xv <- suppressWarnings(as.numeric(unlist(x, use.names = FALSE)))
  xv <- xv[is.finite(xv)]
  if (length(xv) < 1L) return("---")
  format(round(xv[[1]], digits), scientific = FALSE, trim = TRUE)
}

to_named_num <- function(x) {
  if (is.null(x)) return(setNames(numeric(0), character(0)))
  if (is.list(x)) x <- unlist(x, use.names = TRUE)
  v <- suppressWarnings(as.numeric(x))
  nms <- names(x)
  if (is.null(nms)) nms <- rep("", length(v))
  names(v) <- nms
  v
}

build_cm_sem <- function(sem_fit) {
  if (is.null(sem_fit) || is.null(sem_fit$adaptive)) return(NULL)
  lab <- sem_fit$adaptive$adaptive_labelling
  if (is.null(lab)) return(NULL)
  post <- lab[lab$t >= 0, , drop = FALSE]
  if (is.null(post$location_process) || is.null(post$inferred_process)) return(NULL)
  tab <- table(Location = post$location_process, Inferred = post$inferred_process)
  as.data.frame.matrix(tab)
}

write_tex_ok_confusion_F <- function(cm, path, n_post) {
  if (is.null(cm) || nrow(cm) < 1L) {
    lines <- c(
      "% Auto-generated by oklahoma_paper_assets.R — do not edit by hand",
      "\\begin{table}[t]", "\\centering", "\\small",
      sprintf(
        "\\caption{\\label{tab:ok_confusion_F}Model F (SEM KDE, all-free): location vs.\\ SEM-inferred label ($n=%d$ post-treatment events).}",
        as.integer(n_post)
      ),
      "\\begin{tabular}{@{}lcc@{}}", "\\toprule", "\\multicolumn{3}{c}{---}\\\\",
      "\\bottomrule", "\\end{tabular}", "\\end{table}"
    )
    writeLines(lines, path)
    message("Wrote placeholder ", path)
    return(invisible())
  }
  rn <- rownames(cm)
  cn <- colnames(cm)
  ncol_cm <- length(cn)
  colspec <- paste0("l", paste(rep("c", ncol_cm), collapse = ""))
  hdr_inferred <- paste(vapply(cn, function(cj) tex_escape_cfg(as.character(cj)), character(1)), collapse = " & ")
  lines <- c(
    "% Auto-generated by oklahoma_paper_assets.R — do not edit by hand",
    "\\begin{table}[t]",
    "\\centering",
    "\\small",
    sprintf(
      "\\caption{\\label{tab:ok_confusion_F}Model F (SEM KDE, all-free): location vs.\\ SEM-inferred label ($n=%d$ post-treatment events).}",
      as.integer(n_post)
    ),
    sprintf("\\begin{tabular}{@{}%s@{}}", colspec),
    "\\toprule",
    sprintf("Location $\\\\downarrow$ / Inferred $\\\\rightarrow$ & %s\\\\", hdr_inferred),
    "\\midrule"
  )
  for (i in seq_along(rn)) {
    cells <- vapply(cn, function(cj) {
      v <- suppressWarnings(as.numeric(cm[rn[i], cj]))
      if (!is.finite(v)) v <- 0
      pct <- if (n_post > 0) round(100 * v / n_post, 1) else 0
      sprintf("%d (%.1f\\%%)", as.integer(round(v)), pct)
    }, character(1))
    lines <- c(lines, sprintf("%s & %s\\\\", tex_escape_cfg(rn[i]), paste(cells, collapse = " & ")))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, path)
  message("Wrote ", path)
}

ef_param_order <- c(
  "mu_0", "mu_1",
  "A_00", "alpha_m_00",
  "A_11", "alpha_m_11",
  "A_01", "alpha_m_01",
  "A_10", "alpha_m_10",
  "c", "p", "D", "gamma", "q"
)

write_tex_ok_ef_params <- function(E_vec, F_vec, path) {
  lines <- c(
    "% Auto-generated by oklahoma_paper_assets.R — do not edit by hand",
    "\\begin{table}[t]",
    "\\centering",
    "\\footnotesize",
    "\\caption{\\label{tab:ok_ef_params}Full bivariate ETAS parameter estimates: naive KDE (E) vs.\\ SEM KDE (F).}",
    "\\begin{tabular}{@{}lcc@{}}",
    "\\toprule",
    "Parameter & E (naive KDE) & F (SEM KDE)\\\\",
    "\\midrule"
  )
  for (pn in ef_param_order) {
    ev <- E_vec[pn]
    fv <- F_vec[pn]
    lines <- c(
      lines,
      sprintf(
        "%s & %s & %s\\\\",
        tex_escape_cfg(pn),
        fmt_scaled_num(ev, 5L),
        fmt_scaled_num(fv, 5L)
      )
    )
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, path)
  message("Wrote ", path)
}

cm_F <- build_cm_sem(sem_fit_F)
n_post_F <- if (!is.null(cm_F)) sum(cm_F) else 0L
write_tex_ok_confusion_F(cm_F, file.path(tex_dir, "tab_ok_confusion_F.tex"), n_post_F)

E_vec <- to_named_num(if (!is.null(fit_E)) fit_E$params else NULL)
F_vec <- to_named_num(if (!is.null(fit_F_blk)) fit_F_blk$params else NULL)
write_tex_ok_ef_params(E_vec, F_vec, file.path(tex_dir, "tab_ok_ef_params.tex"))

write_tex_ok_sem_config <- function(cfg, path, rds_basename) {
  if (is.null(cfg) || length(cfg) < 1L) {
    cfg <- list()
  }
  pilot_cores <- cfg$SEM_PILOT_CORES
  pc1 <- if (is.null(pilot_cores) || length(pilot_cores) < 1L) {
    NA_integer_
  } else if (is.list(pilot_cores)) {
    suppressWarnings(as.integer(pilot_cores[[1L]]))
  } else {
    suppressWarnings(as.integer(pilot_cores[[1L]]))
  }
  pilot_cores_txt <- if (length(pc1) < 1L || is.na(pc1)) "auto" else as.character(pc1)

  rows <- rbind(
    c("SEM outer rounds ($N_{\\mathrm{iter}}$)", fmt_cfg_val(cfg$SEM_N_ITER)),
    c("Inner iterations (per outer round)", fmt_cfg_val(cfg$SEM_INNER_ITER)),
    c("Inner proposals (per inner iteration)", fmt_cfg_val(cfg$SEM_INNER_PROPS)),
    c("Retained labellings (per outer round)", fmt_cfg_val(cfg$SEM_N_LABELLINGS)),
    c("Outer optim cap (adaptive SEM, uni.~path)", fmt_cfg_val(cfg$SEM_OUTER_MAXIT)),
    c("Outer optim cap (bivariate ETAS path)", fmt_cfg_val(cfg$SEM_OUTER_MAXIT_BIV)),
    c("Warm-start fixed full-parameter step", fmt_cfg_val(cfg$SEM_WARMSTART_FIXED)),
    c("Relabelling selection method", fmt_cfg_val(cfg$SEM_OPTIM_METHOD)),
    c("Relabelling selection temperature", fmt_cfg_val(cfg$SEM_SELECTION_TEMPERATURE)),
    c("Proposal change factor", fmt_cfg_val(cfg$SEM_CHANGE_FACTOR)),
    c("Change-factor min.~multiplier", fmt_cfg_val(cfg$SEM_CHANGE_FACTOR_MIN_MULT)),
    c("Change-factor max.~multiplier", fmt_cfg_val(cfg$SEM_CHANGE_FACTOR_MAX_MULT)),
    c("Max.~relabel step (fraction of points)", fmt_cfg_val(cfg$SEM_MAX_RELABEL_STEP_FRAC)),
    c("Forced param.~update flip fraction", fmt_cfg_val(cfg$SEM_FORCE_PARAM_UPDATE_FLIP_FRAC)),
    c("Temporal relabel weight", fmt_cfg_val(cfg$SEM_TEMPORAL_WEIGHT)),
    c("Temporal relabel scale (days)", fmt_cfg_val(cfg$SEM_TEMPORAL_SCALE_DAYS)),
    c("Pilot tuning enabled", fmt_cfg_val(cfg$RUN_SEM_PILOT)),
    c("Pilot inner iterations", fmt_cfg_val(cfg$SEM_PILOT_INNER_ITER)),
    c("Pilot cores", pilot_cores_txt),
    c("Pilot max.~combinations", fmt_cfg_val(cfg$SEM_PILOT_MAX_COMBOS))
  )

  lines <- c(
    "% Auto-generated by oklahoma_paper_assets.R — do not edit by hand",
    "\\begin{table}[t]",
    "\\centering",
    "\\small",
    sprintf(
      "\\caption{\\label{tab:ok_sem_config}Stochastic SEM configuration for the Oklahoma \\emph{SEM KDE fit}. Values are taken from the archived run \\texttt{%s}.}",
      tex_escape_cfg(rds_basename)
    ),
    "\\begin{tabular}{@{}l@{\\quad}l@{}}",
    "\\toprule",
    "Setting & Value\\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(rows))) {
    lines <- c(lines, sprintf("%s & %s\\\\", rows[i, 1], rows[i, 2]))
  }
  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
  writeLines(lines, path)
  message("Wrote ", path)
}

# ---- LaTeX fragments ----
write_tex_ok_counts(file.path(tex_dir, "tab_ok_counts.tex"))
write_tex_ok_sem_config(
  res$config,
  file.path(tex_dir, "tab_ok_sem_config.tex"),
  basename(input_rds)
)

if (have_boot) {
  summ <- boot_df %>%
    group_by(.data$model) %>%
    summarize(
      n_boot = n(),
      mean_saved = mean(.data$ate_total_mean, na.rm = TRUE),
      sd_saved = sd(.data$ate_total_mean, na.rm = TRUE),
      q025 = as.numeric(quantile(.data$ate_total_mean, 0.025, na.rm = TRUE)),
      q50 = as.numeric(quantile(.data$ate_total_mean, 0.50, na.rm = TRUE)),
      q975 = as.numeric(quantile(.data$ate_total_mean, 0.975, na.rm = TRUE)),
      .groups = "drop"
    )

  lines_boot <- c(
    "% Auto-generated by oklahoma_paper_assets.R",
    "\\begin{table}",
    "\\caption{\\label{tab:ok_bootstrap_summary}Bootstrap summary (bias-corrected replicate means of $\\hat\\Delta_{\\mathrm{AoN}}$, days matching analysis window).}",
    "\\centering",
    "\\begin{tabular}[t]{@{}lrrrrrr@{}}",
    "\\toprule",
    "Model & $n$ & Mean & SD & 2.5\\% & 50\\% & 97.5\\%\\\\",
    "\\midrule"
  )
  fmt1 <- function(x) sprintf("%.1f", as.numeric(x))
  for (i in seq_len(nrow(summ))) {
    lab <- if (grepl("Naive", summ$model[i])) "Naive" else "SEM"
    lines_boot <- c(
      lines_boot,
      sprintf(
        "%s & %d & %s & %s & %s & %s & %s\\\\",
        lab,
        as.integer(summ$n_boot[i]),
        fmt1(summ$mean_saved[i]),
        fmt1(summ$sd_saved[i]),
        fmt1(summ$q025[i]),
        fmt1(summ$q50[i]),
        fmt1(summ$q975[i])
      )
    )
  }
  lines_boot <- c(
    lines_boot,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
  writeLines(lines_boot, file.path(tex_dir, "tab_ok_bootstrap_summary.tex"))
  message("Wrote ", file.path(tex_dir, "tab_ok_bootstrap_summary.tex"))

  utils::write.csv(summ, file.path(tex_dir, "ef_bootstrap_summary.csv"), row.names = FALSE)
} else {
  message("Skipped tab_ok_bootstrap_summary.tex and ef_bootstrap_summary.csv (no bootstrap in RDS).")
}

fig_rel <- function(nm) paste0("inst/oklahoma/paper/generated/figures/", nm)
gen_plots <- c(
  fig_rel("cumulative_count.pdf"),
  fig_rel("ok_partition_county.pdf"),
  fig_rel("ok_point_patterns.pdf"),
  fig_rel("ok_sem_f_confusion_trace.pdf"),
  fig_rel("ok_sem_f_metric_trace.pdf"),
  fig_rel("ok_sem_f_flips.pdf")
)
if (isTRUE(have_boot)) {
  gen_plots <- c(
    fig_rel("ATE_diff.pdf"),
    fig_rel("ok_bootstrap_ecdf.pdf"),
    fig_rel("ok_bootstrap_running_mean.pdf"),
    fig_rel("ok_sim_total_saved_hist.pdf"),
    gen_plots
  )
}
gen_tex <- c(
  "inst/oklahoma/paper/generated/tab_ok_counts.tex",
  "inst/oklahoma/paper/generated/tab_ok_sem_config.tex",
  "inst/oklahoma/paper/generated/tab_ok_confusion_F.tex",
  "inst/oklahoma/paper/generated/tab_ok_ef_params.tex"
)
if (isTRUE(have_boot)) {
  gen_tex <- c(gen_tex, "inst/oklahoma/paper/generated/tab_ok_bootstrap_summary.tex")
}

meta_out <- list(
  input_rds = input_rds,
  ate_days = ate_days,
  plots_dir = plots_dir,
  tex_dir = tex_dir,
  have_bootstrap = have_boot,
  generated = c(gen_plots, gen_tex)
)
saveRDS(meta_out, file.path(tex_dir, "build_manifest.rds"))
message("Done. Manifest: ", file.path(tex_dir, "build_manifest.rds"))
