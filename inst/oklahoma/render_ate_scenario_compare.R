#!/usr/bin/env Rscript
# Render a 4-scenario ATE bootstrap comparison HTML from ate_scenarios/*/for_paper.rds
#
# Usage (from package root):
#   Rscript inst/oklahoma/render_ate_scenario_compare.R
#   Rscript inst/oklahoma/render_ate_scenario_compare.R --output-root /path/to/PPDisentangle-output

args <- commandArgs(trailingOnly = TRUE)
script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- if (length(script_file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[1]), winslash = "/", mustWork = FALSE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
REPO_DIR <- if (basename(SCRIPT_DIR) == "oklahoma" && basename(dirname(SCRIPT_DIR)) == "inst") {
  normalizePath(dirname(dirname(SCRIPT_DIR)), winslash = "/", mustWork = FALSE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) < 1L || idx[[1]] >= length(args)) return(default)
  args[[idx[[1]] + 1L]]
}

output_root <- get_arg("--output-root", Sys.getenv("PPDISENTANGLE_OUTPUT_ROOT", ""))
if (!nzchar(output_root)) {
  output_root <- file.path(dirname(REPO_DIR), "PPDisentangle-output")
}
ok_dir <- file.path(output_root, "oklahoma")
scen_root <- file.path(ok_dir, "ate_scenarios")

scenario_order <- c("univ_aon", "univ_obs", "biv_aon", "biv_obs")
scenario_labels <- c(
  univ_aon = "Univariate · all-or-nothing\n(control everywhere vs treated everywhere)",
  univ_obs = "Univariate · observed vs control\n(control everywhere vs observed allocation)",
  biv_aon  = "Bivariate · all-or-nothing\n(control everywhere vs treated everywhere)",
  biv_obs  = "Bivariate · observed vs control\n(control everywhere vs observed allocation)"
)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

recenter_boot <- function(boot_vals, target_mean) {
  bv <- suppressWarnings(as.numeric(boot_vals))
  if (!is.finite(target_mean) || is.na(target_mean)) return(bv)
  bmean <- mean(bv, na.rm = TRUE)
  if (!is.finite(bmean) || is.na(bmean)) return(bv)
  bv + (target_mean - bmean)
}

point_mean <- function(ate) {
  if (is.null(ate) || is.null(ate$all_nothing_sim)) return(NA_real_)
  v <- suppressWarnings(as.numeric(ate$all_nothing_sim$total_saved))
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  mean(v)
}

load_scenario <- function(scen) {
  path <- file.path(scen_root, scen, "for_paper.rds")
  if (!file.exists(path)) {
    alt <- file.path(ok_dir, paste0("for_paper_", scen, ".rds"))
    if (file.exists(alt)) path <- alt else return(NULL)
  }
  r <- readRDS(path)
  list(
    scenario = scen,
    label = scenario_labels[[scen]] %||% scen,
    path = path,
    method = r$config$ATE_METHOD %||% paste0(
      if (isTRUE(r$config$ATE_BIVARIATE)) "bivariate" else "univariate",
      "_", r$config$ATE_CONTRAST %||% "?"
    ),
    bivariate = isTRUE(r$config$ATE_BIVARIATE),
    contrast = as.character(r$config$ATE_CONTRAST %||% NA),
    point_E = point_mean(r$fits_named$E$ate),
    point_F = point_mean(r$fits_named$F$ate),
    boot_E = r$bootstrap_ate$fit_E$replicate_summary$ate_total_mean,
    boot_F = r$bootstrap_ate$fit_F$replicate_summary$ate_total_mean,
    n_E = length(r$bootstrap_ate$fit_E$replicate_summary$ate_total_mean %||% numeric(0)),
    n_F = length(r$bootstrap_ate$fit_F$replicate_summary$ate_total_mean %||% numeric(0)),
    job_id = r$metadata$job_id %||% NA_character_
  )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

loaded <- lapply(scenario_order, load_scenario)
names(loaded) <- scenario_order
present <- loaded[!vapply(loaded, is.null, logical(1))]
if (!length(present)) {
  stop("No scenario RDS found under ", scen_root,
       "\nExpected ate_scenarios/{univ_aon,univ_obs,biv_aon,biv_obs}/for_paper.rds")
}

summ_rows <- lapply(present, function(s) {
  be <- recenter_boot(s$boot_E, s$point_E)
  bf <- recenter_boot(s$boot_F, s$point_F)
  data.frame(
    scenario = s$scenario,
    label = gsub("\n", " ", s$label, fixed = TRUE),
    method = s$method,
    bivariate = s$bivariate,
    contrast = s$contrast,
    job_id = s$job_id,
    n_E = s$n_E,
    n_F = s$n_F,
    point_E = s$point_E,
    point_F = s$point_F,
    boot_raw_mean_E = mean(s$boot_E, na.rm = TRUE),
    boot_raw_mean_F = mean(s$boot_F, na.rm = TRUE),
    boot_bc_mean_E = mean(be, na.rm = TRUE),
    boot_bc_mean_F = mean(bf, na.rm = TRUE),
    boot_bc_sd_E = stats::sd(be, na.rm = TRUE),
    boot_bc_sd_F = stats::sd(bf, na.rm = TRUE),
    boot_bc_p025_E = as.numeric(stats::quantile(be, 0.025, na.rm = TRUE, names = FALSE)),
    boot_bc_p975_E = as.numeric(stats::quantile(be, 0.975, na.rm = TRUE, names = FALSE)),
    boot_bc_p025_F = as.numeric(stats::quantile(bf, 0.025, na.rm = TRUE, names = FALSE)),
    boot_bc_p975_F = as.numeric(stats::quantile(bf, 0.975, na.rm = TRUE, names = FALSE)),
    stringsAsFactors = FALSE
  )
})
summ <- bind_rows(summ_rows)

boot_df <- bind_rows(lapply(present, function(s) {
  be <- recenter_boot(s$boot_E, s$point_E)
  bf <- recenter_boot(s$boot_F, s$point_F)
  bind_rows(
    data.frame(scenario = s$scenario, label = s$label, model = "Naive (E)",
               ate = be, stringsAsFactors = FALSE),
    data.frame(scenario = s$scenario, label = s$label, model = "SEM (F)",
               ate = bf, stringsAsFactors = FALSE)
  )
})) %>% filter(is.finite(.data$ate))

boot_df$label <- factor(boot_df$label, levels = scenario_labels[scenario_order])
boot_df$model <- factor(boot_df$model, levels = c("Naive (E)", "SEM (F)"))

means_df <- boot_df %>%
  group_by(.data$label, .data$model) %>%
  summarise(ate = mean(.data$ate, na.rm = TRUE), .groups = "drop")

fig_dir <- file.path(ok_dir, "paper", "generated", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
p <- ggplot(boot_df, aes(x = .data$ate, fill = .data$model)) +
  geom_histogram(bins = 30, alpha = 0.55, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_vline(data = means_df, aes(xintercept = .data$ate, colour = .data$model),
             linewidth = 0.9, show.legend = FALSE) +
  facet_wrap(~label, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c("Naive (E)" = "#fdb462", "SEM (F)" = "#80b1d3")) +
  scale_colour_manual(values = c("Naive (E)" = "#e6550d", "SEM (F)" = "#3182bd")) +
  labs(
    x = "Bias-corrected bootstrap replicate mean total saved / 100 days",
    y = "Count",
    fill = "Model",
    title = "Oklahoma ATE scenarios: debiased bootstrap comparison",
    subtitle = "Frozen E/F fits; each panel is a full parametric bootstrap under that estimand"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold", size = 9)
  )
ggsave(file.path(fig_dir, "ATE_scenarios_compare.pdf"), p, width = 11, height = 8.5)
ggsave(file.path(fig_dir, "ATE_scenarios_compare.png"), p, width = 11, height = 8.5, dpi = 150)

utils::write.csv(summ, file.path(ok_dir, "ate_scenarios", "scenario_summary.csv"), row.names = FALSE)

html_path <- file.path(ok_dir, "oklahoma_ate_scenarios.html")
esc <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}
fmt <- function(x, d = 1) ifelse(is.finite(x), sprintf(paste0("%.", d, "f"), x), "NA")

rows_html <- paste0(vapply(seq_len(nrow(summ)), function(i) {
  s <- summ[i, ]
  sprintf(
    paste0(
      "<tr>",
      "<td>%s</td><td><code>%s</code></td><td>%s</td>",
      "<td>%s</td><td>%s</td>",
      "<td>%s</td><td>%s</td>",
      "<td>%s</td><td>%s</td>",
      "<td>[%s, %s]</td><td>[%s, %s]</td>",
      "</tr>"
    ),
    esc(s$label), esc(s$method), esc(s$job_id),
    fmt(s$point_E), fmt(s$point_F),
    fmt(s$boot_bc_mean_E), fmt(s$boot_bc_sd_E),
    fmt(s$boot_bc_mean_F), fmt(s$boot_bc_sd_F),
    fmt(s$boot_bc_p025_E), fmt(s$boot_bc_p975_E),
    fmt(s$boot_bc_p025_F), fmt(s$boot_bc_p975_F)
  )
}, character(1)), collapse = "\n")

html <- paste0(
  "<!DOCTYPE html><html><head><meta charset='utf-8'>",
  "<title>Oklahoma ATE scenario comparison</title>",
  "<style>",
  "body{font-family:Helvetica,Arial,sans-serif;margin:2rem;max-width:1100px;}",
  "table{border-collapse:collapse;width:100%;font-size:14px;}",
  "th,td{border:1px solid #ccc;padding:6px 8px;text-align:left;}",
  "th{background:#f3f3f3;}",
  "img{max-width:100%;height:auto;border:1px solid #ddd;}",
  "code{background:#f6f8fa;padding:1px 4px;}",
  "</style></head><body>",
  "<h1>Oklahoma ATE scenario comparison</h1>",
  "<p>Four clean estimands on <strong>frozen</strong> Naive (E) / SEM (F) fits, ",
  "each with a full parametric bootstrap. Histograms are <em>bias-corrected</em> ",
  "(shape/spread preserved, recentered to the matching point ATE).</p>",
  "<h2>Summary</h2>",
  "<table><thead><tr>",
  "<th>Scenario</th><th>Method</th><th>Job</th>",
  "<th>Point E</th><th>Point F</th>",
  "<th>Boot E mean</th><th>Boot E SD</th>",
  "<th>Boot F mean</th><th>Boot F SD</th>",
  "<th>E 95% (bc)</th><th>F 95% (bc)</th>",
  "</tr></thead><tbody>",
  rows_html,
  "</tbody></table>",
  "<h2>Bias-corrected bootstrap histograms</h2>",
  "<p><img src='paper/generated/figures/ATE_scenarios_compare.png' ",
  "alt='ATE scenario comparison'></p>",
  "<p>Source directory: <code>", esc(scen_root), "</code><br>",
  "Scenarios found: ", esc(paste(names(present), collapse = ", ")), "</p>",
  "</body></html>"
)
writeLines(html, html_path)

# Also mirror figure into docs/paper/plots/oklahoma when present
plots_mirror <- file.path(REPO_DIR, "docs", "paper", "plots", "oklahoma")
if (dir.exists(dirname(plots_mirror))) {
  dir.create(plots_mirror, recursive = TRUE, showWarnings = FALSE)
  file.copy(
    file.path(fig_dir, "ATE_scenarios_compare.pdf"),
    file.path(plots_mirror, "ATE_scenarios_compare.pdf"),
    overwrite = TRUE
  )
}

cat("Wrote:\n")
cat("  ", html_path, "\n", sep = "")
cat("  ", file.path(fig_dir, "ATE_scenarios_compare.pdf"), "\n", sep = "")
cat("  ", file.path(ok_dir, "ate_scenarios", "scenario_summary.csv"), "\n", sep = "")
cat("Scenarios included:", paste(names(present), collapse = ", "), "\n")
invisible(summ)
