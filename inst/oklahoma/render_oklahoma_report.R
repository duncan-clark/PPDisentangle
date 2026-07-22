#!/usr/bin/env Rscript
# Render oklahoma_report.qmd into PPDisentangle-output/oklahoma/.
#
# Usage:
#   Rscript render_oklahoma_report.R
#   Rscript render_oklahoma_report.R --results ../PPDisentangle-output/oklahoma/for_paper.rds
#   OK_RESULTS_FILE=/path/to/results.rds Rscript render_oklahoma_report.R --to html

args <- commandArgs(trailingOnly = TRUE)

script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- if (length(script_file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[1]),
                        winslash = "/", mustWork = FALSE))
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

OUT_DIR <- PPDisentangle::pp_output_path("oklahoma", repo_root = REPO_DIR)
report_out_dir_raw <- trimws(Sys.getenv("OK_REPORT_OUTPUT_DIR", ""))
REPORT_OUT_DIR <- if (nzchar(report_out_dir_raw)) {
  normalizePath(report_out_dir_raw, winslash = "/", mustWork = FALSE)
} else {
  OUT_DIR
}
dir.create(REPORT_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

get_arg_val <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) < 1L || idx[[1]] >= length(args)) return(default)
  args[[idx[[1]] + 1L]]
}

results_arg <- get_arg_val("--results", NULL)
explicit_result <- trimws(Sys.getenv("OK_RESULTS_FILE", ""))
if (is.character(results_arg) && length(results_arg) == 1L && nzchar(results_arg)) {
  results_file <- normalizePath(results_arg, winslash = "/", mustWork = FALSE)
} else if (nzchar(explicit_result)) {
  results_file <- normalizePath(explicit_result, winslash = "/", mustWork = FALSE)
} else {
  candidates <- c(
    file.path(OUT_DIR, "for_paper.rds"),
    file.path(OUT_DIR, "oklahoma_results.rds")
  )
  extra <- list.files(
    OUT_DIR,
    pattern = "^oklahoma_results(_[0-9]+|_job[0-9]+)?\\.rds$",
    full.names = TRUE
  )
  candidates <- unique(c(candidates, extra))
  existing <- candidates[file.exists(candidates)]
  if (length(existing) < 1L) {
    stop(
      "No results RDS found. Pass --results PATH or set OK_RESULTS_FILE.\n",
      "Expected under: ", OUT_DIR
    )
  }
  mtimes <- file.info(existing)$mtime
  results_file <- existing[[which.max(mtimes)]]
}

if (!file.exists(results_file)) {
  stop("Results file not found: ", results_file)
}

formats_raw <- tolower(trimws(Sys.getenv("OK_REPORT_FORMATS", "html")))
if ("--to" %in% args) {
  formats_raw <- get_arg_val("--to", formats_raw)
}
formats <- unique(trimws(unlist(strsplit(formats_raw, ","))))
formats <- formats[formats %in% c("html", "pdf")]
if (length(formats) < 1L) {
  stop("No valid report formats requested (expected html and/or pdf).")
}

quarto_bin <- Sys.which("quarto")
if (!nzchar(quarto_bin)) {
  stop("Quarto not found in PATH.")
}

report_file <- file.path(SCRIPT_DIR, "oklahoma_report.qmd")
if (!file.exists(report_file)) {
  stop("Report source not found: ", report_file)
}

old_wd <- getwd()
old_results_env <- Sys.getenv("OK_RESULTS_FILE", unset = NA_character_)
Sys.setenv(OK_RESULTS_FILE = normalizePath(results_file, winslash = "/", mustWork = FALSE))
setwd(SCRIPT_DIR)
on.exit({
  setwd(old_wd)
  if (is.na(old_results_env)) {
    Sys.unsetenv("OK_RESULTS_FILE")
  } else {
    Sys.setenv(OK_RESULTS_FILE = old_results_env)
  }
}, add = TRUE)

path_rel_to <- function(path, start) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  start <- normalizePath(start, winslash = "/", mustWork = FALSE)
  path_seg <- strsplit(path, "/", fixed = TRUE)[[1L]]
  start_seg <- strsplit(start, "/", fixed = TRUE)[[1L]]
  path_seg <- path_seg[nzchar(path_seg)]
  start_seg <- start_seg[nzchar(start_seg)]
  n <- min(length(path_seg), length(start_seg))
  i <- 1L
  while (i <= n && identical(path_seg[[i]], start_seg[[i]])) i <- i + 1L
  rel_parts <- c(rep("..", length(start_seg) - i + 1L), path_seg[i:length(path_seg)])
  if (length(rel_parts) < 1L) "." else paste(rel_parts, collapse = "/")
}
report_output_dir_arg <- path_rel_to(REPORT_OUT_DIR, SCRIPT_DIR)
if (!nzchar(report_output_dir_arg)) report_output_dir_arg <- "."

cat("Rendering report from: ", results_file, "\n", sep = "")
cat("Output directory: ", REPORT_OUT_DIR, "\n", sep = "")
cat("Quarto output-dir (relative): ", report_output_dir_arg, "\n", sep = "")

errors <- character(0)
for (fmt in formats) {
  report_output_path <- file.path(
    REPORT_OUT_DIR,
    paste0("oklahoma_report.", if (fmt == "html") "html" else "pdf")
  )
  status <- tryCatch(
    system2(
      quarto_bin,
      c(
        "render", basename(report_file), "--to", fmt,
        "--output-dir", report_output_dir_arg
      ),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) e
  )
  if (inherits(status, "error")) {
    errors <- c(errors, sprintf("[%s] %s", fmt, status$message))
    next
  }
  exit_code <- attr(status, "status")
  if (!is.null(exit_code) && !identical(exit_code, 0L)) {
    errors <- c(errors, sprintf("[%s] exit code %s", fmt, exit_code))
    next
  }
  if (!file.exists(report_output_path)) {
    errors <- c(errors, sprintf("[%s] expected output missing: %s", fmt, report_output_path))
    next
  }
  cat("Wrote: ", report_output_path, "\n", sep = "")
}

if (length(errors) > 0L) {
  stop(paste(errors, collapse = "\n"))
}

invisible(NULL)
