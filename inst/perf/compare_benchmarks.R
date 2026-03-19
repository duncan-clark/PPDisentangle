#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) < 1L) return(default)
  sub(paste0("^", name, "="), "", hit[[1]])
}

base_file <- arg_value("--base", NULL)
cand_file <- arg_value("--cand", NULL)
out_file <- arg_value("--out", NULL)
min_speedup <- as.numeric(arg_value("--min-speedup", "1.00"))
if (!is.finite(min_speedup)) min_speedup <- 1.00
if (is.null(base_file) || is.null(cand_file)) {
  stop("Usage: Rscript compare_benchmarks.R --base=<baseline.json> --cand=<candidate.json> [--out=<report.md>] [--min-speedup=1.00]")
}

base <- jsonlite::fromJSON(base_file, simplifyVector = FALSE)
cand <- jsonlite::fromJSON(cand_file, simplifyVector = FALSE)

components <- intersect(names(base$components), names(cand$components))
if (length(components) < 1L) stop("No overlapping components between files")

rows <- lapply(components, function(nm) {
  b <- base$components[[nm]]
  c <- cand$components[[nm]]
  b_time <- as.numeric(b$primary_time_sec)
  c_time <- as.numeric(c$primary_time_sec)
  speedup <- if (is.finite(b_time) && is.finite(c_time) && c_time > 0) b_time / c_time else NA_real_
  list(
    component = nm,
    base_sec = b_time,
    cand_sec = c_time,
    speedup = speedup,
    base_ok = isTRUE(b$correctness_pass),
    cand_ok = isTRUE(c$correctness_pass)
  )
})

df <- as.data.frame(do.call(rbind, lapply(rows, as.data.frame)), stringsAsFactors = FALSE)
df$base_sec <- as.numeric(df$base_sec)
df$cand_sec <- as.numeric(df$cand_sec)
df$speedup <- as.numeric(df$speedup)
df$base_ok <- as.logical(df$base_ok)
df$cand_ok <- as.logical(df$cand_ok)

speedup_pass <- all(is.na(df$speedup) | df$speedup >= min_speedup)
correctness_pass <- all(df$base_ok & df$cand_ok)
overall_pass <- speedup_pass && correctness_pass

cat("A/B benchmark summary\n")
cat("=====================\n")
print(df, row.names = FALSE)
cat(sprintf("\nCorrectness pass: %s\n", correctness_pass))
cat(sprintf("Speedup pass (>= %.3f): %s\n", min_speedup, speedup_pass))
cat(sprintf("Overall pass: %s\n", overall_pass))

if (!is.null(out_file)) {
  lines <- c(
    "# A/B Benchmark Report",
    "",
    sprintf("- Baseline: `%s`", normalizePath(base_file, mustWork = FALSE)),
    sprintf("- Candidate: `%s`", normalizePath(cand_file, mustWork = FALSE)),
    sprintf("- Minimum required speedup: `%.3f`", min_speedup),
    "",
    "## Results",
    "",
    "| Component | Baseline (s) | Candidate (s) | Speedup | Base OK | Cand OK |",
    "|---|---:|---:|---:|:---:|:---:|"
  )
  for (i in seq_len(nrow(df))) {
    lines <- c(
      lines,
      sprintf(
        "| %s | %.4f | %.4f | %s | %s | %s |",
        df$component[i],
        df$base_sec[i],
        df$cand_sec[i],
        if (is.na(df$speedup[i])) "NA" else sprintf("%.3f", df$speedup[i]),
        if (df$base_ok[i]) "yes" else "no",
        if (df$cand_ok[i]) "yes" else "no"
      )
    )
  }
  lines <- c(
    lines,
    "",
    "## Verdict",
    "",
    sprintf("- Correctness pass: **%s**", correctness_pass),
    sprintf("- Speedup pass: **%s**", speedup_pass),
    sprintf("- Overall pass: **%s**", overall_pass)
  )
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, out_file)
  cat("Wrote report:", out_file, "\n")
}

quit(status = if (overall_pass) 0L else 3L)
