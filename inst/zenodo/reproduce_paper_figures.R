#!/usr/bin/env Rscript
#
# Reproduce paper figures/tables from a Zenodo (or local) PPDisentangle-output tree.
#
# Prerequisites:
#   1. Clone PPDisentangle and install package dependencies.
#   2. Download/unpack the Zenodo archive so it looks like:
#        <output-root>/
#          sim_study/paper/main_5228509/time_sweep_5228509_summary_FOR_PAPER.rds
#          sim_study/paper/robustness_merged_tcal/
#          oklahoma/for_paper.rds
#
# Usage (from repo root):
#   Rscript inst/zenodo/reproduce_paper_figures.R
#   Rscript inst/zenodo/reproduce_paper_figures.R --output-root /path/to/PPDisentangle-output
#   Rscript inst/zenodo/reproduce_paper_figures.R --skip-oklahoma
#   Rscript inst/zenodo/reproduce_paper_figures.R --skip-robustness
#
# Canonical paper jobs:
#   main simulation figures  -> paper/main_5228509
#   robustness appendix      -> paper/robustness_merged_tcal
#   Oklahoma application     -> oklahoma/for_paper.rds
#
# Does not rebuild local LaTeX previews such as robustness.pdf.

get_arg_val <- function(args, flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) < 1L || idx[[1]] >= length(args)) return(default)
  args[[idx[[1]] + 1L]]
}

has_flag <- function(args, flag) {
  flag %in% args
}

find_repo_root <- function(start_dir = getwd()) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(cur, "DESCRIPTION")) &&
        dir.exists(file.path(cur, "inst"))) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) {
      stop("Could not find PPDisentangle repo root from: ", start_dir, call. = FALSE)
    }
    cur <- parent
  }
}

args <- commandArgs(trailingOnly = TRUE)
repo_root <- find_repo_root()
source(file.path(repo_root, "R", "paths.R"), local = FALSE)

output_root_arg <- get_arg_val(args, "--output-root", "")
if (nzchar(output_root_arg)) {
  Sys.setenv(PPDISENTANGLE_OUTPUT_ROOT = normalizePath(output_root_arg, winslash = "/", mustWork = TRUE))
}

output_root <- pp_output_root(repo_root = repo_root)
sim_dir <- pp_output_path("sim_study", repo_root = repo_root)
ok_rds <- pp_output_path("oklahoma", "for_paper.rds", repo_root = repo_root)

main_summary <- file.path(
  sim_dir, "paper", "main_5228509", "time_sweep_5228509_summary_FOR_PAPER.rds"
)
robustness_basename <- "robustness_merged_tcal"
skip_main <- has_flag(args, "--skip-main")
skip_robustness <- has_flag(args, "--skip-robustness")
skip_oklahoma <- has_flag(args, "--skip-oklahoma")

message("Repo root:    ", repo_root)
message("Output root:  ", output_root)
message("Sim study:    ", sim_dir)

run_rscript <- function(rel_script, extra_args = character()) {
  script_abs <- file.path(repo_root, rel_script)
  if (!file.exists(script_abs)) stop("Missing script: ", script_abs, call. = FALSE)
  # Run via relative path from repo_root so spaced absolute paths (Google Drive
  # "My Drive") never appear in Rscript --file= / source() resolution.
  message("\n=== Rscript ", paste(c(rel_script, extra_args), collapse = " "), " ===")
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(repo_root)
  status <- system2("Rscript", args = c(rel_script, extra_args))
  if (!identical(as.integer(status), 0L)) {
    stop(
      "Command failed with status ", status, ": Rscript ",
      paste(c(rel_script, extra_args), collapse = " "),
      call. = FALSE
    )
  }
}

if (!skip_main) {
  if (!file.exists(main_summary)) {
    stop(
      "Main simulation FOR_PAPER summary not found:\n  ", main_summary,
      "\nUnpack the Zenodo archive or set --output-root.",
      call. = FALSE
    )
  }
  run_rscript(
    "inst/sim_study/plot_time_sweep_publication.R",
    c(
      "--input", main_summary,
      "--rds-search-dir", file.path(sim_dir, "paper", "main_5228509")
    )
  )

  # Frozen illustrative realisation used by revision.tex fig:pp_realiz.
  # Not regenerated from the time-sweep; copy from the paper archive into generated/.
  hawkes_src <- file.path(sim_dir, "paper", "main_5228509", "simulated_hawkes_hawkes_process.pdf")
  hawkes_dst <- file.path(sim_dir, "generated", "figures", "simulated_hawkes_hawkes_process.pdf")
  if (!file.exists(hawkes_src)) {
    stop(
      "Missing illustrative Hawkes realisation PDF:\n  ", hawkes_src,
      "\nExpected in the Zenodo paper/main_5228509/ archive.",
      call. = FALSE
    )
  }
  dir.create(dirname(hawkes_dst), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(hawkes_src, hawkes_dst, overwrite = TRUE)
  if (!isTRUE(ok)) stop("Failed to copy ", hawkes_src, " -> ", hawkes_dst, call. = FALSE)
  message("Copied illustrative figure: ", hawkes_dst)
}

if (!skip_robustness) {
  manifest <- file.path(sim_dir, "paper", robustness_basename, paste0(robustness_basename, "_manifest.csv"))
  if (!file.exists(manifest)) {
    stop(
      "Robustness manifest not found:\n  ", manifest,
      "\nExpected job family ", robustness_basename, " under sim_study/paper/.",
      call. = FALSE
    )
  }
  run_rscript(
    "inst/sim_study/sim_study_robustness.R",
    c("--replot", robustness_basename)
  )
}

if (!skip_oklahoma) {
  if (!file.exists(ok_rds)) {
    stop(
      "Oklahoma results not found:\n  ", ok_rds,
      "\nUnpack oklahoma/ from Zenodo or pass --skip-oklahoma.",
      call. = FALSE
    )
  }
  run_rscript(
    "inst/oklahoma/paper/oklahoma_paper_assets.R",
    c("--input", ok_rds)
  )
}

message("\nDone. Regenerated assets under:")
message("  ", file.path(sim_dir, "generated/"))
if (!skip_oklahoma) {
  message("  ", pp_output_path("oklahoma", "paper", "generated", repo_root = repo_root))
}

# Lightweight session pin next to regenerated outputs
session_out <- file.path(sim_dir, "generated", "sessionInfo_reproduce.txt")
dir.create(dirname(session_out), recursive = TRUE, showWarnings = FALSE)
pinned <- file.path(repo_root, "inst", "zenodo", "sessionInfo.txt")
header <- c(
  paste0("Captured by reproduce_paper_figures.R: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Pinned reference: ", if (file.exists(pinned)) pinned else "(missing)"),
  ""
)
writeLines(c(header, utils::capture.output(print(utils::sessionInfo()))), session_out)
message("Wrote session snapshot: ", session_out)
if (file.exists(pinned)) {
  message("Compare against repo pin: ", pinned)
}
