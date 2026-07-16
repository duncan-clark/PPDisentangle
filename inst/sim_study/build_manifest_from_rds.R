#!/usr/bin/env Rscript
# Build a robustness manifest CSV from completed per-scenario RDS files.
# Usage:
#   Rscript inst/sim_study/build_manifest_from_rds.R robustness_7671212

suppressPackageStartupMessages(library(dplyr))
`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
run_prefix <- if (length(args) >= 1L) args[[1]] else "robustness_7671212"

fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
repo_dir <- if (length(fa)) {
  normalizePath(file.path(dirname(sub("^--file=", "", fa[[1]])), "..", ".."),
                winslash = "/", mustWork = TRUE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
out_dir <- PPDisentangle::pp_output_path("sim_study", repo_root = repo_dir)
prefix <- paste0(run_prefix, "_")
rds <- list.files(out_dir, pattern = paste0("^", prefix, ".*\\.rds$"), full.names = TRUE)
rds <- rds[base::basename(rds) != paste0(run_prefix, "_summary.rds")]
if (!length(rds)) stop("No RDS files found for basename ", run_prefix, " in ", out_dir)

parse_tag <- function(sid, key) {
  m <- regmatches(sid, regexpr(paste0(key, "[0-9p]+"), sid, perl = TRUE))
  if (!length(m)) return(NA_real_)
  as.numeric(gsub("p", ".", sub(key, "", m)))
}

rows <- lapply(rds, function(p) {
  sid <- sub(paste0("^", prefix), "", sub("\\.rds$", "", base::basename(p)))
  fam <- if (startsWith(sid, "k_separation")) {
    "k_separation"
  } else if (startsWith(sid, "pretreatment_assignment")) {
    "pretreatment_assignment"
  } else if (startsWith(sid, "snr_scale")) {
    "snr_scale"
  } else if (startsWith(sid, "spatiotemporal_kernel_mismatch")) {
    "spatiotemporal_kernel_mismatch"
  } else if (startsWith(sid, "effect_modification")) {
    "effect_modification"
  } else if (startsWith(sid, "geometry_transport")) {
    "geometry_transport"
  } else {
    "unknown"
  }
  result <- readRDS(p)
  cfg <- result$config
  assign <- sub(".*_assign", "", sid)
  simk <- sub(".*_sim", "", sub("_fit.*", "", sid))
  fitk <- sub(".*_fit", "", sub("_ssim.*", "", sid))
  ssim <- sub(".*_ssim", "", sub("_sfit.*", "", sid))
  sfit <- sub(".*_sfit", "", sub("_assign.*", "", sid))
  kc <- parse_tag(sid, "kc")
  kt <- parse_tag(sid, "kt")
  data.frame(
    scenario_id = sid,
    scenario_family = fam,
    control_k = kc,
    treated_k = kt,
    k_delta = abs(kt - kc),
    mu_scale = parse_tag(sid, "mu"),
    sim_kernel = simk,
    fit_kernel = fitk,
    sim_spatial_kernel = ssim,
    fit_spatial_kernel = sfit,
    treatment_assignment = assign,
    h_true = suppressWarnings(as.numeric(cfg$H_TRUE %||% NA_real_)),
    geometry_m = suppressWarnings(as.integer(cfg$GEOMETRY_M %||% NA_integer_)),
    coarseness = suppressWarnings(as.numeric(cfg$COARSENESS %||% NA_real_)),
    target_points = suppressWarnings(as.numeric(cfg$TARGET_POINTS %||% 2500)),
    run_basename = sub("\\.rds$", "", base::basename(p)),
    rds_path = normalizePath(p, winslash = "/", mustWork = TRUE),
    stringsAsFactors = FALSE
  )
})

man <- bind_rows(rows)
paper_dir <- file.path(out_dir, "paper", run_prefix)
dir.create(paper_dir, recursive = TRUE, showWarnings = FALSE)
manifest_flat <- file.path(out_dir, paste0(run_prefix, "_manifest.csv"))
manifest_paper <- file.path(paper_dir, paste0(run_prefix, "_manifest.csv"))
write.csv(man, manifest_flat, row.names = FALSE)
write.csv(man, manifest_paper, row.names = FALSE)
message("Wrote ", nrow(man), " rows to:")
message("  ", manifest_flat)
message("  ", manifest_paper)
print(table(man$scenario_family))
