#!/usr/bin/env Rscript
# Merge per-scenario robustness RDS from a NeSI job into a paper archive folder
# and rewrite the archive manifest.
#
# Usage:
#   Rscript inst/sim_study/merge_robustness_job_into_archive.R \
#     --job 7839871 \
#     --archive robustness_merged_tcal \
#     --family k_spatial_range
#
# Copies matching RDS into paper/<archive>/, appends/replaces those family rows
# in <archive>_manifest.csv (paths point at the archive copies).

suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) < 1L || idx[[1]] >= length(args)) return(default)
  args[[idx[[1]] + 1L]]
}

job_id <- get_arg("--job", Sys.getenv("PP_MERGE_JOB_ID", ""))
archive <- get_arg("--archive", Sys.getenv("PP_MERGE_ARCHIVE", "robustness_merged_tcal"))
family <- get_arg("--family", Sys.getenv("PP_MERGE_FAMILY", "k_spatial_range"))
if (!nzchar(job_id)) stop("Pass --job <SLURM_JOB_ID>")
if (!nzchar(archive)) stop("Pass --archive <basename>")
if (!nzchar(family)) stop("Pass --family <scenario_family>")

fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
repo_dir <- if (length(fa)) {
  normalizePath(file.path(dirname(sub("^--file=", "", fa[[1]])), "..", ".."),
                winslash = "/", mustWork = TRUE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
sim_root <- PPDisentangle::pp_output_path("sim_study", repo_root = repo_dir)
archive_dir <- file.path(sim_root, "paper", archive)
dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)

job_prefix <- paste0("robustness_", job_id, "_")
src_rds <- list.files(
  sim_root,
  pattern = paste0("^", job_prefix, family, "_.*\\.rds$"),
  full.names = TRUE
)
if (!length(src_rds)) {
  stop(sprintf(
    "No RDS found for job=%s family=%s under %s (pattern %s%s_*.rds)",
    job_id, family, sim_root, job_prefix, family
  ))
}

copied <- character()
for (src in src_rds) {
  dest <- file.path(archive_dir, basename(src))
  ok <- file.copy(src, dest, overwrite = TRUE)
  if (!isTRUE(ok)) stop("Failed to copy ", src, " -> ", dest)
  copied <- c(copied, dest)
}
message(sprintf("[merge] copied %d RDS into %s", length(copied), archive_dir))

job_manifest_path <- file.path(sim_root, paste0("robustness_", job_id, "_manifest.csv"))
archive_manifest_path <- file.path(archive_dir, paste0(archive, "_manifest.csv"))
if (!file.exists(archive_manifest_path)) {
  stop("Archive manifest not found: ", archive_manifest_path)
}

archive_man <- read.csv(archive_manifest_path, stringsAsFactors = FALSE)
if (file.exists(job_manifest_path)) {
  job_man <- read.csv(job_manifest_path, stringsAsFactors = FALSE)
  add <- job_man %>% filter(.data$scenario_family == family)
  if (nrow(add) < 1L) {
    stop("Job manifest has no rows for family=", family)
  }
} else {
  warning("Job manifest missing; building rows from copied RDS basenames only")
  add <- NULL
}

# Prefer job manifest rows; rewrite paths/basenames to archive copies.
if (!is.null(add) && nrow(add) > 0L) {
  add$run_basename <- paste0("robustness_", job_id, "_", add$scenario_id)
  add$rds_path <- file.path(archive_dir, paste0(add$run_basename, ".rds"))
  missing <- add$rds_path[!file.exists(add$rds_path)]
  if (length(missing)) {
    stop("Missing archive RDS after copy: ", paste(basename(missing), collapse = ", "))
  }
} else {
  stop("Unable to build merge rows for family=", family)
}

# Align columns: union of archive + add, fill missing with NA.
all_cols <- union(names(archive_man), names(add))
for (nm in setdiff(all_cols, names(archive_man))) archive_man[[nm]] <- NA
for (nm in setdiff(all_cols, names(add))) add[[nm]] <- NA
archive_man <- archive_man[, all_cols, drop = FALSE]
add <- add[, all_cols, drop = FALSE]

keep <- archive_man %>% filter(.data$scenario_family != family)
merged <- bind_rows(keep, add) %>%
  distinct(.data$scenario_id, .keep_all = TRUE) %>%
  arrange(.data$scenario_family, .data$scenario_id)

write.csv(merged, archive_manifest_path, row.names = FALSE)
# Flat copy under sim_root for convenience.
write.csv(merged, file.path(sim_root, paste0(archive, "_manifest.csv")), row.names = FALSE)

message(sprintf(
  "[merge] archive manifest now %d rows (%d %s); wrote %s",
  nrow(merged),
  sum(merged$scenario_family == family, na.rm = TRUE),
  family,
  archive_manifest_path
))
print(table(merged$scenario_family))
