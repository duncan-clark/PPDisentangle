#!/usr/bin/env Rscript
root <- normalizePath(commandArgs(trailingOnly = TRUE)[[1L]], mustWork = TRUE)
git <- Sys.getenv("GIT", "git")
commit <- system2(git, c("rev-parse", "HEAD"), stdout = TRUE)
stopifnot(length(commit) == 1L, grepl("^[a-f0-9]{40}$", commit))
manifest <- list(
  software_version = as.character(read.dcf("DESCRIPTION")[1L, "Version"]),
  source_commit = commit,
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  oklahoma_contrast = "observed",
  oklahoma_source_job = "8804859 + ATE backfill",
  robustness_scenarios = 67L,
  note = "Saved fits are frozen; derived figures and tables regenerated during packaging."
)
jsonlite::write_json(manifest, file.path(root, "release_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
paths <- list.files(root, recursive = TRUE, full.names = FALSE)
paths <- setdiff(paths, "MD5SUMS")
writeLines(paste(unname(tools::md5sum(file.path(root, paths))), paths, sep = "  "),
           file.path(root, "MD5SUMS"))
