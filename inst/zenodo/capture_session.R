#!/usr/bin/env Rscript
#
# Capture a lightweight sessionInfo() pin for Zenodo / figure reproduction.
#
# Usage (from repo root):
#   Rscript inst/zenodo/capture_session.R
#   Rscript inst/zenodo/capture_session.R --out /path/to/sessionInfo.txt

get_arg_val <- function(args, flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) < 1L || idx[[1]] >= length(args)) return(default)
  args[[idx[[1]] + 1L]]
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
out_path <- get_arg_val(
  args,
  "--out",
  file.path(repo_root, "inst", "zenodo", "sessionInfo.txt")
)

# Load package if installed; otherwise still record base + DESCRIPTION Imports.
pkg_ok <- FALSE
tryCatch({
  suppressPackageStartupMessages(library(PPDisentangle))
  pkg_ok <- TRUE
}, error = function(e) {
  message("PPDisentangle not installed in this R; recording base session only.")
})

desc <- read.dcf(file.path(repo_root, "DESCRIPTION"))
header <- c(
  paste0("Captured: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Package: ", desc[, "Package"], " ", desc[, "Version"]),
  paste0("PPDisentangle loaded: ", pkg_ok),
  ""
)

si <- utils::capture.output(print(utils::sessionInfo()))
writeLines(c(header, si), out_path)
message("Wrote ", out_path)
