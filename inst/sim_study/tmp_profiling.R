#!/usr/bin/env Rscript
# Temporary profiling script for the simulation study. Do NOT change sim_study.R.
# Run from package root: Rscript inst/sim_study/tmp_profiling.R
# Or from inst/sim_study: Rscript tmp_profiling.R
# Outputs: tmp_prof.out, tmp_prof_summary.rds, tmp_prof_by_self.csv, tmp_prof_by_total.csv, tmp_profiling_report.txt

if (file.exists("sim_study.R")) {
  # already in inst/sim_study
} else if (file.exists("inst/sim_study/sim_study.R")) {
  setwd("inst/sim_study")
} else {
  stop("Run from package root or inst/sim_study so that sim_study.R is found.")
}

# Force --small so profiling run finishes in reasonable time
override_args <- c("--small")
commandArgs <- function(trailingOnly = TRUE) override_args

out_prof   <- "tmp_prof.out"
out_rds    <- "tmp_prof_summary.rds"
out_self   <- "tmp_prof_by_self.csv"
out_total  <- "tmp_prof_by_total.csv"
out_report <- "tmp_profiling_report.txt"

cat("Profiling sim_study.R with --small ...\n")
Rprof(out_prof, interval = 0.02, memory.profiling = FALSE)
t0 <- proc.time()
source("sim_study.R", local = FALSE)
t1 <- proc.time()
Rprof(NULL)
elapsed <- t1 - t0

s <- summaryRprof(out_prof)
saveRDS(s, out_rds)
if (!is.null(s$by.self) && nrow(s$by.self) > 0) {
  write.csv(s$by.self, out_self)
  write.csv(s$by.total, out_total)
}

# Write a short text report
con <- file(out_report, "w")
writeLines(c(
  "=== Simulation study profiling (--small) ===",
  paste("Total elapsed (s):", round(elapsed[3], 2)),
  "",
  "Top 40 by self time:"
), con)
if (!is.null(s$by.self) && nrow(s$by.self) > 0) {
  write.table(head(s$by.self, 40), file = con, sep = "\t", quote = FALSE)
  writeLines("", con)
  writeLines("Top 40 by total time:", con)
  write.table(head(s$by.total, 40), file = con, sep = "\t", quote = FALSE)
}
close(con)
cat("Done. Summary in", out_report, "\n")
