#!/bin/bash
# First-time setup for PPDisentangle on NeSI Mahuika.
# Run once on the login node before submitting jobs.
#
# Usage:
#   cd /path/to/PPDisentangle
#   bash inst/sim_study/setup_nesi.sh

set -euo pipefail

PKG_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
echo "=== PPDisentangle NeSI Setup ==="
echo "Package root: $PKG_ROOT"
echo ""

# ---- Load modules ----
# R-Geo bundles R + GDAL + GEOS + PROJ + UDUNITS (required for sf, terra, spatstat)
module purge
module load R-Geo/4.3.2-foss-2023a

echo "R version: $(R --version | head -1)"
echo "R library paths:"
Rscript -e 'cat(paste(.libPaths(), collapse="\n"), "\n")'
echo ""

# ---- Install R dependencies ----
echo "Installing R dependencies..."
Rscript -e '
deps <- c("spatstat", "ggplot2", "dplyr", "data.table", "parallel",
          "doParallel", "R.utils", "reshape2", "gridExtra", "scales",
          "Rcpp", "sf", "terra", "numDeriv", "rlang",
          "knitr", "kableExtra", "tidyr", "jsonlite")
installed <- rownames(installed.packages())
to_install <- setdiff(deps, installed)
if (length(to_install) > 0) {
  cat("Installing:", paste(to_install, collapse=", "), "\n")
  install.packages(to_install, repos="https://cloud.r-project.org", Ncpus=4)
} else {
  cat("All dependencies already installed.\n")
}
'
echo ""

# ---- Install PPDisentangle ----
echo "Installing PPDisentangle from source..."
R CMD INSTALL --no-multiarch "$PKG_ROOT" 2>&1
echo ""

# ---- Verify ----
echo "Verifying installation..."
Rscript -e '
library(PPDisentangle)
cat("PPDisentangle loaded successfully.\n")
cat("Functions available:", length(ls("package:PPDisentangle")), "\n")
'

mkdir -p "$PKG_ROOT/output/sim_study"

echo ""
echo "=== Setup complete ==="
echo ""
echo "Next steps:"
echo "  1. Test:   bash inst/sim_study/run_nesi.sh --test --sims 2"
echo "  2. Run:    bash inst/sim_study/run_nesi.sh --sims 50"
echo ""
echo "Monitor:  squeue -u \$USER"
echo "Output:   output/sim_study/<JOB_ID>.rds  output/sim_study/<JOB_ID>.log"
