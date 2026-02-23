#!/usr/bin/env bash
#
# Run the PPDisentangle simulation study (local or after pull/install).
# Creates cluster_output/, saves results list as .rds and timing report there.
#
# Usage:
#   From package root (directory containing R/, inst/, DESCRIPTION):
#     ./inst/sim_study/run_sim_study.sh [OPTIONS]
#
# Options:
#   --sims N     Number of simulated datasets (SIM_SIZE) and overrides N_SIMS;
#                also sets CORES_OVERRIDE to min(N, available_cores).
#   --no-pull    Skip 'git pull' (pull is run by default).
#   --no-install Skip devtools::install() (install is run by default).
#   --pull       Run 'git pull' (default: on).
#   --install    Run devtools::install() (default: on).
#   --small      Pass --small to the R script (reduced iterations for a quick test).
#   --cluster    Force cluster config when not under SLURM (optional; auto-set when run on cluster).
#
# When run on the cluster (module environment detected), cluster config is used automatically
# so you get Mode: CLUSTER and the same settings as sbatch. Locally, no --cluster is added.
#
# Examples:
#   ./inst/sim_study/run_sim_study.sh
#   ./inst/sim_study/run_sim_study.sh --sims 20
#   ./inst/sim_study/run_sim_study.sh --no-pull --no-install   # skip pull and install
#
# Outputs (sim study stdout/stderr go to .out only; console stays free):
#   cluster_output/run_<timestamp>.out                  (full run log; tail -f to follow)
#   cluster_output/sim_study_results_<timestamp>.rds    (rehydrate in R: x <- readRDS(...))
#   cluster_output/sim_study_timing_local.txt
#

set -e

# Defaults: pull and install by default; use --no-pull / --no-install to skip
DO_PULL=true
DO_INSTALL=true
N_SIMS_OVERRIDE=""
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    --sims)
      N_SIMS_OVERRIDE="$2"
      shift 2
      ;;
    --pull)
      DO_PULL=true
      shift
      ;;
    --no-pull)
      DO_PULL=false
      shift
      ;;
    --install)
      DO_INSTALL=true
      shift
      ;;
    --no-install)
      DO_INSTALL=false
      shift
      ;;
    --small)
      EXTRA_ARGS+=("--small")
      shift
      ;;
    --cluster)
      EXTRA_ARGS+=("--cluster")
      shift
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
  esac
done

# Package root = directory containing DESCRIPTION (run from there)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "$SCRIPT_DIR/../../DESCRIPTION" ]]; then
  PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
elif [[ -f "$SCRIPT_DIR/../DESCRIPTION" ]]; then
  PKG_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
else
  PKG_ROOT="$(pwd)"
fi
if [[ ! -f "$PKG_ROOT/DESCRIPTION" ]]; then
  echo "Error: Package root not found (no DESCRIPTION in $PKG_ROOT). Run from package root or ensure inst/sim_study/ is inside the package." >&2
  exit 1
fi

cd "$PKG_ROOT"
mkdir -p cluster_output

# Git pull first so we run the latest code (and get script fixes before loading modules)
if $DO_PULL; then
  echo "=== Git pull ==="
  git pull || { echo "Warning: git pull failed" >&2; }
fi

# Load R and geo modules when on HPC. If 'module' exists (e.g. NeSI), load without
# suppressing errors so libudunits2 etc. are found when R loads the 'units' package.
# Also set UDUNITS2_* so the R package 'units' can find udunits2 at *compile* time (devtools::install).
ON_CLUSTER_ENV=false
if command -v module &>/dev/null; then
  ON_CLUSTER_ENV=true
  module load R
  module load GDAL
  module load PROJ
  module load GEOS
  # Load udunits library (NeSI uses UDUNITS uppercase; try it first to avoid "unknown module" errors)
  if ! module load UDUNITS 2>/dev/null; then
    if ! module load udunits2 2>/dev/null; then
      if ! module load udunits 2>/dev/null; then
        echo "=== Warning: could not load UDUNITS/udunits2/udunits. Run: module spider udunits ===" >&2
        echo "   Then: module load UDUNITS/2.2.28-GCCcore-12.3.0 (or similar) and re-run." >&2
      fi
    fi
  fi
  # Help R's 'units' package configure find udunits2 when building from source
  if [[ -z "$UDUNITS2_INCLUDE" || -z "$UDUNITS2_LIBS" ]]; then
    if [[ -n "$UDUNITS2_DIR" ]]; then
      export UDUNITS2_INCLUDE="${UDUNITS2_INCLUDE:-$UDUNITS2_DIR/include}"
      export UDUNITS2_LIBS="${UDUNITS2_LIBS:--L$UDUNITS2_DIR/lib -ludunits2}"
    else
      for lib_dir in $(echo "${LD_LIBRARY_PATH:-}" | tr ':' ' '); do
        if [[ -n "$lib_dir" ]] && { [[ -f "$lib_dir/libudunits2.so" ]] || [[ -f "$lib_dir/libudunits2.so.0" ]]; }; then
          base="${lib_dir%/lib}"
          [[ -z "$UDUNITS2_LIBS" ]] && export UDUNITS2_LIBS="-L$lib_dir -ludunits2"
          [[ -z "$UDUNITS2_INCLUDE" ]] && export UDUNITS2_INCLUDE="$base/include"
          break
        fi
      done
    fi
  fi
  if [[ -n "$UDUNITS2_INCLUDE" && -n "$UDUNITS2_LIBS" ]]; then
    export UDUNITS2_INCLUDE UDUNITS2_LIBS
    echo "=== UDUNITS2 for R 'units' build: INCLUDE=$UDUNITS2_INCLUDE LIBS=$UDUNITS2_LIBS ==="
  else
    echo "=== Warning: UDUNITS2_INCLUDE/LIBS not set; R package 'units' may fail to build. ===" >&2
    echo "   Run 'module show udunits' to see paths, then set UDUNITS2_DIR or UDUNITS2_INCLUDE and UDUNITS2_LIBS and re-run." >&2
  fi
else
  true
fi

# When running on cluster (module environment), use cluster config so R reports Mode: CLUSTER
if $ON_CLUSTER_ENV && [[ " ${EXTRA_ARGS[*]} " != *" --cluster "* ]]; then
  EXTRA_ARGS+=("--cluster")
  echo "=== Cluster environment detected: using cluster config (Mode: CLUSTER) ==="
fi

# Optional: install package
if $DO_INSTALL; then
  echo "=== Installing package with devtools ==="
  Rscript -e 'if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools", repos = "https://cloud.r-project.org"); devtools::install()' || { echo "Error: devtools::install() failed" >&2; exit 1; }
fi

# Cores: when --sims N is set, use min(N, available_cores)
export SAVE_TO_CLUSTER_OUTPUT=1
if [[ -n "$N_SIMS_OVERRIDE" ]]; then
  export N_SIMS_OVERRIDE="$N_SIMS_OVERRIDE"
  NPROC=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
  CORES=$(( N_SIMS_OVERRIDE < NPROC ? N_SIMS_OVERRIDE : NPROC ))
  export CORES_OVERRIDE=$CORES
  echo "=== Sims: $N_SIMS_OVERRIDE, Cores: $CORES (max $NPROC) ==="
fi

# Run R script with timing; all sim study output goes to .out file (console stays free)
OUT_FILE="cluster_output/run_$(date +%Y%m%d_%H%M%S).out"
{
  echo "=== Starting simulation study at $(date) ==="
  echo "Log file: $OUT_FILE"
} >> "$OUT_FILE"
echo "Sim study output -> $OUT_FILE (tail -f $OUT_FILE to follow)"

START_EPOCH=$(date +%s)
set +e
Rscript inst/sim_study/sim_study.R "${EXTRA_ARGS[@]}" >> "$OUT_FILE" 2>&1
R_EXIT=$?
set -e
END_EPOCH=$(date +%s)
ELAPSED=$(( END_EPOCH - START_EPOCH ))

{
  echo ""
  echo "=== Run finished at $(date) ==="
  echo "Wall time: ${ELAPSED}s ($(( ELAPSED / 60 ))m $(( ELAPSED % 60 ))s)"
  echo "R exit code: $R_EXIT"
} >> "$OUT_FILE"

for timing in cluster_output/sim_study_timing_local.txt cluster_output/sim_study_timing_cluster.txt; do
  if [[ -f "$timing" ]]; then
    { echo ""; echo "=== R timing report ($(basename "$timing")) ==="; cat "$timing"; } >> "$OUT_FILE"
    break
  fi
done

if [[ $R_EXIT -ne 0 ]]; then
  echo "Error: R script exited with code $R_EXIT. Full log: $OUT_FILE" >&2
  exit $R_EXIT
fi

echo "Done. Wall time: ${ELAPSED}s. Exit: $R_EXIT. Full log: $OUT_FILE"
echo "Results in cluster_output/. Load latest .rds in R: x <- readRDS('cluster_output/sim_study_results_<timestamp>.rds')"
