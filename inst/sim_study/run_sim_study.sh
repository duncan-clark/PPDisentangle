#!/usr/bin/env bash
#
# Run the PPDisentangle simulation study (local or after pull/install).
# Creates cluster_output/, saves RData and timing report there.
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
#
# Examples:
#   ./inst/sim_study/run_sim_study.sh
#   ./inst/sim_study/run_sim_study.sh --sims 20
#   ./inst/sim_study/run_sim_study.sh --no-pull --no-install   # skip pull and install
#
# Outputs:
#   cluster_output/sim_study_results_<timestamp>.RData  (rehydrate in R: load(...))
#   cluster_output/sim_study_timing_local.txt
#   cluster_output/run_<timestamp>.log                  (stdout/stderr + timing)
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

# Optional: git pull
if $DO_PULL; then
  echo "=== Git pull ==="
  git pull || { echo "Warning: git pull failed" >&2; }
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

# Run R script with timing; capture stdout/stderr and exit code
LOG_FILE="cluster_output/run_$(date +%Y%m%d_%H%M%S).log"
echo "=== Starting simulation study at $(date) ===" | tee "$LOG_FILE"
echo "Log file: $LOG_FILE" | tee -a "$LOG_FILE"

START_EPOCH=$(date +%s)
set +e
Rscript inst/sim_study/sim_study.R "${EXTRA_ARGS[@]}" 2>&1 | tee -a "$LOG_FILE"
R_EXIT=$?
set -e
END_EPOCH=$(date +%s)
ELAPSED=$(( END_EPOCH - START_EPOCH ))

echo "" | tee -a "$LOG_FILE"
echo "=== Run finished at $(date) ===" | tee -a "$LOG_FILE"
echo "Wall time: ${ELAPSED}s ($(( ELAPSED / 60 ))m $(( ELAPSED % 60 ))s)" | tee -a "$LOG_FILE"
echo "R exit code: $R_EXIT" | tee -a "$LOG_FILE"

# Show timing report if R wrote one
if [[ -f cluster_output/sim_study_timing_local.txt ]]; then
  echo "" | tee -a "$LOG_FILE"
  echo "=== R timing report ===" | tee -a "$LOG_FILE"
  cat cluster_output/sim_study_timing_local.txt | tee -a "$LOG_FILE"
fi

if [[ $R_EXIT -ne 0 ]]; then
  echo "Error: R script exited with code $R_EXIT. Check $LOG_FILE and cluster_output/*.err for details." >&2
  exit $R_EXIT
fi

echo ""
echo "Results and timing are in cluster_output/. Load the latest .RData in R to inspect (e.g. load('cluster_output/sim_study_results_<timestamp>.RData'))."
