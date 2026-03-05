#!/bin/bash
# Gambia radii study: runs IPD and non-IPD analyses in sequence.
# Override cores: sbatch --cpus-per-task=16 run_gambia_radii.sh
#
#SBATCH --job-name=Gambia_Radii
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=64G

set -e

# ---- Resolve paths ----
if [ -z "$PP_SCRIPT_DIR" ]; then
    PP_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
if [ -z "$PP_PKG_ROOT" ]; then
    PP_PKG_ROOT="$(cd "$PP_SCRIPT_DIR/../.." && pwd)"
fi

# ---- If on login node, submit via sbatch and exit ----
if [ -z "$SLURM_JOB_ID" ]; then
    cd "$PP_PKG_ROOT"

    echo "Pulling latest code..."
    git pull origin main

    export PP_SCRIPT_DIR PP_PKG_ROOT

    JOB_ID=$(sbatch --parsable \
        --output="$PP_PKG_ROOT/cluster_output/logs/%j/slurm.out" \
        --error="$PP_PKG_ROOT/cluster_output/logs/%j/slurm.err" \
        --export=PP_SCRIPT_DIR,PP_PKG_ROOT \
        "$PP_SCRIPT_DIR/run_gambia_radii.sh")

    mkdir -p "$PP_PKG_ROOT/cluster_output/logs/$JOB_ID"

    echo "Job $JOB_ID submitted (Gambia IPD + non-IPD)"
    echo "  Logs:    cluster_output/logs/$JOB_ID/"
    echo "  Results: cluster_output/results/gambia_*.rds"
    echo "  tail -f cluster_output/logs/$JOB_ID/slurm.out"
    exit 0
fi

# ---- Inside SLURM job ----
cd "$PP_PKG_ROOT"
JOB_ID="${SLURM_JOB_ID:-local_$(date +%Y%m%d_%H%M%S)}"
LOG_DIR="$PP_PKG_ROOT/cluster_output/logs/$JOB_ID"
RESULTS_DIR="$PP_PKG_ROOT/cluster_output/results"
mkdir -p "$LOG_DIR" "$RESULTS_DIR"

export PP_OUT_DIR="$RESULTS_DIR"

N_CORES="${SLURM_CPUS_PER_TASK:-8}"
export SLURM_CPUS_PER_TASK="$N_CORES"

echo "=== Gambia Radii Study ==="
echo "Job $JOB_ID | $(date)"
echo "Cores: $N_CORES"
echo "Results: $RESULTS_DIR"
echo ""

module load R 2>/dev/null || true
module load UDUNITS 2>/dev/null || module load udunits2 2>/dev/null || true
module load GDAL 2>/dev/null || true
module load GEOS 2>/dev/null || true
module load PROJ 2>/dev/null || true

echo "Rscript: $(which Rscript)"
echo ""

echo "--- IPD radii study ---"
Rscript "$PP_SCRIPT_DIR/gambia_radii_study.R" 2>&1 | tee "$LOG_DIR/gambia_ipd.log"

echo ""
echo "--- Non-IPD radii study ---"
Rscript "$PP_SCRIPT_DIR/gambia_radii_study_nonIPD.R" 2>&1 | tee "$LOG_DIR/gambia_nonipd.log"

echo ""
echo "=== Done $(date) ==="
