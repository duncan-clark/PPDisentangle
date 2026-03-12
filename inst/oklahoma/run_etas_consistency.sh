#!/bin/bash
#SBATCH --job-name=ETAS_Consist
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=15
#SBATCH --time=01:00:00
#SBATCH --mem=32G

set -e

# ---- Resolve paths ----
if [ -z "$PP_SCRIPT_DIR" ]; then
    PP_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
if [ -z "$PP_PKG_ROOT" ]; then
    PP_PKG_ROOT="$(cd "$PP_SCRIPT_DIR/../.." && pwd)"
fi

PP_TEST=""
PP_REPS=""
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --test) PP_TEST="--test"; shift ;;
        --reps) PP_REPS="$2"; shift 2 ;;
        *) shift ;;
    esac
done

if [ -n "$PP_TEST" ]; then
    PP_REPS="${PP_REPS:-3}"
else
    PP_REPS="${PP_REPS:-15}"
fi

# ---- If on login node, submit via sbatch ----
if [ -z "$SLURM_JOB_ID" ]; then
    cd "$PP_PKG_ROOT"

    echo "Pulling latest code..."
    git pull origin main 2>/dev/null || true

    export PP_SCRIPT_DIR PP_PKG_ROOT PP_TEST PP_REPS

    SBATCH_EXPORT="PP_SCRIPT_DIR,PP_PKG_ROOT,PP_TEST,PP_REPS"

    CPUS="${PP_REPS}"
    if [ "$CPUS" -gt 30 ]; then CPUS=30; fi

    JOB_ID=$(sbatch --parsable \
        --cpus-per-task="$CPUS" \
        --output="$PP_PKG_ROOT/cluster_output/logs/%j/slurm.out" \
        --error="$PP_PKG_ROOT/cluster_output/logs/%j/slurm.err" \
        --export="$SBATCH_EXPORT" \
        "$PP_SCRIPT_DIR/run_etas_consistency.sh")

    mkdir -p "$PP_PKG_ROOT/cluster_output/logs/$JOB_ID"

    echo "Job $JOB_ID submitted ($PP_REPS reps${PP_TEST:+, test mode})"
    echo "  Logs:    cluster_output/logs/$JOB_ID/"
    echo "  Results: cluster_output/results/etas_consistency_$JOB_ID.rds"
    echo "  tail -f cluster_output/logs/$JOB_ID/slurm.out"
    exit 0
fi

# ---- Inside SLURM job ----
cd "$PP_PKG_ROOT"
JOB_ID="${SLURM_JOB_ID:-local_$(date +%Y%m%d_%H%M%S)}"
LOG_DIR="$PP_PKG_ROOT/cluster_output/logs/$JOB_ID"
mkdir -p "$LOG_DIR"

echo "=== ETAS Consistency Study ==="
echo "Job $JOB_ID | $(date)"
echo "Reps: $PP_REPS | Cores: ${SLURM_CPUS_PER_TASK:-1}"
echo "Logs: $LOG_DIR"
echo ""

module load R 2>/dev/null || true
module load UDUNITS 2>/dev/null || module load udunits2 2>/dev/null || true
module load GDAL 2>/dev/null || true
module load GEOS 2>/dev/null || true
module load PROJ 2>/dev/null || true

echo "Rscript: $(which Rscript)"
echo ""

Rscript "$PP_SCRIPT_DIR/etas_consistency_study.R" --cluster --reps "$PP_REPS" $PP_TEST

echo "=== Done $(date) ==="
