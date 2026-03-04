#!/bin/bash
#SBATCH --job-name=PPDisentangle_Sim
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=128G

set -e

# ---- Resolve paths ----
# If submitted by the wrapper below, these are set as env vars.
# Otherwise figure them out from BASH_SOURCE.
if [ -z "$PP_SCRIPT_DIR" ]; then
    PP_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
if [ -z "$PP_PKG_ROOT" ]; then
    PP_PKG_ROOT="$(cd "$PP_SCRIPT_DIR/../.." && pwd)"
fi
# Parse --sims and --test from command line (overrides env var)
# Preserve from env when running inside SLURM job (no args)
PP_TEST="${PP_TEST:-}"
PP_SIMS="${PP_SIMS:-}"
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sims) PP_SIMS="$2"; shift 2 ;;
        --test) PP_TEST="--test"; shift ;;
        *) shift ;;
    esac
done
if [ -n "$PP_TEST" ]; then
    PP_SIMS="${PP_SIMS:-2}"
else
    PP_SIMS="${PP_SIMS:-100}"
fi
PP_CPUS="$PP_SIMS"

# ---- If on login node, submit via sbatch and exit ----
if [ -z "$SLURM_JOB_ID" ]; then
    cd "$PP_PKG_ROOT"
    
    # Generate a unique RUN_ID for this run
    RUN_ID=$(date +"%Y-%m-%d_%H-%M-%S")
    RUN_DIR="$PP_PKG_ROOT/cluster_output/logs/$RUN_ID"
    mkdir -p "$RUN_DIR"

    echo "Pulling latest code..."
    git pull origin main

    export PP_SCRIPT_DIR PP_PKG_ROOT PP_SIMS="$PP_SIMS" PP_CPUS="$PP_CPUS" PP_TEST="$PP_TEST" PP_RUN_ID="$RUN_ID"
    # Pass ATE_SEQUENTIAL=1 to avoid OOM during ATE estimation (sequential, lower memory)
    [ -n "$ATE_SEQUENTIAL" ] && export ATE_SEQUENTIAL

    echo "Submitting SLURM job: $PP_SIMS sims, $PP_CPUS cores${PP_TEST:+ (test mode)}"
    echo "Run ID: $RUN_ID"
    echo "Logs: $RUN_DIR"
    
    SBATCH_EXPORT="PP_SCRIPT_DIR,PP_PKG_ROOT,PP_SIMS,PP_CPUS,PP_TEST,PP_RUN_ID"
    [ -n "$ATE_SEQUENTIAL" ] && SBATCH_EXPORT="${SBATCH_EXPORT},ATE_SEQUENTIAL"
    
    JOB_ID=$(sbatch --parsable \
        --cpus-per-task="$PP_CPUS" \
        --output="$RUN_DIR/slurm_%j.out" \
        --error="$RUN_DIR/slurm_%j.err" \
        --export="$SBATCH_EXPORT" \
        "$PP_SCRIPT_DIR/run_sim_study.sh")

    # Create a symlink so the run can be found by Job ID
    ln -sfn "$RUN_ID" "$PP_PKG_ROOT/cluster_output/logs/$JOB_ID"
    # Create a symlink for the latest run
    ln -sfn "$RUN_ID" "$PP_PKG_ROOT/cluster_output/logs/latest"
    # Create a marker file with the Job ID inside the folder
    echo "$JOB_ID" > "$RUN_DIR/slurm_id.txt"

    echo "Job $JOB_ID submitted"
    echo "  tail -f $RUN_DIR/slurm_${JOB_ID}.out"
    exit 0
fi

# ---- Inside SLURM job ----
cd "$PP_PKG_ROOT"
# Use provided RUN_ID or fallback to timestamp
RUN_ID="${PP_RUN_ID:-$(date +"%Y%m%d_%H%M%S")}"
LOG_DIR="$PP_PKG_ROOT/cluster_output/logs/$RUN_ID"
mkdir -p "$LOG_DIR"

echo "=== PPDisentangle Sim Study ==="
echo "Job $SLURM_JOB_ID | $(date)"
echo "Run ID: $RUN_ID"
echo "Sims: $PP_SIMS | Cores: $PP_CPUS"
echo "Logs: $LOG_DIR"
echo "Pkg root: $PP_PKG_ROOT"
echo ""

module load R 2>/dev/null || true
module load UDUNITS 2>/dev/null || module load udunits2 2>/dev/null || true
module load GDAL 2>/dev/null || true
module load GEOS 2>/dev/null || true
module load PROJ 2>/dev/null || true

echo "Rscript: $(which Rscript)"
echo ""

Rscript "$PP_SCRIPT_DIR/sim_study.R" --cluster --sims "$PP_SIMS" $PP_TEST

echo "=== Done $(date) ==="
