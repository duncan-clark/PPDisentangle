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
PP_TEST=""
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
    mkdir -p cluster_output/logs

    echo "Pulling latest code..."
    git pull origin main

    export PP_SCRIPT_DIR PP_PKG_ROOT PP_SIMS="$PP_SIMS" PP_CPUS="$PP_CPUS" PP_TEST="$PP_TEST"
    # Pass ATE_SEQUENTIAL=1 to avoid OOM during ATE estimation (sequential, lower memory)
    [ -n "$ATE_SEQUENTIAL" ] && export ATE_SEQUENTIAL

    echo "Submitting SLURM job: $PP_SIMS sims, $PP_CPUS cores${PP_TEST:+ (test mode)}"
    SBATCH_EXPORT="PP_SCRIPT_DIR,PP_PKG_ROOT,PP_SIMS,PP_CPUS,PP_TEST"
    [ -n "$ATE_SEQUENTIAL" ] && SBATCH_EXPORT="${SBATCH_EXPORT},ATE_SEQUENTIAL"
    JOB_ID=$(sbatch --parsable \
        --cpus-per-task="$PP_CPUS" \
        --output="cluster_output/logs/slurm_%j.out" \
        --error="cluster_output/logs/slurm_%j.err" \
        --export="$SBATCH_EXPORT" \
        "$PP_SCRIPT_DIR/run_sim_study.sh")

    echo "Job $JOB_ID submitted"
    echo "  tail -f cluster_output/logs/slurm_${JOB_ID}.out"
    exit 0
fi

# ---- Inside SLURM job ----
cd "$PP_PKG_ROOT"
mkdir -p cluster_output/logs

echo "=== PPDisentangle Sim Study ==="
echo "Job $SLURM_JOB_ID | $(date)"
echo "Sims: $PP_SIMS | Cores: $PP_CPUS"
echo "Pkg root: $PP_PKG_ROOT"
echo "Script dir: $PP_SCRIPT_DIR"
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
