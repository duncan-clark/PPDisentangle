#!/bin/bash
#SBATCH --job-name=PPDisentangle_SimStudy
#SBATCH --output=cluster_output_profiled/logs/slurm_%j.out
#SBATCH --error=cluster_output_profiled/logs/slurm_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=100
#SBATCH --time=24:00:00
#SBATCH --mem=128G

# Get the directory of this script so paths work regardless of where it's called from
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Parse command line arguments
SIMS=100

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sims) SIMS="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Cap cores to sims so we don't waste resources
if [ "$SIMS" -le 100 ]; then
    CPUS="$SIMS"
else
    CPUS=100
fi

# If not already inside a SLURM job, submit via sbatch and exit
if [ -z "$SLURM_JOB_ID" ]; then
    cd "$PKG_ROOT"
    mkdir -p cluster_output_profiled/logs

    # Pull latest code before submitting
    echo "Pulling latest changes from git..."
    git pull origin main

    echo "Submitting SLURM job: $SIMS sims, $CPUS cores..."
    JOB_ID=$(sbatch --parsable \
        --cpus-per-task="$CPUS" \
        --export=ALL \
        "$SCRIPT_DIR/run_sim_study.sh" --sims "$SIMS")

    echo "Submitted job $JOB_ID"
    echo "  SLURM logs: cluster_output_profiled/logs/slurm_${JOB_ID}.out"
    echo "  Worker logs: cluster_output_profiled/logs/"
    echo "  Check status: squeue -j $JOB_ID"
    echo "  Follow output: tail -f cluster_output_profiled/logs/slurm_${JOB_ID}.out"
    exit 0
fi

# --- Everything below runs inside the SLURM job ---

cd "$PKG_ROOT"
mkdir -p cluster_output_profiled/logs

echo "=== PPDisentangle Simulation Study ==="
echo "Started at: $(date)"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Sims: $SIMS  Cores: $CPUS"
echo "Package root: $PKG_ROOT"
echo ""

# Load necessary modules
echo "--- Loading modules ---"
module load R 2>/dev/null || echo "Warning: 'module load R' failed. Attempting to continue..."
module load UDUNITS 2>/dev/null || module load udunits2 2>/dev/null || module load udunits 2>/dev/null || echo "Warning: could not load UDUNITS module"
module load GDAL 2>/dev/null || true
module load GEOS 2>/dev/null || true
module load PROJ 2>/dev/null || true
echo ""

# Verify Rscript is available
if ! command -v Rscript &>/dev/null; then
    echo "ERROR: Rscript not found after module load. Check 'module avail R' on your cluster."
    exit 1
fi
echo "Using Rscript: $(which Rscript)"
echo ""

# Run the simulation study
echo "--- Starting Rscript ---"
Rscript "$SCRIPT_DIR/sim_study_profile.R" --cluster --sims "$SIMS"
EXIT_CODE=$?

echo ""
echo "=== Job completed at $(date) with exit code $EXIT_CODE ==="
