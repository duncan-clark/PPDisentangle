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
CPUS=100

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sims) SIMS="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [ "$SIMS" -le 16 ]; then
    CPUS="$SIMS"
fi

# 1. Update the code from the repository
echo "Pulling latest changes from git..."
cd "$PKG_ROOT"
git pull origin main

# 2. Load necessary modules
echo "Loading modules..."
module load R 2>/dev/null || echo "Warning: 'module load R' failed. Attempting to continue..."
module load UDUNITS 2>/dev/null || module load udunits2 2>/dev/null || module load udunits 2>/dev/null || echo "Warning: could not load UDUNITS module"
module load GDAL 2>/dev/null || true
module load GEOS 2>/dev/null || true
module load PROJ 2>/dev/null || true

# 3. Verify Rscript is available
if ! command -v Rscript &>/dev/null; then
    echo "ERROR: Rscript not found after module load. Check 'module avail R' on your cluster."
    exit 1
fi
echo "Using Rscript: $(which Rscript)"

# 4. Create necessary directories
mkdir -p cluster_output_profiled/logs

# 5. Run the simulation study
echo "Starting simulation study with $SIMS simulations on $CPUS cores..."
Rscript "$SCRIPT_DIR/sim_study_profile.R" --cluster --sims "$SIMS"

echo "Job completed at $(date)"
