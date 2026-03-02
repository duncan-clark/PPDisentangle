#!/bin/bash
#SBATCH --job-name=PPDisentangle_SimStudy
#SBATCH --output=cluster_output_profiled/logs/slurm_%j.out
#SBATCH --error=cluster_output_profiled/logs/slurm_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=100
#SBATCH --time=24:00:00
#SBATCH --mem=128G

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

# If sims is 16, we also want to reduce the requested CPUs for faster scheduling
if [ "$SIMS" -eq 16 ]; then
    CPUS=16
fi

# 1. Update the code from the repository
echo "Pulling latest changes from git..."
git pull origin main

# 2. Create necessary directories
mkdir -p cluster_output_profiled/logs

# 3. Run the simulation study
echo "Starting simulation study with $SIMS simulations on $CPUS cores..."
# Note: Rscript is called from the directory where the script is located
Rscript sim_study_profile.R --cluster --sims "$SIMS"

echo "Job completed at $(date)"
