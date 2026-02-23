#!/bin/sh
#SBATCH --job-name=PPDisentangle_sim
#SBATCH --output=cluster_output/PPDisentangle_sim_%j.out
#SBATCH --error=cluster_output/PPDisentangle_sim_%j.err
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=200G

# NeSI (Mahuika) notes:
# - Do not set --partition: let Slurm choose (default "large").
# - "milan" is opt-in only; requesting it can cause "configuration not available".
# - large: 72 CPUs/node, 105 GB/node, 3 days max → 128 CPUs = 2 nodes, 210 GB max; we request 200G.
# - For more memory (e.g. 256G) use: --partition=bigmem (if available), same --mem and --cpus-per-task.
# - Very large CPU/RAM requests often wait longer in queue.
#
# Overrides (set before sbatch or in the script):
#   CORES_OVERRIDE=<int>        force N_CORES in the R script (default: SLURM_CPUS_PER_TASK)
#   N_SIMS_OVERRIDE=<int>       number of datasets (SIM_SIZE) and can drive core count when used with shell script
#   SIM_SIZE_OVERRIDE=<int>     number of simulated datasets (overrides cluster default 100)
#
# Submit from the package root (directory containing R/, inst/, DESCRIPTION):
#   mkdir -p cluster_output && sbatch inst/sim_study/run_PPDisentangle_sim.slurm

module load R
module load GDAL
module load PROJ
module load GEOS
module load udunits

mkdir -p cluster_output
echo "Starting job at $(date)"
# Avoid nested BLAS/OpenMP threads stealing CPUs (and helps fork safety).
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Optional: reduce workers if hitting connection limits (e.g. 120 for 128 CPUs)
# export SIM_WORKERS=120

Rscript inst/sim_study/sim_study.R 2>&1
EXIT_CODE=$?
echo "Finished job at $(date)"
exit $EXIT_CODE
