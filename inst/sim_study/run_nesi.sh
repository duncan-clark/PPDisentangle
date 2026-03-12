#!/bin/bash
#SBATCH --job-name=PPDis_sim
#SBATCH --account=FIXME          # <-- your NeSI project code (e.g. nesi00123)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --time=72:00:00
#SBATCH --mem=100G

# ---- NeSI / Mahuika notes ----
# Partition is auto-selected by Slurm based on resources requested.
#   large:  max 3 days,  72 CPUs/node, 1500 MB/CPU (~108 GB/node)
#   long:   max 3 weeks, 72 CPUs/node
#   milan:  max 7 days, 256 CPUs/node (opt-in: --partition=milan)
#   bigmem: max 7 days,  72 CPUs/node, 6300 MB/CPU (~453 GB/node)
#
# Tune --cpus-per-task and --mem for your sim count:
#   50 sims  → --cpus-per-task=50  --mem=64G   --time=24:00:00
#   100 sims → --cpus-per-task=72  --mem=100G  --time=48:00:00  (or use milan)
#
# Submit from the package root:
#   cd /path/to/PPDisentangle
#   sbatch inst/sim_study/run_nesi.sh [--sims 50] [--test]

set -euo pipefail

# ---- Parse arguments ----
PP_SIMS="${PP_SIMS:-}"
PP_TEST=""
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sims) PP_SIMS="$2"; shift 2 ;;
        --test) PP_TEST="--test"; shift ;;
        *) shift ;;
    esac
done
PP_SIMS="${PP_SIMS:-100}"

# ---- Resolve paths ----
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# ---- If on login node (no SLURM_JOB_ID), submit via sbatch ----
if [ -z "${SLURM_JOB_ID:-}" ]; then
    cd "$PKG_ROOT"

    CPUS="$PP_SIMS"
    if [ "$CPUS" -gt 72 ]; then CPUS=72; fi

    echo "Submitting to NeSI: $PP_SIMS sims, $CPUS CPUs${PP_TEST:+, test mode}"

    JOB_ID=$(sbatch --parsable \
        --cpus-per-task="$CPUS" \
        --export=ALL,PP_SIMS="$PP_SIMS",PP_TEST="$PP_TEST" \
        --output="$PKG_ROOT/cluster_output/logs/%j/slurm.out" \
        --error="$PKG_ROOT/cluster_output/logs/%j/slurm.err" \
        "$SCRIPT_DIR/run_nesi.sh")

    mkdir -p "$PKG_ROOT/cluster_output/logs/$JOB_ID"

    echo "Job $JOB_ID submitted"
    echo "  Logs:    cluster_output/logs/$JOB_ID/"
    echo "  Results: cluster_output/results/$JOB_ID.rds"
    echo "  Monitor: tail -f cluster_output/logs/$JOB_ID/slurm.out"
    exit 0
fi

# ---- Inside SLURM job ----
cd "$PKG_ROOT"
JOB_ID="$SLURM_JOB_ID"
LOG_DIR="$PKG_ROOT/cluster_output/logs/$JOB_ID"
mkdir -p "$LOG_DIR"

echo "=== PPDisentangle Sim Study (NeSI) ==="
echo "Job $JOB_ID | $(date)"
echo "Sims: $PP_SIMS | CPUs: $SLURM_CPUS_PER_TASK"
echo "Node: $(hostname) | Partition: ${SLURM_JOB_PARTITION:-unknown}"
echo ""

# ---- Load modules ----
module purge
module load R
module load UDUNITS
module load GDAL
module load GEOS
module load PROJ

echo "R: $(which R) ($(R --version | head -1))"
echo ""

# Prevent nested threading from stealing CPUs
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# ---- Install package if needed ----
if ! Rscript -e 'library(PPDisentangle)' 2>/dev/null; then
    echo "Installing PPDisentangle from source..."
    R CMD INSTALL --no-multiarch "$PKG_ROOT" 2>&1 | tail -5
    echo ""
fi

# ---- Run simulation study ----
Rscript "$PKG_ROOT/inst/sim_study/sim_study.R" --cluster --sims "$PP_SIMS" $PP_TEST 2>&1

echo ""
echo "=== Done $(date) ==="
