#!/bin/bash
# Gambia radii study: runs IPD and non-IPD radii analyses in sequence.
# Both R scripts use parallel::mclapply when SLURM_CPUS_PER_TASK > 1.
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
    mkdir -p inst/gambia/output
    mkdir -p inst/gambia/output/logs

    export PP_SCRIPT_DIR PP_PKG_ROOT

    echo "Submitting Gambia radii job (IPD + non-IPD)..."
    JOB_ID=$(sbatch --parsable \
        --output="inst/gambia/output/logs/slurm_%j.out" \
        --error="inst/gambia/output/logs/slurm_%j.err" \
        --export=PP_SCRIPT_DIR,PP_PKG_ROOT \
        "$PP_SCRIPT_DIR/run_gambia_radii.sh")

    echo "Job $JOB_ID submitted"
    echo "  tail -f inst/gambia/output/logs/slurm_${JOB_ID}.out"
    exit 0
fi

# ---- Inside SLURM job ----
cd "$PP_PKG_ROOT"
mkdir -p inst/gambia/output
mkdir -p inst/gambia/output/logs

N_CORES="${SLURM_CPUS_PER_TASK:-8}"
export SLURM_CPUS_PER_TASK="$N_CORES"

echo "=== Gambia Radii Study ==="
echo "Job $SLURM_JOB_ID | $(date)"
echo "Cores: $N_CORES"
echo "Pkg root: $PP_PKG_ROOT"
echo ""

module load R 2>/dev/null || true
module load UDUNITS 2>/dev/null || module load udunits2 2>/dev/null || true
module load GDAL 2>/dev/null || true
module load GEOS 2>/dev/null || true
module load PROJ 2>/dev/null || true

echo "Rscript: $(which Rscript)"
echo ""

# Run IPD radii study
echo "--- IPD radii study ---"
Rscript "$PP_SCRIPT_DIR/gambia_radii_study.R"

# Run non-IPD radii study
echo ""
echo "--- Non-IPD radii study ---"
Rscript "$PP_SCRIPT_DIR/gambia_radii_study_nonIPD.R"

echo ""
echo "=== Done $(date) ==="
