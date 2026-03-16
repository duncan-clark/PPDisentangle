#!/bin/bash
#SBATCH --job-name=PPDis_sim
#SBATCH --account=uoo04008
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --mem=200G

# ---- NeSI / Mahuika notes ----
# Partition is auto-selected by Slurm based on resources requested.
#   large:  max 3 days,  72 CPUs/node, 1500 MB/CPU (~108 GB/node)
#   long:   max 3 weeks, 72 CPUs/node
#   milan:  max 7 days, 256 CPUs/node (opt-in: --partition=milan)
#   bigmem: max 7 days,  72 CPUs/node, 6300 MB/CPU (~453 GB/node)
#
# Tune --cpus-per-task and --mem for your sim count:
#   50 sims  → --cpus-per-task=50  --mem=64G   --time=24:00:00
#   100 sims → --cpus-per-task=72  --mem=140G  --time=48:00:00  (or use milan)
#
# Submit from the package root:
#   cd /path/to/PPDisentangle
#   sbatch inst/sim_study/run_nesi.sh [--sims 50] [--test]
#
#   --test: uses minimal resources (15 min, 4G) for fast scheduling

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
# When running under SLURM, BASH_SOURCE points to the copied script in /var/spool,
# so we use PKG_ROOT from --export (set at submit time) or SLURM_SUBMIT_DIR.
if [ -n "${SLURM_JOB_ID:-}" ] && [ -n "${PKG_ROOT:-}" ] && [ -d "$PKG_ROOT" ]; then
    :  # PKG_ROOT already set by sbatch --export
elif [ -n "${SLURM_JOB_ID:-}" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -d "$SLURM_SUBMIT_DIR" ]; then
    PKG_ROOT="$SLURM_SUBMIT_DIR"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
fi

# ---- If on login node (no SLURM_JOB_ID), submit via sbatch ----
if [ -z "${SLURM_JOB_ID:-}" ]; then
    cd "$PKG_ROOT"
    git pull origin main 2>/dev/null || true

    CPUS="$PP_SIMS"
    # Robust defaults for the heavy ATE stage.
    if [ -z "${PP_SKIP_CRAZY_PARAMS:-}" ]; then
        PP_SKIP_CRAZY_PARAMS=1
    fi
    if [ -z "${PP_ATE_WORKERS:-}" ]; then
        # Conservative default for memory-heavy ATE fitting.
        PP_ATE_WORKERS=$(( CPUS / 8 ))
        [ "$PP_ATE_WORKERS" -lt 4 ] && PP_ATE_WORKERS=4
        [ "$PP_ATE_WORKERS" -gt 12 ] && PP_ATE_WORKERS=12
    fi
    if [ -z "${PP_ATE_BATCH_SIZE:-}" ]; then
        PP_ATE_BATCH_SIZE=$(( PP_ATE_WORKERS * 2 ))
    fi
    EXTRA_SBATCH=""

    # Test mode: minimal resources for fast scheduling
    if [ -n "$PP_TEST" ]; then
        EXTRA_SBATCH="--time=00:15:00 --mem=4G"
    # large partition: max 72 CPUs/node; milan: up to 256 CPUs/node
    elif [ "$CPUS" -gt 72 ]; then
        if [ "$CPUS" -gt 256 ]; then
            echo "ERROR: max 256 CPUs per node (milan). Reduce --sims or split jobs."
            exit 1
        fi
        EXTRA_SBATCH="--partition=milan"
        echo "Note: >72 CPUs requested, using Milan partition (256 CPUs/node)"
    fi

    echo "Submitting to NeSI: $PP_SIMS sims, $CPUS CPUs (1 sim/core)${PP_TEST:+, test mode (15 min, 4G)}"
    echo "Robust defaults: PP_SKIP_CRAZY_PARAMS=$PP_SKIP_CRAZY_PARAMS  PP_ATE_WORKERS=$PP_ATE_WORKERS  PP_ATE_BATCH_SIZE=$PP_ATE_BATCH_SIZE"

    mkdir -p "$PKG_ROOT/cluster_output"

    SBATCH_EXPORT="ALL,PP_SIMS=$PP_SIMS,PP_TEST=$PP_TEST,PKG_ROOT=$PKG_ROOT"
    [ -n "${ATE_SEQUENTIAL:-}" ] && SBATCH_EXPORT="${SBATCH_EXPORT},ATE_SEQUENTIAL=$ATE_SEQUENTIAL"
    [ -n "${PP_LOG_MEMORY:-}" ] && SBATCH_EXPORT="${SBATCH_EXPORT},PP_LOG_MEMORY=$PP_LOG_MEMORY"
    SBATCH_EXPORT="${SBATCH_EXPORT},PP_SKIP_CRAZY_PARAMS=$PP_SKIP_CRAZY_PARAMS,PP_ATE_WORKERS=$PP_ATE_WORKERS,PP_ATE_BATCH_SIZE=$PP_ATE_BATCH_SIZE"

    JOB_ID=$(sbatch --parsable \
        --cpus-per-task="$CPUS" \
        $EXTRA_SBATCH \
        --export="$SBATCH_EXPORT" \
        --output="$PKG_ROOT/cluster_output/%j_slurm.out" \
        --error="$PKG_ROOT/cluster_output/%j_slurm.err" \
        "$SCRIPT_DIR/run_nesi.sh")

    echo "Job $JOB_ID submitted"
    echo "  Results: cluster_output/$JOB_ID.rds"
    echo "  Log:     cluster_output/$JOB_ID.log"
    echo "  SLURM:   cluster_output/${JOB_ID}_slurm.out"
    exit 0
fi

# ---- Inside SLURM job ----
cd "$PKG_ROOT"
mkdir -p "$PKG_ROOT/cluster_output"

# Re-apply robust defaults inside allocation if not exported explicitly.
if [ -z "${PP_SKIP_CRAZY_PARAMS:-}" ]; then
    PP_SKIP_CRAZY_PARAMS=1
fi
if [ -z "${PP_ATE_WORKERS:-}" ]; then
    PP_ATE_WORKERS=$(( SLURM_CPUS_PER_TASK / 8 ))
    [ "$PP_ATE_WORKERS" -lt 4 ] && PP_ATE_WORKERS=4
    [ "$PP_ATE_WORKERS" -gt 12 ] && PP_ATE_WORKERS=12
fi
if [ -z "${PP_ATE_BATCH_SIZE:-}" ]; then
    PP_ATE_BATCH_SIZE=$(( PP_ATE_WORKERS * 2 ))
fi
export PP_SKIP_CRAZY_PARAMS PP_ATE_WORKERS PP_ATE_BATCH_SIZE

echo "=== PPDisentangle Sim Study (NeSI) ==="
echo "Job $SLURM_JOB_ID | $(date)"
echo "Sims: $PP_SIMS | CPUs: $SLURM_CPUS_PER_TASK"
echo "PP_SKIP_CRAZY_PARAMS=$PP_SKIP_CRAZY_PARAMS | PP_ATE_WORKERS=$PP_ATE_WORKERS | PP_ATE_BATCH_SIZE=$PP_ATE_BATCH_SIZE"
echo "Node: $(hostname) | Partition: ${SLURM_JOB_PARTITION:-unknown}"
echo ""

# ---- Load modules ----
# R-Geo bundles R + GDAL + GEOS + PROJ + UDUNITS (required for sf, terra, spatstat)
module --force purge

# Some NeSI stacks require architecture/compiler parents (e.g., NeSI/zen3)
# before R-Geo can be loaded. Try robust fallbacks in order.
if ! module load R-Geo/4.3.2-foss-2023a; then
    echo "Primary load failed: R-Geo/4.3.2-foss-2023a"
    echo "Trying NeSI/zen3 -> R-Geo/4.3.2-foss-2023a ..."
    if module load NeSI/zen3 && module load R-Geo/4.3.2-foss-2023a; then
        :
    else
        echo "Trying generic R-Geo module ..."
        if ! module load R-Geo; then
            echo "ERROR: Failed to load R-Geo module."
            echo "Run 'module spider R-Geo' and update run_nesi.sh accordingly."
            module spider R-Geo || true
            exit 1
        fi
    fi
fi

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
