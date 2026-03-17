#!/bin/bash
#SBATCH --job-name=PPDis_sim
#SBATCH --account=uoo04008
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --mem=48G

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
PP_MEM="${PP_MEM:-}"
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sims) PP_SIMS="$2"; shift 2 ;;
        --test) PP_TEST="--test"; shift ;;
        --mem) PP_MEM="$2"; shift 2 ;;
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
    MEM_PER_CORE_GB="${PP_MEM_PER_CORE_GB:-1}"
    MEM_MAX_GB="${PP_MEM_MAX_GB:-48}"
    MEM_MIN_GB="${PP_MEM_MIN_GB:-8}"
    if [ -z "$PP_MEM" ]; then
        MEM_GB=$(( CPUS * MEM_PER_CORE_GB ))
        [ "$MEM_GB" -lt "$MEM_MIN_GB" ] && MEM_GB="$MEM_MIN_GB"
        [ "$MEM_GB" -gt "$MEM_MAX_GB" ] && MEM_GB="$MEM_MAX_GB"
        PP_MEM="${MEM_GB}G"
    fi
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
        EXTRA_SBATCH="--time=00:15:00"
        if [ -z "${PP_MEM:-}" ]; then
            PP_MEM="4G"
        fi
    # large partition: max 72 CPUs/node; milan: up to 256 CPUs/node
    elif [ "$CPUS" -gt 72 ]; then
        if [ "$CPUS" -gt 256 ]; then
            echo "ERROR: max 256 CPUs per node (milan). Reduce --sims or split jobs."
            exit 1
        fi
        EXTRA_SBATCH="--partition=milan"
        echo "Note: >72 CPUs requested, using Milan partition (256 CPUs/node)"
    fi

    echo "Submitting to NeSI: $PP_SIMS sims, $CPUS CPUs (1 sim/core), mem=$PP_MEM${PP_TEST:+, test mode (15 min)}"
    echo "Robust defaults: PP_SKIP_CRAZY_PARAMS=$PP_SKIP_CRAZY_PARAMS  PP_ATE_WORKERS=$PP_ATE_WORKERS  PP_ATE_BATCH_SIZE=$PP_ATE_BATCH_SIZE"

    mkdir -p "$PKG_ROOT/cluster_output"

    SBATCH_EXPORT="ALL,PP_SIMS=$PP_SIMS,PP_TEST=$PP_TEST,PKG_ROOT=$PKG_ROOT"
    [ -n "${ATE_SEQUENTIAL:-}" ] && SBATCH_EXPORT="${SBATCH_EXPORT},ATE_SEQUENTIAL=$ATE_SEQUENTIAL"
    [ -n "${PP_LOG_MEMORY:-}" ] && SBATCH_EXPORT="${SBATCH_EXPORT},PP_LOG_MEMORY=$PP_LOG_MEMORY"
    SBATCH_EXPORT="${SBATCH_EXPORT},PP_SKIP_CRAZY_PARAMS=$PP_SKIP_CRAZY_PARAMS,PP_ATE_WORKERS=$PP_ATE_WORKERS,PP_ATE_BATCH_SIZE=$PP_ATE_BATCH_SIZE"

    JOB_ID=$(sbatch --parsable \
        --cpus-per-task="$CPUS" \
        --mem="$PP_MEM" \
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
TARGET_R_GEO="${PP_R_GEO_MODULE:-R-Geo/4.3.2-foss-2023a}"
echo "Requested R-Geo module: $TARGET_R_GEO"

try_load_rgeo() {
    local mod="$1"
    if module load "$mod" >/dev/null 2>&1; then
        echo "Loaded module: $mod"
        return 0
    fi
    return 1
}

try_load_rgeo_with_toolchain() {
    local mod="$1"
    # Example: R-Geo/4.3.2-foss-2023a -> toolchain foss/2023a
    local tail tc_name tc_ver tc_mod
    tail="$(echo "$mod" | awk -F'-' '{print $(NF-1) "-" $NF}')"
    tc_name="${tail%-*}"
    tc_ver="${tail#*-}"
    tc_mod="${tc_name}/${tc_ver}"
    if [ -n "$tc_name" ] && [ -n "$tc_ver" ] && [ "$tc_name" != "$tc_ver" ]; then
        echo "Trying toolchain + R-Geo: $tc_mod -> $mod"
        if module load "$tc_mod" >/dev/null 2>&1 && module load "$mod" >/dev/null 2>&1; then
            echo "Loaded module chain: $tc_mod + $mod"
            return 0
        fi
    fi
    return 1
}

if ! try_load_rgeo "$TARGET_R_GEO"; then
    echo "Primary load failed: $TARGET_R_GEO"
    echo "Trying toolchain-matched load for $TARGET_R_GEO ..."
    module --force purge
    if try_load_rgeo_with_toolchain "$TARGET_R_GEO"; then
        :
    else
    echo "Trying NeSI/zen3 -> $TARGET_R_GEO ..."
    module --force purge
    if module load NeSI/zen3 && try_load_rgeo "$TARGET_R_GEO"; then
        :
    else
        echo "Trying discovered R-Geo/* modules from module avail ..."
        module --force purge
        mapfile -t R_GEO_CANDIDATES < <(module -t avail R-Geo 2>&1 | awk '/^R-Geo\//{print $1}' | sort -Vr | uniq)
        LOADED=0
        for cand in "${R_GEO_CANDIDATES[@]}"; do
            if try_load_rgeo "$cand"; then
                LOADED=1
                break
            fi
            module --force purge
            if try_load_rgeo_with_toolchain "$cand"; then
                LOADED=1
                break
            fi
            module --force purge
            if module load NeSI/zen3 && try_load_rgeo "$cand"; then
                LOADED=1
                break
            fi
            module --force purge
        done
        if [ "$LOADED" -ne 1 ]; then
            echo "Trying generic R-Geo module ..."
            if ! try_load_rgeo R-Geo; then
                module --force purge
                if ! (module load NeSI/zen3 && try_load_rgeo R-Geo); then
                    echo "ERROR: Failed to load any R-Geo module."
                    echo "Set PP_R_GEO_MODULE explicitly, e.g."
                    echo "  sbatch --export=ALL,PP_R_GEO_MODULE=R-Geo/<version> inst/sim_study/run_nesi.sh --sims 32"
                    echo "Diagnostics:"
                    module spider R-Geo || true
                    exit 1
                fi
            fi
        fi
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
