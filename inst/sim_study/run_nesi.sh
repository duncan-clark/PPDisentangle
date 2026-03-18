#!/bin/bash
#SBATCH --job-name=PPDis_sim
#SBATCH --account=uoo04008
#SBATCH --nodes=1
#SBATCH --ntasks=1

set -euo pipefail

# Minimal interface:
#   --sims N    number of simulations / cores
#   --test      quick test profile

PP_SIMS=100
PP_TEST=""
HAVE_SIMS_ARG=0

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --sims) PP_SIMS="$2"; HAVE_SIMS_ARG=1; shift 2 ;;
        --test) PP_TEST="--test"; shift ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

if [ -n "$PP_TEST" ] && [ "$HAVE_SIMS_ARG" -eq 0 ]; then
    PP_SIMS=8
fi

if [ -n "${SLURM_JOB_ID:-}" ] && [ -n "${PKG_ROOT:-}" ] && [ -d "$PKG_ROOT" ]; then
    :
elif [ -n "${SLURM_JOB_ID:-}" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -d "$SLURM_SUBMIT_DIR" ]; then
    PKG_ROOT="$SLURM_SUBMIT_DIR"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
fi

if [ -z "${SLURM_JOB_ID:-}" ]; then
    cd "$PKG_ROOT"
    git pull origin main 2>/dev/null || true
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    mkdir -p "$PKG_ROOT/cluster_output"

    CPUS="$PP_SIMS"
    if [ "$CPUS" -lt 1 ]; then
        echo "ERROR: --sims must be >= 1"
        exit 1
    fi

    # Automatic resources.
    if [ -n "$PP_TEST" ]; then
        SB_TIME="00:20:00"
        SB_MEM="8G"
    else
        MEM_GB=$(( CPUS ))
        [ "$MEM_GB" -lt 8 ] && MEM_GB=8
        [ "$MEM_GB" -gt 48 ] && MEM_GB=48
        SB_MEM="${MEM_GB}G"
        SB_TIME="72:00:00"
    fi

    EXTRA_SBATCH=""
    if [ "$CPUS" -gt 72 ]; then
        if [ "$CPUS" -gt 256 ]; then
            echo "ERROR: max 256 CPUs per node (milan). Reduce --sims or split jobs."
            exit 1
        fi
        EXTRA_SBATCH="--partition=milan"
        echo "Note: >72 CPUs requested, using Milan partition."
    fi

    echo "Submitting to NeSI: sims=$PP_SIMS cpus=$CPUS mem=$SB_MEM time=$SB_TIME ${PP_TEST:+(test)}"

    SBATCH_EXPORT="ALL,PP_SIMS=$PP_SIMS,PP_TEST=$PP_TEST,PKG_ROOT=$PKG_ROOT"
    JOB_ID=$(sbatch --parsable \
        --cpus-per-task="$CPUS" \
        --mem="$SB_MEM" \
        --time="$SB_TIME" \
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

cd "$PKG_ROOT"
mkdir -p "$PKG_ROOT/cluster_output"

echo "=== PPDisentangle Sim Study (NeSI) ==="
echo "Job $SLURM_JOB_ID | $(date)"
echo "Sims: $PP_SIMS | CPUs: $SLURM_CPUS_PER_TASK"
echo "Node: $(hostname) | Partition: ${SLURM_JOB_PARTITION:-unknown}"
echo ""

module --force purge
TARGET_R_GEO="R-Geo/4.3.2-foss-2023a"
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
    local tail tc_name tc_ver tc_mod
    tail="$(echo "$mod" | awk -F'-' '{print $(NF-1) "-" $NF}')"
    tc_name="${tail%-*}"
    tc_ver="${tail#*-}"
    tc_mod="${tc_name}/${tc_ver}"
    if [ -n "$tc_name" ] && [ -n "$tc_ver" ] && [ "$tc_name" != "$tc_ver" ]; then
        if module load "$tc_mod" >/dev/null 2>&1 && module load "$mod" >/dev/null 2>&1; then
            echo "Loaded module chain: $tc_mod + $mod"
            return 0
        fi
    fi
    return 1
}

if ! try_load_rgeo "$TARGET_R_GEO"; then
    module --force purge
    if ! try_load_rgeo_with_toolchain "$TARGET_R_GEO"; then
        module --force purge
        if ! (module load NeSI/zen3 >/dev/null 2>&1 && try_load_rgeo "$TARGET_R_GEO"); then
            module --force purge
            echo "ERROR: Failed to load $TARGET_R_GEO."
            module spider R-Geo || true
            exit 1
        fi
    fi
fi

ensure_r_binaries() {
    if command -v R >/dev/null 2>&1 && command -v Rscript >/dev/null 2>&1; then
        return 0
    fi
    local target_r="R/4.3.2-foss-2023a"
    if module load "$target_r" >/dev/null 2>&1; then
        :
    elif module load NeSI/zen3 >/dev/null 2>&1 && module load "$target_r" >/dev/null 2>&1; then
        :
    elif module load foss/2023a >/dev/null 2>&1 && module load "$target_r" >/dev/null 2>&1; then
        :
    elif module load NeSI/zen3 >/dev/null 2>&1 && module load foss/2023a >/dev/null 2>&1 && module load "$target_r" >/dev/null 2>&1; then
        :
    else
        return 1
    fi
    command -v R >/dev/null 2>&1 && command -v Rscript >/dev/null 2>&1
}

if ! ensure_r_binaries; then
    echo "ERROR: R and/or Rscript not found on PATH after module setup."
    module list 2>&1 || true
    module spider R 2>&1 || true
    exit 1
fi

R_BIN="$(command -v R)"
RSCRIPT_BIN="$(command -v Rscript)"
echo "Resolved R: $R_BIN ($("$R_BIN" --version | head -1))"
echo "Resolved Rscript: $RSCRIPT_BIN"
echo ""

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

echo "Checking/installing sim-study runtime packages if missing..."
"$RSCRIPT_BIN" -e 'pkgs <- c("terra","spatstat","sf","data.table","dplyr","ggplot2","foreach","doParallel","R.utils","reshape2","gridExtra","scales"); miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org", dependencies = TRUE)'

echo "Installing PPDisentangle from source (fresh install)..."
"$R_BIN" CMD INSTALL --preclean --no-multiarch "$PKG_ROOT"
echo ""

"$RSCRIPT_BIN" "$PKG_ROOT/inst/sim_study/sim_study.R" --cluster --sims "$PP_SIMS" $PP_TEST 2>&1

echo ""
echo "=== Done $(date) ==="
