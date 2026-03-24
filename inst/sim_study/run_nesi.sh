#!/bin/bash
#SBATCH --job-name=PPDis_sim
#SBATCH --account=uoo04008
#SBATCH --nodes=1
#SBATCH --ntasks=1

set -euo pipefail

# Minimal interface:
#   --sims N    number of simulations / cores
#   --test      quick test profile

PP_SIMS=32
PP_TEST=""
PP_MODE="${PP_MODE:-}"
PP_POST_TIME_MULTIPLIER="${PP_POST_TIME_MULTIPLIER:-1}"
HAVE_SIMS_ARG=0

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --mode) PP_MODE="$2"; shift 2 ;;
        --sims) PP_SIMS="$2"; HAVE_SIMS_ARG=1; shift 2 ;;
        --post-time-multiplier) PP_POST_TIME_MULTIPLIER="$2"; shift 2 ;;
        --test) PP_TEST="--test"; shift ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

if [ -n "$PP_MODE" ]; then
    mode_norm="$(echo "$PP_MODE" | tr '[:upper:]' '[:lower:]')"
    case "$mode_norm" in
        test)
            PP_TEST="--test"
            if [ "$HAVE_SIMS_ARG" -eq 0 ]; then
                PP_SIMS=2
            fi
            ;;
        quick)
            PP_TEST=""
            if [ "$HAVE_SIMS_ARG" -eq 0 ]; then
                PP_SIMS=32
            fi
            ;;
        long|full|big)
            PP_TEST=""
            if [ "$HAVE_SIMS_ARG" -eq 0 ]; then
                PP_SIMS=32
            fi
            ;;
        *)
            echo "Unknown --mode '$PP_MODE' (expected: test | quick | long)"
            exit 1
            ;;
    esac
fi

if [ -n "$PP_TEST" ] && [ "$HAVE_SIMS_ARG" -eq 0 ]; then
    PP_SIMS=2
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
    OUTPUT_DIR="$PKG_ROOT/output/sim_study"
    LEGACY_OUTPUT_DIR="$PKG_ROOT/inst/sim_study/output"
    mkdir -p "$OUTPUT_DIR" "$LEGACY_OUTPUT_DIR"

    CPUS="$PP_SIMS"
    if [ "$CPUS" -lt 1 ]; then
        echo "ERROR: --sims must be >= 1"
        exit 1
    fi

    # Automatic resources.
    MEM_PER_SIM_GB="${PP_MEM_PER_SIM_GB:-2}"
    MEM_MAX_GB="${PP_MEM_MAX_GB:-200}"
    MEM_MIN_GB="${PP_MEM_MIN_GB:-8}"
    if [ -n "$PP_TEST" ]; then
        SB_TIME="00:20:00"
        SB_MEM="16G"
    else
        MEM_GB=$(( CPUS * MEM_PER_SIM_GB ))
        [ "$MEM_GB" -lt "$MEM_MIN_GB" ] && MEM_GB="$MEM_MIN_GB"
        [ "$MEM_GB" -gt "$MEM_MAX_GB" ] && MEM_GB="$MEM_MAX_GB"
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

    echo "Submitting to NeSI: mode=${PP_MODE:-manual} sims=$PP_SIMS cpus=$CPUS mem=$SB_MEM time=$SB_TIME post_time_multiplier=$PP_POST_TIME_MULTIPLIER ${PP_TEST:+(test)}"

    SBATCH_EXPORT="ALL,PP_SIMS=$PP_SIMS,PP_TEST=$PP_TEST,PP_MODE=$PP_MODE,PP_POST_TIME_MULTIPLIER=$PP_POST_TIME_MULTIPLIER,PKG_ROOT=$PKG_ROOT"
    JOB_ID=$(sbatch --parsable \
        --cpus-per-task="$CPUS" \
        --mem="$SB_MEM" \
        --time="$SB_TIME" \
        $EXTRA_SBATCH \
        --export="$SBATCH_EXPORT" \
        --output="$OUTPUT_DIR/%j_slurm.out" \
        --error="$OUTPUT_DIR/%j_slurm.err" \
        "$SCRIPT_DIR/run_nesi.sh")

    echo "Job $JOB_ID submitted"
    echo "  Results: output/sim_study/$JOB_ID.rds"
    echo "  Log:     output/sim_study/$JOB_ID.log"
    echo "  SLURM:   output/sim_study/${JOB_ID}_slurm.out"
    exit 0
fi

cd "$PKG_ROOT"
mkdir -p "$PKG_ROOT/output/sim_study" "$PKG_ROOT/inst/sim_study/output"

if [ -z "$PP_TEST" ] && [ -n "${SLURM_CPUS_PER_TASK:-}" ] && [ "$PP_SIMS" -ne "$SLURM_CPUS_PER_TASK" ]; then
    echo "Adjusting sims to match allocated CPUs: sims=$PP_SIMS -> ${SLURM_CPUS_PER_TASK}"
    PP_SIMS="$SLURM_CPUS_PER_TASK"
fi

echo "=== PPDisentangle Sim Study (NeSI) ==="
echo "Job $SLURM_JOB_ID | $(date)"
echo "Sims: $PP_SIMS | CPUs: $SLURM_CPUS_PER_TASK"
echo "Mode: ${PP_MODE:-manual}"
echo "Post-time multiplier: ${PP_POST_TIME_MULTIPLIER}"
echo "Node: $(hostname) | Partition: ${SLURM_JOB_PARTITION:-unknown}"
echo ""

export PP_POST_TIME_MULTIPLIER

mode_norm_runtime="$(echo "${PP_MODE:-}" | tr '[:upper:]' '[:lower:]')"
if [ "$mode_norm_runtime" = "quick" ]; then
    if [ -z "${PP_SEM_WORKERS:-}" ]; then
        if [ "$PP_SIMS" -gt 8 ]; then
            export PP_SEM_WORKERS=8
        else
            export PP_SEM_WORKERS="$PP_SIMS"
        fi
    fi
    export PP_SEM_INNER_ITER="${PP_SEM_INNER_ITER:-200}"
    export PP_SEM_OUTER_ITER="${PP_SEM_OUTER_ITER:-3}"
    export PP_SEM_N_PROPS="${PP_SEM_N_PROPS:-20}"
    export PP_SEM_N_LABELLINGS="${PP_SEM_N_LABELLINGS:-10}"
    export PP_ATE_N_SIMS="${PP_ATE_N_SIMS:-1}"
    export PP_ATE_N_TAU_SIMS="${PP_ATE_N_TAU_SIMS:-1}"
    export PP_ATE_N_TAU_I="${PP_ATE_N_TAU_I:-1}"
    echo "Quick profile env:"
    echo "  PP_SEM_INNER_ITER=$PP_SEM_INNER_ITER"
    echo "  PP_SEM_OUTER_ITER=$PP_SEM_OUTER_ITER"
    echo "  PP_SEM_N_PROPS=$PP_SEM_N_PROPS"
    echo "  PP_SEM_N_LABELLINGS=$PP_SEM_N_LABELLINGS"
    echo "  PP_SEM_WORKERS=$PP_SEM_WORKERS"
    echo "  PP_ATE_N_SIMS=$PP_ATE_N_SIMS PP_ATE_N_TAU_SIMS=$PP_ATE_N_TAU_SIMS PP_ATE_N_TAU_I=$PP_ATE_N_TAU_I"
    echo ""
fi

# Shared library path only; guard package install lock collisions.
SHARED_R_LIBS_USER="${R_LIBS_USER:-/nesi/project/uoo04008/Rlibs}"
mkdir -p "$SHARED_R_LIBS_USER"
export R_LIBS_USER="$SHARED_R_LIBS_USER"
# Also expose the same path via R_LIBS so non-interactive calls consistently see it.
export R_LIBS="${R_LIBS_USER}${R_LIBS:+:${R_LIBS}}"
PP_LOCK_DIR="${SHARED_R_LIBS_USER}/00LOCK-PPDisentangle"
echo "R_LIBS_USER=$R_LIBS_USER"

wait_for_pp_lock_clear() {
    local lock_dir="$1"
    local waited_s=0
    local sleep_s=5
    while [ -d "$lock_dir" ]; do
        echo "Waiting for lock release: $lock_dir (waited ${waited_s}s)..."
        sleep "$sleep_s"
        waited_s=$(( waited_s + sleep_s ))
    done
}

cleanup_pp_lock_if_safe() {
    local lock_dir="$1"
    if [ ! -d "$lock_dir" ]; then
        return 0
    fi
    if pgrep -u "${USER:-$(id -un)}" -f "R CMD INSTALL.*PPDisentangle" >/dev/null 2>&1; then
        echo "Lock present but PPDisentangle install still active; leaving lock in place."
        return 0
    fi
    echo "Removing stale lock: $lock_dir"
    rm -rf "$lock_dir" 2>/dev/null || true
}

trap 'cleanup_pp_lock_if_safe "$PP_LOCK_DIR"' EXIT

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
    echo "R binaries not found after initial R-Geo load; retrying module chain fallbacks..."
    local target_rgeo="$TARGET_R_GEO"
    local target_r="R/4.3.2-foss-2023a"
    if [[ "$target_rgeo" =~ ^R-Geo/(.+)$ ]]; then
        target_r="R/${BASH_REMATCH[1]}"
    fi

    # Prefer restoring a clean R-Geo stack first so sf/terra/system deps stay consistent.
    if module --force purge >/dev/null 2>&1 \
        && module load NeSI/zen3 >/dev/null 2>&1 \
        && module load "$target_rgeo" >/dev/null 2>&1; then
        :
    elif module --force purge >/dev/null 2>&1 \
        && module load foss/2023a >/dev/null 2>&1 \
        && module load "$target_rgeo" >/dev/null 2>&1; then
        :
    elif module --force purge >/dev/null 2>&1 \
        && module load NeSI/zen3 >/dev/null 2>&1 \
        && module load foss/2023a >/dev/null 2>&1 \
        && module load "$target_rgeo" >/dev/null 2>&1; then
        :
    # If R-Geo still doesn't expose binaries, fall back to explicit R module.
    elif module load "$target_r" >/dev/null 2>&1; then
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
echo "gdal-config: $(command -v gdal-config || echo 'not found')"
echo "pkg-config: $(command -v pkg-config || echo 'not found')"
echo ""

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

echo "Verifying geospatial toolchain from module stack..."
"$RSCRIPT_BIN" -e 'req <- c("sf","terra","units"); miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]; gdal <- Sys.which("gdal-config"); pkgc <- Sys.which("pkg-config"); cat(sprintf("gdal-config in R: %s\n", ifelse(nzchar(gdal), gdal, "<missing>"))); cat(sprintf("pkg-config in R: %s\n", ifelse(nzchar(pkgc), pkgc, "<missing>"))); if (length(miss) < 1L) { cat("Spatial stack OK (sf/terra/units available).\n"); quit(status = 0L) }; if (!nzchar(gdal) || !nzchar(pkgc)) stop(sprintf("Missing required spatial packages (%s) and missing gdal-config/pkg-config on PATH after module load.", paste(miss, collapse=", "))); stop(sprintf("Missing required spatial packages from R-Geo stack: %s", paste(miss, collapse=", ")))'

echo "Checking/installing non-spatial sim-study runtime packages if missing..."
"$RSCRIPT_BIN" -e 'pkgs <- c("spatstat","data.table","dplyr","ggplot2","foreach","doParallel","R.utils","reshape2","gridExtra","scales"); miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org", dependencies = NA)'

echo "Installing PPDisentangle from source (fresh install)..."
wait_for_pp_lock_clear "$PP_LOCK_DIR"
cleanup_pp_lock_if_safe "$PP_LOCK_DIR"
"$R_BIN" CMD INSTALL --preclean --no-multiarch "$PKG_ROOT"
cleanup_pp_lock_if_safe "$PP_LOCK_DIR"
echo "Verifying PPDisentangle is visible in runtime library paths..."
"$RSCRIPT_BIN" -e 'user_lib <- Sys.getenv("R_LIBS_USER", ""); if (nzchar(user_lib)) { libs <- strsplit(user_lib, .Platform$path.sep, fixed = TRUE)[[1]]; libs <- libs[nzchar(libs)]; if (length(libs) > 0L) .libPaths(c(libs, .libPaths())) }; cat(".libPaths()=", paste(.libPaths(), collapse=" | "), "\n", sep=""); if (!requireNamespace("PPDisentangle", quietly = TRUE)) stop("PPDisentangle not visible after install."); library(PPDisentangle); cat("PPDisentangle load check OK.\n")'
echo ""

"$RSCRIPT_BIN" "$PKG_ROOT/inst/sim_study/sim_study.R" --cluster --sims "$PP_SIMS" $PP_TEST 2>&1

echo ""
echo "=== Done $(date) ==="
