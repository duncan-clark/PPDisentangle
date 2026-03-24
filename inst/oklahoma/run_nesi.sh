#!/bin/bash
#SBATCH --job-name=PPDis_oklahoma
#SBATCH --account=uoo04008
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --mem=48G

set -euo pipefail

# ----------------------------
# Args
# ----------------------------
PP_CORES="${PP_CORES:-32}"
PP_BOOT_REPS="${PP_BOOT_REPS:-}"
PP_SEM_INNER="${PP_SEM_INNER:-100}"
PP_SENS_SEM_INNER="${PP_SENS_SEM_INNER:-}"
PP_BOOT_SEM_INNER="${PP_BOOT_SEM_INNER:-}"
PP_BOOT_TARGETS="${PP_BOOT_TARGETS:-E}"
PP_KDE_VARIANT_MODE="${PP_KDE_VARIANT_MODE:-}"
PP_BOOT_OUTER_CORES="${PP_BOOT_OUTER_CORES:-}"
PP_RUN_SENSITIVITY="${PP_RUN_SENSITIVITY:-auto}"
PP_MEM="${PP_MEM:-}"
PP_TIME="${PP_TIME:-72:00:00}"
PP_SETUP_TEST="${PP_SETUP_TEST:-0}"
PP_MODE="${PP_MODE:-}"
PP_SEED="${PP_SEED:-1}"
CORES_EXPLICIT=0
MEM_EXPLICIT=0
SEM_INNER_EXPLICIT=0
SENS_SEM_INNER_EXPLICIT=0
BOOT_SEM_INNER_EXPLICIT=0
BOOT_TARGETS_EXPLICIT=0
KDE_VARIANT_MODE_EXPLICIT=0
SETUP_TEST_EXPLICIT=0
BOOT_REPS_EXPLICIT=0
BOOT_OUTER_CORES_EXPLICIT=0
RUN_SENS_EXPLICIT=0
if [ -n "$PP_BOOT_REPS" ]; then BOOT_REPS_EXPLICIT=1; fi
if [ -n "$PP_BOOT_OUTER_CORES" ]; then BOOT_OUTER_CORES_EXPLICIT=1; fi
if [ "$PP_RUN_SENSITIVITY" != "auto" ]; then RUN_SENS_EXPLICIT=1; fi

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --mode) PP_MODE="$2"; shift 2 ;;
    --cores) PP_CORES="$2"; CORES_EXPLICIT=1; shift 2 ;;
    --sims) PP_CORES="$2"; CORES_EXPLICIT=1; shift 2 ;;  # alias: match sim_study launcher UX
    --boot-reps) PP_BOOT_REPS="$2"; BOOT_REPS_EXPLICIT=1; shift 2 ;;
    --sem-inner) PP_SEM_INNER="$2"; SEM_INNER_EXPLICIT=1; shift 2 ;;
    --sens-sem-inner) PP_SENS_SEM_INNER="$2"; SENS_SEM_INNER_EXPLICIT=1; shift 2 ;;
    --boot-sem-inner) PP_BOOT_SEM_INNER="$2"; BOOT_SEM_INNER_EXPLICIT=1; shift 2 ;;
    --boot-targets) PP_BOOT_TARGETS="$2"; BOOT_TARGETS_EXPLICIT=1; shift 2 ;;
    --kde-variant-mode) PP_KDE_VARIANT_MODE="$2"; KDE_VARIANT_MODE_EXPLICIT=1; shift 2 ;;
    --boot-outer-cores) PP_BOOT_OUTER_CORES="$2"; BOOT_OUTER_CORES_EXPLICIT=1; shift 2 ;;
    --run-sensitivity) PP_RUN_SENSITIVITY="$2"; RUN_SENS_EXPLICIT=1; shift 2 ;;
    --setup-test) PP_SETUP_TEST=1; SETUP_TEST_EXPLICIT=1; shift ;;
    --mem) PP_MEM="$2"; MEM_EXPLICIT=1; shift 2 ;;
    --time) PP_TIME="$2"; shift 2 ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

if [ -n "$PP_MODE" ]; then
  mode_norm="$(echo "$PP_MODE" | tr '[:upper:]' '[:lower:]')"
  case "$mode_norm" in
    very-quick|veryquick|smoke)
      if [ "$SETUP_TEST_EXPLICIT" -ne 1 ]; then PP_SETUP_TEST=0; fi
      if [ "$CORES_EXPLICIT" -ne 1 ]; then PP_CORES=16; fi
      if [ "$BOOT_REPS_EXPLICIT" -ne 1 ]; then PP_BOOT_REPS=1; fi
      if [ "$SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SEM_INNER=2; fi
      if [ "$SENS_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SENS_SEM_INNER=2; fi
      if [ "$BOOT_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_BOOT_SEM_INNER=2; fi
      if [ "$BOOT_TARGETS_EXPLICIT" -ne 1 ]; then PP_BOOT_TARGETS="E,F"; fi
      if [ "$KDE_VARIANT_MODE_EXPLICIT" -ne 1 ]; then PP_KDE_VARIANT_MODE="triple"; fi
      if [ "$BOOT_OUTER_CORES_EXPLICIT" -ne 1 ]; then PP_BOOT_OUTER_CORES=1; fi
      if [ "$RUN_SENS_EXPLICIT" -ne 1 ]; then PP_RUN_SENSITIVITY=1; fi
      if [ "$MEM_EXPLICIT" -ne 1 ]; then PP_MEM=96G; fi
      ;;
    quick)
      if [ "$SETUP_TEST_EXPLICIT" -ne 1 ]; then PP_SETUP_TEST=0; fi
      if [ "$CORES_EXPLICIT" -ne 1 ]; then PP_CORES=32; fi
      if [ "$BOOT_REPS_EXPLICIT" -ne 1 ]; then PP_BOOT_REPS=0; fi
      if [ "$SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SEM_INNER=200; fi
      if [ "$SENS_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SENS_SEM_INNER=200; fi
      if [ "$BOOT_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_BOOT_SEM_INNER=50; fi
      if [ "$BOOT_TARGETS_EXPLICIT" -ne 1 ]; then PP_BOOT_TARGETS="E"; fi
      if [ "$KDE_VARIANT_MODE_EXPLICIT" -ne 1 ]; then PP_KDE_VARIANT_MODE="triple"; fi
      if [ "$BOOT_OUTER_CORES_EXPLICIT" -ne 1 ]; then PP_BOOT_OUTER_CORES=1; fi
      if [ "$RUN_SENS_EXPLICIT" -ne 1 ]; then PP_RUN_SENSITIVITY=0; fi
      if [ "$MEM_EXPLICIT" -ne 1 ]; then PP_MEM=96G; fi
      ;;
    test|setup-test)
      if [ "$SETUP_TEST_EXPLICIT" -ne 1 ]; then PP_SETUP_TEST=1; fi
      if [ "$CORES_EXPLICIT" -ne 1 ]; then PP_CORES=32; fi
      if [ "$BOOT_REPS_EXPLICIT" -ne 1 ]; then PP_BOOT_REPS=2; fi
      if [ "$SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SEM_INNER=100; fi
      if [ "$SENS_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SENS_SEM_INNER=2; fi
      if [ "$BOOT_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_BOOT_SEM_INNER=2; fi
      if [ "$BOOT_OUTER_CORES_EXPLICIT" -ne 1 ]; then PP_BOOT_OUTER_CORES=1; fi
      if [ "$RUN_SENS_EXPLICIT" -ne 1 ]; then PP_RUN_SENSITIVITY=0; fi
      if [ "$MEM_EXPLICIT" -ne 1 ]; then PP_MEM=64G; fi
      ;;
    default)
      if [ "$SETUP_TEST_EXPLICIT" -ne 1 ]; then PP_SETUP_TEST=0; fi
      if [ "$CORES_EXPLICIT" -ne 1 ]; then PP_CORES=32; fi
      if [ "$BOOT_REPS_EXPLICIT" -ne 1 ]; then PP_BOOT_REPS=12; fi
      if [ "$SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SEM_INNER=1000; fi
      if [ "$SENS_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SENS_SEM_INNER=1000; fi
      if [ "$BOOT_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_BOOT_SEM_INNER=200; fi
      if [ "$BOOT_TARGETS_EXPLICIT" -ne 1 ]; then PP_BOOT_TARGETS="E"; fi
      if [ "$KDE_VARIANT_MODE_EXPLICIT" -ne 1 ]; then PP_KDE_VARIANT_MODE="triple"; fi
      if [ "$BOOT_OUTER_CORES_EXPLICIT" -ne 1 ]; then
        AUTO_BOOT_CORES=$(( PP_CORES / 4 ))
        [ "$AUTO_BOOT_CORES" -lt 2 ] && AUTO_BOOT_CORES=2
        [ "$AUTO_BOOT_CORES" -gt 8 ] && AUTO_BOOT_CORES=8
        PP_BOOT_OUTER_CORES="$AUTO_BOOT_CORES"
      fi
      if [ "$RUN_SENS_EXPLICIT" -ne 1 ]; then PP_RUN_SENSITIVITY=0; fi
      if [ "$MEM_EXPLICIT" -ne 1 ]; then PP_MEM=200G; fi
      ;;
    long|full|big)
      if [ "$SETUP_TEST_EXPLICIT" -ne 1 ]; then PP_SETUP_TEST=0; fi
      if [ "$CORES_EXPLICIT" -ne 1 ]; then PP_CORES=32; fi
      if [ "$BOOT_REPS_EXPLICIT" -ne 1 ]; then PP_BOOT_REPS=12; fi
      if [ "$SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SEM_INNER=1000; fi
      if [ "$SENS_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SENS_SEM_INNER=1000; fi
      if [ "$BOOT_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_BOOT_SEM_INNER=200; fi
      if [ "$BOOT_TARGETS_EXPLICIT" -ne 1 ]; then PP_BOOT_TARGETS="E"; fi
      if [ "$KDE_VARIANT_MODE_EXPLICIT" -ne 1 ]; then PP_KDE_VARIANT_MODE="triple"; fi
      if [ "$BOOT_OUTER_CORES_EXPLICIT" -ne 1 ]; then
        AUTO_BOOT_CORES=$(( PP_CORES / 4 ))
        [ "$AUTO_BOOT_CORES" -lt 2 ] && AUTO_BOOT_CORES=2
        [ "$AUTO_BOOT_CORES" -gt 8 ] && AUTO_BOOT_CORES=8
        PP_BOOT_OUTER_CORES="$AUTO_BOOT_CORES"
      fi
      if [ "$RUN_SENS_EXPLICIT" -ne 1 ]; then PP_RUN_SENSITIVITY=0; fi
      if [ "$MEM_EXPLICIT" -ne 1 ]; then PP_MEM=200G; fi
      ;;
    *)
      echo "Unknown --mode '$PP_MODE' (expected: very-quick | quick | default | test | long)"
      exit 1
      ;;
  esac
fi

MEM_PER_CORE_GB="${PP_MEM_PER_CORE_GB:-2}"
MEM_MAX_GB="${PP_MEM_MAX_GB:-200}"
MEM_MIN_GB="${PP_MEM_MIN_GB:-8}"
if [ -z "$PP_MEM" ]; then
  MEM_GB=$(( PP_CORES * MEM_PER_CORE_GB ))
  [ "$MEM_GB" -lt "$MEM_MIN_GB" ] && MEM_GB="$MEM_MIN_GB"
  [ "$MEM_GB" -gt "$MEM_MAX_GB" ] && MEM_GB="$MEM_MAX_GB"
  PP_MEM="${MEM_GB}G"
fi

if [ -z "$PP_BOOT_REPS" ]; then
  # Bootstrap is the dominant memory consumer; default below total cores.
  if [ "$PP_CORES" -gt 8 ]; then
    PP_BOOT_REPS=8
  else
    PP_BOOT_REPS="$PP_CORES"
  fi
fi
if [ -z "$PP_SENS_SEM_INNER" ]; then
  PP_SENS_SEM_INNER="$PP_SEM_INNER"
fi
if [ -z "$PP_BOOT_SEM_INNER" ]; then
  PP_BOOT_SEM_INNER="$PP_SEM_INNER"
fi
if [ -z "$PP_BOOT_OUTER_CORES" ]; then
  AUTO_BOOT_CORES=$(( PP_CORES / 4 ))
  [ "$AUTO_BOOT_CORES" -lt 2 ] && AUTO_BOOT_CORES=2
  [ "$AUTO_BOOT_CORES" -gt 8 ] && AUTO_BOOT_CORES=8
  PP_BOOT_OUTER_CORES="$AUTO_BOOT_CORES"
fi
if [ -z "$PP_KDE_VARIANT_MODE" ]; then
  PP_KDE_VARIANT_MODE="triple"
fi
case "$(echo "$PP_KDE_VARIANT_MODE" | tr '[:upper:]' '[:lower:]')" in
  single|triple) ;;
  *)
    echo "Invalid --kde-variant-mode '$PP_KDE_VARIANT_MODE' (expected: single | triple)"
    exit 1
    ;;
esac

if [ "$PP_RUN_SENSITIVITY" = "auto" ]; then
  # Prefer memory headroom for bootstrap unless user explicitly enables sensitivity.
  if [ "${PP_BOOT_REPS:-0}" -gt 0 ]; then
    PP_RUN_SENSITIVITY=0
  else
    PP_RUN_SENSITIVITY=1
  fi
fi

# ----------------------------
# Paths
# ----------------------------
if [ -n "${SLURM_JOB_ID:-}" ] && [ -n "${PKG_ROOT:-}" ] && [ -d "$PKG_ROOT" ]; then
  :
elif [ -n "${SLURM_JOB_ID:-}" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -d "$SLURM_SUBMIT_DIR" ]; then
  PKG_ROOT="$SLURM_SUBMIT_DIR"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
fi

# ----------------------------
# Submit mode
# ----------------------------
if [ -z "${SLURM_JOB_ID:-}" ]; then
  cd "$PKG_ROOT"
  git pull origin main 2>/dev/null || true
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  OUTPUT_DIR="$PKG_ROOT/inst/oklahoma/output"
  mkdir -p "$OUTPUT_DIR"

  EXTRA_SBATCH=""
  if [ "$PP_CORES" -gt 72 ]; then
    if [ "$PP_CORES" -gt 256 ]; then
      echo "ERROR: max 256 CPUs per node on milan."
      exit 1
    fi
    EXTRA_SBATCH="--partition=milan"
    echo "Note: using milan partition for >72 cores."
  fi

  SBATCH_EXPORT="ALL,PKG_ROOT=$PKG_ROOT,PP_CORES=$PP_CORES,PP_BOOT_REPS=$PP_BOOT_REPS,PP_SEM_INNER=$PP_SEM_INNER,PP_SENS_SEM_INNER=$PP_SENS_SEM_INNER,PP_BOOT_SEM_INNER=$PP_BOOT_SEM_INNER,PP_BOOT_TARGETS=$PP_BOOT_TARGETS,PP_KDE_VARIANT_MODE=$PP_KDE_VARIANT_MODE,PP_BOOT_OUTER_CORES=$PP_BOOT_OUTER_CORES,PP_RUN_SENSITIVITY=$PP_RUN_SENSITIVITY,PP_MEM=$PP_MEM,PP_TIME=$PP_TIME"
  SBATCH_EXPORT="${SBATCH_EXPORT},PP_SETUP_TEST=$PP_SETUP_TEST"
  [ -n "${PP_R_GEO_MODULE:-}" ] && SBATCH_EXPORT="${SBATCH_EXPORT},PP_R_GEO_MODULE=$PP_R_GEO_MODULE"

  echo "Submitting Oklahoma job: mode=${PP_MODE:-manual} cores=$PP_CORES boot_reps=$PP_BOOT_REPS sem_inner=$PP_SEM_INNER sens_inner=$PP_SENS_SEM_INNER boot_inner=$PP_BOOT_SEM_INNER kde_variants=$PP_KDE_VARIANT_MODE boot_outer_cores=$PP_BOOT_OUTER_CORES setup_test=$PP_SETUP_TEST"

  JOB_ID=$(sbatch --parsable \
    --cpus-per-task="$PP_CORES" \
    --mem="$PP_MEM" \
    --time="$PP_TIME" \
    $EXTRA_SBATCH \
    --export="$SBATCH_EXPORT" \
    --output="$OUTPUT_DIR/%j_oklahoma_slurm.out" \
    --error="$OUTPUT_DIR/%j_oklahoma_slurm.err" \
    "$SCRIPT_DIR/run_nesi.sh")

  echo "Job $JOB_ID submitted"
  echo "SLURM out: inst/oklahoma/output/${JOB_ID}_oklahoma_slurm.out"
  echo "SLURM err: inst/oklahoma/output/${JOB_ID}_oklahoma_slurm.err"
  exit 0
fi

# ----------------------------
# Job mode
# ----------------------------
cd "$PKG_ROOT"
mkdir -p "$PKG_ROOT/inst/oklahoma/output"

echo "=== PPDisentangle Oklahoma (NeSI) ==="
echo "Job: ${SLURM_JOB_ID} | $(date)"
echo "Node: $(hostname) | Partition: ${SLURM_JOB_PARTITION:-unknown}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-$PP_CORES}"
echo "boot_reps=$PP_BOOT_REPS sem_inner=$PP_SEM_INNER sens_inner=$PP_SENS_SEM_INNER boot_inner=$PP_BOOT_SEM_INNER targets=$PP_BOOT_TARGETS kde_variants=$PP_KDE_VARIANT_MODE"
echo "setup_test=$PP_SETUP_TEST mode=${PP_MODE:-manual}"
echo "seed=$PP_SEED (fit jobs RNG de-correlated by model; bootstrap RNG de-correlated by replicate)"
echo ""

# Shared library path only; guard package install lock collisions.
SHARED_R_LIBS_USER="${R_LIBS_USER:-/nesi/project/uoo04008/Rlibs}"
mkdir -p "$SHARED_R_LIBS_USER"
export R_LIBS_USER="$SHARED_R_LIBS_USER"
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
      mapfile -t R_GEO_CANDIDATES < <(module -t avail R-Geo 2>&1 | awk '/^R-Geo\//{print $1}' | sort -Vr | uniq)
      LOADED=0
      for cand in "${R_GEO_CANDIDATES[@]}"; do
        if try_load_rgeo "$cand"; then LOADED=1; break; fi
        module --force purge
        if try_load_rgeo_with_toolchain "$cand"; then LOADED=1; break; fi
        module --force purge
        if module load NeSI/zen3 >/dev/null 2>&1 && try_load_rgeo "$cand"; then LOADED=1; break; fi
        module --force purge
      done
      if [ "$LOADED" -ne 1 ]; then
        echo "ERROR: Failed to load any R-Geo module."
        module spider R-Geo || true
        exit 1
      fi
    fi
  fi
fi

ensure_r_binaries() {
  if command -v R >/dev/null 2>&1 && command -v Rscript >/dev/null 2>&1; then
    return 0
  fi

  echo "R binaries not found after initial R-Geo load; retrying with module chain fallbacks..."
  local target_r=""
  local target_rgeo="$TARGET_R_GEO"
  if [[ "$TARGET_R_GEO" =~ ^R-Geo/(.+)$ ]]; then
    target_r="R/${BASH_REMATCH[1]}"
  else
    target_r="R/4.3.2-foss-2023a"
  fi

  # First, prefer a clean R-Geo stack so geospatial deps (e.g., terra) remain available.
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
  # If that still doesn't expose R binaries, fall back to explicit R module.
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
  echo "Diagnostics:"
  module list 2>&1 || true
  module spider R 2>&1 || true
  exit 1
fi

R_BIN="$(command -v R)"
RSCRIPT_BIN="$(command -v Rscript)"
echo "R: $R_BIN ($("$R_BIN" --version | head -1))"
echo "Rscript: $RSCRIPT_BIN"
echo ""

# Avoid nested threading / OpenMP crashes.
export OMP_NUM_THREADS=1
export OMP_THREAD_LIMIT=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export KMP_INIT_AT_FORK=FALSE

OK_DEPS_STAMP="${SHARED_R_LIBS_USER}/.ppdis_oklahoma_runtime_deps_ok"
if [ "${PP_REFRESH_DEPS:-0}" = "1" ] || [ ! -f "$OK_DEPS_STAMP" ]; then
  echo "Checking/installing Oklahoma runtime packages..."
  "$RSCRIPT_BIN" -e 'pkgs <- c("terra","spatstat","sf","tigris","data.table","dplyr","ggplot2","pkgload","quarto"); miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org", dependencies = TRUE)'
  touch "$OK_DEPS_STAMP"
else
  echo "Skipping runtime dependency bootstrap (set PP_REFRESH_DEPS=1 to recheck)."
fi

echo "Installing PPDisentangle from source (fresh install)..."
wait_for_pp_lock_clear "$PP_LOCK_DIR"
cleanup_pp_lock_if_safe "$PP_LOCK_DIR"
"$R_BIN" CMD INSTALL --preclean --no-multiarch "$PKG_ROOT"
cleanup_pp_lock_if_safe "$PP_LOCK_DIR"
echo ""

# Oklahoma run config.
JOB_CORES="${SLURM_CPUS_PER_TASK:-$PP_CORES}"
SAFE_SHARED_CORES=$(( JOB_CORES < 2 ? JOB_CORES : 2 ))
export OK_MEMORY_SAFE=true
export OK_PARALLEL_BACKEND=psock
export OK_CORES="${JOB_CORES}"
export OK_SENS_CORES="${SAFE_SHARED_CORES}"
export OK_ATE_SIM_CORES="${SAFE_SHARED_CORES}"
export OK_RUN_DECODE=false
export OK_VERBOSE=false
export OK_SEM_INNER_ITER="$PP_SEM_INNER"
export OK_SENS_SEM_INNER_ITER="$PP_SENS_SEM_INNER"
if [ "$PP_RUN_SENSITIVITY" = "1" ] || [ "$PP_RUN_SENSITIVITY" = "true" ] || [ "$PP_RUN_SENSITIVITY" = "yes" ]; then
  export OK_RUN_SENSITIVITY=true
else
  export OK_RUN_SENSITIVITY=false
fi
export OK_RUN_BOOTSTRAP_ATE=true
export OK_BOOT_N_REPS="$PP_BOOT_REPS"
export OK_BOOT_TARGETS="$PP_BOOT_TARGETS"
export OK_KDE_VARIANT_MODE="$PP_KDE_VARIANT_MODE"
export OK_BOOT_SEM_INNER_ITER="$PP_BOOT_SEM_INNER"
export OK_BOOT_OUTER_CORES="$PP_BOOT_OUTER_CORES"
export OK_GLOBAL_SEED="$PP_SEED"
export OK_IDENTICAL_RANDOMNESS=false
export OK_BOOT_IDENTICAL_RANDOMNESS=false
export OK_BOOT_GUARD_DEGENERATE=true
export OK_REPORT_FORMATS=html

if [ "${PP_BOOT_REPS:-0}" -le 0 ]; then
  export OK_RUN_BOOTSTRAP_ATE=false
fi

if [ -n "${PP_MODE:-}" ]; then
  mode_norm_runtime="$(echo "$PP_MODE" | tr '[:upper:]' '[:lower:]')"
  if [ "$mode_norm_runtime" = "very-quick" ] || [ "$mode_norm_runtime" = "veryquick" ] || [ "$mode_norm_runtime" = "smoke" ]; then
    echo "Applying very-quick profile: decode+sens+bootstrap enabled with 2 inner iterations."
    export OK_RUN_DECODE=true
    export OK_DECODE_ITER=2
    export OK_RUN_SENSITIVITY=true
    export OK_RUN_BOOTSTRAP_ATE=true
    export OK_SEM_INNER_ITER=2
    export OK_SENS_SEM_INNER_ITER=2
    export OK_BOOT_SEM_INNER_ITER=2
    export OK_BOOT_OUTER_CORES=1
    export OK_BOOT_N_REPS="${PP_BOOT_REPS:-1}"
    export OK_BOOT_TARGETS="${PP_BOOT_TARGETS:-E,F}"
  fi
fi

if [ "$PP_SETUP_TEST" = "1" ]; then
  echo "Applying setup-test profile: main SEM inner=100, decode inner=2, sensitivity inner=2, bootstrap inner=2, sequential bootstrap."
  export OK_SEM_INNER_ITER=100
  export OK_DECODE_ITER=2
  export OK_SENS_SEM_INNER_ITER=2
  export OK_BOOT_SEM_INNER_ITER=2
  export OK_SENS_CORES=1
  export OK_ATE_SIM_CORES=1
  export OK_BOOT_OUTER_CORES=1
  if [ "$BOOT_REPS_EXPLICIT" -ne 1 ]; then
    export OK_BOOT_N_REPS=2
  fi
  if [ "$RUN_SENS_EXPLICIT" -ne 1 ]; then
    export OK_RUN_SENSITIVITY=false
  fi
  export OK_BOOT_TARGETS="E,F"
  export OK_RUN_DECODE=true
  export OK_RUN_BOOTSTRAP_ATE=true
fi

"$RSCRIPT_BIN" "$PKG_ROOT/inst/oklahoma/oklahoma_analysis.R" 2>&1

echo ""
echo "=== Done $(date) ==="
