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
PP_BOOT_TARGETS="${PP_BOOT_TARGETS:-E,F}"
PP_MEM="${PP_MEM:-}"
PP_TIME="${PP_TIME:-72:00:00}"
PP_SETUP_TEST="${PP_SETUP_TEST:-0}"

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --cores) PP_CORES="$2"; shift 2 ;;
    --sims) PP_CORES="$2"; shift 2 ;;  # alias: match sim_study launcher UX
    --boot-reps) PP_BOOT_REPS="$2"; shift 2 ;;
    --sem-inner) PP_SEM_INNER="$2"; shift 2 ;;
    --sens-sem-inner) PP_SENS_SEM_INNER="$2"; shift 2 ;;
    --boot-sem-inner) PP_BOOT_SEM_INNER="$2"; shift 2 ;;
    --boot-targets) PP_BOOT_TARGETS="$2"; shift 2 ;;
    --setup-test) PP_SETUP_TEST=1; shift ;;
    --mem) PP_MEM="$2"; shift 2 ;;
    --time) PP_TIME="$2"; shift 2 ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

MEM_PER_CORE_GB="${PP_MEM_PER_CORE_GB:-1}"
MEM_MAX_GB="${PP_MEM_MAX_GB:-48}"
MEM_MIN_GB="${PP_MEM_MIN_GB:-8}"
if [ -z "$PP_MEM" ]; then
  MEM_GB=$(( PP_CORES * MEM_PER_CORE_GB ))
  [ "$MEM_GB" -lt "$MEM_MIN_GB" ] && MEM_GB="$MEM_MIN_GB"
  [ "$MEM_GB" -gt "$MEM_MAX_GB" ] && MEM_GB="$MEM_MAX_GB"
  PP_MEM="${MEM_GB}G"
fi

if [ -z "$PP_BOOT_REPS" ]; then
  PP_BOOT_REPS="$PP_CORES"
fi
if [ -z "$PP_SENS_SEM_INNER" ]; then
  PP_SENS_SEM_INNER="$PP_SEM_INNER"
fi
if [ -z "$PP_BOOT_SEM_INNER" ]; then
  PP_BOOT_SEM_INNER="$PP_SEM_INNER"
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
  mkdir -p "$PKG_ROOT/cluster_output"

  EXTRA_SBATCH=""
  if [ "$PP_CORES" -gt 72 ]; then
    if [ "$PP_CORES" -gt 256 ]; then
      echo "ERROR: max 256 CPUs per node on milan."
      exit 1
    fi
    EXTRA_SBATCH="--partition=milan"
    echo "Note: using milan partition for >72 cores."
  fi

  SBATCH_EXPORT="ALL,PKG_ROOT=$PKG_ROOT,PP_CORES=$PP_CORES,PP_BOOT_REPS=$PP_BOOT_REPS,PP_SEM_INNER=$PP_SEM_INNER,PP_SENS_SEM_INNER=$PP_SENS_SEM_INNER,PP_BOOT_SEM_INNER=$PP_BOOT_SEM_INNER,PP_BOOT_TARGETS=$PP_BOOT_TARGETS,PP_MEM=$PP_MEM,PP_TIME=$PP_TIME"
  SBATCH_EXPORT="${SBATCH_EXPORT},PP_SETUP_TEST=$PP_SETUP_TEST"
  [ -n "${PP_R_GEO_MODULE:-}" ] && SBATCH_EXPORT="${SBATCH_EXPORT},PP_R_GEO_MODULE=$PP_R_GEO_MODULE"

  echo "Submitting Oklahoma job: cores=$PP_CORES boot_reps=$PP_BOOT_REPS sem_inner=$PP_SEM_INNER sens_inner=$PP_SENS_SEM_INNER boot_inner=$PP_BOOT_SEM_INNER setup_test=$PP_SETUP_TEST"

  JOB_ID=$(sbatch --parsable \
    --cpus-per-task="$PP_CORES" \
    --mem="$PP_MEM" \
    --time="$PP_TIME" \
    $EXTRA_SBATCH \
    --export="$SBATCH_EXPORT" \
    --output="$PKG_ROOT/cluster_output/%j_oklahoma_slurm.out" \
    --error="$PKG_ROOT/cluster_output/%j_oklahoma_slurm.err" \
    "$SCRIPT_DIR/run_nesi.sh")

  echo "Job $JOB_ID submitted"
  echo "SLURM out: cluster_output/${JOB_ID}_oklahoma_slurm.out"
  echo "SLURM err: cluster_output/${JOB_ID}_oklahoma_slurm.err"
  exit 0
fi

# ----------------------------
# Job mode
# ----------------------------
cd "$PKG_ROOT"
mkdir -p "$PKG_ROOT/cluster_output"

echo "=== PPDisentangle Oklahoma (NeSI) ==="
echo "Job: ${SLURM_JOB_ID} | $(date)"
echo "Node: $(hostname) | Partition: ${SLURM_JOB_PARTITION:-unknown}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-$PP_CORES}"
echo "boot_reps=$PP_BOOT_REPS sem_inner=$PP_SEM_INNER sens_inner=$PP_SENS_SEM_INNER boot_inner=$PP_BOOT_SEM_INNER targets=$PP_BOOT_TARGETS"
echo "setup_test=$PP_SETUP_TEST"
echo ""

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

  echo "R binaries not found after R-Geo load; trying explicit R module fallbacks..."
  local target_r=""
  if [[ "$TARGET_R_GEO" =~ ^R-Geo/(.+)$ ]]; then
    target_r="R/${BASH_REMATCH[1]}"
  else
    target_r="R/4.3.2-foss-2023a"
  fi

  # Keep this order identical across launchers.
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

if ! "$RSCRIPT_BIN" -e 'library(PPDisentangle)' 2>/dev/null; then
  echo "Installing PPDisentangle from source..."
  "$R_BIN" CMD INSTALL --no-multiarch "$PKG_ROOT" 2>&1 | tail -5
  echo ""
fi

# Oklahoma run config.
export OK_MEMORY_SAFE=false
export OK_PARALLEL_BACKEND=psock
export OK_CORES="${SLURM_CPUS_PER_TASK:-$PP_CORES}"
export OK_SENS_CORES="${SLURM_CPUS_PER_TASK:-$PP_CORES}"
export OK_ATE_SIM_CORES="${SLURM_CPUS_PER_TASK:-$PP_CORES}"
export OK_RUN_DECODE=false
export OK_SEM_INNER_ITER="$PP_SEM_INNER"
export OK_SENS_SEM_INNER_ITER="$PP_SENS_SEM_INNER"
export OK_RUN_SENSITIVITY=true
export OK_RUN_BOOTSTRAP_ATE=true
export OK_BOOT_N_REPS="$PP_BOOT_REPS"
export OK_BOOT_TARGETS="$PP_BOOT_TARGETS"
export OK_BOOT_SEM_INNER_ITER="$PP_BOOT_SEM_INNER"
export OK_BOOT_OUTER_CORES="${SLURM_CPUS_PER_TASK:-$PP_CORES}"
export OK_REPORT_FORMATS=html

if [ "$PP_SETUP_TEST" = "1" ]; then
  echo "Applying setup-test profile: main SEM inner=100, decode inner=2, sensitivity inner=10, bootstrap inner=10."
  export OK_SEM_INNER_ITER=100
  export OK_DECODE_ITER=2
  export OK_SENS_SEM_INNER_ITER=10
  export OK_BOOT_SEM_INNER_ITER=10
  export OK_RUN_DECODE=true
  export OK_RUN_SENSITIVITY=true
  export OK_RUN_BOOTSTRAP_ATE=true
fi

"$RSCRIPT_BIN" "$PKG_ROOT/inst/oklahoma/oklahoma_analysis.R" 2>&1

echo ""
echo "=== Done $(date) ==="
