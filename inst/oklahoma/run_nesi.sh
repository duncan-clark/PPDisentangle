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
PP_SEM_WARMSTART_FIXED="${PP_SEM_WARMSTART_FIXED:-0}"
PP_SEM_N_ITER="${PP_SEM_N_ITER:-}"
PP_SEM_N_LABELLINGS="${PP_SEM_N_LABELLINGS:-}"
PP_SEM_OUTER_MAXIT="${PP_SEM_OUTER_MAXIT:-}"
PP_SEM_OUTER_MAXIT_BIV="${PP_SEM_OUTER_MAXIT_BIV:-}"
PP_SEM_T_TRUNC_DAYS="${PP_SEM_T_TRUNC_DAYS:-}"
PP_SEM_T_TRUNC_REL="${PP_SEM_T_TRUNC_REL:-0.05}"
PP_SEM_TEMPORAL_WEIGHT="${PP_SEM_TEMPORAL_WEIGHT:-0}"
PP_SEM_OPTIM_METHOD="${PP_SEM_OPTIM_METHOD:-sample_weighted}"
PP_SEM_SELECTION_TEMPERATURE="${PP_SEM_SELECTION_TEMPERATURE:-0.08}"
PP_SEM_CHANGE_FACTOR_MIN_MULT="${PP_SEM_CHANGE_FACTOR_MIN_MULT:-0.2}"
PP_SEM_CHANGE_FACTOR_MAX_MULT="${PP_SEM_CHANGE_FACTOR_MAX_MULT:-2.0}"
PP_SEM_MAX_RELABEL_STEP_FRAC="${PP_SEM_MAX_RELABEL_STEP_FRAC:-0.05}"
PP_SEM_FORCE_PARAM_UPDATE_FLIP_FRAC="${PP_SEM_FORCE_PARAM_UPDATE_FLIP_FRAC:-0.05}"
PP_RUN_SEM_PILOT="${PP_RUN_SEM_PILOT:-0}"
PP_SEM_PILOT_INNER="${PP_SEM_PILOT_INNER:-100}"
PP_SEM_PILOT_CORES="${PP_SEM_PILOT_CORES:-}"
PP_SEM_PILOT_MAX_COMBOS="${PP_SEM_PILOT_MAX_COMBOS:-24}"
PP_SEM_PILOT_CHANGE_FACTORS="${PP_SEM_PILOT_CHANGE_FACTORS:-}"
PP_SEM_PILOT_MIN_MULTS="${PP_SEM_PILOT_MIN_MULTS:-}"
PP_SEM_PILOT_MAX_MULTS="${PP_SEM_PILOT_MAX_MULTS:-}"
PP_SEM_PILOT_TEMPS="${PP_SEM_PILOT_TEMPS:-}"
PP_SEM_WORKER_LOGS="${PP_SEM_WORKER_LOGS:-1}"
PP_SEM_WORKER_LOG_VERBOSE="${PP_SEM_WORKER_LOG_VERBOSE:-1}"
PP_SEM_WORKER_LOG_SPLIT="${PP_SEM_WORKER_LOG_SPLIT:-0}"
PP_SEM_TIMING_VERBOSE="${PP_SEM_TIMING_VERBOSE:-1}"
PP_SEM_PROPOSAL_VERBOSE="${PP_SEM_PROPOSAL_VERBOSE:-1}"
PP_SIM_PROGRESS_EVERY="${PP_SIM_PROGRESS_EVERY:-10000}"
PP_SENS_SEM_INNER="${PP_SENS_SEM_INNER:-}"
PP_BOOT_SEM_INNER="${PP_BOOT_SEM_INNER:-}"
PP_BOOT_REFIT_SCOPE="${PP_BOOT_REFIT_SCOPE:-}"
PP_BOOT_TARGETS="${PP_BOOT_TARGETS:-C,D}"
PP_KDE_VARIANT_MODE="${PP_KDE_VARIANT_MODE:-}"
PP_BOOT_OUTER_CORES="${PP_BOOT_OUTER_CORES:-}"
PP_RUN_SENSITIVITY="${PP_RUN_SENSITIVITY:-auto}"
PP_MEM="${PP_MEM:-}"
PP_TIME="${PP_TIME:-72:00:00}"
PP_SETUP_TEST="${PP_SETUP_TEST:-0}"
PP_MODE="${PP_MODE:-}"
PP_SEED="${PP_SEED:-1}"
PP_ATE_N_SIMS="${PP_ATE_N_SIMS:-}"
CORES_EXPLICIT=0
MEM_EXPLICIT=0
SEM_INNER_EXPLICIT=0
SEM_WARMSTART_EXPLICIT=0
SEM_N_ITER_EXPLICIT=0
SEM_N_LABELLINGS_EXPLICIT=0
SEM_OUTER_MAXIT_EXPLICIT=0
SEM_OUTER_MAXIT_BIV_EXPLICIT=0
SEM_OPTIM_METHOD_EXPLICIT=0
SEM_SELECTION_TEMP_EXPLICIT=0
SEM_MIN_MULT_EXPLICIT=0
SEM_MAX_MULT_EXPLICIT=0
RUN_SEM_PILOT_EXPLICIT=0
SEM_PILOT_INNER_EXPLICIT=0
SEM_PILOT_CORES_EXPLICIT=0
SEM_PILOT_MAX_COMBOS_EXPLICIT=0
SENS_SEM_INNER_EXPLICIT=0
BOOT_SEM_INNER_EXPLICIT=0
BOOT_REFIT_SCOPE_EXPLICIT=0
BOOT_TARGETS_EXPLICIT=0
KDE_VARIANT_MODE_EXPLICIT=0
SETUP_TEST_EXPLICIT=0
BOOT_REPS_EXPLICIT=0
BOOT_OUTER_CORES_EXPLICIT=0
RUN_SENS_EXPLICIT=0
ATE_N_SIMS_EXPLICIT=0
if [ -n "$PP_BOOT_REPS" ]; then BOOT_REPS_EXPLICIT=1; fi
if [ -n "$PP_BOOT_OUTER_CORES" ]; then BOOT_OUTER_CORES_EXPLICIT=1; fi
if [ "$PP_RUN_SENSITIVITY" != "auto" ]; then RUN_SENS_EXPLICIT=1; fi
if [ -n "$PP_ATE_N_SIMS" ]; then ATE_N_SIMS_EXPLICIT=1; fi

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --mode) PP_MODE="$2"; shift 2 ;;
    --cores) PP_CORES="$2"; CORES_EXPLICIT=1; shift 2 ;;
    --sims) PP_CORES="$2"; CORES_EXPLICIT=1; shift 2 ;;  # alias: match sim_study launcher UX
    --boot-reps) PP_BOOT_REPS="$2"; BOOT_REPS_EXPLICIT=1; shift 2 ;;
    --sem-inner) PP_SEM_INNER="$2"; SEM_INNER_EXPLICIT=1; shift 2 ;;
    --sem-warmstart-fixed) PP_SEM_WARMSTART_FIXED="$2"; SEM_WARMSTART_EXPLICIT=1; shift 2 ;;
    --sem-n-iter) PP_SEM_N_ITER="$2"; SEM_N_ITER_EXPLICIT=1; shift 2 ;;
    --sem-n-labellings) PP_SEM_N_LABELLINGS="$2"; SEM_N_LABELLINGS_EXPLICIT=1; shift 2 ;;
    --sem-outer-maxit) PP_SEM_OUTER_MAXIT="$2"; SEM_OUTER_MAXIT_EXPLICIT=1; shift 2 ;;
    --sem-outer-maxit-biv) PP_SEM_OUTER_MAXIT_BIV="$2"; SEM_OUTER_MAXIT_BIV_EXPLICIT=1; shift 2 ;;
    --sem-t-trunc-days) PP_SEM_T_TRUNC_DAYS="$2"; shift 2 ;;
    --sem-t-trunc-rel) PP_SEM_T_TRUNC_REL="$2"; shift 2 ;;
    --sem-temporal-weight) PP_SEM_TEMPORAL_WEIGHT="$2"; shift 2 ;;
    --sem-optim-method) PP_SEM_OPTIM_METHOD="$2"; SEM_OPTIM_METHOD_EXPLICIT=1; shift 2 ;;
    --sem-selection-temperature) PP_SEM_SELECTION_TEMPERATURE="$2"; SEM_SELECTION_TEMP_EXPLICIT=1; shift 2 ;;
    --sem-change-factor-min-mult) PP_SEM_CHANGE_FACTOR_MIN_MULT="$2"; SEM_MIN_MULT_EXPLICIT=1; shift 2 ;;
    --sem-change-factor-max-mult) PP_SEM_CHANGE_FACTOR_MAX_MULT="$2"; SEM_MAX_MULT_EXPLICIT=1; shift 2 ;;
    --sem-max-relabel-step-frac) PP_SEM_MAX_RELABEL_STEP_FRAC="$2"; shift 2 ;;
    --sem-force-param-update-flip-frac) PP_SEM_FORCE_PARAM_UPDATE_FLIP_FRAC="$2"; shift 2 ;;
    --run-sem-pilot) PP_RUN_SEM_PILOT="$2"; RUN_SEM_PILOT_EXPLICIT=1; shift 2 ;;
    --sem-pilot-inner) PP_SEM_PILOT_INNER="$2"; SEM_PILOT_INNER_EXPLICIT=1; shift 2 ;;
    --sem-pilot-cores) PP_SEM_PILOT_CORES="$2"; SEM_PILOT_CORES_EXPLICIT=1; shift 2 ;;
    --sem-pilot-max-combos) PP_SEM_PILOT_MAX_COMBOS="$2"; SEM_PILOT_MAX_COMBOS_EXPLICIT=1; shift 2 ;;
    --sem-pilot-change-factors) PP_SEM_PILOT_CHANGE_FACTORS="$2"; shift 2 ;;
    --sem-pilot-min-mults) PP_SEM_PILOT_MIN_MULTS="$2"; shift 2 ;;
    --sem-pilot-max-mults) PP_SEM_PILOT_MAX_MULTS="$2"; shift 2 ;;
    --sem-pilot-temps) PP_SEM_PILOT_TEMPS="$2"; shift 2 ;;
    --sem-worker-logs) PP_SEM_WORKER_LOGS="$2"; shift 2 ;;
    --sem-worker-log-verbose) PP_SEM_WORKER_LOG_VERBOSE="$2"; shift 2 ;;
    --sem-worker-log-split) PP_SEM_WORKER_LOG_SPLIT="$2"; shift 2 ;;
    --sem-timing-verbose) PP_SEM_TIMING_VERBOSE="$2"; shift 2 ;;
    --sem-proposal-verbose) PP_SEM_PROPOSAL_VERBOSE="$2"; shift 2 ;;
    --sim-progress-every) PP_SIM_PROGRESS_EVERY="$2"; shift 2 ;;
    --sens-sem-inner) PP_SENS_SEM_INNER="$2"; SENS_SEM_INNER_EXPLICIT=1; shift 2 ;;
    --boot-sem-inner) PP_BOOT_SEM_INNER="$2"; BOOT_SEM_INNER_EXPLICIT=1; shift 2 ;;
    --boot-refit-scope) PP_BOOT_REFIT_SCOPE="$2"; BOOT_REFIT_SCOPE_EXPLICIT=1; shift 2 ;;
    --boot-targets) PP_BOOT_TARGETS="$2"; BOOT_TARGETS_EXPLICIT=1; shift 2 ;;
    --kde-variant-mode) PP_KDE_VARIANT_MODE="$2"; KDE_VARIANT_MODE_EXPLICIT=1; shift 2 ;;
    --boot-outer-cores) PP_BOOT_OUTER_CORES="$2"; BOOT_OUTER_CORES_EXPLICIT=1; shift 2 ;;
    --run-sensitivity) PP_RUN_SENSITIVITY="$2"; RUN_SENS_EXPLICIT=1; shift 2 ;;
    --ate-n-sims) PP_ATE_N_SIMS="$2"; ATE_N_SIMS_EXPLICIT=1; shift 2 ;;
    --setup-test) PP_SETUP_TEST=1; SETUP_TEST_EXPLICIT=1; shift ;;
    --mem) PP_MEM="$2"; MEM_EXPLICIT=1; shift 2 ;;
    --time) PP_TIME="$2"; shift 2 ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

if [ -n "$PP_MODE" ]; then
  mode_norm="$(echo "$PP_MODE" | tr '[:upper:]' '[:lower:]')"
  case "$mode_norm" in
    quick)
      if [ "$SETUP_TEST_EXPLICIT" -ne 1 ]; then PP_SETUP_TEST=0; fi
      if [ "$CORES_EXPLICIT" -ne 1 ]; then PP_CORES=32; fi
      # Quick profile default: disable bootstrap to minimize turnaround.
      if [ "$BOOT_REPS_EXPLICIT" -ne 1 ]; then PP_BOOT_REPS=0; fi
      if [ "$SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SEM_INNER=200; fi
      if [ "$SEM_N_ITER_EXPLICIT" -ne 1 ]; then PP_SEM_N_ITER=1; fi
      if [ "$SEM_N_LABELLINGS_EXPLICIT" -ne 1 ]; then PP_SEM_N_LABELLINGS=20; fi
      if [ "$SEM_OUTER_MAXIT_EXPLICIT" -ne 1 ]; then PP_SEM_OUTER_MAXIT=120; fi
      if [ "$SEM_OUTER_MAXIT_BIV_EXPLICIT" -ne 1 ]; then PP_SEM_OUTER_MAXIT_BIV=1000; fi
      if [ "$SENS_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SENS_SEM_INNER=200; fi
      if [ "$BOOT_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_BOOT_SEM_INNER=50; fi
      if [ "$BOOT_REFIT_SCOPE_EXPLICIT" -ne 1 ]; then PP_BOOT_REFIT_SCOPE="none"; fi
      if [ "$BOOT_TARGETS_EXPLICIT" -ne 1 ]; then PP_BOOT_TARGETS="C,D"; fi
      if [ "$KDE_VARIANT_MODE_EXPLICIT" -ne 1 ]; then PP_KDE_VARIANT_MODE="triple"; fi
      if [ "$BOOT_OUTER_CORES_EXPLICIT" -ne 1 ]; then PP_BOOT_OUTER_CORES="$PP_CORES"; fi
      if [ "$RUN_SENS_EXPLICIT" -ne 1 ]; then PP_RUN_SENSITIVITY=1; fi
      if [ "$ATE_N_SIMS_EXPLICIT" -ne 1 ]; then PP_ATE_N_SIMS=100; fi
      if [ "$MEM_EXPLICIT" -ne 1 ]; then PP_MEM=96G; fi
      ;;
    test|setup-test|very-quick|veryquick|smoke)
      if [ "$SETUP_TEST_EXPLICIT" -ne 1 ]; then PP_SETUP_TEST=1; fi
      if [ "$CORES_EXPLICIT" -ne 1 ]; then PP_CORES=32; fi
      if [ "$BOOT_REPS_EXPLICIT" -ne 1 ]; then PP_BOOT_REPS="$PP_CORES"; fi
      if [ "$SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SEM_INNER=100; fi
      if [ "$SEM_N_ITER_EXPLICIT" -ne 1 ]; then PP_SEM_N_ITER=1; fi
      if [ "$SEM_N_LABELLINGS_EXPLICIT" -ne 1 ]; then PP_SEM_N_LABELLINGS=5; fi
      if [ "$SEM_OUTER_MAXIT_EXPLICIT" -ne 1 ]; then PP_SEM_OUTER_MAXIT=40; fi
      if [ "$SEM_OUTER_MAXIT_BIV_EXPLICIT" -ne 1 ]; then PP_SEM_OUTER_MAXIT_BIV=1000; fi
      if [ "$SENS_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SENS_SEM_INNER=2; fi
      if [ "$BOOT_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_BOOT_SEM_INNER=2; fi
      if [ "$BOOT_REFIT_SCOPE_EXPLICIT" -ne 1 ]; then PP_BOOT_REFIT_SCOPE="none"; fi
      if [ "$BOOT_OUTER_CORES_EXPLICIT" -ne 1 ]; then PP_BOOT_OUTER_CORES=1; fi
      if [ "$RUN_SENS_EXPLICIT" -ne 1 ]; then PP_RUN_SENSITIVITY=0; fi
      if [ "$ATE_N_SIMS_EXPLICIT" -ne 1 ]; then PP_ATE_N_SIMS=20; fi
      if [ "$MEM_EXPLICIT" -ne 1 ]; then PP_MEM=64G; fi
      ;;
    full|default|long|big)
      if [ "$SETUP_TEST_EXPLICIT" -ne 1 ]; then PP_SETUP_TEST=0; fi
      if [ "$CORES_EXPLICIT" -ne 1 ]; then PP_CORES=32; fi
      # Full/default production profile: enable enough refits for stable ATE SD.
      if [ "$BOOT_REPS_EXPLICIT" -ne 1 ]; then PP_BOOT_REPS="$PP_CORES"; fi
      if [ "$SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SEM_INNER=1000; fi
      if [ "$SEM_N_ITER_EXPLICIT" -ne 1 ]; then PP_SEM_N_ITER=1; fi
      if [ "$SEM_N_LABELLINGS_EXPLICIT" -ne 1 ]; then PP_SEM_N_LABELLINGS=20; fi
      if [ "$SEM_OUTER_MAXIT_EXPLICIT" -ne 1 ]; then PP_SEM_OUTER_MAXIT=220; fi
      if [ "$SEM_OUTER_MAXIT_BIV_EXPLICIT" -ne 1 ]; then PP_SEM_OUTER_MAXIT_BIV=1000; fi
      if [ "$SENS_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SENS_SEM_INNER=1000; fi
      if [ "$BOOT_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_BOOT_SEM_INNER=1000; fi
      # Full profile default: full parametric bootstrap (simulate + refit + ATE).
      if [ "$BOOT_REFIT_SCOPE_EXPLICIT" -ne 1 ]; then PP_BOOT_REFIT_SCOPE="full"; fi
      if [ "$BOOT_TARGETS_EXPLICIT" -ne 1 ]; then PP_BOOT_TARGETS="C,D"; fi
      if [ "$KDE_VARIANT_MODE_EXPLICIT" -ne 1 ]; then PP_KDE_VARIANT_MODE="triple"; fi
      if [ "$BOOT_OUTER_CORES_EXPLICIT" -ne 1 ]; then PP_BOOT_OUTER_CORES=$(( PP_CORES < 6 ? PP_CORES : 6 )); fi
      if [ "$RUN_SENS_EXPLICIT" -ne 1 ]; then PP_RUN_SENSITIVITY=0; fi
      if [ "$ATE_N_SIMS_EXPLICIT" -ne 1 ]; then PP_ATE_N_SIMS=500; fi
      if [ "$MEM_EXPLICIT" -ne 1 ]; then PP_MEM=200G; fi
      ;;
    *)
      echo "Unknown --mode '$PP_MODE' (expected: test | quick | full)"
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

if [ -z "$PP_SENS_SEM_INNER" ]; then
  PP_SENS_SEM_INNER="$PP_SEM_INNER"
fi
if [ -z "$PP_SEM_N_ITER" ]; then
  PP_SEM_N_ITER=10
fi
if [ -z "$PP_SEM_N_LABELLINGS" ]; then
  PP_SEM_N_LABELLINGS=20
fi
if [ -z "$PP_SEM_OUTER_MAXIT" ]; then
  PP_SEM_OUTER_MAXIT=220
fi
if [ -z "$PP_SEM_OUTER_MAXIT_BIV" ]; then
  PP_SEM_OUTER_MAXIT_BIV=1000
fi
if [ -z "$PP_BOOT_SEM_INNER" ]; then
  PP_BOOT_SEM_INNER="$PP_SEM_INNER"
fi
if [ -z "$PP_BOOT_REFIT_SCOPE" ]; then
  PP_BOOT_REFIT_SCOPE="none"
fi
boot_refit_norm="$(echo "$PP_BOOT_REFIT_SCOPE" | tr '[:upper:]' '[:lower:]')"
if [ "$boot_refit_norm" = "none" ] || [ "$boot_refit_norm" = "partial" ] || [ "$boot_refit_norm" = "full" ]; then
  PP_BOOT_REFIT_SCOPE="$boot_refit_norm"
else
  echo "Invalid --boot-refit-scope '$PP_BOOT_REFIT_SCOPE' (expected: none | partial | full)"
  exit 1
fi
if [ -z "$PP_BOOT_OUTER_CORES" ]; then
  PP_BOOT_OUTER_CORES="$PP_CORES"
fi
# Default policy (manual mode): one bootstrap replicate per available outer bootstrap core.
# Do not override mode-specific defaults already assigned above.
if [ "$BOOT_REPS_EXPLICIT" -ne 1 ] && [ -z "${PP_BOOT_REPS:-}" ]; then
  PP_BOOT_REPS="$PP_BOOT_OUTER_CORES"
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
  OUTPUT_DIR="$PKG_ROOT/output/oklahoma"
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

  # Export vars in parent shell, then forward all env with --export=ALL.
  # This is more robust than a giant --export=... CSV and avoids truncation/parsing edge cases.
  export PKG_ROOT PP_MODE PP_CORES PP_BOOT_REPS PP_SEM_INNER PP_SEM_WARMSTART_FIXED PP_SEM_N_ITER PP_SEM_N_LABELLINGS
  export PP_SEM_OUTER_MAXIT PP_SEM_OUTER_MAXIT_BIV PP_SEM_T_TRUNC_DAYS PP_SEM_T_TRUNC_REL PP_SEM_TEMPORAL_WEIGHT
  export PP_SEM_WORKER_LOGS PP_SEM_WORKER_LOG_VERBOSE PP_SEM_WORKER_LOG_SPLIT PP_SEM_TIMING_VERBOSE PP_SEM_PROPOSAL_VERBOSE
  export PP_SIM_PROGRESS_EVERY PP_SENS_SEM_INNER PP_BOOT_SEM_INNER PP_BOOT_REFIT_SCOPE PP_BOOT_TARGETS PP_KDE_VARIANT_MODE
  export PP_BOOT_OUTER_CORES PP_ATE_N_SIMS PP_RUN_SENSITIVITY PP_MEM PP_TIME PP_SEM_OPTIM_METHOD PP_SEM_SELECTION_TEMPERATURE
  export PP_SEM_CHANGE_FACTOR_MIN_MULT PP_SEM_CHANGE_FACTOR_MAX_MULT PP_SEM_MAX_RELABEL_STEP_FRAC PP_SEM_FORCE_PARAM_UPDATE_FLIP_FRAC
  export PP_RUN_SEM_PILOT PP_SEM_PILOT_INNER PP_SEM_PILOT_CORES PP_SEM_PILOT_MAX_COMBOS PP_SEM_PILOT_CHANGE_FACTORS
  export PP_SEM_PILOT_MIN_MULTS PP_SEM_PILOT_MAX_MULTS PP_SEM_PILOT_TEMPS PP_SETUP_TEST
  [ -n "${PP_R_GEO_MODULE:-}" ] && export PP_R_GEO_MODULE

  echo "Submitting Oklahoma job: mode=${PP_MODE:-manual} cores=$PP_CORES sem_n_iter=$PP_SEM_N_ITER sem_n_labellings=$PP_SEM_N_LABELLINGS sem_outer_maxit=$PP_SEM_OUTER_MAXIT sem_outer_maxit_biv=$PP_SEM_OUTER_MAXIT_BIV sem_t_trunc_days=${PP_SEM_T_TRUNC_DAYS:-auto} sem_t_trunc_rel=$PP_SEM_T_TRUNC_REL sem_temporal_weight=$PP_SEM_TEMPORAL_WEIGHT sem_warmstart_fixed=$PP_SEM_WARMSTART_FIXED sem_optim=$PP_SEM_OPTIM_METHOD sem_temp=$PP_SEM_SELECTION_TEMPERATURE sem_cf_min=$PP_SEM_CHANGE_FACTOR_MIN_MULT sem_cf_max=$PP_SEM_CHANGE_FACTOR_MAX_MULT sem_step_cap=$PP_SEM_MAX_RELABEL_STEP_FRAC sem_force_refit_frac=$PP_SEM_FORCE_PARAM_UPDATE_FLIP_FRAC run_sem_pilot=$PP_RUN_SEM_PILOT sem_pilot_inner=$PP_SEM_PILOT_INNER sem_pilot_cores=${PP_SEM_PILOT_CORES:-auto} sem_pilot_max_combos=$PP_SEM_PILOT_MAX_COMBOS sem_worker_logs=$PP_SEM_WORKER_LOGS sem_worker_log_verbose=$PP_SEM_WORKER_LOG_VERBOSE sem_worker_log_split=$PP_SEM_WORKER_LOG_SPLIT sem_timing_verbose=$PP_SEM_TIMING_VERBOSE sem_proposal_verbose=$PP_SEM_PROPOSAL_VERBOSE sim_progress_every=$PP_SIM_PROGRESS_EVERY ate_n_sims=$PP_ATE_N_SIMS boot_reps=$PP_BOOT_REPS sem_inner=$PP_SEM_INNER sens_inner=$PP_SENS_SEM_INNER boot_inner=$PP_BOOT_SEM_INNER boot_refit_scope=$PP_BOOT_REFIT_SCOPE kde_variants=$PP_KDE_VARIANT_MODE boot_outer_cores=$PP_BOOT_OUTER_CORES setup_test=$PP_SETUP_TEST"

  JOB_ID=$(sbatch --parsable \
    --cpus-per-task="$PP_CORES" \
    --mem="$PP_MEM" \
    --time="$PP_TIME" \
    $EXTRA_SBATCH \
    --export=ALL \
    --output="$OUTPUT_DIR/%j_oklahoma_slurm.out" \
    --error="$OUTPUT_DIR/%j_oklahoma_slurm.err" \
    "$SCRIPT_DIR/run_nesi.sh")

  echo "Job $JOB_ID submitted"
  echo "SLURM out: output/oklahoma/${JOB_ID}_oklahoma_slurm.out"
  echo "SLURM err: output/oklahoma/${JOB_ID}_oklahoma_slurm.err"
  exit 0
fi

# ----------------------------
# Job mode
# ----------------------------
cd "$PKG_ROOT"
mkdir -p "$PKG_ROOT/output/oklahoma"

echo "=== PPDisentangle Oklahoma (NeSI) ==="
echo "Job: ${SLURM_JOB_ID} | $(date)"
echo "Node: $(hostname) | Partition: ${SLURM_JOB_PARTITION:-unknown}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-$PP_CORES}"
echo "boot_reps=$PP_BOOT_REPS sem_n_iter=$PP_SEM_N_ITER sem_outer_maxit=$PP_SEM_OUTER_MAXIT sem_outer_maxit_biv=$PP_SEM_OUTER_MAXIT_BIV sem_t_trunc_days=${PP_SEM_T_TRUNC_DAYS:-auto} sem_t_trunc_rel=$PP_SEM_T_TRUNC_REL sem_temporal_weight=$PP_SEM_TEMPORAL_WEIGHT sem_warmstart_fixed=$PP_SEM_WARMSTART_FIXED sem_optim=$PP_SEM_OPTIM_METHOD sem_temp=$PP_SEM_SELECTION_TEMPERATURE sem_cf_min=$PP_SEM_CHANGE_FACTOR_MIN_MULT sem_cf_max=$PP_SEM_CHANGE_FACTOR_MAX_MULT run_sem_pilot=$PP_RUN_SEM_PILOT sem_pilot_inner=$PP_SEM_PILOT_INNER sem_pilot_cores=${PP_SEM_PILOT_CORES:-auto} sem_pilot_max_combos=$PP_SEM_PILOT_MAX_COMBOS sem_worker_logs=$PP_SEM_WORKER_LOGS sem_worker_log_verbose=$PP_SEM_WORKER_LOG_VERBOSE sem_worker_log_split=$PP_SEM_WORKER_LOG_SPLIT sem_timing_verbose=$PP_SEM_TIMING_VERBOSE sem_proposal_verbose=$PP_SEM_PROPOSAL_VERBOSE sim_progress_every=$PP_SIM_PROGRESS_EVERY sem_inner=$PP_SEM_INNER sens_inner=$PP_SENS_SEM_INNER boot_inner=$PP_BOOT_SEM_INNER boot_refit_scope=$PP_BOOT_REFIT_SCOPE targets=$PP_BOOT_TARGETS kde_variants=$PP_KDE_VARIANT_MODE"
echo "setup_test=$PP_SETUP_TEST mode=${PP_MODE:-manual}"
echo "seed=$PP_SEED (fit jobs RNG de-correlated by model; bootstrap RNG de-correlated by replicate)"
echo "ENV CHECK: PP_SEM_INNER=$PP_SEM_INNER | PP_SENS_SEM_INNER=$PP_SENS_SEM_INNER | PP_BOOT_SEM_INNER=$PP_BOOT_SEM_INNER"
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
export OK_MEMORY_SAFE=true
export OK_PARALLEL_BACKEND=psock
export OK_CORES="${JOB_CORES}"
export OK_SENS_CORES="${JOB_CORES}"
# Default ATE sims to single-core execution unless explicitly overridden.
export OK_ATE_SIM_CORES="${OK_ATE_SIM_CORES:-${PP_ATE_SIM_CORES:-1}}"
export OK_BOOT_OUTER_CAP_MEMSAFE="${JOB_CORES}"
export OK_VERBOSE=false
export OK_SEM_INNER_ITER="$PP_SEM_INNER"
export OK_SEM_WARMSTART_FIXED="$PP_SEM_WARMSTART_FIXED"
export OK_SEM_N_ITER="$PP_SEM_N_ITER"
export OK_SEM_N_LABELLINGS="$PP_SEM_N_LABELLINGS"
export OK_SEM_OUTER_MAXIT="$PP_SEM_OUTER_MAXIT"
export OK_SEM_OUTER_MAXIT_BIV="$PP_SEM_OUTER_MAXIT_BIV"
export OK_SEM_T_TRUNC_DAYS="$PP_SEM_T_TRUNC_DAYS"
export OK_SEM_T_TRUNC_REL="$PP_SEM_T_TRUNC_REL"
export OK_SEM_TEMPORAL_WEIGHT="$PP_SEM_TEMPORAL_WEIGHT"
export OK_SEM_OPTIM_METHOD="$PP_SEM_OPTIM_METHOD"
export OK_SEM_SELECTION_TEMPERATURE="$PP_SEM_SELECTION_TEMPERATURE"
export OK_SEM_CHANGE_FACTOR_MIN_MULT="$PP_SEM_CHANGE_FACTOR_MIN_MULT"
export OK_SEM_CHANGE_FACTOR_MAX_MULT="$PP_SEM_CHANGE_FACTOR_MAX_MULT"
export OK_SEM_MAX_RELABEL_STEP_FRAC="$PP_SEM_MAX_RELABEL_STEP_FRAC"
export OK_SEM_FORCE_PARAM_UPDATE_FLIP_FRAC="$PP_SEM_FORCE_PARAM_UPDATE_FLIP_FRAC"
if [ "$PP_RUN_SEM_PILOT" = "1" ] || [ "$PP_RUN_SEM_PILOT" = "true" ] || [ "$PP_RUN_SEM_PILOT" = "yes" ]; then
  export OK_RUN_SEM_PILOT=true
else
  export OK_RUN_SEM_PILOT=false
fi
export OK_SEM_PILOT_INNER_ITER="$PP_SEM_PILOT_INNER"
if [ -n "${PP_SEM_PILOT_CORES:-}" ]; then export OK_SEM_PILOT_CORES="$PP_SEM_PILOT_CORES"; fi
export OK_SEM_PILOT_MAX_COMBOS="$PP_SEM_PILOT_MAX_COMBOS"
if [ -n "${PP_SEM_PILOT_CHANGE_FACTORS:-}" ]; then export OK_SEM_PILOT_CHANGE_FACTORS="$PP_SEM_PILOT_CHANGE_FACTORS"; fi
if [ -n "${PP_SEM_PILOT_MIN_MULTS:-}" ]; then export OK_SEM_PILOT_MIN_MULTS="$PP_SEM_PILOT_MIN_MULTS"; fi
if [ -n "${PP_SEM_PILOT_MAX_MULTS:-}" ]; then export OK_SEM_PILOT_MAX_MULTS="$PP_SEM_PILOT_MAX_MULTS"; fi
if [ -n "${PP_SEM_PILOT_TEMPS:-}" ]; then export OK_SEM_PILOT_TEMPS="$PP_SEM_PILOT_TEMPS"; fi
export OK_SEM_WORKER_LOGS="$PP_SEM_WORKER_LOGS"
export OK_SEM_WORKER_LOG_VERBOSE="$PP_SEM_WORKER_LOG_VERBOSE"
export OK_SEM_WORKER_LOG_SPLIT="$PP_SEM_WORKER_LOG_SPLIT"
export OK_SEM_TIMING_VERBOSE="$PP_SEM_TIMING_VERBOSE"
export OK_SEM_PROPOSAL_VERBOSE="$PP_SEM_PROPOSAL_VERBOSE"
export OK_SIM_PROGRESS_VERBOSE="$PP_SEM_PROPOSAL_VERBOSE"
export OK_SIM_PROGRESS_EVERY="$PP_SIM_PROGRESS_EVERY"
export OK_SENS_SEM_INNER_ITER="$PP_SENS_SEM_INNER"
if [ "$PP_RUN_SENSITIVITY" = "1" ] || [ "$PP_RUN_SENSITIVITY" = "true" ] || [ "$PP_RUN_SENSITIVITY" = "yes" ]; then
  export OK_RUN_SENSITIVITY=true
else
  export OK_RUN_SENSITIVITY=false
fi
export OK_RUN_BOOTSTRAP_ATE=true
export OK_BOOT_N_REPS="$PP_BOOT_REPS"
export OK_BOOT_REFIT_SCOPE="$PP_BOOT_REFIT_SCOPE"
export OK_BOOT_TARGETS="$PP_BOOT_TARGETS"
export OK_KDE_VARIANT_MODE="$PP_KDE_VARIANT_MODE"
export OK_BOOT_SEM_INNER_ITER="$PP_BOOT_SEM_INNER"
export OK_BOOT_OUTER_CORES="$PP_BOOT_OUTER_CORES"
export OK_ATE_N_SIMS="$PP_ATE_N_SIMS"
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
    echo "Applying legacy very-quick alias runtime overrides (equivalent to test profile intent)."
    export OK_RUN_SENSITIVITY=true
    export OK_RUN_BOOTSTRAP_ATE=true
    export OK_SEM_INNER_ITER=2
    export OK_SEM_N_ITER=1
    export OK_SEM_OUTER_MAXIT=20
    export OK_SEM_OUTER_MAXIT_BIV=20
    export OK_SENS_SEM_INNER_ITER=2
    export OK_BOOT_SEM_INNER_ITER=2
    export OK_BOOT_OUTER_CORES="${JOB_CORES}"
    export OK_BOOT_N_REPS="${PP_BOOT_REPS:-1}"
    export OK_BOOT_TARGETS="${PP_BOOT_TARGETS:-C,D}"
  fi
fi

if [ "$PP_SETUP_TEST" = "1" ]; then
  echo "Applying setup-test profile: main SEM inner=100, sensitivity inner=2, bootstrap inner=2, sequential bootstrap."
  export OK_SEM_INNER_ITER=100
  export OK_SEM_N_ITER=1
  export OK_SENS_SEM_INNER_ITER=2
  export OK_BOOT_SEM_INNER_ITER=2
  export OK_SENS_CORES=1
  export OK_ATE_SIM_CORES=1
  export OK_BOOT_OUTER_CORES=1
  if [ "$BOOT_REPS_EXPLICIT" -ne 1 ]; then
    export OK_BOOT_N_REPS="$OK_BOOT_OUTER_CORES"
  fi
  if [ "$RUN_SENS_EXPLICIT" -ne 1 ]; then
    export OK_RUN_SENSITIVITY=false
  fi
  export OK_BOOT_TARGETS="C,D"
  export OK_RUN_BOOTSTRAP_ATE=true
fi

"$RSCRIPT_BIN" "$PKG_ROOT/inst/oklahoma/oklahoma_analysis.R" 2>&1

echo ""
echo "=== Done $(date) ==="
