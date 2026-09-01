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
PP_SEM_T_TRUNC_DAYS="${PP_SEM_T_TRUNC_DAYS:-90}"
PP_SEM_T_TRUNC_REL="${PP_SEM_T_TRUNC_REL:-0.05}"
PP_T_TRUNC_SENS_DAYS="${PP_T_TRUNC_SENS_DAYS:-1,5,7,10,14,21}"
PP_RUN_T_TRUNC_SENSITIVITY="${PP_RUN_T_TRUNC_SENSITIVITY:-1}"
PP_SEM_TEMPORAL_WEIGHT="${PP_SEM_TEMPORAL_WEIGHT:-0}"
PP_SEM_OPTIM_METHOD="${PP_SEM_OPTIM_METHOD:-sample_weighted}"
PP_SEM_SELECTION_TEMPERATURE="${PP_SEM_SELECTION_TEMPERATURE:-0.2}"
PP_SEM_CHANGE_FACTOR_MIN_MULT="${PP_SEM_CHANGE_FACTOR_MIN_MULT:-1.0}"
PP_SEM_CHANGE_FACTOR_MAX_MULT="${PP_SEM_CHANGE_FACTOR_MAX_MULT:-5.0}"
PP_SEM_MAX_RELABEL_STEP_FRAC="${PP_SEM_MAX_RELABEL_STEP_FRAC:-0.05}"
PP_SEM_FORCE_PARAM_UPDATE_FLIP_FRAC="${PP_SEM_FORCE_PARAM_UPDATE_FLIP_FRAC:-0.05}"
PP_SEM_MONOTONE_COMPLETE_LL="${PP_SEM_MONOTONE_COMPLETE_LL:-0}"
PP_SEM_START_FROM_C="${PP_SEM_START_FROM_C:-0}"
PP_SEM_BIV_N_THREADS="${PP_SEM_BIV_N_THREADS:-1}"
PP_SEM_SINGLE_FLIP_FROM_ITER="${PP_SEM_SINGLE_FLIP_FROM_ITER:-}"
PP_RUN_SEM_PILOT="${PP_RUN_SEM_PILOT:-0}"
PP_SEM_PILOT_INNER="${PP_SEM_PILOT_INNER:-100}"
PP_SEM_PILOT_CORES="${PP_SEM_PILOT_CORES:-}"
PP_SEM_PILOT_MAX_COMBOS="${PP_SEM_PILOT_MAX_COMBOS:-24}"
PP_SEM_PILOT_CHANGE_FACTORS="${PP_SEM_PILOT_CHANGE_FACTORS:-}"
PP_SEM_PILOT_MIN_MULTS="${PP_SEM_PILOT_MIN_MULTS:-}"
PP_SEM_PILOT_MAX_MULTS="${PP_SEM_PILOT_MAX_MULTS:-}"
PP_SEM_PILOT_TEMPS="${PP_SEM_PILOT_TEMPS:-}"
PP_SEM_WORKER_LOGS="${PP_SEM_WORKER_LOGS:-1}"
PP_SEM_WORKER_LOG_VERBOSE="${PP_SEM_WORKER_LOG_VERBOSE:-0}"
PP_SEM_WORKER_LOG_SPLIT="${PP_SEM_WORKER_LOG_SPLIT:-0}"
PP_SEM_TIMING_VERBOSE="${PP_SEM_TIMING_VERBOSE:-0}"
PP_SEM_PROPOSAL_VERBOSE="${PP_SEM_PROPOSAL_VERBOSE:-0}"
PP_SIM_PROGRESS_EVERY="${PP_SIM_PROGRESS_EVERY:-10000}"
PP_SENS_SEM_INNER="${PP_SENS_SEM_INNER:-}"
PP_BOOT_SEM_INNER="${PP_BOOT_SEM_INNER:-}"
PP_BOOT_REFIT_SCOPE="${PP_BOOT_REFIT_SCOPE:-}"
PP_BOOT_TARGETS="${PP_BOOT_TARGETS:-C,D}"
PP_KDE_VARIANT_MODE="${PP_KDE_VARIANT_MODE:-}"
PP_KDE_BW_METHOD="${PP_KDE_BW_METHOD:-}"
PP_KDE_BW_SENS_KM="${PP_KDE_BW_SENS_KM:-}"
PP_PRIMARY_PARTITION="${PP_PRIMARY_PARTITION:-}"
PP_RUN_PARTITION_SENSITIVITY="${PP_RUN_PARTITION_SENSITIVITY:-}"
PP_BOOT_OUTER_CORES="${PP_BOOT_OUTER_CORES:-}"
PP_RUN_SENSITIVITY="${PP_RUN_SENSITIVITY:-auto}"
PP_RUN_FIT_VARIABILITY="${PP_RUN_FIT_VARIABILITY:-0}"
PP_FIT_VARIABILITY_REPS="${PP_FIT_VARIABILITY_REPS:-}"
PP_FIT_VARIABILITY_CORES="${PP_FIT_VARIABILITY_CORES:-}"
PP_FIT_VARIABILITY_PATCH_FILE="${PP_FIT_VARIABILITY_PATCH_FILE:-}"
PP_FIT_VARIABILITY_ONLY="${PP_FIT_VARIABILITY_ONLY:-0}"
PP_BOOTSTRAP_PATCH_FILE="${PP_BOOTSTRAP_PATCH_FILE:-}"
PP_BOOTSTRAP_ONLY="${PP_BOOTSTRAP_ONLY:-0}"
PP_T_TRUNC_SENS_PATCH_FILE="${PP_T_TRUNC_SENS_PATCH_FILE:-}"
PP_T_TRUNC_SENS_ONLY="${PP_T_TRUNC_SENS_ONLY:-0}"
PP_SMOKE_SEM_D_SEEDS="${PP_SMOKE_SEM_D_SEEDS:-0}"
PP_SMOKE_SEM_D_TRUNC="${PP_SMOKE_SEM_D_TRUNC:-3}"
PP_SKIP_FULL_REPORT="${PP_SKIP_FULL_REPORT:-0}"
PP_CD_ONLY="${PP_CD_ONLY:-0}"
PP_UNIV_KDE_ONLY="${PP_UNIV_KDE_ONLY:-0}"
PP_SKIP_CONTROL_SNAPSHOTS="${PP_SKIP_CONTROL_SNAPSHOTS:-}"
PP_MEM="${PP_MEM:-}"
PP_TIME="${PP_TIME:-72:00:00}"
PP_SETUP_TEST="${PP_SETUP_TEST:-0}"
PP_MODE="${PP_MODE:-}"
PP_SEED="${PP_SEED:-1}"
PP_ATE_N_SIMS="${PP_ATE_N_SIMS:-}"
PP_ATE_BIVARIATE="${PP_ATE_BIVARIATE:-}"
PP_ATE_CONTRAST="${PP_ATE_CONTRAST:-}"
PP_ATE_SCENARIO="${PP_ATE_SCENARIO:-}"
CORES_EXPLICIT=0
MEM_EXPLICIT=0
SEM_INNER_EXPLICIT=0
SEM_WARMSTART_EXPLICIT=0
SEM_N_ITER_EXPLICIT=0
SEM_N_LABELLINGS_EXPLICIT=0
SEM_T_TRUNC_EXPLICIT=0
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
RUN_FIT_VARIABILITY_EXPLICIT=0
FIT_VARIABILITY_REPS_EXPLICIT=0
FIT_VARIABILITY_CORES_EXPLICIT=0
ATE_N_SIMS_EXPLICIT=0
if [ -n "$PP_BOOT_REPS" ]; then BOOT_REPS_EXPLICIT=1; fi
if [ -n "$PP_BOOT_OUTER_CORES" ]; then BOOT_OUTER_CORES_EXPLICIT=1; fi
if [ "$PP_RUN_SENSITIVITY" != "auto" ]; then RUN_SENS_EXPLICIT=1; fi
if [ "$PP_RUN_FIT_VARIABILITY" = "1" ] || [ "$PP_RUN_FIT_VARIABILITY" = "true" ] || [ "$PP_RUN_FIT_VARIABILITY" = "yes" ]; then RUN_FIT_VARIABILITY_EXPLICIT=1; fi
if [ -n "$PP_FIT_VARIABILITY_REPS" ]; then FIT_VARIABILITY_REPS_EXPLICIT=1; fi
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
    --sem-t-trunc-days) PP_SEM_T_TRUNC_DAYS="$2"; SEM_T_TRUNC_EXPLICIT=1; shift 2 ;;
    --sem-t-trunc-rel) PP_SEM_T_TRUNC_REL="$2"; shift 2 ;;
    --t-trunc-sens-days) PP_T_TRUNC_SENS_DAYS="$2"; shift 2 ;;
    --run-t-trunc-sensitivity) PP_RUN_T_TRUNC_SENSITIVITY="$2"; shift 2 ;;
    --t-trunc-sens-patch-file) PP_T_TRUNC_SENS_PATCH_FILE="$2"; shift 2 ;;
    --t-trunc-sens-only) PP_T_TRUNC_SENS_ONLY=1; shift ;;
    --smoke-sem-d-seeds) PP_SMOKE_SEM_D_SEEDS="$2"; shift 2 ;;
    --smoke-sem-d-trunc) PP_SMOKE_SEM_D_TRUNC="$2"; shift 2 ;;
    --sem-temporal-weight) PP_SEM_TEMPORAL_WEIGHT="$2"; shift 2 ;;
    --sem-optim-method) PP_SEM_OPTIM_METHOD="$2"; SEM_OPTIM_METHOD_EXPLICIT=1; shift 2 ;;
    --sem-selection-temperature) PP_SEM_SELECTION_TEMPERATURE="$2"; SEM_SELECTION_TEMP_EXPLICIT=1; shift 2 ;;
    --sem-change-factor-min-mult) PP_SEM_CHANGE_FACTOR_MIN_MULT="$2"; SEM_MIN_MULT_EXPLICIT=1; shift 2 ;;
    --sem-change-factor-max-mult) PP_SEM_CHANGE_FACTOR_MAX_MULT="$2"; SEM_MAX_MULT_EXPLICIT=1; shift 2 ;;
    --sem-max-relabel-step-frac) PP_SEM_MAX_RELABEL_STEP_FRAC="$2"; shift 2 ;;
    --sem-force-param-update-flip-frac) PP_SEM_FORCE_PARAM_UPDATE_FLIP_FRAC="$2"; shift 2 ;;
    --sem-monotone-complete-ll) PP_SEM_MONOTONE_COMPLETE_LL="$2"; shift 2 ;;
    --sem-start-from-c) PP_SEM_START_FROM_C="$2"; shift 2 ;;
    --sem-biv-n-threads) PP_SEM_BIV_N_THREADS="$2"; shift 2 ;;
    --sem-single-flip-from-iter) PP_SEM_SINGLE_FLIP_FROM_ITER="$2"; shift 2 ;;
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
    --kde-bw-method) PP_KDE_BW_METHOD="$2"; shift 2 ;;
    --kde-bw-sens-km) PP_KDE_BW_SENS_KM="$2"; shift 2 ;;
    --primary-partition) PP_PRIMARY_PARTITION="$2"; shift 2 ;;
    --run-partition-sensitivity) PP_RUN_PARTITION_SENSITIVITY="$2"; shift 2 ;;
    --boot-outer-cores) PP_BOOT_OUTER_CORES="$2"; BOOT_OUTER_CORES_EXPLICIT=1; shift 2 ;;
    --run-sensitivity) PP_RUN_SENSITIVITY="$2"; RUN_SENS_EXPLICIT=1; shift 2 ;;
    --run-fit-variability) PP_RUN_FIT_VARIABILITY="$2"; RUN_FIT_VARIABILITY_EXPLICIT=1; shift 2 ;;
    --fit-variability-reps) PP_FIT_VARIABILITY_REPS="$2"; FIT_VARIABILITY_REPS_EXPLICIT=1; shift 2 ;;
    --fit-variability-cores) PP_FIT_VARIABILITY_CORES="$2"; FIT_VARIABILITY_CORES_EXPLICIT=1; shift 2 ;;
    --fit-variability-patch-file) PP_FIT_VARIABILITY_PATCH_FILE="$2"; shift 2 ;;
    --fit-variability-only) PP_FIT_VARIABILITY_ONLY=1; shift ;;
    --bootstrap-patch-file) PP_BOOTSTRAP_PATCH_FILE="$2"; shift 2 ;;
    --bootstrap-only) PP_BOOTSTRAP_ONLY=1; shift ;;
    --ate-n-sims) PP_ATE_N_SIMS="$2"; ATE_N_SIMS_EXPLICIT=1; shift 2 ;;
    --skip-control-snapshots) PP_SKIP_CONTROL_SNAPSHOTS="$2"; shift 2 ;;
    --ate-bivariate) PP_ATE_BIVARIATE="$2"; shift 2 ;;
    --ate-contrast) PP_ATE_CONTRAST="$2"; shift 2 ;;
    --ate-scenario) PP_ATE_SCENARIO="$2"; shift 2 ;;
    --seed) PP_SEED="$2"; shift 2 ;;
    --setup-test) PP_SETUP_TEST=1; SETUP_TEST_EXPLICIT=1; shift ;;
    --mem) PP_MEM="$2"; MEM_EXPLICIT=1; shift 2 ;;
    --time) PP_TIME="$2"; shift 2 ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

if [ -n "$PP_MODE" ] && [ -z "${SLURM_JOB_ID:-}" ]; then
  mode_norm="$(echo "$PP_MODE" | tr '[:upper:]' '[:lower:]')"
  apply_paper_science_suite() {
    # Publication science: county primary, t_trunc=250, gated D from that C,
    # A–J (not C/D-only), scott-iso, bandwidth + partition C/D, no t-trunc grid.
    if [ "$SETUP_TEST_EXPLICIT" -ne 1 ]; then PP_SETUP_TEST=0; fi
    if [ "$CORES_EXPLICIT" -ne 1 ]; then PP_CORES=64; fi
    if [ "$SEM_N_ITER_EXPLICIT" -ne 1 ]; then PP_SEM_N_ITER=1; fi
    if [ "$SEM_N_LABELLINGS_EXPLICIT" -ne 1 ]; then PP_SEM_N_LABELLINGS=0; fi
    if [ "$SEM_OUTER_MAXIT_EXPLICIT" -ne 1 ]; then PP_SEM_OUTER_MAXIT=5000; fi
    if [ "$SEM_OUTER_MAXIT_BIV_EXPLICIT" -ne 1 ]; then PP_SEM_OUTER_MAXIT_BIV=5000; fi
    if [ "$SEM_OPTIM_METHOD_EXPLICIT" -ne 1 ]; then PP_SEM_OPTIM_METHOD="max"; fi
    if [ "$KDE_VARIANT_MODE_EXPLICIT" -ne 1 ]; then PP_KDE_VARIANT_MODE="single"; fi
    if [ "$ATE_N_SIMS_EXPLICIT" -ne 1 ]; then PP_ATE_N_SIMS=500; fi
    if [ "$MEM_EXPLICIT" -ne 1 ]; then PP_MEM=200G; fi
    if [ "$BOOT_TARGETS_EXPLICIT" -ne 1 ]; then PP_BOOT_TARGETS="C,D"; fi
    if [ "$RUN_FIT_VARIABILITY_EXPLICIT" -ne 1 ]; then PP_RUN_FIT_VARIABILITY=0; fi
    PP_RUN_SENSITIVITY=1
    RUN_SENS_EXPLICIT=1
    PP_RUN_PARTITION_SENSITIVITY=1
    PP_RUN_T_TRUNC_SENSITIVITY=0
    if [ -z "$PP_KDE_BW_METHOD" ]; then PP_KDE_BW_METHOD="scott-iso"; fi
    if [ -z "$PP_PRIMARY_PARTITION" ]; then PP_PRIMARY_PARTITION="county"; fi
    if [ "${SEM_T_TRUNC_EXPLICIT:-0}" -ne 1 ]; then PP_SEM_T_TRUNC_DAYS=250; fi
    if [ -z "$PP_ATE_BIVARIATE" ]; then PP_ATE_BIVARIATE=true; fi
    if [ -z "$PP_ATE_CONTRAST" ]; then PP_ATE_CONTRAST=all_or_nothing; fi
    PP_SEM_MONOTONE_COMPLETE_LL=1
    PP_SEM_START_FROM_C=1
    PP_SEM_BIV_N_THREADS=1
    if [ -z "${PP_SEM_SINGLE_FLIP_FROM_ITER:-}" ]; then
      PP_SEM_SINGLE_FLIP_FROM_ITER=1001
    fi
    PP_SKIP_FULL_REPORT=1
    PP_CD_ONLY=0
    PP_UNIV_KDE_ONLY=0
    if [ -z "${PP_SKIP_CONTROL_SNAPSHOTS}" ]; then PP_SKIP_CONTROL_SNAPSHOTS=1; fi
  }
  case "$mode_norm" in
    paper|paper-final|publication)
      apply_paper_science_suite
      if [ "$SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SEM_INNER=1250; fi
      if [ "$SENS_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SENS_SEM_INNER=1250; fi
      if [ "$BOOT_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_BOOT_SEM_INNER=1000; fi
      if [ "$BOOT_REPS_EXPLICIT" -ne 1 ]; then PP_BOOT_REPS=256; fi
      if [ "$BOOT_REFIT_SCOPE_EXPLICIT" -ne 1 ]; then PP_BOOT_REFIT_SCOPE="full"; fi
      if [ "$BOOT_OUTER_CORES_EXPLICIT" -ne 1 ]; then PP_BOOT_OUTER_CORES=64; fi
      if [ -z "${PP_TIME:-}" ] || [ "$PP_TIME" = "72:00:00" ]; then PP_TIME="24:00:00"; fi
      ;;
    paper-quick|paper_quick|preview)
      apply_paper_science_suite
      if [ "$SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SEM_INNER=1250; fi
      if [ "$SENS_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SENS_SEM_INNER=1250; fi
      if [ "$BOOT_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_BOOT_SEM_INNER=1000; fi
      if [ "$BOOT_REPS_EXPLICIT" -ne 1 ]; then PP_BOOT_REPS=0; fi
      if [ "$BOOT_REFIT_SCOPE_EXPLICIT" -ne 1 ]; then PP_BOOT_REFIT_SCOPE="none"; fi
      if [ -z "${PP_TIME:-}" ] || [ "$PP_TIME" = "72:00:00" ]; then PP_TIME="12:00:00"; fi
      ;;
    quick)
      if [ "$SETUP_TEST_EXPLICIT" -ne 1 ]; then PP_SETUP_TEST=0; fi
      if [ "$CORES_EXPLICIT" -ne 1 ]; then PP_CORES=32; fi
      # Quick profile default: disable bootstrap to minimize turnaround.
      if [ "$BOOT_REPS_EXPLICIT" -ne 1 ]; then PP_BOOT_REPS=0; fi
      if [ "$SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SEM_INNER=200; fi
      if [ "$SEM_N_ITER_EXPLICIT" -ne 1 ]; then PP_SEM_N_ITER=1; fi
      if [ "$SEM_N_LABELLINGS_EXPLICIT" -ne 1 ]; then PP_SEM_N_LABELLINGS=20; fi
      if [ "$SEM_OUTER_MAXIT_EXPLICIT" -ne 1 ]; then PP_SEM_OUTER_MAXIT=120; fi
      if [ "$SEM_OUTER_MAXIT_BIV_EXPLICIT" -ne 1 ]; then PP_SEM_OUTER_MAXIT_BIV=5000; fi
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
      if [ "$SEM_OUTER_MAXIT_BIV_EXPLICIT" -ne 1 ]; then PP_SEM_OUTER_MAXIT_BIV=5000; fi
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
      if [ "$SEM_OUTER_MAXIT_EXPLICIT" -ne 1 ]; then PP_SEM_OUTER_MAXIT=5000; fi
      if [ "$SEM_OUTER_MAXIT_BIV_EXPLICIT" -ne 1 ]; then PP_SEM_OUTER_MAXIT_BIV=5000; fi
      if [ "$SENS_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SENS_SEM_INNER=1000; fi
      if [ "$BOOT_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_BOOT_SEM_INNER=1000; fi
      # Full profile default: full parametric bootstrap (simulate + refit + ATE).
      if [ "$BOOT_REFIT_SCOPE_EXPLICIT" -ne 1 ]; then PP_BOOT_REFIT_SCOPE="full"; fi
      if [ "$BOOT_TARGETS_EXPLICIT" -ne 1 ]; then PP_BOOT_TARGETS="C,D"; fi
      if [ "$KDE_VARIANT_MODE_EXPLICIT" -ne 1 ]; then PP_KDE_VARIANT_MODE="triple"; fi
      # One mostly-single-threaded bootstrap worker per allocated core.
      if [ "$BOOT_OUTER_CORES_EXPLICIT" -ne 1 ]; then PP_BOOT_OUTER_CORES="$PP_CORES"; fi
      if [ "$RUN_SENS_EXPLICIT" -ne 1 ]; then PP_RUN_SENSITIVITY=0; fi
      if [ "$ATE_N_SIMS_EXPLICIT" -ne 1 ]; then PP_ATE_N_SIMS=500; fi
      if [ "$MEM_EXPLICIT" -ne 1 ]; then PP_MEM=200G; fi
      ;;
    fitvar|fit-variability|variability)
      if [ "$SETUP_TEST_EXPLICIT" -ne 1 ]; then PP_SETUP_TEST=0; fi
      if [ "$CORES_EXPLICIT" -ne 1 ]; then PP_CORES=32; fi
      if [ "$BOOT_REPS_EXPLICIT" -ne 1 ]; then PP_BOOT_REPS=0; fi
      if [ "$RUN_SENS_EXPLICIT" -ne 1 ]; then PP_RUN_SENSITIVITY=0; fi
      if [ "$RUN_FIT_VARIABILITY_EXPLICIT" -ne 1 ]; then PP_RUN_FIT_VARIABILITY=1; fi
      if [ "$FIT_VARIABILITY_REPS_EXPLICIT" -ne 1 ]; then PP_FIT_VARIABILITY_REPS="$PP_CORES"; fi
      if [ "$FIT_VARIABILITY_CORES_EXPLICIT" -ne 1 ]; then PP_FIT_VARIABILITY_CORES="$PP_CORES"; fi
      if [ "$ATE_N_SIMS_EXPLICIT" -ne 1 ]; then PP_ATE_N_SIMS=200; fi
      if [ "$MEM_EXPLICIT" -ne 1 ]; then PP_MEM=160G; fi
      ;;
    fit-variability-only|fitvar-only|variability-only)
      PP_FIT_VARIABILITY_ONLY=1
      if [ "$SETUP_TEST_EXPLICIT" -ne 1 ]; then PP_SETUP_TEST=0; fi
      if [ "$CORES_EXPLICIT" -ne 1 ]; then PP_CORES=32; fi
      if [ "$BOOT_REPS_EXPLICIT" -ne 1 ]; then PP_BOOT_REPS=0; fi
      if [ "$RUN_SENS_EXPLICIT" -ne 1 ]; then PP_RUN_SENSITIVITY=0; fi
      if [ "$RUN_FIT_VARIABILITY_EXPLICIT" -ne 1 ]; then PP_RUN_FIT_VARIABILITY=1; fi
      if [ "$FIT_VARIABILITY_REPS_EXPLICIT" -ne 1 ]; then PP_FIT_VARIABILITY_REPS="$PP_CORES"; fi
      if [ "$FIT_VARIABILITY_CORES_EXPLICIT" -ne 1 ]; then PP_FIT_VARIABILITY_CORES="$PP_CORES"; fi
      if [ "$ATE_N_SIMS_EXPLICIT" -ne 1 ]; then PP_ATE_N_SIMS=200; fi
      if [ "$MEM_EXPLICIT" -ne 1 ]; then PP_MEM=160G; fi
      ;;
    bootstrap-only|boot-only)
      PP_BOOTSTRAP_ONLY=1
      if [ "$SETUP_TEST_EXPLICIT" -ne 1 ]; then PP_SETUP_TEST=0; fi
      if [ "$CORES_EXPLICIT" -ne 1 ]; then PP_CORES=32; fi
      if [ "$BOOT_REPS_EXPLICIT" -ne 1 ]; then PP_BOOT_REPS=64; fi
      if [ "$RUN_SENS_EXPLICIT" -ne 1 ]; then PP_RUN_SENSITIVITY=0; fi
      if [ "$RUN_FIT_VARIABILITY_EXPLICIT" -ne 1 ]; then PP_RUN_FIT_VARIABILITY=0; fi
      if [ "$BOOT_SEM_INNER_EXPLICIT" -ne 1 ]; then PP_BOOT_SEM_INNER=2000; fi
      if [ "$BOOT_REFIT_SCOPE_EXPLICIT" -ne 1 ]; then PP_BOOT_REFIT_SCOPE="full"; fi
      if [ "$BOOT_TARGETS_EXPLICIT" -ne 1 ]; then PP_BOOT_TARGETS="C,D"; fi
      if [ "$KDE_VARIANT_MODE_EXPLICIT" -ne 1 ]; then PP_KDE_VARIANT_MODE="single"; fi
      # One mostly-single-threaded bootstrap worker per allocated core.
      if [ "$BOOT_OUTER_CORES_EXPLICIT" -ne 1 ]; then PP_BOOT_OUTER_CORES="$PP_CORES"; fi
      if [ "$ATE_N_SIMS_EXPLICIT" -ne 1 ]; then PP_ATE_N_SIMS=500; fi
      if [ "$MEM_EXPLICIT" -ne 1 ]; then PP_MEM=200G; fi
      if [ -z "${PP_TIME:-}" ] || [ "$PP_TIME" = "72:00:00" ]; then PP_TIME="24:00:00"; fi
      ;;
    t-trunc-sens-only|trunc-sens-only|t_trunc_sens_only)
      PP_T_TRUNC_SENS_ONLY=1
      PP_RUN_T_TRUNC_SENSITIVITY=1
      if [ "$SETUP_TEST_EXPLICIT" -ne 1 ]; then PP_SETUP_TEST=0; fi
      if [ "$CORES_EXPLICIT" -ne 1 ]; then PP_CORES=16; fi
      if [ "$BOOT_REPS_EXPLICIT" -ne 1 ]; then PP_BOOT_REPS=0; fi
      if [ "$RUN_SENS_EXPLICIT" -ne 1 ]; then PP_RUN_SENSITIVITY=0; fi
      if [ "$RUN_FIT_VARIABILITY_EXPLICIT" -ne 1 ]; then PP_RUN_FIT_VARIABILITY=0; fi
      if [ "$SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SEM_INNER=2000; fi
      if [ "$SEM_N_ITER_EXPLICIT" -ne 1 ]; then PP_SEM_N_ITER=1; fi
      if [ "$SEM_OUTER_MAXIT_BIV_EXPLICIT" -ne 1 ]; then PP_SEM_OUTER_MAXIT_BIV=5000; fi
      if [ "$KDE_VARIANT_MODE_EXPLICIT" -ne 1 ]; then PP_KDE_VARIANT_MODE="single"; fi
      if [ "$ATE_N_SIMS_EXPLICIT" -ne 1 ]; then PP_ATE_N_SIMS=500; fi
      if [ "$MEM_EXPLICIT" -ne 1 ]; then PP_MEM=120G; fi
      if [ -z "${PP_TIME:-}" ] || [ "$PP_TIME" = "72:00:00" ]; then PP_TIME="06:00:00"; fi
      ;;
    smoke-sem-d|smoke-d-trunc|smoke_sem_d)
      # Multi-seed Fit D SEM at fixed t_trunc (no ATE); reuses trunc-only skip path.
      PP_T_TRUNC_SENS_ONLY=1
      PP_RUN_T_TRUNC_SENSITIVITY=0
      if [ "${PP_SMOKE_SEM_D_SEEDS:-0}" = "0" ]; then PP_SMOKE_SEM_D_SEEDS=8; fi
      if [ -z "${PP_SMOKE_SEM_D_TRUNC:-}" ]; then PP_SMOKE_SEM_D_TRUNC=3; fi
      if [ "$SETUP_TEST_EXPLICIT" -ne 1 ]; then PP_SETUP_TEST=0; fi
      if [ "$CORES_EXPLICIT" -ne 1 ]; then PP_CORES="$PP_SMOKE_SEM_D_SEEDS"; fi
      if [ "$BOOT_REPS_EXPLICIT" -ne 1 ]; then PP_BOOT_REPS=0; fi
      if [ "$RUN_SENS_EXPLICIT" -ne 1 ]; then PP_RUN_SENSITIVITY=0; fi
      if [ "$RUN_FIT_VARIABILITY_EXPLICIT" -ne 1 ]; then PP_RUN_FIT_VARIABILITY=0; fi
      if [ "$SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SEM_INNER=1000; fi
      if [ "$SEM_N_ITER_EXPLICIT" -ne 1 ]; then PP_SEM_N_ITER=1; fi
      if [ "$SEM_OUTER_MAXIT_BIV_EXPLICIT" -ne 1 ]; then PP_SEM_OUTER_MAXIT_BIV=5000; fi
      if [ "$KDE_VARIANT_MODE_EXPLICIT" -ne 1 ]; then PP_KDE_VARIANT_MODE="single"; fi
      if [ "$ATE_N_SIMS_EXPLICIT" -ne 1 ]; then PP_ATE_N_SIMS=1; fi
      if [ "$MEM_EXPLICIT" -ne 1 ]; then PP_MEM=80G; fi
      if [ -z "${PP_TIME:-}" ] || [ "$PP_TIME" = "72:00:00" ]; then PP_TIME="02:00:00"; fi
      ;;
    cd-primary|cd|cd_primary|cd-slim)
      # C/D-focused refresh (still runs A/B + univ unless --mode cd-only):
      # large fixed t_trunc, no bootstrap/sens, no full HTML overwrite.
      if [ "$SETUP_TEST_EXPLICIT" -ne 1 ]; then PP_SETUP_TEST=0; fi
      if [ "$CORES_EXPLICIT" -ne 1 ]; then PP_CORES=64; fi
      if [ "$BOOT_REPS_EXPLICIT" -ne 1 ]; then PP_BOOT_REPS=0; fi
      if [ "$SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SEM_INNER=2000; fi
      if [ "$SEM_N_ITER_EXPLICIT" -ne 1 ]; then PP_SEM_N_ITER=1; fi
      if [ "$SEM_N_LABELLINGS_EXPLICIT" -ne 1 ]; then PP_SEM_N_LABELLINGS=20; fi
      if [ "$SEM_OUTER_MAXIT_EXPLICIT" -ne 1 ]; then PP_SEM_OUTER_MAXIT=5000; fi
      if [ "$SEM_OUTER_MAXIT_BIV_EXPLICIT" -ne 1 ]; then PP_SEM_OUTER_MAXIT_BIV=5000; fi
      if [ -z "${PP_SEM_T_TRUNC_DAYS}" ]; then PP_SEM_T_TRUNC_DAYS=90; fi
      PP_RUN_T_TRUNC_SENSITIVITY=0
      if [ "$RUN_SENS_EXPLICIT" -ne 1 ]; then PP_RUN_SENSITIVITY=0; fi
      if [ "$RUN_FIT_VARIABILITY_EXPLICIT" -ne 1 ]; then PP_RUN_FIT_VARIABILITY=0; fi
      if [ "$BOOT_REFIT_SCOPE_EXPLICIT" -ne 1 ]; then PP_BOOT_REFIT_SCOPE="none"; fi
      if [ "$BOOT_TARGETS_EXPLICIT" -ne 1 ]; then PP_BOOT_TARGETS="C,D"; fi
      if [ "$KDE_VARIANT_MODE_EXPLICIT" -ne 1 ]; then PP_KDE_VARIANT_MODE="single"; fi
      if [ "$ATE_N_SIMS_EXPLICIT" -ne 1 ]; then PP_ATE_N_SIMS=500; fi
      if [ "$MEM_EXPLICIT" -ne 1 ]; then PP_MEM=200G; fi
      if [ -z "${PP_TIME:-}" ] || [ "$PP_TIME" = "72:00:00" ]; then PP_TIME="06:00:00"; fi
      # Keep the full diagnostic HTML untouched; slim report is rendered locally.
      PP_SKIP_FULL_REPORT=1
      ;;
    cd-only|cd_only|just-cd|just_cd)
      # Strict C/D-only: skip homogeneous A/B and univariate fits/ATEs.
      if [ "$SETUP_TEST_EXPLICIT" -ne 1 ]; then PP_SETUP_TEST=0; fi
      if [ "$CORES_EXPLICIT" -ne 1 ]; then PP_CORES=64; fi
      if [ "$BOOT_REPS_EXPLICIT" -ne 1 ]; then PP_BOOT_REPS=0; fi
      if [ "$SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SEM_INNER=2000; fi
      if [ "$SEM_N_ITER_EXPLICIT" -ne 1 ]; then PP_SEM_N_ITER=1; fi
      if [ "$SEM_N_LABELLINGS_EXPLICIT" -ne 1 ]; then PP_SEM_N_LABELLINGS=20; fi
      if [ "$SEM_OUTER_MAXIT_EXPLICIT" -ne 1 ]; then PP_SEM_OUTER_MAXIT=5000; fi
      if [ "$SEM_OUTER_MAXIT_BIV_EXPLICIT" -ne 1 ]; then PP_SEM_OUTER_MAXIT_BIV=5000; fi
      if [ -z "${PP_SEM_T_TRUNC_DAYS}" ]; then PP_SEM_T_TRUNC_DAYS=90; fi
      PP_RUN_T_TRUNC_SENSITIVITY=0
      if [ "$RUN_SENS_EXPLICIT" -ne 1 ]; then PP_RUN_SENSITIVITY=0; fi
      if [ "$RUN_FIT_VARIABILITY_EXPLICIT" -ne 1 ]; then PP_RUN_FIT_VARIABILITY=0; fi
      if [ "$BOOT_REFIT_SCOPE_EXPLICIT" -ne 1 ]; then PP_BOOT_REFIT_SCOPE="none"; fi
      if [ "$BOOT_TARGETS_EXPLICIT" -ne 1 ]; then PP_BOOT_TARGETS="C,D"; fi
      if [ "$KDE_VARIANT_MODE_EXPLICIT" -ne 1 ]; then PP_KDE_VARIANT_MODE="single"; fi
      if [ "$ATE_N_SIMS_EXPLICIT" -ne 1 ]; then PP_ATE_N_SIMS=500; fi
      if [ "$MEM_EXPLICIT" -ne 1 ]; then PP_MEM=160G; fi
      if [ -z "${PP_TIME:-}" ] || [ "$PP_TIME" = "72:00:00" ]; then PP_TIME="04:00:00"; fi
      PP_SKIP_FULL_REPORT=1
      PP_CD_ONLY=1
      if [ -z "${PP_SKIP_CONTROL_SNAPSHOTS}" ]; then PP_SKIP_CONTROL_SNAPSHOTS=1; fi
      ;;
    univ-kde-only|univ_kde_only|just-gh|just_gh|gh-only|gh_only)
      # Lean univariate+KDE (public G/H): skip A–F and homogeneous I/J.
      if [ "$SETUP_TEST_EXPLICIT" -ne 1 ]; then PP_SETUP_TEST=0; fi
      if [ "$CORES_EXPLICIT" -ne 1 ]; then PP_CORES=64; fi
      if [ "$BOOT_REPS_EXPLICIT" -ne 1 ]; then PP_BOOT_REPS=0; fi
      if [ "$SEM_INNER_EXPLICIT" -ne 1 ]; then PP_SEM_INNER=1500; fi
      if [ "$SEM_N_ITER_EXPLICIT" -ne 1 ]; then PP_SEM_N_ITER=1; fi
      if [ "$SEM_N_LABELLINGS_EXPLICIT" -ne 1 ]; then PP_SEM_N_LABELLINGS=20; fi
      if [ "$SEM_OUTER_MAXIT_EXPLICIT" -ne 1 ]; then PP_SEM_OUTER_MAXIT=5000; fi
      if [ "$SEM_OUTER_MAXIT_BIV_EXPLICIT" -ne 1 ]; then PP_SEM_OUTER_MAXIT_BIV=5000; fi
      if [ -z "${PP_SEM_T_TRUNC_DAYS}" ]; then PP_SEM_T_TRUNC_DAYS=90; fi
      PP_RUN_T_TRUNC_SENSITIVITY=0
      if [ "$RUN_SENS_EXPLICIT" -ne 1 ]; then PP_RUN_SENSITIVITY=0; fi
      if [ "$RUN_FIT_VARIABILITY_EXPLICIT" -ne 1 ]; then PP_RUN_FIT_VARIABILITY=0; fi
      if [ "$BOOT_REFIT_SCOPE_EXPLICIT" -ne 1 ]; then PP_BOOT_REFIT_SCOPE="none"; fi
      if [ "$KDE_VARIANT_MODE_EXPLICIT" -ne 1 ]; then PP_KDE_VARIANT_MODE="single"; fi
      if [ "$ATE_N_SIMS_EXPLICIT" -ne 1 ]; then PP_ATE_N_SIMS=20; fi
      if [ "$MEM_EXPLICIT" -ne 1 ]; then PP_MEM=96G; fi
      if [ -z "${PP_TIME:-}" ] || [ "$PP_TIME" = "72:00:00" ]; then PP_TIME="01:00:00"; fi
      PP_SKIP_FULL_REPORT=1
      PP_UNIV_KDE_ONLY=1
      if [ -z "${PP_SKIP_CONTROL_SNAPSHOTS}" ]; then PP_SKIP_CONTROL_SNAPSHOTS=1; fi
      ;;
    *)
      echo "Unknown --mode '$PP_MODE' (expected: paper | paper-quick | test | quick | full | fit-variability | fit-variability-only | bootstrap-only | t-trunc-sens-only | smoke-sem-d | cd-primary | cd-only | univ-kde-only)"
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
  PP_SEM_OUTER_MAXIT=5000
fi
if [ -z "$PP_SEM_OUTER_MAXIT_BIV" ]; then
  PP_SEM_OUTER_MAXIT_BIV=5000
fi
if [ -z "$PP_BOOT_SEM_INNER" ]; then
  PP_BOOT_SEM_INNER="$PP_SEM_INNER"
fi
if [ -z "$PP_BOOT_REFIT_SCOPE" ]; then
  # Keep "none" when bootstrap is off; require full refits when reps > 0.
  if [ "${PP_BOOT_REPS:-0}" -gt 0 ]; then
    PP_BOOT_REFIT_SCOPE="full"
  else
    PP_BOOT_REFIT_SCOPE="none"
  fi
fi
boot_refit_norm="$(echo "$PP_BOOT_REFIT_SCOPE" | tr '[:upper:]' '[:lower:]')"
if [ "$boot_refit_norm" = "none" ] || [ "$boot_refit_norm" = "partial" ] || [ "$boot_refit_norm" = "full" ]; then
  PP_BOOT_REFIT_SCOPE="$boot_refit_norm"
else
  echo "Invalid --boot-refit-scope '$PP_BOOT_REFIT_SCOPE' (expected: none | partial | full)"
  exit 1
fi
# Parametric bootstrap requires per-replicate refits; upgrade none -> full.
if [ "${PP_BOOT_REPS:-0}" -gt 0 ] && [ "$PP_BOOT_REFIT_SCOPE" = "none" ]; then
  echo "NOTE: boot_reps=${PP_BOOT_REPS} with boot-refit-scope=none is invalid; using full."
  PP_BOOT_REFIT_SCOPE="full"
fi
if [ -z "$PP_BOOT_OUTER_CORES" ]; then
  PP_BOOT_OUTER_CORES="$PP_CORES"
fi
# Keep sensitivity/bootstrap targets aligned on the all-free pair.
boot_targets_norm="$(echo "$PP_BOOT_TARGETS" | tr '[:lower:]' '[:upper:]' | tr ';| ' ',,,')"
boot_targets_norm="$(echo "$boot_targets_norm" | sed 's/,,*/,/g; s/^,//; s/,$//')"
has_e=0
has_f=0
IFS=',' read -r -a _boot_target_tokens <<< "$boot_targets_norm"
for tok in "${_boot_target_tokens[@]}"; do
  [ "$tok" = "C" ] && has_e=1
  [ "$tok" = "D" ] && has_f=1
done
if [ "$has_e" -ne 1 ] || [ "$has_f" -ne 1 ]; then
  echo "Invalid --boot-targets '$PP_BOOT_TARGETS': sensitivity/bootstrap alignment requires C,D"
  exit 1
fi
PP_BOOT_TARGETS="C,D"
# Default policy (manual mode): one bootstrap replicate per available outer bootstrap core.
# Do not override mode-specific defaults already assigned above.
if [ "$BOOT_REPS_EXPLICIT" -ne 1 ] && [ -z "${PP_BOOT_REPS:-}" ]; then
  PP_BOOT_REPS="$PP_BOOT_OUTER_CORES"
fi
if [ -z "$PP_KDE_VARIANT_MODE" ]; then
  PP_KDE_VARIANT_MODE="triple"
fi
if [ -z "$PP_KDE_BW_METHOD" ]; then
  PP_KDE_BW_METHOD="scott-iso"
fi
case "$(echo "$PP_KDE_BW_METHOD" | tr '[:upper:]' '[:lower:]')" in
  scott|bw.scott|scott_aniso) PP_KDE_BW_METHOD="scott" ;;
  scott-iso|scott_iso|bw.scott.iso) PP_KDE_BW_METHOD="scott-iso" ;;
  diggle|bw.diggle) PP_KDE_BW_METHOD="diggle" ;;
  diggle2|2diggle|"2*diggle"|digglex2) PP_KDE_BW_METHOD="diggle2" ;;
  *)
    echo "Invalid --kde-bw-method '$PP_KDE_BW_METHOD' (expected: scott-iso | scott | diggle | diggle2)"
    exit 1
    ;;
esac
if [ -z "$PP_PRIMARY_PARTITION" ]; then
  PP_PRIMARY_PARTITION="grid_1.0R"
fi
case "$(echo "$PP_PRIMARY_PARTITION" | tr '[:upper:]' '[:lower:]' | tr -d '_-' | tr -d '[:space:]')" in
  county|counties|admin) PP_PRIMARY_PARTITION="county" ;;
  gridcoarse|coarse|quickgrid) PP_PRIMARY_PARTITION="grid_coarse" ;;
  grid10r|grid1r|grid1.0r|1r|1.0r|1) PP_PRIMARY_PARTITION="grid_1.0R" ;;
  *)
    echo "Invalid --primary-partition '$PP_PRIMARY_PARTITION' (expected: grid_1.0R | county)"
    exit 1
    ;;
esac
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
if [ -z "$PP_FIT_VARIABILITY_REPS" ]; then
  PP_FIT_VARIABILITY_REPS="$PP_CORES"
fi
if [ -z "$PP_FIT_VARIABILITY_CORES" ]; then
  PP_FIT_VARIABILITY_CORES="$PP_CORES"
fi
if [ "${PP_BOOTSTRAP_ONLY:-0}" = "1" ] || [ "${PP_BOOTSTRAP_ONLY:-0}" = "true" ]; then
  if [ -z "${PP_BOOTSTRAP_PATCH_FILE:-}" ]; then
    echo "ERROR: --bootstrap-only (or mode bootstrap-only) requires --bootstrap-patch-file"
    exit 1
  fi
fi
if [ "${PP_FIT_VARIABILITY_ONLY:-0}" = "1" ] || [ "${PP_FIT_VARIABILITY_ONLY:-0}" = "true" ]; then
  if [ -z "${PP_FIT_VARIABILITY_PATCH_FILE:-}" ]; then
    echo "ERROR: --fit-variability-only (or mode fit-variability-only) requires --fit-variability-patch-file"
    exit 1
  fi
fi
if [ "${PP_T_TRUNC_SENS_ONLY:-0}" = "1" ] || [ "${PP_T_TRUNC_SENS_ONLY:-0}" = "true" ]; then
  if [ -z "${PP_T_TRUNC_SENS_PATCH_FILE:-}" ]; then
    echo "ERROR: --t-trunc-sens-only / --mode smoke-sem-d requires --t-trunc-sens-patch-file"
    exit 1
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
source "$PKG_ROOT/inst/include/output_root.sh"
# shellcheck source=../include/nesi_mem.sh
source "$PKG_ROOT/inst/include/nesi_mem.sh"

# Cap bootstrap outer workers so they fit in PP_MEM. 2 GB/worker plus 8 GB
# for the master/OS: 64G → 28 workers (covers ~1.9 GB full-mode spikes;
# bootstrap-only jobs measured ~0.85 GB). Set PP_BOOT_RAM_CAP=0 to disable.
if [ "${PP_BOOT_RAM_CAP:-1}" != "0" ] && [ "${PP_BOOT_RAM_CAP:-1}" != "false" ]; then
  PP_BOOT_WORKER_GB="${PP_BOOT_WORKER_GB:-2}"
  PP_BOOT_MEM_RESERVE_GB="${PP_BOOT_MEM_RESERVE_GB:-8}"
  BOOT_RAM_CAP="$(pp_boot_outer_from_mem "$PP_CORES" "$PP_MEM" "$PP_BOOT_WORKER_GB" "$PP_BOOT_MEM_RESERVE_GB")"
  if [ "$PP_BOOT_OUTER_CORES" -gt "$BOOT_RAM_CAP" ]; then
    echo "RAM-aware bootstrap outer cap: ${PP_BOOT_OUTER_CORES} -> ${BOOT_RAM_CAP} (${PP_MEM}, ${PP_BOOT_WORKER_GB}G/worker, ${PP_BOOT_MEM_RESERVE_GB}G reserve)"
    PP_BOOT_OUTER_CORES="$BOOT_RAM_CAP"
  fi
fi

# ----------------------------
# Submit mode
# ----------------------------
if [ -z "${SLURM_JOB_ID:-}" ]; then
  cd "$PKG_ROOT"
  # shellcheck source=../include/git_sync.sh
  source "$PKG_ROOT/inst/include/git_sync.sh"
  pp_git_sync_repo "$PKG_ROOT"
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  OUTPUT_DIR="$(pp_disentangle_output_path oklahoma)"
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
  export PP_SEM_OUTER_MAXIT PP_SEM_OUTER_MAXIT_BIV PP_SEM_T_TRUNC_DAYS PP_SEM_T_TRUNC_REL PP_T_TRUNC_SENS_DAYS PP_RUN_T_TRUNC_SENSITIVITY PP_SEM_TEMPORAL_WEIGHT
  export PP_SEM_WORKER_LOGS PP_SEM_WORKER_LOG_VERBOSE PP_SEM_WORKER_LOG_SPLIT PP_SEM_TIMING_VERBOSE PP_SEM_PROPOSAL_VERBOSE
  export PP_SIM_PROGRESS_EVERY PP_SENS_SEM_INNER PP_BOOT_SEM_INNER PP_BOOT_REFIT_SCOPE PP_BOOT_TARGETS PP_KDE_VARIANT_MODE PP_KDE_BW_METHOD PP_KDE_BW_SENS_KM PP_PRIMARY_PARTITION PP_RUN_PARTITION_SENSITIVITY
  export PP_BOOT_OUTER_CORES PP_ATE_N_SIMS PP_ATE_BIVARIATE PP_ATE_CONTRAST PP_ATE_SCENARIO
  export PP_RUN_SENSITIVITY PP_RUN_FIT_VARIABILITY PP_FIT_VARIABILITY_REPS PP_FIT_VARIABILITY_CORES PP_FIT_VARIABILITY_PATCH_FILE PP_FIT_VARIABILITY_ONLY PP_BOOTSTRAP_PATCH_FILE PP_BOOTSTRAP_ONLY PP_T_TRUNC_SENS_PATCH_FILE PP_T_TRUNC_SENS_ONLY PP_SMOKE_SEM_D_SEEDS PP_SMOKE_SEM_D_TRUNC PP_SKIP_FULL_REPORT PP_CD_ONLY PP_UNIV_KDE_ONLY PP_MEM PP_TIME PP_SEM_OPTIM_METHOD PP_SEM_SELECTION_TEMPERATURE
  export PP_SEM_CHANGE_FACTOR_MIN_MULT PP_SEM_CHANGE_FACTOR_MAX_MULT PP_SEM_MAX_RELABEL_STEP_FRAC PP_SEM_FORCE_PARAM_UPDATE_FLIP_FRAC
  export PP_SEM_MONOTONE_COMPLETE_LL PP_SEM_START_FROM_C PP_SEM_BIV_N_THREADS PP_SEM_SINGLE_FLIP_FROM_ITER
  export PP_RUN_SEM_PILOT PP_SEM_PILOT_INNER PP_SEM_PILOT_CORES PP_SEM_PILOT_MAX_COMBOS PP_SEM_PILOT_CHANGE_FACTORS
  export PP_SEM_PILOT_MIN_MULTS PP_SEM_PILOT_MAX_MULTS PP_SEM_PILOT_TEMPS PP_SETUP_TEST
  [ -n "${PP_R_GEO_MODULE:-}" ] && export PP_R_GEO_MODULE

  echo "Submitting Oklahoma job: mode=${PP_MODE:-manual} cores=$PP_CORES ate_bivariate=${PP_ATE_BIVARIATE:-default} ate_contrast=${PP_ATE_CONTRAST:-default} ate_scenario=${PP_ATE_SCENARIO:-} ate_n_sims=$PP_ATE_N_SIMS boot_reps=$PP_BOOT_REPS boot_refit_scope=$PP_BOOT_REFIT_SCOPE boot_outer_cores=$PP_BOOT_OUTER_CORES setup_test=$PP_SETUP_TEST skip_full_report=$PP_SKIP_FULL_REPORT"

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
  echo "SLURM out: PPDisentangle-output/oklahoma/${JOB_ID}_oklahoma_slurm.out"
  echo "SLURM err: PPDisentangle-output/oklahoma/${JOB_ID}_oklahoma_slurm.err"
  exit 0
fi

# ----------------------------
# Job mode
# ----------------------------
cd "$PKG_ROOT"
mkdir -p "$(pp_disentangle_output_path oklahoma)"

echo "=== PPDisentangle Oklahoma (NeSI) ==="
echo "Job: ${SLURM_JOB_ID} | $(date)"
echo "Node: $(hostname) | Partition: ${SLURM_JOB_PARTITION:-unknown}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-$PP_CORES}"
echo "boot_reps=$PP_BOOT_REPS sem_n_iter=$PP_SEM_N_ITER sem_outer_maxit=$PP_SEM_OUTER_MAXIT sem_outer_maxit_biv=$PP_SEM_OUTER_MAXIT_BIV sem_t_trunc_days=$PP_SEM_T_TRUNC_DAYS sem_t_trunc_rel=$PP_SEM_T_TRUNC_REL sem_temporal_weight=$PP_SEM_TEMPORAL_WEIGHT sem_warmstart_fixed=$PP_SEM_WARMSTART_FIXED sem_optim=$PP_SEM_OPTIM_METHOD sem_temp=$PP_SEM_SELECTION_TEMPERATURE sem_cf_min=$PP_SEM_CHANGE_FACTOR_MIN_MULT sem_cf_max=$PP_SEM_CHANGE_FACTOR_MAX_MULT run_sem_pilot=$PP_RUN_SEM_PILOT sem_pilot_inner=$PP_SEM_PILOT_INNER sem_pilot_cores=${PP_SEM_PILOT_CORES:-auto} sem_pilot_max_combos=$PP_SEM_PILOT_MAX_COMBOS sem_worker_logs=$PP_SEM_WORKER_LOGS sem_worker_log_verbose=$PP_SEM_WORKER_LOG_VERBOSE sem_worker_log_split=$PP_SEM_WORKER_LOG_SPLIT sem_timing_verbose=$PP_SEM_TIMING_VERBOSE sem_proposal_verbose=$PP_SEM_PROPOSAL_VERBOSE sim_progress_every=$PP_SIM_PROGRESS_EVERY sem_inner=$PP_SEM_INNER sens_inner=$PP_SENS_SEM_INNER boot_inner=$PP_BOOT_SEM_INNER boot_refit_scope=$PP_BOOT_REFIT_SCOPE targets=$PP_BOOT_TARGETS kde_variants=$PP_KDE_VARIANT_MODE kde_bw=$PP_KDE_BW_METHOD primary_partition=$PP_PRIMARY_PARTITION kde_bw_sens_km=${PP_KDE_BW_SENS_KM:-} run_partition_sens=${PP_RUN_PARTITION_SENSITIVITY:-auto} run_fit_variability=$PP_RUN_FIT_VARIABILITY fit_variability_reps=$PP_FIT_VARIABILITY_REPS fit_variability_cores=$PP_FIT_VARIABILITY_CORES fit_variability_only=$PP_FIT_VARIABILITY_ONLY"
echo "setup_test=$PP_SETUP_TEST mode=${PP_MODE:-manual}"
echo "seed=$PP_SEED (fit jobs RNG de-correlated by model; bootstrap RNG de-correlated by replicate)"
echo "ENV CHECK: PP_SEM_INNER=$PP_SEM_INNER | PP_SENS_SEM_INNER=$PP_SENS_SEM_INNER | PP_BOOT_SEM_INNER=$PP_BOOT_SEM_INNER | PP_SEM_MONOTONE_COMPLETE_LL=$PP_SEM_MONOTONE_COMPLETE_LL | PP_SEM_START_FROM_C=$PP_SEM_START_FROM_C | PP_SEM_BIV_N_THREADS=$PP_SEM_BIV_N_THREADS | PP_SEM_SINGLE_FLIP_FROM_ITER=${PP_SEM_SINGLE_FLIP_FROM_ITER:-}"
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

# Skip a full --preclean rebuild when this git SHA is already installed
# for this R version. Dirty trees and PP_FORCE_INSTALL=1 always rebuild.
PP_FORCE_INSTALL="${PP_FORCE_INSTALL:-0}"
PKG_SHA=""
PKG_DIRTY=0
if git -C "$PKG_ROOT" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  PKG_SHA="$(git -C "$PKG_ROOT" rev-parse HEAD 2>/dev/null || true)"
  if ! git -C "$PKG_ROOT" diff --quiet --ignore-submodules HEAD 2>/dev/null ||
     ! git -C "$PKG_ROOT" diff --cached --quiet --ignore-submodules 2>/dev/null; then
    PKG_DIRTY=1
  fi
fi
R_VER_STAMP="$("$RSCRIPT_BIN" -e 'cat(R.version$major, ".", R.version$minor, sep="")' 2>/dev/null || echo unknown)"
INSTALL_STAMP=""
if [ -n "$PKG_SHA" ] && [ "$PKG_DIRTY" -eq 0 ]; then
  INSTALL_STAMP="${SHARED_R_LIBS_USER}/.ppdis_pkg_install_${PKG_SHA}_${R_VER_STAMP}"
fi
NEED_INSTALL=1
if [ "$PP_FORCE_INSTALL" = "1" ] || [ "$PP_FORCE_INSTALL" = "true" ] || [ "$PP_FORCE_INSTALL" = "yes" ]; then
  echo "Forcing PPDisentangle rebuild (PP_FORCE_INSTALL=1)."
elif [ -n "$INSTALL_STAMP" ] && [ -f "$INSTALL_STAMP" ] && [ -d "${SHARED_R_LIBS_USER}/PPDisentangle" ]; then
  NEED_INSTALL=0
  echo "Skipping PPDisentangle install; SHA ${PKG_SHA} already built for R ${R_VER_STAMP}."
  echo "  stamp: $INSTALL_STAMP"
  echo "  set PP_FORCE_INSTALL=1 to rebuild."
elif [ "$PKG_DIRTY" -eq 1 ]; then
  echo "Installing PPDisentangle from source (dirty tree; SHA ${PKG_SHA:-unknown})..."
else
  echo "Installing PPDisentangle from source (fresh install)..."
fi
if [ "$NEED_INSTALL" -eq 1 ]; then
  wait_for_pp_lock_clear "$PP_LOCK_DIR"
  cleanup_pp_lock_if_safe "$PP_LOCK_DIR"
  if [ -n "$INSTALL_STAMP" ] && [ -f "$INSTALL_STAMP" ] && [ -d "${SHARED_R_LIBS_USER}/PPDisentangle" ]; then
    echo "Skipping PPDisentangle install after lock wait; SHA ${PKG_SHA} already built."
  else
    "$R_BIN" CMD INSTALL --preclean --no-multiarch "$PKG_ROOT"
    rm -f "${SHARED_R_LIBS_USER}"/.ppdis_pkg_install_*
    if [ -n "$INSTALL_STAMP" ]; then
      touch "$INSTALL_STAMP"
      echo "Recorded install stamp: $INSTALL_STAMP"
    fi
  fi
  cleanup_pp_lock_if_safe "$PP_LOCK_DIR"
fi
echo ""

# Oklahoma run config.
JOB_CORES="${SLURM_CPUS_PER_TASK:-$PP_CORES}"
export OK_MEMORY_SAFE=true
export OK_PARALLEL_BACKEND=psock
export OK_CORES="${JOB_CORES}"
export OK_SENS_CORES="${JOB_CORES}"
# ATE forward sims are independent (CRN seeds are per-replicate); use the
# full job allocation unless the caller overrides OK_ATE_SIM_CORES / PP_ATE_SIM_CORES.
export OK_ATE_SIM_CORES="${OK_ATE_SIM_CORES:-${PP_ATE_SIM_CORES:-$JOB_CORES}}"
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
export OK_T_TRUNC_SENS_DAYS="$PP_T_TRUNC_SENS_DAYS"
if [ "$PP_RUN_T_TRUNC_SENSITIVITY" = "1" ] || [ "$PP_RUN_T_TRUNC_SENSITIVITY" = "true" ] || [ "$PP_RUN_T_TRUNC_SENSITIVITY" = "yes" ]; then
  export OK_RUN_T_TRUNC_SENSITIVITY=true
else
  export OK_RUN_T_TRUNC_SENSITIVITY=false
fi
export OK_SEM_TEMPORAL_WEIGHT="$PP_SEM_TEMPORAL_WEIGHT"
export OK_SEM_OPTIM_METHOD="$PP_SEM_OPTIM_METHOD"
export OK_SEM_SELECTION_TEMPERATURE="$PP_SEM_SELECTION_TEMPERATURE"
export OK_SEM_CHANGE_FACTOR_MIN_MULT="$PP_SEM_CHANGE_FACTOR_MIN_MULT"
export OK_SEM_CHANGE_FACTOR_MAX_MULT="$PP_SEM_CHANGE_FACTOR_MAX_MULT"
export OK_SEM_MAX_RELABEL_STEP_FRAC="$PP_SEM_MAX_RELABEL_STEP_FRAC"
export OK_SEM_FORCE_PARAM_UPDATE_FLIP_FRAC="$PP_SEM_FORCE_PARAM_UPDATE_FLIP_FRAC"
if [ "$PP_SEM_MONOTONE_COMPLETE_LL" = "1" ] || [ "$PP_SEM_MONOTONE_COMPLETE_LL" = "true" ] || [ "$PP_SEM_MONOTONE_COMPLETE_LL" = "yes" ]; then
  export OK_SEM_MONOTONE_COMPLETE_LL=true
else
  export OK_SEM_MONOTONE_COMPLETE_LL=false
fi
if [ "$PP_SEM_START_FROM_C" = "1" ] || [ "$PP_SEM_START_FROM_C" = "true" ] || [ "$PP_SEM_START_FROM_C" = "yes" ]; then
  export OK_SEM_START_FROM_C=true
else
  export OK_SEM_START_FROM_C=false
fi
export OK_SEM_BIV_N_THREADS="$PP_SEM_BIV_N_THREADS"
if [ -n "${PP_SEM_SINGLE_FLIP_FROM_ITER:-}" ]; then
  export OK_SEM_SINGLE_FLIP_FROM_ITER="$PP_SEM_SINGLE_FLIP_FROM_ITER"
fi
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
export OK_RUN_FIT_VARIABILITY="$PP_RUN_FIT_VARIABILITY"
export OK_FIT_VARIABILITY_REPS="$PP_FIT_VARIABILITY_REPS"
export OK_FIT_VARIABILITY_CORES="$PP_FIT_VARIABILITY_CORES"
if [ -n "${PP_FIT_VARIABILITY_PATCH_FILE:-}" ]; then
  export OK_FIT_VARIABILITY_PATCH_FILE="$PP_FIT_VARIABILITY_PATCH_FILE"
fi
if [ "${PP_FIT_VARIABILITY_ONLY:-0}" = "1" ] || [ "${PP_FIT_VARIABILITY_ONLY:-0}" = "true" ]; then
  export OK_FIT_VARIABILITY_ONLY=true
else
  export OK_FIT_VARIABILITY_ONLY=false
fi
if [ -n "${PP_BOOTSTRAP_PATCH_FILE:-}" ]; then
  export OK_BOOTSTRAP_PATCH_FILE="$PP_BOOTSTRAP_PATCH_FILE"
fi
if [ "${PP_BOOTSTRAP_ONLY:-0}" = "1" ] || [ "${PP_BOOTSTRAP_ONLY:-0}" = "true" ]; then
  export OK_BOOTSTRAP_ONLY=true
else
  export OK_BOOTSTRAP_ONLY=false
fi
if [ -n "${PP_T_TRUNC_SENS_PATCH_FILE:-}" ]; then
  export OK_T_TRUNC_SENS_PATCH_FILE="$PP_T_TRUNC_SENS_PATCH_FILE"
fi
if [ "${PP_T_TRUNC_SENS_ONLY:-0}" = "1" ] || [ "${PP_T_TRUNC_SENS_ONLY:-0}" = "true" ]; then
  export OK_T_TRUNC_SENS_ONLY=true
  export OK_RUN_T_TRUNC_SENSITIVITY=true
  export OK_RUN_BOOTSTRAP_ATE=false
  export OK_BOOT_N_REPS=0
else
  export OK_T_TRUNC_SENS_ONLY=false
fi
export OK_SMOKE_SEM_D_SEEDS="${PP_SMOKE_SEM_D_SEEDS:-0}"
export OK_SMOKE_SEM_D_TRUNC="${PP_SMOKE_SEM_D_TRUNC:-3}"
if [ "${PP_CD_ONLY:-0}" = "1" ] || [ "${PP_CD_ONLY:-0}" = "true" ] || [ "${PP_CD_ONLY:-0}" = "yes" ]; then
  export OK_CD_ONLY=true
else
  export OK_CD_ONLY=false
fi
if [ "${PP_UNIV_KDE_ONLY:-0}" = "1" ] || [ "${PP_UNIV_KDE_ONLY:-0}" = "true" ] || [ "${PP_UNIV_KDE_ONLY:-0}" = "yes" ]; then
  export OK_UNIV_KDE_ONLY=true
else
  export OK_UNIV_KDE_ONLY=false
fi
if [ -n "${PP_SKIP_CONTROL_SNAPSHOTS:-}" ]; then
  export OK_SKIP_CONTROL_SNAPSHOTS="$PP_SKIP_CONTROL_SNAPSHOTS"
fi
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
export OK_KDE_BW_METHOD="$PP_KDE_BW_METHOD"
export OK_PRIMARY_PARTITION="$PP_PRIMARY_PARTITION"
if [ -n "${PP_KDE_BW_SENS_KM:-}" ]; then
  export OK_KDE_BW_SENS_KM="$PP_KDE_BW_SENS_KM"
fi
if [ -n "${PP_RUN_PARTITION_SENSITIVITY:-}" ]; then
  case "$(echo "$PP_RUN_PARTITION_SENSITIVITY" | tr '[:upper:]' '[:lower:]')" in
    1|true|yes|y|on) export OK_RUN_PARTITION_SENSITIVITY=true ;;
    *) export OK_RUN_PARTITION_SENSITIVITY=false ;;
  esac
fi
export OK_BOOT_SEM_INNER_ITER="$PP_BOOT_SEM_INNER"
export OK_BOOT_OUTER_CORES="$PP_BOOT_OUTER_CORES"
export OK_ATE_N_SIMS="$PP_ATE_N_SIMS"
# ATE evaluation: bivariate|univariate × all_or_nothing|observed.
# CLI --ate-bivariate / --ate-contrast override env defaults.
if [ -n "${PP_ATE_BIVARIATE}" ]; then
  case "$(echo "$PP_ATE_BIVARIATE" | tr '[:upper:]' '[:lower:]')" in
    0|false|no|n|off|marginal|univariate|univ) export OK_ATE_BIVARIATE=false ;;
    *) export OK_ATE_BIVARIATE=true ;;
  esac
else
  export OK_ATE_BIVARIATE="${OK_ATE_BIVARIATE:-true}"
fi
if [ -n "${PP_ATE_CONTRAST}" ]; then
  export OK_ATE_CONTRAST="$(echo "$PP_ATE_CONTRAST" | tr '[:upper:]' '[:lower:]')"
else
  export OK_ATE_CONTRAST="${OK_ATE_CONTRAST:-all_or_nothing}"
fi
if [ -n "${PP_ATE_SCENARIO}" ]; then
  export OK_ATE_SCENARIO="$PP_ATE_SCENARIO"
fi
export OK_GLOBAL_SEED="$PP_SEED"
export OK_BOOT_SEED="$PP_SEED"
export OK_IDENTICAL_RANDOMNESS=false
export OK_BOOT_IDENTICAL_RANDOMNESS=false
export OK_BOOT_GUARD_DEGENERATE=true
if [ "${PP_SKIP_FULL_REPORT:-0}" = "1" ] || [ "${PP_SKIP_FULL_REPORT:-0}" = "true" ] || [ "${PP_SKIP_FULL_REPORT:-0}" = "yes" ]; then
  # Empty formats skips oklahoma_report.html so diagnostic HTML stays intact.
  export OK_REPORT_FORMATS=""
  echo "Skipping full oklahoma_report.html render (PP_SKIP_FULL_REPORT=1); use slim C/D report locally."
else
  export OK_REPORT_FORMATS=html
fi
echo "ATE settings: OK_ATE_BIVARIATE=$OK_ATE_BIVARIATE OK_ATE_CONTRAST=$OK_ATE_CONTRAST OK_ATE_SCENARIO=${OK_ATE_SCENARIO:-} OK_ATE_N_SIMS=$OK_ATE_N_SIMS OK_ATE_SIM_CORES=$OK_ATE_SIM_CORES"

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
