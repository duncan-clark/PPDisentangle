#!/usr/bin/env bash
# Launch a C/D-focused Oklahoma refresh on NeSI with large t_trunc.
#
# Defaults:
#   --mode cd-primary
#   --sem-t-trunc-days 90
#   boot_reps=0, KDE single, no BW/trunc sensitivity
#   SEM inner 2000, outer maxit 5000, ATE 500 sims
#
# Usage (from package root):
#   bash inst/oklahoma/launch_oklahoma_cd_nesi.sh
#   PP_SEM_T_TRUNC_DAYS=120 bash inst/oklahoma/launch_oklahoma_cd_nesi.sh
#   PP_NESI_PUSH=1 bash inst/oklahoma/launch_oklahoma_cd_nesi.sh
#
# After the job completes:
#   bash inst/oklahoma/watch_and_publish_cd.sh JOB_ID

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

PP_SEM_T_TRUNC_DAYS="${PP_SEM_T_TRUNC_DAYS:-90}"
PP_CORES="${PP_CORES:-64}"
PP_MEM="${PP_MEM:-200G}"
PP_TIME="${PP_TIME:-06:00:00}"
PP_SEM_INNER="${PP_SEM_INNER:-2000}"
PP_ATE_N_SIMS="${PP_ATE_N_SIMS:-500}"
PP_NESI_PUSH="${PP_NESI_PUSH:-1}"

# shellcheck source=../include/nesi_config.sh
source "$PKG_ROOT/inst/include/nesi_config.sh"
pp_nesi_init_config "$PKG_ROOT"
LOCAL_OK="${PP_NESI_LOCAL_OUTPUT}/oklahoma"
mkdir -p "$LOCAL_OK"
PREV_RDS="${LOCAL_OK}/oklahoma_cd_previous.rds"
# Seed previous from today's diagnostic run so the first 90d slim report can show diffs.
if [ ! -f "$PREV_RDS" ] && [ -f "${LOCAL_OK}/oklahoma_results_job8430360.rds" ]; then
  cp -p "${LOCAL_OK}/oklahoma_results_job8430360.rds" "$PREV_RDS"
  echo "Seeded oklahoma_cd_previous.rds from oklahoma_results_job8430360.rds (diagnostic baseline)."
fi

echo "=== Launch Oklahoma C/D primary (t_trunc=${PP_SEM_T_TRUNC_DAYS}d) ==="
echo "cores=${PP_CORES} mem=${PP_MEM} time=${PP_TIME} sem_inner=${PP_SEM_INNER} ate_n_sims=${PP_ATE_N_SIMS}"

PP_NESI_PUSH="$PP_NESI_PUSH" bash "$PKG_ROOT/inst/nesi/submit.sh" oklahoma \
  --mode cd-primary \
  --sem-t-trunc-days "$PP_SEM_T_TRUNC_DAYS" \
  --cores "$PP_CORES" \
  --mem "$PP_MEM" \
  --time "$PP_TIME" \
  --sem-inner "$PP_SEM_INNER" \
  --ate-n-sims "$PP_ATE_N_SIMS" \
  --ate-bivariate true \
  --ate-contrast all_or_nothing \
  "$@"

echo ""
echo "When the job finishes, publish the slim C/D HTML (keeps diagnostic HTML untouched):"
echo "  bash inst/oklahoma/watch_and_publish_cd.sh JOB_ID"
