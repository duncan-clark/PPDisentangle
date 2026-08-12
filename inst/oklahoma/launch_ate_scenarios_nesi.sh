#!/usr/bin/env bash
# Launch 4 parallel Oklahoma bootstrap-only ATE scenarios on NeSI.
#
# Scenarios (frozen fits from seed RDS; full debiased parametric bootstrap):
#   univ_aon  — univariate margins × all-or-nothing
#   univ_obs  — univariate margins × observed vs control-everywhere
#   biv_aon   — bivariate law × all-or-nothing
#   biv_obs   — bivariate law × observed vs control-everywhere
#
# Usage (from package root, with working `ssh mahuika`):
#   bash inst/oklahoma/launch_ate_scenarios_nesi.sh
#   PP_BOOT_REPS=512 PP_SCENARIO_SUFFIX=_boot512 bash inst/oklahoma/launch_ate_scenarios_nesi.sh
#   PP_NESI_PUSH=1 bash inst/oklahoma/launch_ate_scenarios_nesi.sh
#
# Writes job ids to:
#   PPDisentangle-output/oklahoma/ate_scenarios/launch_jobs${PP_SCENARIO_SUFFIX}.txt
#
# When jobs finish:
#   bash inst/oklahoma/watch_and_publish_ate_scenarios.sh launch_jobs_boot512.txt

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
# shellcheck source=../include/nesi_config.sh
source "$PKG_ROOT/inst/include/nesi_config.sh"
# shellcheck source=../include/git_sync.sh
source "$PKG_ROOT/inst/include/git_sync.sh"

pp_nesi_init_config "$PKG_ROOT"
pp_nesi_print_config

# Tunables (env overrides)
PP_BOOT_REPS="${PP_BOOT_REPS:-512}"
PP_CORES="${PP_CORES:-64}"
PP_ATE_N_SIMS="${PP_ATE_N_SIMS:-500}"
PP_BOOT_SEM_INNER="${PP_BOOT_SEM_INNER:-2000}"
PP_TIME="${PP_TIME:-36:00:00}"
PP_MEM="${PP_MEM:-64G}"
# RAM-aware outer workers: 2 GB/worker + 8 GB reserve. On the default 64G
# allocation this is 28 workers (job 8353220 was ~0.85 GB RSS/worker).
# Override with PP_BOOT_OUTER_CORES; raise PP_MEM to get more workers.
# shellcheck source=../include/nesi_mem.sh
source "$PKG_ROOT/inst/include/nesi_mem.sh"
if [ -z "${PP_BOOT_OUTER_CORES:-}" ]; then
  PP_BOOT_OUTER_CORES="$(pp_boot_outer_from_mem "$PP_CORES" "$PP_MEM" "${PP_BOOT_WORKER_GB:-2}" "${PP_BOOT_MEM_RESERVE_GB:-8}")"
fi
PP_SCENARIO_SUFFIX="${PP_SCENARIO_SUFFIX:-_boot512}"
# Refresh immutable seed from current remote for_paper.rds (fit params only matter).
PP_REFRESH_SEED="${PP_REFRESH_SEED:-1}"

LOCAL_OK="${PP_NESI_LOCAL_OUTPUT}/oklahoma"
LOCAL_SEED="${LOCAL_OK}/for_paper.rds"
if [ ! -f "$LOCAL_SEED" ]; then
  echo "ERROR: missing local seed RDS: $LOCAL_SEED" >&2
  exit 1
fi

REMOTE_OK="${PP_NESI_REMOTE_OUTPUT}/oklahoma"
LAUNCH_TAG="launch_jobs${PP_SCENARIO_SUFFIX}"
REMOTE_SEED_IMMUTABLE="${REMOTE_OK}/ate_scenarios/seed_for_paper${PP_SCENARIO_SUFFIX}.rds"
mkdir -p "${LOCAL_OK}/ate_scenarios"

pp_git_push_branch "$PKG_ROOT" || exit 1
REMOTE_BRANCH="${PP_GIT_BRANCH:-}"
if [ -z "$REMOTE_BRANCH" ]; then
  REMOTE_BRANCH="$(pp_git_current_branch "$PKG_ROOT")"
fi

echo "Syncing seed + submitting 4 scenario jobs (boot_reps=${PP_BOOT_REPS} cores=${PP_CORES} outer=${PP_BOOT_OUTER_CORES} mem=${PP_MEM} suffix=${PP_SCENARIO_SUFFIX}) via ${PP_NESI_SSH} ..."
ssh "$PP_NESI_SSH" bash -s -- \
  "$PP_NESI_REMOTE_PKG" \
  "$PP_NESI_REMOTE_OUTPUT" \
  "$REMOTE_BRANCH" \
  "$REMOTE_SEED_IMMUTABLE" \
  "$PP_BOOT_REPS" \
  "$PP_CORES" \
  "$PP_BOOT_OUTER_CORES" \
  "$PP_ATE_N_SIMS" \
  "$PP_BOOT_SEM_INNER" \
  "$PP_TIME" \
  "$PP_MEM" \
  "$PP_SCENARIO_SUFFIX" \
  "$PP_REFRESH_SEED" \
  "$LAUNCH_TAG" \
  <<'REMOTE'
set -euo pipefail
REMOTE_PKG="$1"
REMOTE_OUTPUT="$2"
REMOTE_BRANCH="$3"
REMOTE_SEED_IMMUTABLE="$4"
BOOT_REPS="$5"
CORES="$6"
OUTER="$7"
ATE_N_SIMS="$8"
BOOT_INNER="$9"
WALLTIME="${10}"
MEM="${11}"
SUFFIX="${12}"
REFRESH_SEED="${13}"
LAUNCH_TAG="${14}"

cd "$REMOTE_PKG"
export PPDISENTANGLE_OUTPUT_ROOT="$REMOTE_OUTPUT"
export PP_GIT_BRANCH="$REMOTE_BRANCH"
source inst/include/git_sync.sh
pp_git_sync_repo "$REMOTE_PKG"

OK_DIR="${REMOTE_OUTPUT}/oklahoma"
SCEN_DIR="${OK_DIR}/ate_scenarios"
mkdir -p "$SCEN_DIR"

if [ "$REFRESH_SEED" = "1" ] || [ ! -f "$REMOTE_SEED_IMMUTABLE" ]; then
  if [ ! -f "${OK_DIR}/for_paper.rds" ]; then
    echo "ERROR: no for_paper.rds on NeSI to seed scenarios" >&2
    exit 2
  fi
  cp -p "${OK_DIR}/for_paper.rds" "$REMOTE_SEED_IMMUTABLE"
  echo "Wrote seed from current for_paper.rds: $REMOTE_SEED_IMMUTABLE"
else
  echo "Using existing seed: $REMOTE_SEED_IMMUTABLE"
fi

COMMON=(
  --mode bootstrap-only
  --cores "$CORES"
  --boot-reps "$BOOT_REPS"
  --ate-n-sims "$ATE_N_SIMS"
  --boot-refit-scope full
  --boot-outer-cores "$OUTER"
  --boot-sem-inner "$BOOT_INNER"
  --time "$WALLTIME"
  --mem "$MEM"
)

submit_one() {
  local base="$1" biv="$2" contrast="$3"
  local scen="${base}${SUFFIX}"
  local patch="${SCEN_DIR}/seed_${scen}.rds"
  cp -p "$REMOTE_SEED_IMMUTABLE" "$patch"
  echo "=== Submitting scenario=${scen} bivariate=${biv} contrast=${contrast} boot_reps=${BOOT_REPS} ==="
  local out
  out="$(
    OK_ATE_SCENARIO="$scen" \
    bash inst/oklahoma/run_nesi.sh \
      "${COMMON[@]}" \
      --bootstrap-patch-file "$patch" \
      --ate-bivariate "$biv" \
      --ate-contrast "$contrast" \
      --ate-scenario "$scen" 2>&1
  )"
  echo "$out"
  local jid
  jid="$(printf '%s\n' "$out" | sed -n 's/^Job \([0-9][0-9]*\) submitted$/\1/p' | tail -1)"
  if [ -z "$jid" ]; then
    jid="$(printf '%s\n' "$out" | sed -n 's/.*Submitted batch job \([0-9][0-9]*\).*/\1/p' | tail -1)"
  fi
  if [ -z "$jid" ]; then
    echo "ERROR: could not parse job id for scenario=${scen}" >&2
    exit 3
  fi
  # base \t job \t biv \t contrast \t full_scenario_id \t boot_reps
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$base" "$jid" "$biv" "$contrast" "$scen" "$BOOT_REPS" \
    >> "${SCEN_DIR}/${LAUNCH_TAG}.txt"
  echo "scenario=${scen} job_id=${jid} base=${base} boot_reps=${BOOT_REPS}" > "${SCEN_DIR}/job_${scen}.txt"
  echo "RECORD ${scen} ${jid}"
}

: > "${SCEN_DIR}/${LAUNCH_TAG}.txt"
submit_one univ_aon false all_or_nothing
submit_one univ_obs false observed
submit_one biv_aon  true  all_or_nothing
submit_one biv_obs  true  observed

echo "=== All 4 scenarios submitted (boot_reps=${BOOT_REPS}) ==="
cat "${SCEN_DIR}/${LAUNCH_TAG}.txt"
REMOTE

# Pull the launch record locally
mkdir -p "${LOCAL_OK}/ate_scenarios"
rsync -avz "${PP_NESI_SSH}:${REMOTE_OK}/ate_scenarios/${LAUNCH_TAG}.txt" \
  "${LOCAL_OK}/ate_scenarios/${LAUNCH_TAG}.txt" || true
rsync -avz "${PP_NESI_SSH}:${REMOTE_OK}/ate_scenarios/job_"*"${PP_SCENARIO_SUFFIX}".txt \
  "${LOCAL_OK}/ate_scenarios/" || true

echo ""
echo "Launch record (${LAUNCH_TAG}.txt):"
cat "${LOCAL_OK}/ate_scenarios/${LAUNCH_TAG}.txt" 2>/dev/null || echo "(not fetched yet)"
echo ""
echo "Watch / merge / regen when done:"
echo "  bash inst/oklahoma/watch_and_publish_ate_scenarios.sh ${LAUNCH_TAG}.txt"
echo "Or fetch-only compare:"
echo "  bash inst/oklahoma/fetch_and_render_ate_scenarios.sh"
