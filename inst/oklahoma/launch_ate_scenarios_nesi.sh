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
#   PP_NESI_PUSH=1 bash inst/oklahoma/launch_ate_scenarios_nesi.sh
#
# Writes job ids to PPDisentangle-output/oklahoma/ate_scenarios/launch_jobs.txt

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
# shellcheck source=../include/nesi_config.sh
source "$PKG_ROOT/inst/include/nesi_config.sh"
# shellcheck source=../include/git_sync.sh
source "$PKG_ROOT/inst/include/git_sync.sh"

pp_nesi_init_config "$PKG_ROOT"
pp_nesi_print_config

LOCAL_OK="${PP_NESI_LOCAL_OUTPUT}/oklahoma"
LOCAL_SEED="${LOCAL_OK}/for_paper.rds"
if [ ! -f "$LOCAL_SEED" ]; then
  echo "ERROR: missing local seed RDS: $LOCAL_SEED" >&2
  exit 1
fi

REMOTE_OK="${PP_NESI_REMOTE_OUTPUT}/oklahoma"
REMOTE_SEED_IMMUTABLE="${REMOTE_OK}/ate_scenarios/seed_for_paper.rds"
mkdir -p "${LOCAL_OK}/ate_scenarios"

pp_git_push_branch "$PKG_ROOT" || exit 1
REMOTE_BRANCH="${PP_GIT_BRANCH:-}"
if [ -z "$REMOTE_BRANCH" ]; then
  REMOTE_BRANCH="$(pp_git_current_branch "$PKG_ROOT")"
fi

echo "Syncing immutable seed + submitting 4 scenario jobs via ${PP_NESI_SSH} ..."
ssh "$PP_NESI_SSH" bash -s -- \
  "$PP_NESI_REMOTE_PKG" \
  "$PP_NESI_REMOTE_OUTPUT" \
  "$REMOTE_BRANCH" \
  "$REMOTE_SEED_IMMUTABLE" \
  <<'REMOTE'
set -euo pipefail
REMOTE_PKG="$1"
REMOTE_OUTPUT="$2"
REMOTE_BRANCH="$3"
REMOTE_SEED_IMMUTABLE="$4"
cd "$REMOTE_PKG"
export PPDISENTANGLE_OUTPUT_ROOT="$REMOTE_OUTPUT"
export PP_GIT_BRANCH="$REMOTE_BRANCH"
source inst/include/git_sync.sh
pp_git_sync_repo "$REMOTE_PKG"

OK_DIR="${REMOTE_OUTPUT}/oklahoma"
SCEN_DIR="${OK_DIR}/ate_scenarios"
mkdir -p "$SCEN_DIR"
# Prefer existing immutable seed; else copy current for_paper.rds once.
if [ ! -f "$REMOTE_SEED_IMMUTABLE" ]; then
  if [ ! -f "${OK_DIR}/for_paper.rds" ]; then
    echo "ERROR: no for_paper.rds on NeSI to seed scenarios" >&2
    exit 2
  fi
  cp -p "${OK_DIR}/for_paper.rds" "$REMOTE_SEED_IMMUTABLE"
  echo "Created immutable seed: $REMOTE_SEED_IMMUTABLE"
else
  echo "Using existing immutable seed: $REMOTE_SEED_IMMUTABLE"
fi

COMMON=(
  --mode bootstrap-only
  --cores 32
  --boot-reps 64
  --ate-n-sims 500
  --boot-refit-scope full
  --boot-outer-cores 6
  --boot-sem-inner 2000
  --time 24:00:00
  --mem 200G
)

submit_one() {
  local scen="$1" biv="$2" contrast="$3"
  local patch="${SCEN_DIR}/seed_${scen}.rds"
  cp -p "$REMOTE_SEED_IMMUTABLE" "$patch"
  echo "=== Submitting scenario=${scen} bivariate=${biv} contrast=${contrast} ==="
  # Capture sbatch job id from run_nesi.sh output
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
  printf '%s\t%s\t%s\t%s\n' "$scen" "$jid" "$biv" "$contrast" >> "${SCEN_DIR}/launch_jobs.txt"
  echo "scenario=${scen} job_id=${jid}" > "${SCEN_DIR}/job_${scen}.txt"
  echo "RECORD ${scen} ${jid}"
}

: > "${SCEN_DIR}/launch_jobs.txt"
submit_one univ_aon false all_or_nothing
submit_one univ_obs false observed
submit_one biv_aon  true  all_or_nothing
submit_one biv_obs  true  observed

echo "=== All 4 scenarios submitted ==="
cat "${SCEN_DIR}/launch_jobs.txt"
REMOTE

# Pull the launch record locally
mkdir -p "${LOCAL_OK}/ate_scenarios"
rsync -avz "${PP_NESI_SSH}:${REMOTE_OK}/ate_scenarios/launch_jobs.txt" \
  "${LOCAL_OK}/ate_scenarios/launch_jobs.txt" || true
rsync -avz "${PP_NESI_SSH}:${REMOTE_OK}/ate_scenarios/job_"*.txt \
  "${LOCAL_OK}/ate_scenarios/" || true

echo ""
echo "Launch record:"
cat "${LOCAL_OK}/ate_scenarios/launch_jobs.txt" 2>/dev/null || echo "(not fetched yet)"
echo ""
echo "When all jobs finish, fetch + render compare HTML with:"
echo "  bash inst/oklahoma/fetch_and_render_ate_scenarios.sh"
