#!/usr/bin/env bash
# Wait for a C/D-primary Oklahoma NeSI job, fetch results, promote
# oklahoma_cd_current.rds, and render the slim C/D HTML.
#
# Does NOT overwrite oklahoma_report.html / oklahoma_report_diagnostic.html
# and does NOT promote for_paper.rds (diagnostic / paper paths stay intact).
#
# Usage (from package root):
#   bash inst/oklahoma/watch_and_publish_cd.sh JOB_ID
#   PP_NESI_POLL_SECS=900 bash inst/oklahoma/watch_and_publish_cd.sh JOB_ID

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
# shellcheck source=../include/nesi_config.sh
source "$PKG_ROOT/inst/include/nesi_config.sh"

JOB_ID="${1:-}"
POLL_SECS="${PP_NESI_POLL_SECS:-1800}"

if [ -z "$JOB_ID" ]; then
  echo "Usage: bash inst/oklahoma/watch_and_publish_cd.sh JOB_ID" >&2
  exit 1
fi

pp_nesi_init_config "$PKG_ROOT"
LOCAL_OK="${PP_NESI_LOCAL_OUTPUT}/oklahoma"
LOG="${LOCAL_OK}/watch_and_publish_cd_${JOB_ID}.log"
mkdir -p "$LOCAL_OK" "${HOME}/.ssh/sockets"

SSH_OPTS=(
  -o BatchMode=yes
  -o ControlMaster=auto
  -o ControlPersist=24h
  -o "ControlPath=${HOME}/.ssh/sockets/ssh_mux_%h_%p_%r"
  -o ServerAliveInterval=60
  -o ServerAliveCountMax=3
  -o ConnectTimeout=30
)
ssh_nesi() { ssh "${SSH_OPTS[@]}" "$PP_NESI_SSH" "$@"; }
RSYNC_RSH="ssh ${SSH_OPTS[*]}"
export RSYNC_RSH

exec > >(tee -a "$LOG") 2>&1

echo "=== Oklahoma C/D watch_and_publish job=${JOB_ID} started $(date) ==="
echo "Polling ${PP_NESI_SSH} every ${POLL_SECS}s; log=${LOG}"

ensure_ssh() {
  local why="${1:-ssh}"
  while ! ssh_nesi "echo SSH_OK" >/dev/null 2>&1; do
    rm -f "${HOME}/.ssh/sockets/ssh_mux_"* 2>/dev/null || true
    echo "OKLAHOMA_CD_WATCH_NEED_SSH {\"job_id\":\"${JOB_ID}\",\"why\":\"${why}\",\"hint\":\"run: ssh mahuika\"}"
    sleep "$POLL_SECS"
  done
  echo "SSH ready via ControlMaster (${why}) ($(date))"
}

ensure_ssh "startup"

while true; do
  state="$(ssh_nesi "squeue -h -j ${JOB_ID} -o '%T' 2>/dev/null | head -1" || true)"
  ssh_rc=$?
  if [ "$ssh_rc" -ne 0 ]; then
    echo "SSH poll failed (rc=${ssh_rc}); will retry after ${POLL_SECS}s ($(date))"
    ensure_ssh "poll-recover"
    continue
  fi
  if [ -z "$state" ]; then
    echo "Job ${JOB_ID} left the queue ($(date))."
    break
  fi
  echo "Job ${JOB_ID} state=${state} ($(date)); next poll in ${POLL_SECS}s"
  sleep "$POLL_SECS"
done

ensure_ssh "pre-fetch"

sacct_line="$(ssh_nesi "sacct -n -X -j ${JOB_ID} -o State,ExitCode --parsable2 2>/dev/null | head -1" || true)"
echo "sacct: ${sacct_line:-unavailable}"
state_main="${sacct_line%%|*}"
if [ -n "$state_main" ] && ! echo "$state_main" | grep -Eqi 'COMPLETED|COMPLETING'; then
  echo "WARNING: job did not report COMPLETED (state=${state_main}). Continuing fetch anyway."
fi

echo "Fetching Oklahoma artifacts for job ${JOB_ID} ..."
export PP_NESI_SSH
bash "$PKG_ROOT/inst/nesi/fetch.sh" oklahoma "$JOB_ID"

rsync -avz --human-readable \
  "${PP_NESI_SSH}:${PP_NESI_REMOTE_OUTPUT}/oklahoma/" \
  "${LOCAL_OK}/" \
  --include="*job${JOB_ID}*" \
  --include="${JOB_ID}_oklahoma_slurm.*" \
  --include="oklahoma_results.rds" \
  --exclude="*" || true

SRC_RDS=""
for c in \
  "${LOCAL_OK}/oklahoma_results_job${JOB_ID}.rds" \
  "${LOCAL_OK}/oklahoma_results.rds"
do
  if [ -f "$c" ]; then SRC_RDS="$c"; break; fi
done
if [ -z "$SRC_RDS" ]; then
  SRC_RDS="$(ls -t "${LOCAL_OK}"/oklahoma_results_job${JOB_ID}*.rds 2>/dev/null | head -1 || true)"
fi
if [ -z "$SRC_RDS" ] || [ ! -f "$SRC_RDS" ]; then
  echo "ERROR: no oklahoma_results*.rds found for job ${JOB_ID} under ${LOCAL_OK}" >&2
  echo "OKLAHOMA_CD_PUBLISH_FAILED {\"job_id\":\"${JOB_ID}\",\"reason\":\"missing_rds\"}"
  exit 2
fi

CURRENT_RDS="${LOCAL_OK}/oklahoma_cd_current.rds"
PREV_RDS="${LOCAL_OK}/oklahoma_cd_previous.rds"

# If an update launcher already archived previous, keep it. Otherwise, on first
# promote of a new current when one already exists, archive automatically.
if [ -f "$CURRENT_RDS" ] && [ ! -f "$PREV_RDS" ]; then
  cp -p "$CURRENT_RDS" "$PREV_RDS"
  echo "Auto-archived prior oklahoma_cd_current.rds -> ${PREV_RDS}"
fi

cp -p "$SRC_RDS" "$CURRENT_RDS"
echo "Promoted ${SRC_RDS} -> ${CURRENT_RDS}"

echo "Rendering slim C/D report (diagnostic HTML untouched) ..."
(
  cd "$PKG_ROOT/inst/oklahoma"
  if [ -f "$PREV_RDS" ]; then
    Rscript render_oklahoma_report_cd.R --results "$CURRENT_RDS" --prev "$PREV_RDS"
  else
    Rscript render_oklahoma_report_cd.R --results "$CURRENT_RDS"
  fi
)

HTML="${LOCAL_OK}/oklahoma_report_cd.html"
if [ -f "$HTML" ]; then
  cp -p "$HTML" "${LOCAL_OK}/oklahoma_report_cd_job${JOB_ID}.html"
fi

echo "=== Oklahoma C/D watch_and_publish finished $(date) ==="
echo "OKLAHOMA_CD_PUBLISH_DONE {\"job_id\":\"${JOB_ID}\",\"rds\":\"${CURRENT_RDS}\",\"prev\":\"${PREV_RDS}\",\"html\":\"${HTML}\"}"
