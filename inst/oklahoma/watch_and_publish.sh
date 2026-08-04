#!/usr/bin/env bash
# Wait for an Oklahoma NeSI job, fetch results, promote for_paper.rds,
# regenerate paper figures/tables, and render the HTML report.
#
# Usage (from package root):
#   bash inst/oklahoma/watch_and_publish.sh 8212755
#   PP_NESI_POLL_SECS=120 bash inst/oklahoma/watch_and_publish.sh 8212755
#
# On success prints:
#   OKLAHOMA_PUBLISH_DONE {"job_id":"...","rds":"...","html":"...","ate_diff":"..."}

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
# shellcheck source=../include/nesi_config.sh
source "$PKG_ROOT/inst/include/nesi_config.sh"

JOB_ID="${1:-}"
POLL_SECS="${PP_NESI_POLL_SECS:-120}"

if [ -z "$JOB_ID" ]; then
  echo "Usage: bash inst/oklahoma/watch_and_publish.sh JOB_ID" >&2
  exit 1
fi

pp_nesi_init_config "$PKG_ROOT"
LOCAL_OK="${PP_NESI_LOCAL_OUTPUT}/oklahoma"
LOG="${LOCAL_OK}/watch_and_publish_${JOB_ID}.log"
mkdir -p "$LOCAL_OK" "${HOME}/.ssh/sockets"

# Prefer a long-lived multiplexed session (NeSI lander is MFA / keyboard-interactive).
# Open once in another terminal:  ssh mahuika
# then this watcher reuses ControlMaster for the rest of the run.
# Options must come before the host; anything after is treated as the remote command.
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
# Quote-safe rsync transport (ControlPath may contain spaces in theory).
RSYNC_RSH="ssh ${SSH_OPTS[*]}"
export RSYNC_RSH
rsync_nesi() { rsync "$@"; }

exec > >(tee -a "$LOG") 2>&1

echo "=== Oklahoma watch_and_publish job=${JOB_ID} started $(date) ==="
echo "Polling ${PP_NESI_SSH} every ${POLL_SECS}s; log=${LOG}"

ensure_ssh() {
  local why="${1:-ssh}"
  while ! ssh_nesi "echo SSH_OK" >/dev/null 2>&1; do
    # Drop stale mux sockets that break with "Broken pipe"
    rm -f "${HOME}/.ssh/sockets/ssh_mux_"* 2>/dev/null || true
    echo "OKLAHOMA_WATCH_NEED_SSH {\"job_id\":\"${JOB_ID}\",\"why\":\"${why}\",\"hint\":\"run: ssh mahuika\"}"
    sleep "$POLL_SECS"
  done
  echo "SSH ready via ControlMaster (${why}) ($(date))"
}

# --- wait for working SSH (MFA mux) ---
ensure_ssh "startup"

# --- wait for job leave queue ---
while true; do
  if ! ssh_nesi "echo SSH_OK" >/dev/null 2>&1; then
    ensure_ssh "poll"
    continue
  fi
  state="$(ssh_nesi "squeue -h -j ${JOB_ID} -o '%T' 2>/dev/null | head -1" || true)"
  if [ -z "$state" ]; then
    echo "Job ${JOB_ID} left the queue ($(date))."
    break
  fi
  echo "Job ${JOB_ID} state=${state} ($(date))"
  sleep "$POLL_SECS"
done

# Re-auth if the long wait killed ControlPersist / mux
ensure_ssh "pre-fetch"

# Confirm Slurm terminal state if sacct available
sacct_line="$(ssh_nesi "sacct -n -X -j ${JOB_ID} -o State,ExitCode --parsable2 2>/dev/null | head -1" || true)"
echo "sacct: ${sacct_line:-unavailable}"
state_main="${sacct_line%%|*}"
if [ -n "$state_main" ] && ! echo "$state_main" | grep -Eqi 'COMPLETED|COMPLETING'; then
  echo "WARNING: job did not report COMPLETED (state=${state_main}). Continuing fetch anyway."
fi

# --- fetch ---
echo "Fetching Oklahoma artifacts for job ${JOB_ID} ..."
# Ensure child scripts inherit a usable PP_NESI_SSH alias (same host).
export PP_NESI_SSH
bash "$PKG_ROOT/inst/nesi/fetch.sh" oklahoma "$JOB_ID"

# Broader sync of oklahoma tree so untagged helpers (html, paper/) arrive too
rsync_nesi -avz --human-readable \
  "${PP_NESI_SSH}:${PP_NESI_REMOTE_OUTPUT}/oklahoma/" \
  "${LOCAL_OK}/" \
  --include="*job${JOB_ID}*" \
  --include="${JOB_ID}_oklahoma_slurm.*" \
  --include="oklahoma_results.rds" \
  --include="oklahoma_report.html" \
  --include="oklahoma_report_files/***" \
  --include="plots/***" \
  --include="paper/***" \
  --include="last_run_sync_stamp.txt" \
  --exclude="*" || true

# --- promote for_paper.rds ---
CANDIDATES=(
  "${LOCAL_OK}/oklahoma_results_job${JOB_ID}.rds"
  "${LOCAL_OK}/oklahoma_results.rds"
)
SRC_RDS=""
for c in "${CANDIDATES[@]}"; do
  if [ -f "$c" ]; then SRC_RDS="$c"; break; fi
done
# Fallback: newest job-tagged results
if [ -z "$SRC_RDS" ]; then
  SRC_RDS="$(ls -t "${LOCAL_OK}"/oklahoma_results_job${JOB_ID}*.rds 2>/dev/null | head -1 || true)"
fi
if [ -z "$SRC_RDS" ] || [ ! -f "$SRC_RDS" ]; then
  echo "ERROR: no oklahoma_results*.rds found for job ${JOB_ID} under ${LOCAL_OK}" >&2
  echo "OKLAHOMA_PUBLISH_FAILED {\"job_id\":\"${JOB_ID}\",\"reason\":\"missing_rds\"}"
  exit 2
fi

FOR_PAPER="${LOCAL_OK}/for_paper.rds"
if [ -f "$FOR_PAPER" ]; then
  bak="${LOCAL_OK}/for_paper_backup_before_${JOB_ID}_$(date +%Y%m%d%H%M%S).rds"
  cp -p "$FOR_PAPER" "$bak"
  echo "Backed up previous for_paper.rds -> ${bak}"
fi
cp -p "$SRC_RDS" "$FOR_PAPER"
echo "Promoted ${SRC_RDS} -> ${FOR_PAPER}"

# Quick sanity: ATE method/contrast if readable
Rscript -e "
`%||%` <- function(a, b) if (!is.null(a)) a else b
r <- readRDS('${FOR_PAPER}')
cfg <- r\$config
cat(sprintf('ATE_BIVARIATE=%s ATE_CONTRAST=%s ATE_N_SIMS=%s BOOT_N=%s\n',
  as.character(cfg\$ATE_BIVARIATE %||% NA),
  as.character(cfg\$ATE_CONTRAST %||% NA),
  as.character(cfg\$ATE_N_SIMS %||% NA),
  as.character(cfg\$BOOT_N_REPS %||% NA)))
for (lab in c('E', 'F')) {
  ate <- r\$fits_named[[lab]]\$ate
  if (is.null(ate)) next
  m <- mean(ate\$all_nothing_sim\$total_saved, na.rm = TRUE)
  cat(sprintf('Fit %s mean saved=%.2f method=%s\n', lab, m,
              as.character(ate\$ate_method %||% 'NA')))
}
" || echo "WARNING: could not summarize ATE from RDS"

# --- paper assets ---
echo "Regenerating oklahoma paper assets ..."
(
  cd "$PKG_ROOT"
  Rscript inst/oklahoma/paper/oklahoma_paper_assets.R --input "$FOR_PAPER"
)

ATE_DIFF="${LOCAL_OK}/paper/generated/figures/ATE_diff.pdf"
# Local Overleaf-style mirror used by revised.tex paths
PLOTS_MIRROR="${PKG_ROOT}/docs/paper/plots/oklahoma"
mkdir -p "$PLOTS_MIRROR"
if [ -d "${LOCAL_OK}/paper/generated/figures" ]; then
  rsync -a "${LOCAL_OK}/paper/generated/figures/" "$PLOTS_MIRROR/"
  echo "Synced figures -> ${PLOTS_MIRROR}"
fi

# --- HTML report ---
echo "Rendering oklahoma_report.html ..."
(
  cd "$PKG_ROOT/inst/oklahoma"
  Rscript render_oklahoma_report.R --results "$FOR_PAPER" --to html
)
HTML="${LOCAL_OK}/oklahoma_report.html"
if [ -f "$HTML" ]; then
  cp -p "$HTML" "${LOCAL_OK}/oklahoma_report_job${JOB_ID}.html"
fi

echo "=== Oklahoma watch_and_publish finished $(date) ==="
echo "OKLAHOMA_PUBLISH_DONE {\"job_id\":\"${JOB_ID}\",\"rds\":\"${FOR_PAPER}\",\"html\":\"${HTML}\",\"ate_diff\":\"${ATE_DIFF}\",\"plots_mirror\":\"${PLOTS_MIRROR}\"}"
