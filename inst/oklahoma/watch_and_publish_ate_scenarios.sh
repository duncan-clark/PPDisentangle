#!/usr/bin/env bash
# Wait for 4 Oklahoma ATE-scenario NeSI jobs, fetch, promote canonical RDS,
# then regenerate scenario compare + paper assets + HTML report.
#
# Usage (from package root):
#   bash inst/oklahoma/watch_and_publish_ate_scenarios.sh
#   bash inst/oklahoma/watch_and_publish_ate_scenarios.sh launch_jobs_boot512.txt
#   PP_NESI_POLL_SECS=1800 bash inst/oklahoma/watch_and_publish_ate_scenarios.sh
#
# Launch-record columns (tab-separated):
#   base_scen  job_id  bivariate  contrast  full_scen_id  [boot_reps]
#
# Promotion rules:
#   full_scen → for_paper_<base>.rds
#   full_scen → ate_scenarios/<base>/for_paper.rds
#   full_scen → for_paper_<full_scen>.rds (provenance alias)
#   biv_aon   → for_paper.rds (canonical paper RDS; prior backed up)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
# shellcheck source=../include/nesi_config.sh
source "$PKG_ROOT/inst/include/nesi_config.sh"

LAUNCH_FILE_ARG="${1:-launch_jobs_boot512.txt}"
POLL_SECS="${PP_NESI_POLL_SECS:-1800}"
PROMOTE_CANONICAL="${PP_PROMOTE_CANONICAL:-1}"

pp_nesi_init_config "$PKG_ROOT"
LOCAL_OK="${PP_NESI_LOCAL_OUTPUT}/oklahoma"
SCEN_DIR="${LOCAL_OK}/ate_scenarios"
mkdir -p "$SCEN_DIR" "${HOME}/.ssh/sockets"

if [[ "$LAUNCH_FILE_ARG" = /* ]]; then
  LAUNCH_FILE="$LAUNCH_FILE_ARG"
else
  LAUNCH_FILE="${SCEN_DIR}/${LAUNCH_FILE_ARG}"
fi

LOG="${LOCAL_OK}/watch_ate_scenarios_$(basename "$LAUNCH_FILE" .txt).log"

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

# Append-only logging (avoid process-substitution tee; more robust under nohup/disown).
exec >>"$LOG" 2>&1

echo "=== watch_and_publish_ate_scenarios started $(date) ==="
echo "Launch file: ${LAUNCH_FILE}"
echo "Polling every ${POLL_SECS}s; log=${LOG}"

ensure_ssh() {
  local why="${1:-ssh}"
  while ! ssh_nesi "echo SSH_OK" >/dev/null 2>&1; do
    rm -f "${HOME}/.ssh/sockets/ssh_mux_"* 2>/dev/null || true
    echo "OKLAHOMA_WATCH_NEED_SSH {\"why\":\"${why}\",\"hint\":\"run: ssh mahuika\"}"
    sleep "$POLL_SECS"
  done
  echo "SSH ready (${why}) ($(date))"
}

# Prefer local launch file; else pull from NeSI.
if [ ! -f "$LAUNCH_FILE" ]; then
  ensure_ssh "fetch-launch"
  rsync -avz "${PP_NESI_SSH}:${PP_NESI_REMOTE_OUTPUT}/oklahoma/ate_scenarios/$(basename "$LAUNCH_FILE")" \
    "$LAUNCH_FILE" || true
fi
if [ ! -f "$LAUNCH_FILE" ]; then
  echo "ERROR: launch record not found: $LAUNCH_FILE" >&2
  exit 1
fi

# Bash 3.2 compatible (macOS /bin/bash): no mapfile.
BASES=()
JOBS=()
FULLS=()
while IFS= read -r row || [ -n "$row" ]; do
  case "$row" in
    ''|\#*) continue ;;
  esac
  # New: base \t jid \t biv \t contrast \t full \t boot_reps
  # Old: scen \t jid \t biv \t contrast
  nfields="$(printf '%s\n' "$row" | awk -F'\t' '{print NF}')"
  if [ "$nfields" -ge 5 ]; then
    base="$(printf '%s\n' "$row" | cut -f1)"
    jid="$(printf '%s\n' "$row" | cut -f2)"
    full="$(printf '%s\n' "$row" | cut -f5)"
  else
    base="$(printf '%s\n' "$row" | cut -f1)"
    jid="$(printf '%s\n' "$row" | cut -f2)"
    full="$base"
  fi
  BASES+=("$base")
  JOBS+=("$jid")
  FULLS+=("$full")
  echo "Tracked: base=${base} full=${full} job=${jid}"
done < "$LAUNCH_FILE"

if [ "${#JOBS[@]}" -lt 1 ]; then
  echo "ERROR: empty launch record: $LAUNCH_FILE" >&2
  exit 1
fi

ensure_ssh "startup"

# --- wait until all jobs leave the queue ---
while true; do
  pending=0
  status_line=""
  for i in "${!JOBS[@]}"; do
    jid="${JOBS[$i]}"
    state="$(ssh_nesi "squeue -h -j ${jid} -o '%T' 2>/dev/null | head -1" || true)"
    if [ -n "$state" ]; then
      pending=$((pending + 1))
      status_line="${status_line}${FULLS[$i]}=${state} "
    else
      status_line="${status_line}${FULLS[$i]}=done "
    fi
  done
  echo "Queue status ($(date)): ${status_line}"
  if [ "$pending" -eq 0 ]; then
    echo "All ${#JOBS[@]} jobs left the queue ($(date))."
    break
  fi
  sleep "$POLL_SECS"
  ensure_ssh "poll-recover" || true
done

ensure_ssh "pre-fetch"

# Confirm sacct states
for i in "${!JOBS[@]}"; do
  jid="${JOBS[$i]}"
  sacct_line="$(ssh_nesi "sacct -n -X -j ${jid} -o JobID,State,ExitCode,Elapsed --parsable2 2>/dev/null | head -1" || true)"
  echo "sacct ${FULLS[$i]} (${jid}): ${sacct_line:-unavailable}"
  state_main="$(printf '%s' "$sacct_line" | cut -d'|' -f2)"
  if [ -n "$state_main" ] && ! echo "$state_main" | grep -Eqi 'COMPLETED|COMPLETING'; then
    echo "WARNING: ${FULLS[$i]} job ${jid} state=${state_main}"
  fi
done

# --- fetch scenario artifacts ---
echo "Fetching ate_scenarios + for_paper_* aliases ..."
rsync -avz --human-readable \
  "${PP_NESI_SSH}:${PP_NESI_REMOTE_OUTPUT}/oklahoma/ate_scenarios/" \
  "${LOCAL_OK}/ate_scenarios/"

# Job-tagged results + scenario aliases
INCLUDE_ARGS=()
for i in "${!JOBS[@]}"; do
  INCLUDE_ARGS+=(--include="oklahoma_results_job${JOBS[$i]}.rds")
  INCLUDE_ARGS+=(--include="for_paper_${FULLS[$i]}.rds")
  INCLUDE_ARGS+=(--include="${JOBS[$i]}_oklahoma_slurm.out")
  INCLUDE_ARGS+=(--include="${JOBS[$i]}_oklahoma_slurm.err")
done
rsync -avz --human-readable \
  "${PP_NESI_SSH}:${PP_NESI_REMOTE_OUTPUT}/oklahoma/" \
  "${LOCAL_OK}/" \
  "${INCLUDE_ARGS[@]}" \
  --exclude='*' || true

# --- promote to canonical names ---
promote_one() {
  local base="$1" full="$2" jid="$3"
  local src=""
  local cands=(
    "${LOCAL_OK}/ate_scenarios/${full}/for_paper.rds"
    "${LOCAL_OK}/for_paper_${full}.rds"
    "${LOCAL_OK}/oklahoma_results_job${jid}.rds"
  )
  for c in "${cands[@]}"; do
    if [ -f "$c" ]; then src="$c"; break; fi
  done
  if [ -z "$src" ]; then
    echo "ERROR: missing RDS for ${full} (job ${jid})" >&2
    return 1
  fi
  mkdir -p "${LOCAL_OK}/ate_scenarios/${base}"
  cp -p "$src" "${LOCAL_OK}/for_paper_${full}.rds"
  cp -p "$src" "${LOCAL_OK}/for_paper_${base}.rds"
  cp -p "$src" "${LOCAL_OK}/ate_scenarios/${base}/for_paper.rds"
  cp -p "$src" "${LOCAL_OK}/ate_scenarios/${full}/for_paper.rds" 2>/dev/null || \
    (mkdir -p "${LOCAL_OK}/ate_scenarios/${full}" && cp -p "$src" "${LOCAL_OK}/ate_scenarios/${full}/for_paper.rds")
  echo "Promoted ${full} <- ${src}"
  echo "  -> for_paper_${base}.rds"
  echo "  -> ate_scenarios/${base}/for_paper.rds"
  if [ "$base" = "biv_aon" ] && [ "$PROMOTE_CANONICAL" = "1" ]; then
    local for_paper="${LOCAL_OK}/for_paper.rds"
    if [ -f "$for_paper" ]; then
      local bak="${LOCAL_OK}/for_paper_backup_before_${jid}_$(date +%Y%m%d%H%M%S).rds"
      cp -p "$for_paper" "$bak"
      echo "  backed up previous for_paper.rds -> ${bak}"
    fi
    cp -p "$src" "$for_paper"
    echo "  -> for_paper.rds (canonical bivariate AoN)"
  fi
}

fail=0
for i in "${!JOBS[@]}"; do
  if ! promote_one "${BASES[$i]}" "${FULLS[$i]}" "${JOBS[$i]}"; then
    fail=1
  fi
done
if [ "$fail" -ne 0 ]; then
  echo "OKLAHOMA_SCENARIO_PUBLISH_FAILED {\"reason\":\"missing_rds\"}"
  exit 2
fi

# --- regenerate compare HTML / figure ---
echo "Rendering four-scenario compare ..."
(
  cd "$PKG_ROOT"
  Rscript inst/oklahoma/render_ate_scenario_compare.R
)

# --- paper assets from canonical for_paper.rds ---
FOR_PAPER="${LOCAL_OK}/for_paper.rds"
if [ -f "$FOR_PAPER" ]; then
  echo "Regenerating oklahoma paper assets from for_paper.rds ..."
  (
    cd "$PKG_ROOT"
    Rscript inst/oklahoma/paper/oklahoma_paper_assets.R --input "$FOR_PAPER"
  )
  PLOTS_MIRROR="${PKG_ROOT}/docs/paper/plots/oklahoma"
  mkdir -p "$PLOTS_MIRROR"
  if [ -d "${LOCAL_OK}/paper/generated/figures" ]; then
    rsync -a "${LOCAL_OK}/paper/generated/figures/" "$PLOTS_MIRROR/"
    echo "Synced figures -> ${PLOTS_MIRROR}"
  fi
fi

# --- HTML report ---
echo "Rendering oklahoma_report.html ..."
(
  cd "$PKG_ROOT/inst/oklahoma"
  Rscript render_oklahoma_report.R --results "$FOR_PAPER" --to html
)

echo "=== watch_and_publish_ate_scenarios finished $(date) ==="
echo "OKLAHOMA_SCENARIO_PUBLISH_DONE {\"launch_file\":\"${LAUNCH_FILE}\",\"for_paper\":\"${FOR_PAPER}\",\"html\":\"${LOCAL_OK}/oklahoma_report.html\",\"compare\":\"${LOCAL_OK}/oklahoma_ate_scenarios.html\"}"
