#!/usr/bin/env bash
# Infrastructure smoke test: trivial SLURM job + rsync back to local output.
# Does NOT run SEM, sim_study, or Oklahoma analysis.
#
# Usage:
#   bash inst/nesi/smoke_test.sh
#   bash inst/nesi/smoke_test.sh --no-fetch    # submit only
#   bash inst/nesi/smoke_test.sh --no-wait     # submit + fetch without polling
#
# Requires: SSH to NeSI (PP_NESI_SSH in nesi.env), rsync.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
# shellcheck source=../include/nesi_config.sh
source "$PKG_ROOT/inst/include/nesi_config.sh"

DO_WAIT=1
DO_FETCH=1
POLL_SECS="${PP_NESI_POLL_SECS:-10}"
SMOKE_SUBDIR="infra_smoke"

while [ "$#" -gt 0 ]; do
  case "$1" in
    --no-wait) DO_WAIT=0; shift ;;
    --no-fetch) DO_FETCH=0; shift ;;
    -h|--help)
      echo "Usage: bash inst/nesi/smoke_test.sh [--no-wait] [--no-fetch]"
      exit 0
      ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

pp_nesi_init_config "$PKG_ROOT"
: "${PP_NESI_ACCOUNT:=uoo04008}"

REMOTE_SMOKE="${PP_NESI_REMOTE_OUTPUT}/${SMOKE_SUBDIR}"
LOCAL_SMOKE="${PP_NESI_LOCAL_OUTPUT}/${SMOKE_SUBDIR}"

pp_nesi_print_config
echo "  PP_NESI_ACCOUNT          = ${PP_NESI_ACCOUNT}"
echo "  Remote smoke dir         = ${REMOTE_SMOKE}"
echo "  Local smoke dir          = ${LOCAL_SMOKE}"
echo ""

echo "=== Step 1: submit trivial SLURM job on NeSI ==="
SUBMIT_OUT="$(ssh "$PP_NESI_SSH" bash -s -- \
  "$PP_NESI_ACCOUNT" \
  "$REMOTE_SMOKE" <<'REMOTE'
set -euo pipefail
ACCOUNT="$1"
SMOKE_DIR="$2"
mkdir -p "$SMOKE_DIR"
WRAP=$(cat <<EOF
set -euo pipefail
MARKER="${SMOKE_DIR}/\${SLURM_JOB_ID}_marker.txt"
{
  echo "PPDisentangle infra smoke test"
  echo "job_id=\${SLURM_JOB_ID}"
  echo "host=\$(hostname)"
  echo "date=\$(date -Is)"
  echo "output_root=${SMOKE_DIR}"
} > "\$MARKER"
echo "Wrote \$MARKER"
EOF
)
JOB_ID=$(sbatch --parsable \
  --account="$ACCOUNT" \
  --job-name=PPDis_smoke \
  --time=00:02:00 \
  --mem=512M \
  --cpus-per-task=1 \
  --output="${SMOKE_DIR}/%j_slurm.out" \
  --error="${SMOKE_DIR}/%j_slurm.err" \
  --wrap="$WRAP")
echo "JOB_ID=${JOB_ID}"
REMOTE
)"

echo "$SUBMIT_OUT"
JOB_ID="$(printf '%s\n' "$SUBMIT_OUT" | sed -n 's/^JOB_ID=//p' | tail -1)"
if [ -z "$JOB_ID" ]; then
  echo "ERROR: could not parse JOB_ID from remote submit output." >&2
  exit 1
fi
echo "Submitted job ${JOB_ID}"

if [ "$DO_WAIT" -eq 1 ]; then
  echo ""
  echo "=== Step 2: wait for job ${JOB_ID} ==="
  while true; do
    STATE="$(ssh "$PP_NESI_SSH" "sacct -n -j ${JOB_ID} --format=State --noheader 2>/dev/null | head -1 | awk '{print \$1}'" || true)"
    STATE="${STATE:-PENDING}"
    echo "$(date +%H:%M:%S) state=${STATE}"
    case "$STATE" in
      COMPLETED) break ;;
      FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL)
        echo "ERROR: job ended with state=${STATE}" >&2
        ssh "$PP_NESI_SSH" "tail -20 '${REMOTE_SMOKE}/${JOB_ID}_slurm.err' 2>/dev/null || true"
        exit 1
        ;;
    esac
    sleep "$POLL_SECS"
  done
fi

if [ "$DO_FETCH" -eq 1 ]; then
  echo ""
  echo "=== Step 3: rsync results to local machine ==="
  mkdir -p "$LOCAL_SMOKE"
  rsync -avz --progress --human-readable \
    "${PP_NESI_SSH}:${REMOTE_SMOKE}/" \
    "${LOCAL_SMOKE}/"

  MARKER="${LOCAL_SMOKE}/${JOB_ID}_marker.txt"
  if [ -f "$MARKER" ]; then
    echo ""
    echo "=== SUCCESS: local marker file found ==="
    cat "$MARKER"
    echo ""
    echo "Infra smoke test passed (SLURM submit + rsync)."
  else
    echo "ERROR: expected marker not found locally: $MARKER" >&2
    echo "Remote listing:" >&2
    ssh "$PP_NESI_SSH" "ls -la '${REMOTE_SMOKE}/'" >&2 || true
    exit 1
  fi
else
  echo ""
  echo "Fetch skipped. When ready:"
  echo "  rsync -avz ${PP_NESI_SSH}:${REMOTE_SMOKE}/ ${LOCAL_SMOKE}/"
fi
