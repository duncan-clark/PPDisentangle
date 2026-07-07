#!/usr/bin/env bash
# Poll NeSI queue and optionally fetch results when a job finishes.
#
# Usage:
#   bash inst/nesi/wait_and_fetch.sh sim_study 1234567
#   bash inst/nesi/wait_and_fetch.sh oklahoma 1234567
#   bash inst/nesi/wait_and_fetch.sh sim_study 1234567 --no-fetch
#
# Environment:
#   PP_NESI_POLL_SECS=60     poll interval (default 60)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
# shellcheck source=../include/nesi_config.sh
source "$PKG_ROOT/inst/include/nesi_config.sh"

TARGET="${1:-}"
JOB_ID="${2:-}"
DO_FETCH=1
POLL_SECS="${PP_NESI_POLL_SECS:-60}"

if [ -z "$TARGET" ] || [ -z "$JOB_ID" ]; then
  cat <<EOF
Usage: bash inst/nesi/wait_and_fetch.sh <sim_study|oklahoma> JOB_ID [--no-fetch]
EOF
  exit 1
fi

shift 2 || true
while [ "$#" -gt 0 ]; do
  case "$1" in
    --no-fetch) DO_FETCH=0 ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
  shift
done

pp_nesi_init_config "$PKG_ROOT"
echo "Watching job $JOB_ID on ${PP_NESI_SSH} (poll every ${POLL_SECS}s) ..."

while true; do
  state="$(ssh "$PP_NESI_SSH" "squeue -h -j ${JOB_ID} 2>/dev/null | awk '{print \$5}'" || true)"
  if [ -z "$state" ]; then
    echo "Job $JOB_ID no longer in queue ($(date))."
    break
  fi
  echo "Job $JOB_ID state=$state ($(date))"
  sleep "$POLL_SECS"
done

if [ "$DO_FETCH" -eq 1 ]; then
  bash "$SCRIPT_DIR/fetch.sh" "$TARGET" "$JOB_ID"
else
  echo "Skipping fetch (--no-fetch)."
fi
