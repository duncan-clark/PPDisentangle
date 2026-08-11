#!/usr/bin/env bash
# Re-run the C/D primary profile and archive the previous current RDS so the
# slim report can show deltas (settings / params / ATEs).
#
# Steps:
#   1. If oklahoma_cd_current.rds exists, copy it -> oklahoma_cd_previous.rds
#   2. Launch the same cd-primary NeSI job (default t_trunc=90)
#   3. Print the watch command that fetches + renders oklahoma_report_cd.html
#      with previous vs current diffs
#
# Usage (from package root):
#   bash inst/oklahoma/launch_oklahoma_cd_update_nesi.sh
#   PP_SEM_T_TRUNC_DAYS=90 bash inst/oklahoma/launch_oklahoma_cd_update_nesi.sh
#   PP_SKIP_ARCHIVE=1 bash ...   # do not overwrite oklahoma_cd_previous.rds

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
# shellcheck source=../include/nesi_config.sh
source "$PKG_ROOT/inst/include/nesi_config.sh"

pp_nesi_init_config "$PKG_ROOT"
LOCAL_OK="${PP_NESI_LOCAL_OUTPUT}/oklahoma"
mkdir -p "$LOCAL_OK"

CURRENT_RDS="${LOCAL_OK}/oklahoma_cd_current.rds"
PREV_RDS="${LOCAL_OK}/oklahoma_cd_previous.rds"
PP_SKIP_ARCHIVE="${PP_SKIP_ARCHIVE:-0}"

if [ "$PP_SKIP_ARCHIVE" != "1" ] && [ -f "$CURRENT_RDS" ]; then
  cp -p "$CURRENT_RDS" "$PREV_RDS"
  echo "Archived previous current -> $PREV_RDS"
elif [ -f "$PREV_RDS" ]; then
  echo "Keeping existing previous RDS: $PREV_RDS"
else
  echo "No oklahoma_cd_current.rds yet; first update will have no diff baseline until this job is published."
fi

# Record intent for the watcher / humans
stamp="$(date +%Y%m%d_%H%M%S)"
cat > "${LOCAL_OK}/oklahoma_cd_update_${stamp}.txt" <<EOF
update_requested_at=$(date -Iseconds)
prev_rds=${PREV_RDS}
current_rds=${CURRENT_RDS}
t_trunc_days=${PP_SEM_T_TRUNC_DAYS:-90}
EOF

echo "=== Launching C/D update run ==="
bash "$SCRIPT_DIR/launch_oklahoma_cd_nesi.sh" "$@"
echo ""
echo "After submit, watch with:"
echo "  bash inst/oklahoma/watch_and_publish_cd.sh JOB_ID"
echo "That will promote oklahoma_cd_current.rds and render oklahoma_report_cd.html"
echo "with diffs against oklahoma_cd_previous.rds when present."
