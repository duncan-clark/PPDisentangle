#!/usr/bin/env bash
# Submit a PPDisentangle NeSI job from your laptop via SSH.
#
# Usage:
#   bash inst/nesi/submit.sh sim_study --sims 100
#   bash inst/nesi/submit.sh sim_study --test --sims 2
#   bash inst/nesi/submit.sh oklahoma --cores 32 --mode long
#
# Optional env:
#   PP_NESI_PUSH=1          push local git branch before remote pull/submit
#   PP_GIT_BRANCH=name      branch to pull on NeSI (default: local current branch)
#
# Configure SSH/paths in inst/nesi/nesi.env (see config.example.env).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
# shellcheck source=../include/nesi_config.sh
source "$PKG_ROOT/inst/include/nesi_config.sh"
# shellcheck source=../include/git_sync.sh
source "$PKG_ROOT/inst/include/git_sync.sh"

usage() {
  cat <<EOF
Usage: bash inst/nesi/submit.sh <sim_study|oklahoma> [run_nesi.sh args...]

Examples:
  bash inst/nesi/submit.sh sim_study --test --sims 2
  bash inst/nesi/submit.sh sim_study --sims 100
  bash inst/nesi/submit.sh oklahoma --cores 32

Run 'bash inst/nesi/submit.sh --show-config' to print resolved settings.
EOF
}

if [ "${1:-}" = "--show-config" ]; then
  pp_nesi_init_config "$PKG_ROOT"
  pp_nesi_print_config
  exit 0
fi

if [ "$#" -lt 1 ]; then
  usage
  exit 1
fi

TARGET="$1"
shift

case "$TARGET" in
  sim_study) RUNNER="inst/sim_study/run_nesi.sh" ;;
  oklahoma) RUNNER="inst/oklahoma/run_nesi.sh" ;;
  -h|--help|help) usage; exit 0 ;;
  *) echo "Unknown target: $TARGET (expected sim_study or oklahoma)" >&2; usage; exit 1 ;;
esac

pp_nesi_init_config "$PKG_ROOT"
pp_nesi_print_config
echo ""

pp_git_push_branch "$PKG_ROOT" || exit 1

REMOTE_BRANCH="${PP_GIT_BRANCH:-}"
if [ -z "$REMOTE_BRANCH" ]; then
  REMOTE_BRANCH="$(pp_git_current_branch "$PKG_ROOT")"
fi

echo "Submitting via SSH to ${PP_NESI_SSH} ..."
ssh "$PP_NESI_SSH" bash -s -- \
  "$PP_NESI_REMOTE_PKG" \
  "$PP_NESI_REMOTE_OUTPUT" \
  "$REMOTE_BRANCH" \
  "$RUNNER" \
  "$@" <<'REMOTE'
set -euo pipefail
REMOTE_PKG="$1"
REMOTE_OUTPUT="$2"
REMOTE_BRANCH="$3"
RUNNER="$4"
shift 4
cd "$REMOTE_PKG"
export PPDISENTANGLE_OUTPUT_ROOT="$REMOTE_OUTPUT"
export PP_GIT_BRANCH="$REMOTE_BRANCH"
source inst/include/git_sync.sh
pp_git_sync_repo "$REMOTE_PKG"
exec bash "$RUNNER" "$@"
REMOTE
