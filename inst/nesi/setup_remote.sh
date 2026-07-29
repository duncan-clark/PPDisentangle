#!/usr/bin/env bash
# Run first-time NeSI setup on the remote login node via SSH.
#
# Usage: bash inst/nesi/setup_remote.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
# shellcheck source=../include/nesi_config.sh
source "$PKG_ROOT/inst/include/nesi_config.sh"

pp_nesi_init_config "$PKG_ROOT"
pp_nesi_print_config
echo ""

echo "Running remote setup on ${PP_NESI_SSH} ..."
ssh "$PP_NESI_SSH" bash -s -- \
  "$PP_NESI_REMOTE_PKG" \
  "$PP_NESI_REMOTE_OUTPUT" <<'REMOTE'
set -euo pipefail
REMOTE_PKG="$1"
REMOTE_OUTPUT="$2"
cd "$REMOTE_PKG"
export PPDISENTANGLE_OUTPUT_ROOT="$REMOTE_OUTPUT"
bash inst/sim_study/setup_nesi.sh
REMOTE
