#!/usr/bin/env bash
# Show NeSI queue status for your user via SSH.
#
# Usage: bash inst/nesi/status.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
# shellcheck source=../include/nesi_config.sh
source "$PKG_ROOT/inst/include/nesi_config.sh"

pp_nesi_init_config "$PKG_ROOT"
echo "Queue for ${PP_NESI_SSH}:"
ssh "$PP_NESI_SSH" "squeue -u \${USER:-$(printf '%q' "$PP_NESI_USER")}"
