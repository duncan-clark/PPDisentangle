#!/usr/bin/env bash
# Fetch 4 ATE scenario RDS files from NeSI and render comparison HTML.
# For full promote + report regen after a boot512 campaign, prefer:
#   bash inst/oklahoma/watch_and_publish_ate_scenarios.sh launch_jobs_boot512.txt
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
# shellcheck source=../include/nesi_config.sh
source "$PKG_ROOT/inst/include/nesi_config.sh"
pp_nesi_init_config "$PKG_ROOT"

LOCAL_OK="${PP_NESI_LOCAL_OUTPUT}/oklahoma"
REMOTE_OK="${PP_NESI_REMOTE_OUTPUT}/oklahoma"
mkdir -p "${LOCAL_OK}/ate_scenarios"

echo "Fetching ate_scenarios/ and for_paper_*.rds ..."
rsync -avz --human-readable \
  "${PP_NESI_SSH}:${REMOTE_OK}/ate_scenarios/" \
  "${LOCAL_OK}/ate_scenarios/"

rsync -avz --human-readable \
  "${PP_NESI_SSH}:${REMOTE_OK}/" \
  "${LOCAL_OK}/" \
  --include='for_paper_univ_aon.rds' \
  --include='for_paper_univ_obs.rds' \
  --include='for_paper_biv_aon.rds' \
  --include='for_paper_biv_obs.rds' \
  --include='for_paper_*_boot512.rds' \
  --include='for_paper_*_boot256.rds' \
  --include='oklahoma_results_job*.rds' \
  --exclude='*' || true

echo "Rendering comparison HTML ..."
(
  cd "$PKG_ROOT"
  Rscript inst/oklahoma/render_ate_scenario_compare.R
)

echo "Done."
echo "  HTML: ${LOCAL_OK}/oklahoma_ate_scenarios.html"
echo "  Fig:  ${LOCAL_OK}/paper/generated/figures/ATE_scenarios_compare.pdf"
