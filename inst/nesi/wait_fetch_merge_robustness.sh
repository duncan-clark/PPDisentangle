#!/usr/bin/env bash
# Wait for a NeSI robustness job, fetch results, merge a scenario family into a
# local paper archive, replot, and rebuild robustness.pdf.
#
# Usage:
#   bash inst/nesi/wait_fetch_merge_robustness.sh 7839871
#   bash inst/nesi/wait_fetch_merge_robustness.sh 7839871 \
#     --archive robustness_merged_tcal --family k_spatial_range
#
# Environment:
#   PP_NESI_POLL_SECS=120   poll interval (default 120)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
# shellcheck source=../include/nesi_config.sh
source "$PKG_ROOT/inst/include/nesi_config.sh"

JOB_ID="${1:-}"
ARCHIVE="robustness_merged_tcal"
FAMILY="k_spatial_range"
POLL_SECS="${PP_NESI_POLL_SECS:-120}"
DO_PDF=1

if [ -z "$JOB_ID" ]; then
  cat <<EOF
Usage: bash inst/nesi/wait_fetch_merge_robustness.sh JOB_ID [options]

Options:
  --archive NAME   paper archive basename (default: robustness_merged_tcal)
  --family NAME    scenario family to merge (default: k_spatial_range)
  --no-pdf         skip pdflatex preview rebuild
EOF
  exit 1
fi
shift

while [ "$#" -gt 0 ]; do
  case "$1" in
    --archive) ARCHIVE="$2"; shift 2 ;;
    --family) FAMILY="$2"; shift 2 ;;
    --no-pdf) DO_PDF=0; shift ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

pp_nesi_init_config "$PKG_ROOT"
cd "$PKG_ROOT"

echo "[worker] Watching job $JOB_ID on ${PP_NESI_SSH} (poll every ${POLL_SECS}s)"
echo "[worker] On completion: fetch → merge ${FAMILY} into ${ARCHIVE} → replot → pdf"

while true; do
  state="$(ssh "$PP_NESI_SSH" "squeue -h -j ${JOB_ID} -o '%T' 2>/dev/null | head -1" || true)"
  if [ -z "$state" ]; then
    echo "[worker] Job $JOB_ID no longer in queue ($(date))"
    break
  fi
  echo "[worker] Job $JOB_ID state=$state ($(date))"
  sleep "$POLL_SECS"
done

# Confirm exit via sacct when available (non-fatal if sacct unavailable).
exit_code="$(ssh "$PP_NESI_SSH" "sacct -n -X -j ${JOB_ID} -o State,ExitCode --parsable2 2>/dev/null | head -1" || true)"
echo "[worker] sacct: ${exit_code:-unknown}"
if echo "${exit_code}" | grep -Eqi 'FAILED|CANCELLED|TIMEOUT|NODE_FAIL|OUT_OF_MEMORY'; then
  echo "[worker] ERROR: job did not complete successfully; refusing merge." >&2
  printf 'AGENT_MERGE_FAILED {"job":"%s","sacct":"%s"}\n' "$JOB_ID" "$exit_code"
  exit 2
fi

echo "[worker] Fetching sim_study job $JOB_ID ..."
bash "$SCRIPT_DIR/fetch.sh" sim_study "$JOB_ID"

n_rds="$(ls -1 "${PP_NESI_LOCAL_OUTPUT}/sim_study/robustness_${JOB_ID}_${FAMILY}_"*.rds 2>/dev/null | wc -l | tr -d ' ')"
echo "[worker] Found ${n_rds} local RDS for family=${FAMILY}"
if [ "${n_rds}" -lt 1 ]; then
  echo "[worker] ERROR: no fetched RDS for ${FAMILY}; refusing merge." >&2
  printf 'AGENT_MERGE_FAILED {"job":"%s","reason":"no_rds"}\n' "$JOB_ID"
  exit 3
fi

echo "[worker] Merging into paper/${ARCHIVE} ..."
Rscript inst/sim_study/merge_robustness_job_into_archive.R \
  --job "$JOB_ID" \
  --archive "$ARCHIVE" \
  --family "$FAMILY"

echo "[worker] Replotting ${ARCHIVE} ..."
Rscript inst/sim_study/sim_study_robustness.R --replot "$ARCHIVE"

if [ "$DO_PDF" -eq 1 ]; then
  gen_dir="${PP_NESI_LOCAL_OUTPUT}/sim_study/generated/robustness"
  if [ -d "$gen_dir" ]; then
    echo "[worker] Building robustness.pdf ..."
    (
      cd "$gen_dir"
      pdflatex -interaction=nonstopmode -jobname=robustness \
        '\def\robustnessstandalone{}\input{simulation_robustness_appendix.tex}' \
        >/dev/null
    )
    echo "[worker] Wrote ${gen_dir}/robustness.pdf"
  else
    echo "[worker] WARN: generated/robustness missing; skipped PDF"
  fi
fi

echo "[worker] Done ($(date)): job=${JOB_ID} family=${FAMILY} archive=${ARCHIVE} n_rds=${n_rds}"
echo "AGENT_MERGE_DONE {\"job\":\"${JOB_ID}\",\"family\":\"${FAMILY}\",\"archive\":\"${ARCHIVE}\",\"n_rds\":${n_rds}}"
