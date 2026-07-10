#!/usr/bin/env bash
# Fetch NeSI results into the local PPDisentangle-output directory.
#
# Usage:
#   bash inst/nesi/fetch.sh sim_study              # all sim_study outputs
#   bash inst/nesi/fetch.sh sim_study 1234567      # one SLURM job id
#   bash inst/nesi/fetch.sh oklahoma               # all oklahoma outputs
#   bash inst/nesi/fetch.sh oklahoma 1234567       # job-tagged oklahoma files
#   bash inst/nesi/fetch.sh all                    # sim_study + oklahoma
#
# Uses rsync over SSH. Configure paths in inst/nesi/nesi.env.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
# shellcheck source=../include/nesi_config.sh
source "$PKG_ROOT/inst/include/nesi_config.sh"

usage() {
  cat <<EOF
Usage: bash inst/nesi/fetch.sh <sim_study|oklahoma|all> [JOB_ID]

Examples:
  bash inst/nesi/fetch.sh sim_study
  bash inst/nesi/fetch.sh sim_study 1234567
  bash inst/nesi/fetch.sh oklahoma 1234567
  bash inst/nesi/fetch.sh all
EOF
}

pp_nesi_rsync_dir() {
  local remote_subdir="$1"
  local local_dir="$PP_NESI_LOCAL_OUTPUT/$remote_subdir"
  mkdir -p "$local_dir"
  echo "rsync ${PP_NESI_SSH}:${PP_NESI_REMOTE_OUTPUT}/${remote_subdir}/ -> ${local_dir}/"
  rsync -avz --progress --human-readable \
    "${PP_NESI_SSH}:${PP_NESI_REMOTE_OUTPUT}/${remote_subdir}/" \
    "${local_dir}/"
}

pp_nesi_fetch_sim_study() {
  local job_id="${1:-}"
  if [ -z "$job_id" ]; then
    pp_nesi_rsync_dir "sim_study"
    return 0
  fi

  local local_dir="$PP_NESI_LOCAL_OUTPUT/sim_study"
  mkdir -p "$local_dir" "$local_dir/generated"
  local remote="$PP_NESI_REMOTE_OUTPUT/sim_study"
  local patterns=(
    "${job_id}.rds"
    "${job_id}.log"
    "${job_id}_slurm.out"
    "${job_id}_slurm.err"
    "time_sweep_${job_id}_summary.rds"
  )
  local pat
  for pat in "${patterns[@]}"; do
    rsync -avz --progress --human-readable \
      "${PP_NESI_SSH}:${remote}/${pat}" \
      "${local_dir}/" 2>/dev/null || true
  done
  rsync -avz --progress --human-readable \
    "${PP_NESI_SSH}:${remote}/time_sweep_${job_id}_"* \
    "${local_dir}/" 2>/dev/null || true
  # Robustness suite: scenario RDS/CSV/logs + shared summaries
  rsync -avz --progress --human-readable \
    "${PP_NESI_SSH}:${remote}/robustness_${job_id}_"* \
    "${local_dir}/" 2>/dev/null || true
  rsync -avz --progress --human-readable \
    "${PP_NESI_SSH}:${remote}/generated/" \
    "${local_dir}/generated/" 2>/dev/null || true
}

pp_nesi_fetch_oklahoma() {
  local job_id="${1:-}"
  if [ -z "$job_id" ]; then
    pp_nesi_rsync_dir "oklahoma"
    return 0
  fi

  local local_dir="$PP_NESI_LOCAL_OUTPUT/oklahoma"
  mkdir -p "$local_dir" "$local_dir/plots" "$local_dir/paper/generated"
  local remote="$PP_NESI_REMOTE_OUTPUT/oklahoma"

  rsync -avz --progress --human-readable \
    "${PP_NESI_SSH}:${remote}/" \
    "${local_dir}/" \
    --include="*job${job_id}*" \
    --include="${job_id}_oklahoma_slurm.out" \
    --include="${job_id}_oklahoma_slurm.err" \
    --include="plots/***" \
    --include="paper/generated/***" \
    --exclude="*" 2>/dev/null || true
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
JOB_ID="${2:-}"

case "$TARGET" in
  -h|--help|help) usage; exit 0 ;;
esac

pp_nesi_init_config "$PKG_ROOT"
pp_nesi_print_config
echo ""

case "$TARGET" in
  sim_study) pp_nesi_fetch_sim_study "$JOB_ID" ;;
  oklahoma) pp_nesi_fetch_oklahoma "$JOB_ID" ;;
  all)
    pp_nesi_fetch_sim_study "$JOB_ID"
    pp_nesi_fetch_oklahoma "$JOB_ID"
    ;;
  *)
    echo "Unknown target: $TARGET" >&2
    usage
    exit 1
    ;;
esac

echo ""
echo "Local output: $PP_NESI_LOCAL_OUTPUT"
