#!/bin/bash
set -euo pipefail

if [ "$#" -lt 2 ]; then
  echo "Usage: bash inst/perf/run_ab_perf.sh <baseline-ref> <candidate-ref> [--quick] [--run-decode] [--min-speedup 1.00]"
  exit 1
fi

BASE_REF="$1"
CAND_REF="$2"
shift 2

EXTRA_ARGS=()
MIN_SPEEDUP="1.00"
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --quick|--run-decode) EXTRA_ARGS+=("$1"); shift ;;
    --min-speedup) MIN_SPEEDUP="$2"; shift 2 ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
WT_DIR="$REPO_ROOT/.perf-worktrees"
BASE_WT="$WT_DIR/base"
CAND_WT="$WT_DIR/cand"
OUT_DIR="$REPO_ROOT/inst/perf/results"
BASE_JSON="$OUT_DIR/base.json"
CAND_JSON="$OUT_DIR/cand.json"
REPORT_MD="$OUT_DIR/ab_report.md"

mkdir -p "$WT_DIR" "$OUT_DIR"

cleanup() {
  git -C "$REPO_ROOT" worktree remove --force "$BASE_WT" >/dev/null 2>&1 || true
  git -C "$REPO_ROOT" worktree remove --force "$CAND_WT" >/dev/null 2>&1 || true
}
trap cleanup EXIT

echo "Creating baseline worktree: $BASE_REF"
git -C "$REPO_ROOT" worktree add --force "$BASE_WT" "$BASE_REF" >/dev/null
echo "Creating candidate worktree: $CAND_REF"
git -C "$REPO_ROOT" worktree add --force "$CAND_WT" "$CAND_REF" >/dev/null

echo "Running baseline benchmark..."
Rscript "$BASE_WT/inst/perf/benchmark_hotspots.R" --label=baseline --out="$BASE_JSON" "${EXTRA_ARGS[@]}"

echo "Running candidate benchmark..."
Rscript "$CAND_WT/inst/perf/benchmark_hotspots.R" --label=candidate --out="$CAND_JSON" "${EXTRA_ARGS[@]}"

echo "Comparing baseline vs candidate..."
Rscript "$REPO_ROOT/inst/perf/compare_benchmarks.R" \
  --base="$BASE_JSON" \
  --cand="$CAND_JSON" \
  --out="$REPORT_MD" \
  --min-speedup="$MIN_SPEEDUP"

echo "Done."
echo "Report: $REPORT_MD"
