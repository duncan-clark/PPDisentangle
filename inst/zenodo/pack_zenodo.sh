#!/usr/bin/env bash
# Pack a Zenodo-ready tarball from the local PPDisentangle-output tree.
#
# Includes paper-canonical runs only:
#   - sim_study/paper/main_5228509/   (FOR_PAPER summary + raw horizon RDS + pub assets)
#   - sim_study/paper/robustness_7568933/  (scenario RDS + summaries)
#   - sim_study/generated/{figures,robustness figures+tex fragment,tab_*.tex}
#   - oklahoma/ (application RDS + paper assets)
#
# Excludes:
#   - sim_study/local/ (e.g. 7568879 SEM refresh, standalone PDF previews)
#   - duplicate time_sweep_5228509_summary.rds (FOR_PAPER only)
#   - slurm/logs (optional clutter; kept on disk but not archived)
#
# Usage (from repo root):
#   bash inst/zenodo/pack_zenodo.sh
#   bash inst/zenodo/pack_zenodo.sh --output-root /path/to/PPDisentangle-output
#   bash inst/zenodo/pack_zenodo.sh --out /tmp/PPDisentangle-sim-outputs.tar.gz
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

OUTPUT_ROOT=""
OUT_TGZ=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --output-root)
      OUTPUT_ROOT="${2:?}"
      shift 2
      ;;
    --out)
      OUT_TGZ="${2:?}"
      shift 2
      ;;
    -h|--help)
      sed -n '1,25p' "$0"
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "$OUTPUT_ROOT" ]]; then
  if [[ -n "${PPDISENTANGLE_OUTPUT_ROOT:-}" ]]; then
    OUTPUT_ROOT="$PPDISENTANGLE_OUTPUT_ROOT"
  else
    OUTPUT_ROOT="$(cd "$PKG_ROOT/.." && pwd)/PPDisentangle-output"
  fi
fi
OUTPUT_ROOT="$(cd "$OUTPUT_ROOT" && pwd)"

if [[ -z "$OUT_TGZ" ]]; then
  stamp="$(date +%Y%m%d)"
  OUT_TGZ="$(cd "$PKG_ROOT/.." && pwd)/PPDisentangle-zenodo-outputs_${stamp}.tar.gz"
fi

STAGE="$(mktemp -d "${TMPDIR:-/tmp}/ppdisentangle-zenodo.XXXXXX")"
trap 'rm -rf "$STAGE"' EXIT

DEST="$STAGE/PPDisentangle-output"
mkdir -p "$DEST/sim_study/paper" "$DEST/sim_study/generated" "$DEST/oklahoma"

echo "Staging from: $OUTPUT_ROOT"
echo "Tarball:      $OUT_TGZ"

MAIN_DIR="$OUTPUT_ROOT/sim_study/paper/main_5228509"
ROB_DIR="$OUTPUT_ROOT/sim_study/paper/robustness_7568933"

if [[ ! -d "$MAIN_DIR" ]]; then
  echo "ERROR: missing $MAIN_DIR" >&2
  exit 1
fi
if [[ ! -d "$ROB_DIR" ]]; then
  echo "ERROR: missing $ROB_DIR" >&2
  exit 1
fi

# --- main paper run: FOR_PAPER + raws + frozen pub copies; skip duplicate summary + logs ---
rsync -a \
  --include='time_sweep_5228509_summary_FOR_PAPER.rds' \
  --include='time_sweep_5228509_tmp*.rds' \
  --include='time_sweep_5228509_all_nothing_grouped_data*' \
  --include='time_sweep_5228509_core_true_control_pub*' \
  --include='simulated_hawkes_hawkes_process.pdf' \
  --exclude='*' \
  "$MAIN_DIR/" "$DEST/sim_study/paper/main_5228509/"

# --- robustness paper run: RDS + CSV summaries; skip logs/slurm ---
rsync -a \
  --include='robustness_7568933_*.rds' \
  --include='robustness_7568933_*.csv' \
  --exclude='*' \
  "$ROB_DIR/" "$DEST/sim_study/paper/robustness_7568933/"

# --- generated paper figures ---
if [[ -d "$OUTPUT_ROOT/sim_study/generated" ]]; then
  rsync -a \
    --exclude='figures_sem5000_*/' \
    --exclude='*.aux' \
    --exclude='*.log' \
    --exclude='*.out' \
    --exclude='*.fls' \
    --exclude='*.fdb_latexmk' \
    --exclude='robustness_standalone*' \
    --exclude='tex/' \
    "$OUTPUT_ROOT/sim_study/generated/" "$DEST/sim_study/generated/"
fi

# --- oklahoma ---
if [[ -d "$OUTPUT_ROOT/oklahoma" ]]; then
  rsync -a \
    --exclude='.DS_Store' \
    "$OUTPUT_ROOT/oklahoma/" "$DEST/oklahoma/"
fi

cat > "$DEST/README.md" <<'EOF'
# PPDisentangle paper outputs (Zenodo)

Companion data archive for the PPDisentangle software repository.

## Canonical jobs

| Component | Path | Role |
|-----------|------|------|
| Main simulation | `sim_study/paper/main_5228509/` | Paper simulation figures |
| Illustrative realisation | `.../simulated_hawkes_hawkes_process.pdf` | `fig:pp_realiz` |
| Robustness appendix | `sim_study/paper/robustness_7568933/` | Appendix robustness figures |
| Oklahoma application | `oklahoma/for_paper.rds` | Application figures/tables |

## Layout

```text
PPDisentangle-output/
  README.md
  sim_study/
    paper/
      main_5228509/
        time_sweep_5228509_summary_FOR_PAPER.rds
        time_sweep_5228509_tmp*.rds
        simulated_hawkes_hawkes_process.pdf
      robustness_7568933/
        robustness_7568933_*.rds
        robustness_7568933_*_summary.csv
        robustness_7568933_manifest.csv
    generated/
      figures/   # results_* + simulated_hawkes_*.pdf
      robustness/figures/
      robustness/simulation_robustness_appendix.tex
  oklahoma/
    for_paper.rds
    paper/generated/
```

## Reproduce figures

```bash
export PPDISENTANGLE_OUTPUT_ROOT=/path/to/PPDisentangle-output
Rscript inst/zenodo/reproduce_paper_figures.R
```
EOF

echo "Staged size:"
du -sh "$DEST" "$DEST/sim_study" "$DEST/oklahoma" 2>/dev/null || true

mkdir -p "$(dirname "$OUT_TGZ")"
tar -C "$STAGE" -czf "$OUT_TGZ" PPDisentangle-output
echo "Wrote $OUT_TGZ"
ls -lh "$OUT_TGZ"
