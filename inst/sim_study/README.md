# Simulation study

## Entry points

- **R script:** `sim_study.R` — adapts to environment via `ON_CLUSTER` and flags. Supports `--sims N`, `--test`.
- **NeSI cluster:** `run_nesi.sh` — submits via `sbatch`; when in SLURM, runs `sim_study.R --cluster --sims N`. Use `--mode long` (default long profile) or `--mode test` (lightweight path check), and optionally override with `--sims N`.
- **Laptop → NeSI → laptop:** `inst/nesi/submit.sh`, `fetch.sh`, `wait_and_fetch.sh` — see [`../nesi/README.md`](../nesi/README.md).

```bash
cd /path/to/PPDisentangle
bash inst/sim_study/run_nesi.sh --sims 100
bash inst/sim_study/run_nesi.sh --test --sims 2   # quick smoke test
bash inst/sim_study/run_nesi.sh --mode test       # lightweight mode preset
bash inst/sim_study/run_nesi.sh --mode long       # long-run mode preset
```

For robustness runs, `--sims` controls the number of independent replications.
Use `--cpus` only when you want to request a larger SLURM allocation without
increasing replications. The robustness launcher uses extra CPUs by running
multiple scenarios concurrently. By default it uses
`floor(SLURM_CPUS_PER_TASK / sims)` scenario workers; override with
`--scenario-workers`. SEM worker parallelism inside each scenario is naturally
capped by the number of replications. Extra CPUs within a scenario can also help
ATE/off-support work via `--ate-workers`.

## Output layout

Each run is identified by its SLURM job ID (or `local_YYYYMMDD_HHMMSS` for local runs).
All paths are under **`../PPDisentangle-output/sim_study/`** (sibling to this repo):

```
PPDisentangle-output/sim_study/
  <JOB_ID>.rds           # full results and plots
  <JOB_ID>.log           # R log
  <JOB_ID>_slurm.out     # SLURM stdout
  <JOB_ID>_slurm.err     # SLURM stderr
```

Override the output root with `PPDISENTANGLE_OUTPUT_ROOT`. See [`../OUTPUT.md`](../OUTPUT.md).

## Publication outputs

Paper-canonical main run: **`time_sweep_5228509`** (`*_summary_FOR_PAPER.rds` +
raw horizon `tmp*.rds`). Robustness appendix: **`robustness_7568933`**.

To regenerate all paper figures from a Zenodo unpack, see
[`../zenodo/README.md`](../zenodo/README.md):

```bash
Rscript inst/zenodo/reproduce_paper_figures.R
```

(`robustness.pdf` is a local preview only and is not part of Zenodo.)

The publication plotting helper reads a time-sweep summary from
`PPDisentangle-output/sim_study/` and writes paper figures/tables to
`PPDisentangle-output/sim_study/generated/`:

```bash
Rscript inst/sim_study/plot_time_sweep_publication.R
```

Default outputs:

```text
PPDisentangle-output/sim_study/generated/figures/results_true_control.{pdf,png}
PPDisentangle-output/sim_study/generated/figures/results_estimated_control.{pdf,png}
PPDisentangle-output/sim_study/generated/tab_sim_time_sweep_param_tables.tex
```

## Robustness suite

The robustness launcher runs the complete robustness suite in one job: grids
over K separation, signal-to-noise,
a $4\times 4$ spatiotemporal kernel misspecification matrix
(temporal × spatial, each exponential or power-law, under sim and fit),
pretreatment-informed treatment assignment
(highest count, lowest count, count propensity, contiguous AOI, Voronoi random),
off-support allocation contrasts (including the
all-or-nothing DTAITE as the global contrast), label recovery,
the forward-simulation decay diagnostics, binary-covariate effect modification,
and transport across balanced allocation geometries. It writes per-scenario
result objects plus structured-study results, summary CSV/RDS files, and
paper-ready figures/LaTeX fragments.

**K-separation grid (default):** control `K = 0.8`, treated `K ∈ {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7}` (treatment lowers `K`) with μ calibrated to `target_points`. Fixed-K scenarios use `(0.8, 0.2)`. Override with `PP_K_VALUES` or `--k-values`.

**SNR grid (default):** `μ_scale ∈ {0.25, 0.5, 1, 1.5, 2}` at fixed `K = (0.8, 0.2)`. Override with `PP_MU_SCALES` or `--mu-scales`.

**Pretreatment assignment family (`pretreatment_assignment`):** five rules at
`K = (0.8, 0.2)`, each treating 50% of cells from a reference pre-treatment
catalogue — `highest_count_50pct`, `lowest_count_50pct`,
`count_propensity_50pct`, `contiguous_aoi_50pct`, `voronoi_random_50pct`
(`PP_TREATMENT_ASSIGNMENT` / scenario column `treatment_assignment`).
Legacy `--scenario-set high_count_assignment` is aliased to this family.

```bash
Rscript inst/sim_study/sim_study_robustness.R --sims 32 --target-points 2500
```

Default generated assets:

```text
PPDisentangle-output/sim_study/generated/robustness/figures/
  robustness_all_scenarios_label_recovery.{pdf,png}
  robustness_parameter_sweep_support_contrasts.{pdf,png}
  robustness_spatiotemporal_kernel_mismatch_heatmap.{pdf,png}
  robustness_pretreatment_assignment_support_contrasts.{pdf,png}
  robustness_decay_validation.{pdf,png}
  robustness_temporal_decay_validation.{pdf,png}
PPDisentangle-output/sim_study/generated/robustness/simulation_robustness_appendix.tex
```

Effect modification and geometry transport are automatically run after the
standard scenario grid by `--robustness`. They can still be rerun independently
for inspection or debugging:

```bash
Rscript inst/sim_study/sim_study_structured_robustness.R --study both --sims 32
Rscript inst/sim_study/sim_study_structured_robustness.R --study both --pilot-only
bash inst/sim_study/run_nesi.sh --structured-inspect
bash inst/sim_study/run_nesi.sh --structured-robustness --structured-study both --sims 32
```

Use `--skip-structured` only for a deliberate grid-only timing probe.

They add two composite figures:

```text
generated/robustness/figures/robustness_effect_modification.{pdf,png}
generated/robustness/figures/robustness_geometry_transport.{pdf,png}
generated/robustness/structured/*.rds
```

Label recovery is summarised once across all scenarios (naive vs SEM). Bias support
figures report naive vs SEM only (oracle omitted) for the all-or-nothing DTAITE.
Standalone absolute ATE error figures are not generated.

Copy `simulation_robustness_appendix.tex` and `figures/*.pdf` to
`plots/sim_study/robustness/` on Overleaf, then
`\input{plots/sim_study/robustness/simulation_robustness_appendix.tex}` inside
`\section{Additional simulation results}`.

Local PDF preview:

```bash
cd PPDisentangle-output/sim_study/generated/robustness
pdflatex -jobname=robustness '\def\robustnessstandalone{}\input{simulation_robustness_appendix.tex}'
```

Small NeSI inspection run, intended to finish quickly enough to inspect before
an overnight run:

```bash
bash inst/sim_study/run_nesi.sh --robustness-inspect
```

Quick-mode timing probe (one scenario from each of three standard families;
use this only for grid timing, not as the complete appendix):

```bash
bash inst/sim_study/run_nesi.sh --robustness-quick-probe
```

Full NeSI robustness run (33 standard scenarios followed by both structured
studies, then a combined replot). Recommended resources: **64 CPUs / 48h /
128G** (better backfill than 100-core requests):

```bash
PP_NESI_PUSH=1 bash inst/nesi/submit.sh sim_study \
  --robustness --mode long \
  --sims 32 --cpus 64 --time 48:00:00 \
  --target-points 2500 --sem-inner 2000 --skip-ate-tau
```

Larger-node alternative (100 CPUs / default 72h wall):

```bash
bash inst/sim_study/run_nesi.sh --robustness --mode long \
  --sims 32 --cpus 100 --target-points 2500 --skip-ate-tau
```

Quick NeSI run for the full `robustness.pdf` appendix (33 standard scenarios,
both structured studies with reduced forward Monte Carlo, no estimated
`tau_i`):

```bash
bash inst/nesi/submit.sh sim_study \
  --robustness --mode quick \
  --sims 32 --cpus 100 --target-points 2500
```

After fetch, replot and compile locally:

```bash
Rscript inst/sim_study/sim_study_robustness.R --replot robustness_<JOBID>
cd PPDisentangle-output/sim_study/generated/robustness
pdflatex -jobname=robustness '\def\robustnessstandalone{}\input{simulation_robustness_appendix.tex}'
```

Partial timing probe (subset only — not sufficient for the full PDF):

```bash
bash inst/sim_study/run_nesi.sh --robustness --mode quick \
  --scenario-set pretreatment_assignment,snr_scale --target-points 2500 \
  --sims 32 --cpus 100 --ate-workers 32 --skip-ate-tau --skip-structured
```

## OOM during ATE estimation

1. **Sequential ATE** (lower memory, slower):
   ```bash
   ATE_SEQUENTIAL=1 bash inst/sim_study/run_nesi.sh --sims 50
   ```

2. **Request more memory** in `run_nesi.sh` (`#SBATCH --mem=...`).

3. **Skip crazy-param tasks** (K≥0.95, mu>1e5):
   ```bash
   PP_SKIP_CRAZY_PARAMS=1 bash inst/sim_study/run_nesi.sh --test
   ```

## Debugging

**Memory logging:**
```bash
PP_LOG_MEMORY=1 bash inst/sim_study/run_nesi.sh --test
```

**Skip explosive params:**
```bash
PP_SKIP_CRAZY_PARAMS=1 bash inst/sim_study/run_nesi.sh --test
```

Check `PPDisentangle-output/sim_study/<JOB_ID>.log` for `[CRAZY PARAMS]` and `[MEM ...]` lines.
