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
increasing replications. Extra CPUs help most in ATE/off-support work via
`--ate-workers`; SEM worker parallelism is naturally capped by the number of
replications.

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

The robustness launcher runs grids over K separation, signal-to-noise,
kernel misspecification, off-support allocation contrasts, label recovery,
and the forward-simulation decay diagnostic. It writes per-scenario result
objects plus summary CSV/RDS files and paper-ready figures/LaTeX fragments.

**K-separation grid (default):** control `K = 0.8`, treated `K ∈ {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7}` (treatment lowers `K`) with μ calibrated to `target_points`. Fixed-K scenarios use `(0.8, 0.2)`. Override with `PP_K_VALUES` or `--k-values`.

```bash
Rscript inst/sim_study/sim_study_robustness.R --sims 32 --target-points 2500
```

Default generated assets:

```text
PPDisentangle-output/sim_study/generated/robustness/figures/
  robustness_k_separation_label_recovery.{pdf,png}
  robustness_k_separation_ate_error.{pdf,png}
  robustness_k_separation_support_contrasts.{pdf,png}
  robustness_kernel_mismatch_label_recovery.{pdf,png}
  robustness_kernel_mismatch_ate_error.{pdf,png}
  robustness_kernel_mismatch_support_contrasts.{pdf,png}
  robustness_decay_validation.{pdf,png}
PPDisentangle-output/sim_study/generated/robustness/simulation_robustness_appendix.tex
```

Copy `simulation_robustness_appendix.tex` and `figures/*.pdf` to
`plots/sim_study/robustness/` on Overleaf, then
`\input{plots/sim_study/robustness/simulation_robustness_appendix.tex}` inside
`\section{Additional simulation results}`.

Local PDF preview:

```bash
cd PPDisentangle-output/sim_study/generated/robustness
pdflatex -jobname=robustness_standalone '\def\robustnessstandalone{}\input{simulation_robustness_appendix.tex}'
```

Small NeSI inspection run, intended to finish quickly enough to inspect before
an overnight run:

```bash
bash inst/sim_study/run_nesi.sh --robustness-inspect
```

Full NeSI robustness run:

```bash
bash inst/sim_study/run_nesi.sh --robustness --mode long --sims 32 --target-points 2500
```

Timing probe with 32 replications and a larger CPU allocation:

```bash
bash inst/sim_study/run_nesi.sh --robustness --mode quick \
  --scenario-set high_count_assignment --target-points 2500 \
  --sims 32 --cpus 100 --ate-workers 32 --skip-ate-tau
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
