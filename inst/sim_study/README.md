# Simulation study

## Entry points

- **R script:** `sim_study.R` — adapts to environment via `ON_CLUSTER` and flags. Supports `--sims N`, `--test`.
- **NeSI cluster:** `run_nesi.sh` — submits via `sbatch`; when in SLURM, runs `sim_study.R --cluster --sims N`. Use `--mode long` (default long profile) or `--mode test` (lightweight path check), and optionally override with `--sims N`.

```bash
cd /path/to/PPDisentangle
bash inst/sim_study/run_nesi.sh --sims 100
bash inst/sim_study/run_nesi.sh --test --sims 2   # quick smoke test
bash inst/sim_study/run_nesi.sh --mode test       # lightweight mode preset
bash inst/sim_study/run_nesi.sh --mode long       # long-run mode preset
```

## Output layout

Each run is identified by its SLURM job ID (or `local_YYYYMMDD_HHMMSS` for local runs):

```
output/sim_study/
  <JOB_ID>.rds           # full results and plots
  <JOB_ID>.log           # R log
  <JOB_ID>_slurm.out     # SLURM stdout
  <JOB_ID>_slurm.err     # SLURM stderr
```

Large local or cluster run outputs are ignored by git. Historical local/HPC
results should stay in `output/sim_study/` or `cluster_output/`, not in the
package source tree.

## Publication artifacts

The publication plotting helper reads a time-sweep summary from
`output/sim_study/` and writes versioned paper artifacts to
`artifacts/sim_study/generated/`:

```bash
Rscript inst/sim_study/plot_time_sweep_publication.R
```

Default outputs:

```text
artifacts/sim_study/generated/figures/results_true_control.{pdf,png}
artifacts/sim_study/generated/figures/results_estimated_control.{pdf,png}
artifacts/sim_study/generated/tab_sim_time_sweep_param_tables.tex
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

Check `output/sim_study/<JOB_ID>.log` for `[CRAZY PARAMS]` and `[MEM ...]` lines.
