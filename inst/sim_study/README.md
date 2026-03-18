# Simulation study

## Entry points

- **R script:** `sim_study.R` — adapts to environment via `ON_CLUSTER` and flags. Supports `--sims N`, `--test`.
- **NeSI cluster:** `run_nesi.sh` — submits via `sbatch`; when in SLURM, runs `sim_study.R --cluster --sims N`. Use `--sims 100` (default), `--sims 50`, or `--test --sims 2` for a quick check.

```bash
cd /path/to/PPDisentangle
bash inst/sim_study/run_nesi.sh --sims 100
bash inst/sim_study/run_nesi.sh --test --sims 2   # quick smoke test
```

## Output layout

Each run is identified by its SLURM job ID (or `local_YYYYMMDD_HHMMSS` for local runs):

```
inst/sim_study/output/
  <JOB_ID>.rds           # full results and plots
  <JOB_ID>.log           # R log
  <JOB_ID>_slurm.out     # SLURM stdout
  <JOB_ID>_slurm.err     # SLURM stderr
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

Check `inst/sim_study/output/<JOB_ID>.log` for `[CRAZY PARAMS]` and `[MEM ...]` lines.
