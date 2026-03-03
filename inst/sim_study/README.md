# Simulation study

## Entry points (use these only)

- **Unified R script:** `sim_study.R` — adapts to environment via `ON_CLUSTER` and flags. Supports `--sims N` for number of simulations.
- **Shell (local or login node):** `run_sim_study.sh` — when not in SLURM, submits itself via `sbatch`; when in SLURM, runs `sim_study.R --cluster --sims N`. Use `--sims 100` (default) or `--sims 50` etc.
- **Consistency study:** `consistency_study.R` — sanity check for Hawkes fit (run locally).

## Why "Mode: LOCAL" on the cluster?

Mode is **LOCAL** when the script is run **without** SLURM (e.g. you ran `run_sim_study.sh` or `Rscript sim_study.R` on the login node). In that case `SLURM_JOB_ID` is unset, so the script uses local defaults (fewer sims/iterations, different save dir logic).

To get **cluster config** when not using sbatch (e.g. interactive or `run_sim_study.sh` on the login node), pass `--cluster`:

```bash
./inst/sim_study/run_sim_study.sh --cluster
# or
Rscript inst/sim_study/sim_study.R --cluster
```

Then you’ll see `Mode: CLUSTER` and the same settings as in the SLURM job.

## Logging and output

- **Log file:** `cluster_output/logs/sim_study_YYYYMMDD_HHMMSS.log` — timestamped messages, phase timings, and ETA estimates.
- **Results:** `cluster_output/sim_study_results_*.rds` — full results and plots.
- **SLURM logs:** `cluster_output/logs/slurm_*.out` and `slurm_*.err` when run via sbatch.

## Old files in deprecated/

Legacy and profiling scripts live in `deprecated/`. Do **not** run those. Use `sim_study.R` and `run_sim_study.sh` only.
- **On the cluster:** after `git pull`, if the old files still appear at top level, remove or move them so you only use `sim_study.R` and the run scripts above.

## Why is adaptive SEM slower than EM-style (same iter/n_props)?

With the same `iter=100` and `n_props=10`, **adaptive SEM** does **more work per dataset**:

- **EM-style:** one run of `em_style_labelling(..., iter=100, n_props=10)` per dataset → 100 iterations with 10 proposals.
- **Adaptive SEM:** an **outer** SEM loop (`N_iter` × `N_labellings`); inside it repeatedly calls the **same** `em_style_labelling` (with `adaptive_control$iter=100`, `n_props=10`) when weights are refreshed, plus it generates multiple labellings, computes importance weights, and runs parameter optimisation. So you get at least one full 100-iter labelling run per dataset, plus extra labelling generation, weight calculation, and `optim()` in the SEM loop.

So SEM is expected to take roughly 2–3× longer than EM-style for the same inner iter/n_props; that’s due to the extra SEM structure (extra adaptive steps when weights change, labelling simulation, and M-step optimisation), not a bug.

## OOM during ATE estimation

If you see `oom_kill` or `error reading from connection` during "Estimating ATEs", the parallel workers may be running out of memory. Try:

1. **Sequential ATE fallback** (lower memory, slower):
   ```bash
   ATE_SEQUENTIAL=1 ./run_sim_study.sh
   # or: export ATE_SEQUENTIAL=1 before Rscript in your SLURM script
   ```

2. **Request more memory** in your SLURM job (e.g. `#SBATCH --mem=32G` or higher).

3. **Reduce simulations** with `--sims 16` or `--sims 8` to lower memory use.
