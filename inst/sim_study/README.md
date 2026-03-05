# Simulation study

## Entry points (use these only)

- **Unified R script:** `sim_study.R` — adapts to environment via `ON_CLUSTER` and flags. Supports `--sims N` for number of simulations, `--test` for quick-test mode.
- **Shell (local or login node):** `run_sim_study.sh` — when not in SLURM, submits itself via `sbatch`; when in SLURM, runs `sim_study.R --cluster --sims N`. Use `--sims 100` (default), `--sims 50`, or `--test --sims 2` for a quick machinery check (10 iters, 2 sims/cores).
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

## Output layout

Each run is identified by its SLURM job ID (or `local_YYYYMMDD_HHMMSS` for local runs):

```
cluster_output/
  logs/<JOB_ID>/
    slurm.out        # SLURM stdout
    slurm.err        # SLURM stderr
    sim.log          # R log (timestamped messages, phase timings, ETA)
  results/<JOB_ID>.rds   # full results and plots
```

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

4. **Skip crazy-param tasks** — when SEM estimates K≥0.95 or mu>1e5, Hawkes simulations can explode. Skip those tasks:
   ```bash
   PP_SKIP_CRAZY_PARAMS=1 ./run_sim_study.sh --test
   ```

## Debugging memory and crazy params

Use **test mode** for fast iteration on the cluster:

```bash
./run_sim_study.sh --test --sims 2
```

Test mode uses reduced iters (5 EM, 5 SEM adaptive, 3 SEM outer), sequential ATE, and fewer tau sims.

**Memory logging** (log R heap at each phase):
```bash
PP_LOG_MEMORY=1 ./run_sim_study.sh --test
```

**Skip tasks with explosive params** (K≥0.95, mu>1e5) to avoid OOM:
```bash
PP_SKIP_CRAZY_PARAMS=1 ./run_sim_study.sh --test
```

Check `cluster_output/logs/<JOB_ID>/sim.log` for `[CRAZY PARAMS]` and `[MEM ...]` lines.
