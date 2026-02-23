# Simulation study

## Entry points (use these only)

- **Unified R script:** `sim_study.R` — adapts to environment via `ON_CLUSTER` and flags.
- **Shell (local or login node):** `run_sim_study.sh` — pull, install, then runs `sim_study.R`; use `--cluster` on the cluster so config matches batch jobs.
- **Batch (NeSI):** `run_PPDisentangle_sim.slurm` — submits a SLURM job that runs the study in cluster mode.

## Why "Mode: LOCAL" on the cluster?

Mode is **LOCAL** when the script is run **without** SLURM (e.g. you ran `run_sim_study.sh` or `Rscript sim_study.R` on the login node). In that case `SLURM_JOB_ID` is unset, so the script uses local defaults (fewer sims/iterations, different save dir logic).

To get **cluster config** when not using sbatch (e.g. interactive or `run_sim_study.sh` on the login node), pass `--cluster`:

```bash
./inst/sim_study/run_sim_study.sh --cluster
# or
Rscript inst/sim_study/sim_study.R --cluster
```

Then you’ll see `Mode: CLUSTER` and the same settings as in the SLURM job.

## Old files in this directory

If you still see **`sim_study_cluster.R`**, **`sim_study_local.R`**, or **`sim_study.slurm`** in the main `inst/sim_study/` folder (not in `deprecated/`), your clone is out of date or the move into `deprecated/` wasn’t applied there. Do **not** run those scripts.

- **Correct layout:** only `sim_study.R`, `run_sim_study.sh`, and `run_PPDisentangle_sim.slurm` (and docs) at top level; legacy scripts live in `deprecated/`.
- **On the cluster:** after `git pull`, if the old files still appear at top level, remove or move them so you only use `sim_study.R` and the run scripts above.

## Why is adaptive SEM slower than EM-style (same iter/n_props)?

With the same `iter=100` and `n_props=10`, **adaptive SEM** does **more work per dataset**:

- **EM-style:** one run of `em_style_labelling(..., iter=100, n_props=10)` per dataset → 100 iterations with 10 proposals.
- **Adaptive SEM:** an **outer** SEM loop (`N_iter` × `N_labellings`); inside it repeatedly calls the **same** `em_style_labelling` (with `adaptive_control$iter=100`, `n_props=10`) when weights are refreshed, plus it generates multiple labellings, computes importance weights, and runs parameter optimisation. So you get at least one full 100-iter labelling run per dataset, plus extra labelling generation, weight calculation, and `optim()` in the SEM loop.

So SEM is expected to take roughly 2–3× longer than EM-style for the same inner iter/n_props; that’s due to the extra SEM structure (extra adaptive steps when weights change, labelling simulation, and M-step optimisation), not a bug.
