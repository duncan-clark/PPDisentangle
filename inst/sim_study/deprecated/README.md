# Deprecated: old sim study scripts and profiling tools

These scripts are superseded by the unified `../sim_study.R` and `../run_sim_study.sh`.

**Use instead:**
- `../sim_study.R` – main simulation study (supports `--small`, `--cluster`, `--sims N`)
- `../run_sim_study.sh` – SLURM wrapper (submits via sbatch, runs sim_study.R)
- `../consistency_study.R` – sanity check for Hawkes fit

**Original local/cluster split:**
- `sim_study_local.R` – old local run
- `sim_study_cluster.R` – old cluster run

**Profiling (no longer needed):**
- `sim_study_profile.R` – profiled version of sim study
- `profile_sim_hawkes.R`, `profile_em_proposals.R`, `tmp_profiling.R` – profiling helpers
- `sim_study.slurm`, `run_PPDisentangle_sim.slurm` – standalone SLURM scripts
- `tmp_sim_study_speedups.md`, `EM_PROPOSAL_SPEEDUPS.md`, `SEM_TIMING.md` – profiling notes

Kept here for reference only.
