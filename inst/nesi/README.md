# NeSI workflow (laptop → cluster → laptop)

Run heavy jobs on [NeSI Mahuika](https://docs.nesi.org.nz) from your laptop via
SSH. Results are written to project-space on the cluster, then **rsync'd** into
your local `PPDisentangle-output/` for review.

## One-time setup

### 1. Configure SSH

Ensure you can log in:

```bash
ssh your_username@mahuika.nesi.org.nz
```

Optional: add a `Host` alias in `~/.ssh/config` and set `PP_NESI_SSH` in your
env file.

### 2. Create local config

```bash
cp inst/nesi/config.example.env inst/nesi/nesi.env
# edit inst/nesi/nesi.env — set PP_NESI_USER, paths, and any local secrets
```

**Credentials:** all `*.env` files under `inst/nesi/` are gitignored (only
`config.example.env` is tracked). Prefer SSH keys; never put secrets in the
example template or commit `nesi.env`.

Config is loaded from (first match):

1. `$PPDISENTANGLE_NESI_CONFIG`
2. `~/.config/ppdisentangle/nesi.env`
3. `inst/nesi/nesi.env` (gitignored)

Print resolved settings:

```bash
bash inst/nesi/submit.sh --show-config
```

### 3. Clone on NeSI and install dependencies

On NeSI (once):

```bash
git clone git@github.com:duncan-clark/PPDisentangle.git ~/PPDisentangle
cd ~/PPDisentangle
git checkout major-revisions-cleanup   # or your branch
```

From your laptop:

```bash
bash inst/nesi/setup_remote.sh
```

This runs `inst/sim_study/setup_nesi.sh` on the login node (R-Geo, deps,
package install).

## Infrastructure smoke test (no SEM / no sim study)

Verify SSH, SLURM submit, and rsync **without** running analysis code:

```bash
bash inst/nesi/smoke_test.sh
```

This submits a ~seconds-long batch job that writes a marker file under
`PPDisentangle-output/infra_smoke/` on NeSI, waits for completion, and rsyncs
it to your local `PPDisentangle-output/infra_smoke/`.

Use this before any real `sim_study` or Oklahoma job.

## Typical workflow

```text
Edit locally → (optional push) → submit via SSH → SLURM runs → fetch rsync → review locally
```

### Submit a job

```bash
# Quick smoke test
bash inst/nesi/submit.sh sim_study --test --sims 2

# Full sim study
bash inst/nesi/submit.sh sim_study --sims 100

# Oklahoma application
bash inst/nesi/submit.sh oklahoma --cores 32 --mode long
```

The remote launcher:

1. `git pull` your branch (see `PP_GIT_BRANCH`, default: local current branch)
2. Sets `PPDISENTANGLE_OUTPUT_ROOT` to `PP_NESI_REMOTE_OUTPUT` on the cluster
3. Runs `inst/*/run_nesi.sh`, which `sbatch`s the compute job

Push local commits before submit:

```bash
PP_NESI_PUSH=1 bash inst/nesi/submit.sh sim_study --sims 50
```

### Monitor

```bash
bash inst/nesi/status.sh
```

Or on NeSI: `squeue -u $USER`

### Fetch results

After the job finishes:

```bash
# Entire sim_study output tree
bash inst/nesi/fetch.sh sim_study

# One job id (RDS, logs, slurm files, sweep artifacts)
bash inst/nesi/fetch.sh sim_study 1234567

# Oklahoma outputs
bash inst/nesi/fetch.sh oklahoma 1234567

# Both applications
bash inst/nesi/fetch.sh all
```

Results land in your **local** `../PPDisentangle-output/` (same layout as on
the cluster).

### Wait and fetch automatically

```bash
bash inst/nesi/wait_and_fetch.sh sim_study 1234567
```

Polls `squeue` every 60s (`PP_NESI_POLL_SECS`) then runs `fetch.sh`.

## Output locations

| Location | Role |
|----------|------|
| `PP_NESI_REMOTE_OUTPUT` on NeSI | Canonical cluster outputs (default: `/nesi/project/uoo04008/PPDisentangle-output`) |
| Local `PPDisentangle-output/` | Your review copy after `fetch.sh` |

Override locally or on cluster with `PPDISENTANGLE_OUTPUT_ROOT`.

## Scripts reference

| Script | Runs where | Purpose |
|--------|------------|---------|
| `inst/nesi/submit.sh` | Laptop | SSH → remote `run_nesi.sh` → sbatch |
| `inst/nesi/fetch.sh` | Laptop | rsync cluster → local output |
| `inst/nesi/wait_and_fetch.sh` | Laptop | Poll queue, then fetch |
| `inst/nesi/status.sh` | Laptop | Remote `squeue` |
| `inst/nesi/setup_remote.sh` | Laptop | One-time remote setup |
| `inst/sim_study/run_nesi.sh` | NeSI | Submit + worker for sim study |
| `inst/oklahoma/run_nesi.sh` | NeSI | Submit + worker for Oklahoma |

You can still SSH manually and run `run_nesi.sh` on the login node; the
`inst/nesi/*` scripts are convenience wrappers for the laptop → cluster →
laptop loop.

## Git branch sync

Remote submit no longer hard-codes `git pull origin main`. It uses
`inst/include/git_sync.sh`:

- `PP_GIT_BRANCH` — explicit branch to pull on NeSI
- default — same branch name as your **local** checkout when you run `submit.sh`
- `PP_GIT_SYNC=0` — skip pull entirely
