# Why adaptive SEM is slow (and how to speed it up)

## What the 3320 s is

The "Running adaptive SEM ..." block runs `adaptive_SEM()` **once per simulated dataset** (e.g. 100 times on cluster, or `SIM_SIZE` times locally). The total time you see (e.g. 3320.7 s) is the wall-clock time for all those runs (possibly in parallel). So:

- **Per-dataset cost** ≈ total time ÷ number of datasets (or ÷ effective parallel cores).
- Example: 3320 s for 100 datasets on 16 cores → ~33 s per run if parallelized, or 3320 s for 8 datasets on 1 core → ~415 s per run.

Either way, each single `adaptive_SEM()` call is doing a lot of work.

---

## Why each `adaptive_SEM()` call is expensive

`adaptive_SEM()` has an **outer loop** (`N_iter = 10`) and, when it decides to refresh, an **adaptive step** plus **labelling generation**.

### 1. Adaptive step = full EM with 1000 iterations

When `check_weights(weights) == TRUE` or on the first iteration, the code runs **one full** `em_style_labelling()` with:

- **iter** = `SEM_EM_ADAPTIVE_ITER` = **1000**
- **n_props** = **10**

So each adaptive step does:

- 1000 EM iterations × 10 proposals = **10,000** calls to `simulation_labeling_hawkes_hawkes_fast` and many `loglik_hawk_fast` evaluations.

That is comparable to the whole of **Section 6 (EM-style labelling)** for one dataset, which runs a single `em_style_labelling(..., iter = EM_ITER = 1000, n_props = 10)`. So **one adaptive step ≈ one full Section-6 run per dataset**.

### 2. Generating N_labellings after each adaptive step

Right after each adaptive step, the code generates **N_labellings** new labellings:

- `N_labellings = max(10, N_PROPOSALS / 10)` → with `N_PROPOSALS = 1000` this is **100** labellings.
- Each is one call to `simulation_labeling_hawkes_hawkes_fast()`.

So per adaptive step: **100** expensive labelling simulations.

### 3. How often does the adaptive step run?

- **First** outer iteration: always runs (adaptive_counter == 0).
- **Later** outer iterations: runs again only when `check_weights(weights) == TRUE` (top 5% of importance weights sum to > 0.95).

So you typically get **1–2 adaptive steps** per `adaptive_SEM()` run, each step being (1000-iter EM + 100 labellings).

### 4. Every outer iteration: weights + M-step optim

On **every** of the 10 outer iterations the code:

1. **calculate_weights(labellings)**  
   For each of the ~100 labellings: 2 × `loglik_hawk_fast` (control + treated) → **~200 loglik** calls per outer iteration.

2. **optim(...)** for the M-step  
   - Objective: weighted sum of loglik over all “keeper” labellings (~100).
   - Nelder-Mead with **maxit = 1000**.
   - Each evaluation of the objective = **~100** `loglik_hawk_fast` calls (one per keeper).
   - So in the worst case: 1000 × 100 = **100,000** loglik evaluations in a single M-step.

So the **dominant costs** per `adaptive_SEM()` run are:

- **Adaptive steps**: 1–2 × (em_style_labelling with 1000 iter × 10 props + 100 labellings).
- **M-step**: 10 outer iterations × optim with many evaluations × ~100 loglik per evaluation.

---

## Are the sim_study settings “bad”?

They are **very heavy** by design:

- **SEM_EM_ADAPTIVE_ITER = 1000** and **EM_ITER = 1000**: full-length EM inside SEM.
- **N_labellings = 100**: many labellings per adaptive step.
- **N_iter = 10**: 10 outer SEM iterations, each with weight computation and a long optim.
- **optim(..., maxit = 1000)**: up to 1000 M-step evaluations per outer iteration.

So the 3320 s is **expected** for the current settings: they are chosen for quality/convergence, not speed.

---

## How to make adaptive SEM faster

### 1. Reduce SEM-specific settings (recommended for testing / development)

In `sim_study.R`, you can lower the SEM workload without changing the rest of the study:

- **SEM_EM_ADAPTIVE_ITER**: e.g. **100** or **200** (fewer EM iterations inside each adaptive step).
- **N_labellings**: e.g. **20** (e.g. set `N_PROPOSALS` lower for SEM only, or add a separate `SEM_N_LABELLINGS` and use `N_labellings = SEM_N_LABELLINGS` in `run_sem`).
- **N_iter** (outer SEM iterations): e.g. **5** instead of 10.
- In `adaptive_SEM`, **optim maxit**: reduce from 1000 to e.g. **200** (in `em_algorithm.R`).

Using **--small** already sets `SEM_EM_ADAPTIVE_ITER <- 10` and `N_PROPOSALS <- 5` (so `N_labellings = 10`), which makes Section 7 much faster for quick checks.

### 2. Keep full settings for production

For final simulation results, keeping **SEM_EM_ADAPTIVE_ITER = 1000** and the current N_labellings/N_iter is reasonable; then long runtimes are expected and parallelization (e.g. cluster with many cores and SIM_SIZE = 100) is the way to get results in reasonable wall-clock time.

### 3. Optional: reduce verbose output

In `R/em_algorithm.R`, inside `adaptive_step`, `em_style_labelling` is called with `verbose = TRUE`. You can set that to `adaptive_control$verbose` (or `FALSE`) to reduce I/O and clutter when running many SEM runs.

---

## Summary

- **Why it’s slow**: Each `adaptive_SEM()` does 1–2 full EM runs (1000 iter × 10 props), generates 100 labellings, then 10 outer iterations each with ~200 loglik calls for weights and an M-step optim that can do up to ~100k loglik evaluations.
- **Settings**: They are not “wrong”; they are high-quality and therefore expensive.
- **To speed up**: Lower `SEM_EM_ADAPTIVE_ITER`, `N_labellings`, `N_iter`, and/or M-step `maxit` for development; use **--small** for quick runs; keep current settings for production and rely on parallel runs on a cluster.
