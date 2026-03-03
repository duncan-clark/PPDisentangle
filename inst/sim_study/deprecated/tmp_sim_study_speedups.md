# Simulation study: deep profiling and 10 speedups (by potential impact)

**Temporary analysis file.** Do not change `sim_study.R`; use this only to guide future optimisations.

---

## 1. How to run the profiler

From the **package root** (e.g. `PPDisentangle/`):

```bash
Rscript inst/sim_study/tmp_profiling.R
```

Or from `inst/sim_study/`:

```bash
Rscript tmp_profiling.R
```

This runs `sim_study.R` with `--small` under `Rprof()` and writes:

- `tmp_prof.out` – raw Rprof log  
- `tmp_prof_summary.rds` – `summaryRprof()` result  
- `tmp_prof_by_self.csv` – time by function (self)  
- `tmp_prof_by_total.csv` – time by function (total)  
- `tmp_profiling_report.txt` – short human-readable summary  

Interpret “self” as time spent inside that function; “total” as self + callees. Focus on “by.self” for bottlenecks.

---

## 2. Where time goes (high level)

From the current design of `sim_study.R`, time is spent roughly as follows.

| Phase | What it does | Likely share of total |
|-------|----------------|-------------------------|
| **2. True tau_i** | `parSapply` over cells: each calls `tau_i()` → `generate_inhomogeneous_hawkes` × 2 × N_TAU_I_TRUE | Medium–high |
| **3. Obs data** | `lapply`: for each sim, `sim_hawkes` (pre) + `generate_inhomogeneous_hawkes` (post) | High |
| **5. Labelling proposals** | For each sim, N_PROPOSALS × `simulation_labeling_hawkes_hawkes` (each does `generate_inhomogeneous_hawkes` + relabel) | **Very high** |
| **5b. Best proposal** | Per sim, N_PROPOSALS × `caret::confusionMatrix` | Medium (if N_PROPOSALS large) |
| **6. EM-style** | `parLapply(run_em)`: each `em_style_labelling(..., n_props=10)` – inner loop is `loglik_hawk_fast` + `simulation_labeling_hawkes_hawkes_fast` | **Very high** |
| **7. Adaptive SEM** | `parLapply(run_sem)`: `adaptive_SEM` with inner `em_style_labelling` + weighted M-step | **Very high** |
| **10. ATE estimation** | `parLapply(task_function)`: each `ATE_estim_hawkes` (optional `fit_hawkes` + `tau_i` × n_tau_i + all-or-nothing sims × n_sims) | High |

So the main cost drivers are:

1. **Hawkes simulation** – `generate_inhomogeneous_hawkes` / `sim_hawkes_fast` (used in data gen, proposals, tau_i, ATE sims).  
2. **Likelihood** – `loglik_hawk_fast` (and C++ `hawkes_loglik_inhom_cpp`) in EM/SEM.  
3. **Proposal generation** – `simulation_labeling_hawkes_hawkes` (or `_fast`) and how many proposals per sim / per EM iter.  
4. **ATE_estim_hawkes** – `fit_hawkes`, repeated `tau_i`, and all-or-nothing Monte Carlo.

---

## 3. Already done (checked against combined sim study + package)

These were implemented earlier; the **combined** sim study and package already use them where applicable.

| What | Where it's done |
|------|------------------|
| **Precomp for loglik in EM** | `em_style_labelling` uses `precompute_loglik_args(ref_post, ...)` and passes `precomp` into every `loglik_hawk_fast` call over proposals (`R/labelings.R` ~683–706). So Sections 6 and 7 (run_em / run_sem) already get the fast likelihood path. |
| **Accuracy without caret inside EM** | Inside `em_style_labelling`, accuracy is `mean(y_post$inferred_process == y_post$process)` (labelings.R ~724, 842, 850). No `confusionMatrix` in the inner EM loop. |
| **Proposals _fast inside EM** | `em_style_labelling` calls `simulation_labeling_hawkes_hawkes_fast` for its own proposals (labelings.R ~657). So Section 6 (EM-style) uses the fast proposal path. |
| **normalize_weights log-sum-exp** | `normalize_weights()` in `R/utils.R` uses log-sum-exp; used by `adaptive_SEM` for importance weights. |
| **C++ loglik + binary search** | `hawkes_loglik_inhom_cpp` uses temporal cutoff and binary search; used by `loglik_hawk_fast`. |
| **generate_inhomogeneous_hawkes** | Already refactored to reduce `as.data.frame`/list overhead in hot paths (conversation history). |

**Still not done:**

- **sim_study.R Section 5** still uses `simulation_labeling_hawkes_hawkes` (slow) and `caret::confusionMatrix` in two places (best proposal + Section 9 class_metrics).
- **adaptive_SEM** (`R/em_algorithm.R` ~175) still uses `simulation_labeling_hawkes_hawkes` (slow) when generating labellings after the adaptive step—only `em_style_labelling` uses `_fast` for its inner proposals.

---

## 4. Updated list: remaining speedups (by potential impact)

These are the ones still open after the “already done” check.

### 1. **[Sim study]** Use `simulation_labeling_hawkes_hawkes_fast` in Section 5

**Where:** `sim_study.R` Section 5: `tmp <- simulation_labeling_hawkes_hawkes(...)`.  
**Change:** Call `simulation_labeling_hawkes_hawkes_fast` with the same arguments (arg name is `partition_process` in both).  
**Why big:** Section 5 does SIM_SIZE × N_PROPOSALS proposal steps; each call runs a full Hawkes sim + relabel.  
**Risk:** Low – API is compatible; `_fast` is already used inside `em_style_labelling`.

---

### 2. **[Sim study]** Replace `caret::confusionMatrix` in Section 5 (best proposal) and Section 9 (class_metrics)

**Where:** (a) `pp_labeled_best_proposal`: `caret::confusionMatrix(...)$overall[["Accuracy"]]`; (b) `class_metrics`: same pattern.  
**Change:** Use `mean(y$inferred_process[keep] == y$process[keep])` in both places.  
**Why big:** N_PROPOSALS × SIM_SIZE accuracy evals in (a), and 5 methods × SIM_SIZE in (b).  
**Risk:** None – same quantity, much cheaper.

---

### 3. **[Package]** Use `simulation_labeling_hawkes_hawkes_fast` in `adaptive_SEM`

**Where:** `R/em_algorithm.R` ~175: `labellings <- lapply(1:N_labellings, function(i) { simulation_labeling_hawkes_hawkes(...) })`.  
**Change:** Call `simulation_labeling_hawkes_hawkes_fast` instead (same args; note `partition_process`).  
**Why big:** Section 7 runs adaptive_SEM; each outer iteration generates N_labellings proposals with the slow function.  
**Risk:** Low – same API; _fast is the intended path.

---

### 4. Reduce Hawkes simulation cost in `generate_inhomogeneous_hawkes` / `sim_hawkes_fast`

**Where:** `R/hawkes.R`: `generate_inhomogeneous_hawkes`, `sim_hawkes_fast`; C++ offspring if used.  
**Changes:** (a) More C++ for offspring generation; (b) further reduce data.frame/list churn; (c) batch `inside.owin` where possible.  
**Why big:** Simulation is used in data gen, Section 5, tau_i, and ATE; 20–30% gain here helps everywhere.  
**Risk:** Medium – preserve behaviour and boundaries.

**Done (partial):** Background phase in `generate_inhomogeneous_hawkes` now builds a single list of vectors (no big `rbind` data.frame) and passes list subsets to `sim_hawkes_fast` as `background_realization` instead of data.frame subsets, reducing allocations. Offspring generation remains in C++; `inside.owin` is already batched per call.

---

### 5. Fewer N_PROPOSALS in Section 5 or parallelise the proposal loop

**Where:** `sim_study.R` Section 5: N_PROPOSALS (e.g. 1000) per sim, single `lapply`.  
**Changes:** (a) Reduce N_PROPOSALS when acceptable (e.g. 100–200); (b) or parallelise over proposals (e.g. parLapply over 1:N_PROPOSALS with chunking).  
**Why big:** Section 5 is a major fraction of total time.  
**Risk:** (a) Slight change in “best” proposal; (b) more complex, more memory on cluster.

---

### 6. Reduce work in `ATE_estim_hawkes` per task (Section 10)

**Where:** `task_function` calls `ATE_estim_hawkes` with `n_sims`, `n_tau_sims`, `n_tau_i`; some tasks refit via `fit_hawkes`.  
**Changes:** (a) Lower `n_tau_i` / `n_tau_sims` for exploratory runs; (b) skip or shorten refit when `hawkes_params` provided (EM_full_params); (c) fewer `optim` iterations for fit.  
**Why big:** Many tasks; each does fits and many sims.  
**Risk:** Low–medium; keep a “full” config for final runs.

---

### 7. True one-flip ATE (Section 2): fewer cells or sims

**Where:** Section 2: `parSapply(..., tau_i, ..., n_sim = N_TAU_I_TRUE)` over all cells.  
**Change:** Reduce `N_TAU_I_TRUE` or sample a subset of cells for quick runs.  
**Why big:** 100 cells × 100 sims × 2 = 20k sims.  
**Risk:** Low if used as tuning/exploratory only.

---

### 8. Cluster export and eval overhead

**Where:** `make_cluster`, `export_globals`, and parLapply/parSapply over large objects.  
**Changes:** Export once; chunk tasks so each worker does fewer, larger jobs; avoid re-exporting unchanged data.  
**Why big:** Can dominate if tasks are short or data is large.  
**Risk:** Low; refactor only.

---

### 9. Data generation (Section 3): lighter or shared pre

**Where:** Section 3: per sim, `sim_hawkes` (pre) + `generate_inhomogeneous_hawkes` (post).  
**Change:** (a) For experiments only: reuse one pre-treatment, vary post; (b) or reduce SIM_SIZE for profiling.  
**Why big:** Section 3 is heavy.  
**Risk:** (a) Design change; (b) profiling-only.

---

### 10. Plotting and results assembly (Sections 11–13)

**Where:** Building `results_df`, `sim_study_plots`, ggplot/grid.  
**Change:** Skip plot building in batch/cluster or when a flag is set (e.g. `if (!SKIP_PLOTS) { ... }`).  
**Why big:** Smaller than 5–7 and 10 but non-negligible; keeps timing runs clean.  
**Risk:** None if gated by a flag.

---

## 5. Suggested order of implementation

1. **Quick wins in sim study:** #1 (Section 5 → _fast), #2 (confusionMatrix → mean), then #3 (adaptive_SEM → _fast).  
2. **Package/core:** #4 (simulation speed).  
3. **Tuning/parallel:** #5 (N_PROPOSALS or parallel), #6, #7, #9.  
4. **Infrastructure:** #8 (cluster), #10 (optional plotting).

Run `tmp_profiling.R` before and after each change to confirm improvements in `by.self` / `by.total` and total elapsed time.
