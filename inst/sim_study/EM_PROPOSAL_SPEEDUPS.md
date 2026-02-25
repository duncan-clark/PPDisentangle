# EM-style labelling proposal speedups

## Profiling

Use `profile_em_proposals.R` for detailed timing and Rprof:

```bash
Rscript inst/sim_study/profile_em_proposals.R
```

It reports:

1. **Phase timing** for one EM iteration: proposal generation (A), precompute_loglik_args (B), rbind (C), metric/loglik (D).
2. **Per-call breakdown** for a single `simulation_labeling_hawkes_hawkes_fast`: tileindex, two `generate_inhomogeneous_hawkes` calls, full call.
3. **Rprof** on a short `em_style_labelling` run and on 30× `simulation_labeling_hawkes_hawkes_fast` to see where time is spent (e.g. `generate_inhomogeneous_hawkes`, `sim_hawkes_fast`, `tileindex`, `inside.owin`, `loglik_hawk_fast`).
4. **Full EM benchmark** (iter=5, n_props=10, 2 runs) for before/after comparison.

## Implemented speedups

### 1. **Tile index reuse** (`R/labelings.R`, `em_style_labelling` + `simulation_labeling_hawkes_hawkes_fast`)

- **Problem:** Every proposal in an EM iteration uses the same `post`; `simulation_labeling_hawkes_hawkes_fast` was calling `tileindex(dat$x, dat$y, partition)` once per call (n_props times per iteration).
- **Change:** `em_style_labelling` computes `post_inds <- as.numeric(tileindex(post$x, post$y, partition))` once when updating proposals and passes `points_tile_index = post_inds` to `simulation_labeling_hawkes_hawkes_fast`. The fast labeling function accepts optional `points_tile_index` and uses it when length matches `nrow(dat)`, skipping `tileindex()`.
- **Effect:** Saves (n_props − 1) `tileindex` calls per proposal-update step (Rprof had tileindex ~10–22% of labeling time).

### 2. **Filtration split reuse** (`R/hawkes.R` `generate_inhomogeneous_hawkes` + `R/labelings.R` `em_style_labelling`)

- **Problem:** Each `generate_inhomogeneous_hawkes` call (2 per proposal: control and treated) was doing `split(filtration, filtration$location_process)` and, when needed, `tileindex(filtration$x, filtration$y, partition)`. So 2×n_props splits per EM iteration.
- **Change:** `generate_inhomogeneous_hawkes` accepts optional `filt_by_proc` in `...`. When provided, it skips building `filt_by_proc` from `filtration`. `em_style_labelling` computes `filt_by_proc <- split(pre, pre$location_process)` once per proposal-update and passes `filt_by_proc = filt_by_proc` into `simulation_labeling_hawkes_hawkes_fast`, which forwards it via `...` to `generate_inhomogeneous_hawkes`.
- **Effect:** Removes 2×n_props split (and any filtration tileindex) per iteration.

### 3. **Statespace as owin once** (`R/labelings.R` `em_style_labelling`)

- **Problem:** `statespace` (e.g. numeric `OMEGA`) was passed through to `precompute_loglik_args` and `generate_inhomogeneous_hawkes`, which call `as.owin(statespace)` repeatedly.
- **Change:** At the start of `em_style_labelling`, `if (!inherits(statespace, "owin")) statespace <- as.owin(statespace)`. All later uses (precomp, labeling calls) then receive an owin.
- **Effect:** Fewer redundant `as.owin`/owin-related work in the hot path.

## Benchmark (typical)

On the script’s setup (591 post points, 80 pre, 5×5 grid, iter=5, n_props=10):

- **Before:** Full `em_style_labelling` mean ~2.3 s (2 runs).
- **After:** Full `em_style_labelling` mean ~1.4 s (2 runs).

So about **1.5–1.7× faster** for the EM proposal path, with no change to results (all tests pass).

## Optional further work

- **adaptive_SEM:** The same `points_tile_index` and `filt_by_proc` pattern could be used when generating the N_labellings in the adaptive step (one tile index for the current labelling, one split of pre/filtration per step).
- **Rprof** still shows notable time in `inside.owin`, `sim_hawkes_fast`, and `loglik_hawk_fast`; any further gains would likely need changes there (e.g. fewer point-in-polygon checks or faster likelihood code).
