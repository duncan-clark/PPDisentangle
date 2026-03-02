# Gambia Analysis — Deprecated Scripts

This folder contains the iterative exploratory scripts that led to the final
radii study design in `inst/gambia/gambia_radii_study.R`.  They are preserved
here as a record of the analytical reasoning process.

---

## The Gambia Vaccine Study

The dataset comes from a pneumococcal conjugate vaccine (PCV) trial in The
Gambia.  Health centres administered vaccines starting around May 2011.
We model IPD (invasive pneumococcal disease) cases as a spatio-temporal Hawkes
process and use the PPDisentangle SEM framework to disentangle background
(endemic) and triggered (clustering) cases, then estimate vaccine savings by
comparing treated (near health centres) and control (far from health centres)
regions.

---

## Evolution of the Analysis

### Phase 1: First exploratory fit (`gambia_fit.R`)

The initial script used **all cases** (IPD and non-IPD together) with
a population raster (`GMB_ppp_v2b_2020.tif`) for the inhomogeneous background
rate.  Key choices:

- **Treatment radius**: 500 m (very small — only the immediate vicinity of
  each health centre).
- **Beta**: free (estimated from data).
- **Background**: population density raster, normalized per-region.
- **Spatial jitter**: sd = 1 m added to break ties.
- **Fits**: pre/post vanilla Hawkes on each partition, then adaptive SEM with
  simulation-based proposals.

**Problems identified**:
- The population raster is an external dependency and a proxy, not directly
  related to disease incidence patterns.
- The 500 m treatment radius is too narrow — PCV rollout covered much larger
  catchment areas.
- Mixing IPD and non-IPD cases conflates the disease process we want to model
  with general healthcare utilization.

### Phase 2: Inhomogeneous fit with non-IPD background (`gambia_fit_inhom.R`)

Addressed the background rate problem by using **non-IPD cases as a spatial
proxy for healthcare access** — the idea being that non-IPD presentations
reflect the spatial distribution of the at-risk population without being
contaminated by the vaccine effect we want to measure.

Key changes:
- **Data**: IPD cases only for the Hawkes model.
- **Background**: KDE from sampled non-IPD cases (`N_BG_SAMPLE = 1000`) using
  `bw.diggle`, regionally normalized so the mark W integrates to the
  partition area.
- **Treatment radius**: 10 km (much more realistic catchment).
- **Beta**: fixed at 0.05 (~20 day mean trigger time) as a profile likelihood
  — avoids identifiability issues between mu and beta.
- **Fits**: vanilla Hawkes, SEM with `single_flip` (greedy), and SEM with
  `simulation` (EM-style) proposals — compared both to understand sensitivity.
- Also ran a **non-IPD fit** (reversing the roles: IPD background for non-IPD
  estimation) as a falsification check.

This was the first clean, well-structured analysis. But a single 10 km radius
is arbitrary — the question remained: *how sensitive are results to the
definition of "treated"?*

### Phase 3: Sensitivity analysis (`gambia_sensitivity.R`)

Before varying radii, we needed to understand how sensitive K estimates are
to the **background specification itself**.  This script ran vanilla Hawkes
fits (no SEM) across five background settings:

1. Homogeneous (no spatial mark)
2. IPD pre-treatment KDE, `bw.diggle` bandwidth
3. IPD pre-treatment KDE, `bw * 20` (over-smoothed)
4. Non-IPD sample KDE, `bw.diggle`
5. Non-IPD sample KDE, `bw * 20`

**Key finding** (see `gambia_sensitivity_output.txt`): K estimates are fairly
robust across specifications (range ~0.25–0.32), with the non-IPD `bw.diggle`
giving the most stable results across pre/post and treated/control.  The IPD
`bw.diggle` showed some instability (very narrow bandwidth → noisy KDE).
This confirmed non-IPD Diggle as the right background choice.

### Phase 4: Tuning the SEM inner loop (`gambia_pilot_tuning.R`, `three_min_pilot.R`)

The SEM inner loop has two key tuning parameters: `change_factor` (proposal
perturbation size) and `n_props` (number of proposals per iteration).

`gambia_pilot_tuning.R` ran a grid search over:
- `change_factor` ∈ {0.01, 0.05, 0.1, 0.2}
- `n_props` ∈ {10, 50, 100}

measuring average flips (label changes accepted) and convergence speed.  The
sweet spot was `change_factor = 0.05`, `n_props = 20`, which balanced
exploration with computation time.

`three_min_pilot.R` was a quick sanity-check pilot — a stripped-down version
of the full fit to verify the pipeline runs end-to-end in ~3 minutes.

### Phase 5: Radii study — the final design

With the background specification and SEM tuning settled, the remaining
question was the treatment radius.  Two radii study scripts were developed:

**`gambia_radii_study.R`** (original):
- IPD cases, non-IPD Diggle background.
- Beta **fixed** at 0.05.
- Radii: 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20 km.
- Vanilla + SEM at each radius.

**`gambia_radii_study_trunc.R`** (truncated variant):
- Same as above but with **temporal kernel truncation at 30 days** and
  **beta free** (estimated from data).
- Rationale: pneumococcal transmission is unlikely beyond ~30 days, so
  truncation provides a biologically motivated regularization that replaces
  the need to fix beta.

These two scripts were nearly identical (same data, background, radii grid,
SEM settings) — differing only in the temporal kernel specification.

---

## Final Consolidated Script

The two radii studies have been merged into a single
`inst/gambia/gambia_radii_study.R` with a `TEMPORAL_MODE` switch:

- `"fixed_beta"` — beta = 0.05, no truncation (original study)
- `"truncated"` — beta free, kernel truncated at 30 days

This is the only active analysis script.  Change `TEMPORAL_MODE` at the top
to run either variant.

---

## File Inventory

| File | Purpose | Key Decisions |
|------|---------|---------------|
| `gambia_fit.R` | First exploratory fit | All cases, population raster background, 500m radius |
| `gambia_fit_inhom.R` | Clean inhomogeneous fit | IPD only, non-IPD background, 10km, fixed beta, single_flip + simulation SEM |
| `gambia_sensitivity.R` | Background specification sensitivity | 5 background variants, vanilla only, confirmed non-IPD Diggle |
| `gambia_sensitivity_output.txt` | Output from sensitivity analysis | K estimates across all specifications |
| `gambia_pilot_tuning.R` | SEM inner loop tuning | Grid over change_factor × n_props |
| `three_min_pilot.R` | Quick end-to-end sanity check | Minimal pilot run |
| `gambia_radii_study_trunc.R` | Truncated radii study (superseded) | Now merged into `../gambia_radii_study.R` with `TEMPORAL_MODE = "truncated"` |
