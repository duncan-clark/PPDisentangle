# Oklahoma Induced Seismicity Analysis

Causal analysis of the Oklahoma Corporation Commission (OCC) directive
AOI_20150318 on earthquake rates using the PPDisentangle bivariate ETAS
framework. The target estimand is the **all-or-nothing Average Treatment
Effect (ATE)**: the expected per-county difference in earthquake counts
between a world where every county receives the directive and one where
none do.

## Study Design

**Treatment:** On 18 March 2015, the OCC issued directive AOI_20150318
requiring significant reductions in wastewater disposal volumes within
a defined Area of Interest (AOI) in central Oklahoma. This region overlaps
the most seismically active zone linked to injection-induced seismicity.

**Outcome:** Earthquake rate and clustering, measured via a
spatio-temporal Epidemic-Type Aftershock Sequence (ETAS) model.

## Partitioning and Treatment Assignment

**Spatial tessellation:** Oklahoma’s **77 counties** (US Census Bureau
boundaries via `tigris`, year 2022), projected to EPSG:5070 (NAD83
Conus Albers). The analysis builds a `spatstat` tessellation from these
county polygons; events are assigned to counties by spatial location.

**Treatment assignment:** A county is **treated** if its **centroid** lies
inside the OCC AOI polygon; otherwise it is **control**. The AOI geometry
is taken from `occ_aoi_layer_2.geojson` (OCC Area of Interest layer 2).
This rule yields **9 treated counties**: Pawnee, Love, Alfalfa, Oklahoma,
Logan, Payne, Noble, Grant, Lincoln. The remaining 68 counties are
control.

**Rationale:** The centroid rule aligns treatment with the OCC’s
regulatory boundary: counties whose geographic center falls within the
AOI are subject to the directive. Events are then labelled by the
county they occur in (naive labels) or by the SEM-inferred process
(SEM labels).

**Pre-treatment window:** 1 Jan 2014 to 18 Mar 2015  
**Post-treatment window:** 18 Mar 2015 to 24 Jun 2015  
(ends before the next regional directive AOI_20150624 takes effect)

**Magnitude threshold:** $m_0 = 2.5$

## Four-Way Model Comparison

| Label | Model | Labels | Cross-excitation |
|:---:|:---|:---|:---|
| **A** | Naive independent ETAS | Location-based | None |
| **B** | Naive bivariate ETAS | Location-based | Estimated |
| **C** | SEM independent ETAS | SEM-corrected | None |
| **D** | SEM bivariate ETAS | SEM-corrected | Estimated |

Each fit produces an estimate of the all-or-nothing ATE via forward
simulation over a 1-year (365-day) horizon.

## Model

The bivariate ETAS model allows cross-excitation between the treated and
control processes:

```
λ_k(t,x,y) = (μ_k / |S_k|) + Σ_{j: t_j < t} κ_{kl}(m_j) g(t - t_j) f(x-x_j, y-y_j | m_j)
```

where `k, l ∈ {control, treated}` and the 2×2 productivity matrix
`κ_{kl}(m) = A_{kl} exp(α_{kl}(m - m₀))` captures self-excitation
(`A_00`, `A_11`) and spillover (`A_01`, `A_10`).

Structural parameters `(c, p, D, γ, q)` are shared and held fixed during
estimation. The free parameters are:

| Parameter | Description |
|-----------|-------------|
| `mu_0`, `mu_1` | Background rates (control, treated) |
| `A_00`, `alpha_m_00` | Control self-excitation |
| `A_11`, `alpha_m_11` | Treated self-excitation |
| `A_01`, `alpha_m_01` | Treated → control cross-excitation |
| `A_10`, `alpha_m_10` | Control → treated cross-excitation |

## Estimation

1. **Naive fits (A, B):** MLE on location-labeled events. A fits
   independent ETAS to control and treated separately; B fits a joint
   bivariate ETAS with cross-excitation.

2. **SEM fits (C, D):** The Stochastic EM algorithm jointly estimates
   ETAS parameters and latent process labels. The inner adaptive step
   proposes labellings via simulated discrepancy; the outer loop
   re-optimizes the likelihood using importance-weighted proposals.
   Runs 1 outer iteration with 10 labelling proposals (full mode).

## Files

| File | Purpose |
|------|---------|
| `Oklahoma_data_and_viz.R` | Downloads USGS earthquake catalog and OCC AOI geometry, builds regional grid, saves CSVs and GeoJSON |
| `oklahoma_analysis.R` | Main analysis: county tessellation, E/F fits, ATE simulation, plots, results |
| `ate_bivariate.R` | Bivariate vs univariate × all-or-nothing / observed ATE evaluation |
| `oklahoma_report.qmd` | Quarto report source (HTML); defaults to bivariate all-or-nothing `for_paper.rds` |
| `render_oklahoma_report.R` | Render report into `PPDisentangle-output/oklahoma/` |
| `paper/oklahoma_paper_assets.R` | Builds publication figures and LaTeX table fragments from saved analysis results |
| `launch_ate_scenarios_nesi.sh` | Submit four ATE-scenario bootstrap-only jobs on NeSI |
| `launch_oklahoma_cd_nesi.sh` | C/D-primary refresh (default `t_trunc=90`, no bootstrap) |
| `launch_oklahoma_cd_update_nesi.sh` | Re-run C/D primary and archive previous for slim-report diffs |
| `watch_and_publish_cd.sh` | Fetch C/D job; write `oklahoma_cd_current.rds` + slim HTML |
| `oklahoma_report_cd.qmd` / `render_oklahoma_report_cd.R` | Slim HTML: C/D settings, params, ATEs (+ optional prev diff) |
| `render_ate_scenario_compare.R` | Build the four-way scenario-compare HTML/PDF |
| `../sim_study/consistency_study.R` | Simulation study verifying point-process parameter recovery |

## Usage

```bash
# 1. Prepare data (downloads from USGS and OCC, writes to
#    oklahoma_induced_seismicity_data_regional20150318/)
Rscript Oklahoma_data_and_viz.R

# 2. Run analysis (test mode for quick check)
Rscript oklahoma_analysis.R --test

# 3. Run full analysis
Rscript oklahoma_analysis.R

# 4. Render HTML report from saved results (no re-analysis)
#    Prefers PPDisentangle-output/oklahoma/for_paper.rds (bivariate AoN)
Rscript render_oklahoma_report.R
Rscript render_oklahoma_report.R --results ../../../PPDisentangle-output/oklahoma/for_paper.rds
```

### NeSI (from your laptop)

```bash
# Production bivariate all-or-nothing bootstrap (hydrate fits from for_paper.rds):
bash inst/oklahoma/run_nesi.sh --mode bootstrap-only \
  --ate-bivariate true --ate-contrast all_or_nothing \
  --boot-reps 64 --cores 64 --boot-outer-cores 12 \
  --boot-sem-inner 2000 --ate-n-sims 500 \
  --bootstrap-patch-file /path/to/for_paper.rds

# Four-way estimand compare (univ/biv × AoN/observed):
bash inst/oklahoma/launch_ate_scenarios_nesi.sh
bash inst/oklahoma/fetch_and_render_ate_scenarios.sh
```

Bootstrap-only inputs must come from a fit run produced with the current
`beta_gr` stability constraints. The job rejects a source DGP outside
`OK_ETAS_BRANCHING_MAX` (default `0.98`); it does not silently reuse an
explosive legacy fit. Magnitude productivity is also constrained by default:
all `alpha_m` / `alpha_m_*` satisfy `0 < alpha < beta_gr - gap` (larger events
more triggering). Override with `alpha_m_lower_bound = -Inf` only if needed.
Spatial scale is magnitude-independent: `gamma` is fixed at `0` in every fit
(`d(m) = D`). Univariate comparisons include I/J (KDE) and K/L (homogeneous). Bivariate scenarios simulate and refit the full
cross-exciting bivariate law, while univariate scenarios simulate and refit
independent univariate laws. The requested replicate count is the number
attempted. Failed or explosive refits (`eta >= 1` or `rho >= 1`) are excluded,
and the retained count is recorded before bootstrap recentering.

See [`../nesi/README.md`](../nesi/README.md).

## Output

Running `oklahoma_analysis.R` now does the following automatically on success:

1. Saves results and plots to `../PPDisentangle-output/oklahoma/`
2. Renders `oklahoma_report.qmd` (HTML) into the same output folder
3. Writes `last_run_sync_stamp.txt` in canonical output to trigger cloud sync tools

Primary artifacts:

- `PPDisentangle-output/oklahoma/for_paper.rds` — canonical bivariate all-or-nothing results (E/F + bootstrap)
- `PPDisentangle-output/oklahoma/oklahoma_report.html` — interactive review report (Quarto)
- `PPDisentangle-output/oklahoma/oklahoma_ate_scenarios.html` — four-way estimand compare
- `PPDisentangle-output/oklahoma/ate_scenarios/` — per-scenario RDS + summary CSV

Publication-ready paper figures and tables are written under
`PPDisentangle-output/oklahoma/paper/generated/` (outside this repo). Regenerate from
the repository root with:

```bash
Rscript inst/oklahoma/paper/oklahoma_paper_assets.R
```

## Data

The prepared data lives in
`oklahoma_induced_seismicity_data_regional20150318/`. This is a small
versioned input snapshot for reproducibility; rerun `Oklahoma_data_and_viz.R`
to refresh it from USGS/OCC sources.

| File | Contents |
|------|----------|
| `events_all.csv` | All earthquakes (pre + post) with projected coords |
| `events_pre.csv` | Pre-treatment events |
| `events_post.csv` | Post-treatment events |
| `events_raw_usgs.csv` | Raw downloaded USGS events before analysis filtering |
| `cells.csv`, `grid_cells.gpkg` | Regional grid outputs from the data-prep script |
| `analysis_window.gpkg`, `aoi.gpkg` | Analysis and AOI geometries |
| `metadata.json` | Design parameters and counts |
| `occ_aoi_layer_2.geojson` | OCC AOI boundary |

The analysis uses county boundaries from `tigris` (not the grid from
the data prep script); treatment is assigned by county centroid inside
the AOI.

## References

- Ogata, Y. (1988). Statistical models for earthquake occurrences and
  residual analysis for point processes. *JASA*, 83(401), 9–27.
- Zhuang, J., Ogata, Y., & Vere-Jones, D. (2002). Stochastic declustering
  of space-time earthquake occurrences. *JASA*, 97(458), 369–380.
