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
| `oklahoma_analysis.R` | Main analysis: county tessellation, four-way fits, ATE simulation, plots, results |
| `oklahoma_report.qmd` | Quarto report (HTML + PDF) |
| `etas_consistency_study.R` | Simulation study verifying `fit_etas` recovers known ETAS parameters |

## Usage

```bash
# 1. Prepare data (downloads from USGS and OCC, writes to
#    oklahoma_induced_seismicity_data_regional20150318/)
Rscript Oklahoma_data_and_viz.R

# 2. Run analysis (test mode for quick check)
Rscript oklahoma_analysis.R --test

# 3. Run full analysis
Rscript oklahoma_analysis.R

# 4. Optional manual render (analysis auto-renders by default)
quarto render oklahoma_report.qmd
```

## Output

Running `oklahoma_analysis.R` now does the following automatically on success:

1. Saves results and plots to `output/oklahoma/`
2. Mirrors core artifacts to legacy paths (`output/` and `inst/oklahoma/output/`) for compatibility
3. Renders `oklahoma_report.qmd` (HTML + PDF + TeX)
4. Writes `last_run_sync_stamp.txt` to output folders to trigger cloud sync tools

Primary artifacts:

- `output/oklahoma/oklahoma_results.rds` — full results (fits A–D, ATE, config, counties)
- `output/oklahoma/plots/partition_map.png` — county partition (treated vs control)
- `output/oklahoma/plots/pp_pre_treatment.png` — pre-treatment point pattern
- `output/oklahoma/plots/pp_post_treatment.png` — post-treatment (location labels)
- `output/oklahoma/plots/pp_post_sem_indep.png` — post-treatment (SEM independent labels)
- `output/oklahoma/plots/pp_post_sem_biv.png` — post-treatment (SEM bivariate labels)
- `output/oklahoma/plots/sem_flips_indep.png` — SEM convergence (independent)
- `output/oklahoma/plots/sem_flips_biv.png` — SEM convergence (bivariate)

## Data

The prepared data lives in
`oklahoma_induced_seismicity_data_regional20150318/`:

| File | Contents |
|------|----------|
| `events_all.csv` | All earthquakes (pre + post) with projected coords |
| `events_pre.csv` | Pre-treatment events |
| `events_post.csv` | Post-treatment events |
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
