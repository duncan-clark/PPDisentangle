# Oklahoma Induced Seismicity Analysis

Causal analysis of the Oklahoma Corporation Commission (OCC) directive
AOI_20150318 on earthquake rates using the PPDisentangle bivariate ETAS
framework.

## Study Design

**Treatment:** On 18 March 2015, the OCC issued a directive requiring
significant reductions in wastewater disposal volumes within a defined Area
of Interest (AOI) in central Oklahoma. This region overlaps the most
seismically active zone linked to injection-induced seismicity.

**Outcome:** Earthquake rate and clustering, measured via a spatio-temporal
Epidemic-Type Aftershock Sequence (ETAS) model.

**Control region:** A 300 km buffer around the AOI, tessellated into 20 km
grid cells. Cells whose centroid falls inside the AOI are designated
"treated"; all others are "control."

**Pre-treatment window:** 1 Jan 2014 to 18 Mar 2015  
**Post-treatment window:** 18 Mar 2015 to 24 Jun 2015  
(ends before the next regional directive AOI_20150624 takes effect)

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

1. **Naive fit:** Independent ETAS MLE on location-labeled control and
   treated events (no label correction).

2. **Bivariate SEM:** The Stochastic EM algorithm jointly estimates the
   bivariate ETAS parameters and the latent process labels. The inner
   adaptive step proposes labellings via simulated discrepancy; the outer
   loop re-optimizes the 15-parameter bivariate likelihood using
   importance-weighted proposals. Runs 100 outer iterations.

## Files

| File | Purpose |
|------|---------|
| `Oklahoma_data_and_viz.R` | Downloads USGS earthquake catalog and OCC AOI geometry, builds the grid, assigns treatment, saves CSVs and GeoJSON |
| `oklahoma_analysis.R` | Main analysis: loads data, fits naive + bivariate SEM, produces plots, saves results |
| `etas_consistency_study.R` | Simulation study verifying `fit_etas` recovers known ETAS parameters |

## Usage

```bash
# 1. Prepare data (downloads from USGS and OCC, writes to
#    oklahoma_induced_seismicity_data_regional20150318/)
Rscript Oklahoma_data_and_viz.R

# 2. Run analysis (test mode for quick check)
Rscript oklahoma_analysis.R --test

# 3. Run full analysis (100 SEM iterations)
Rscript oklahoma_analysis.R
```

## Output

Results are saved to `output/`:

- `oklahoma_results.rds` — full results list (naive params, SEM result,
  bivariate params, config)
- `output/plots/partition_map.png` — 2D map of the treatment partition
- `output/plots/pp_pre_treatment.png` — pre-treatment point pattern
- `output/plots/pp_post_treatment.png` — post-treatment (location labels)
- `output/plots/pp_post_sem_labels.png` — post-treatment (SEM-inferred labels)
- `output/plots/sem_flips.png` — SEM convergence (label flips per iteration)

## Data

The prepared data lives in
`oklahoma_induced_seismicity_data_regional20150318/`:

| File | Contents |
|------|----------|
| `events_all.csv` | All earthquakes (pre + post) with projected coords |
| `events_pre.csv` | Pre-treatment events |
| `events_post.csv` | Post-treatment events |
| `cells.csv` | Grid cells with treatment assignment (Z) |
| `metadata.json` | Design parameters and counts |
| `occ_aoi_layer_2.geojson` | OCC AOI boundary |
| `grid_cells.gpkg` | Grid geometry (GeoPackage) |

## References

- Ogata, Y. (1988). Statistical models for earthquake occurrences and
  residual analysis for point processes. *JASA*, 83(401), 9–27.
- Zhuang, J., Ogata, Y., & Vere-Jones, D. (2002). Stochastic declustering
  of space-time earthquake occurrences. *JASA*, 97(458), 369–380.
