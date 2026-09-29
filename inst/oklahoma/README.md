# Oklahoma induced-seismicity application

The paper uses two bivariate ETAS fits with a background KDE: a naive fit
with location labels and a SEM fit with inferred component labels. The
canonical saved results are `PPDisentangle-output/oklahoma/for_paper.rds`
(job 8804859 plus the ATE backfill), distributed in the companion Zenodo
results archive. Use software v0.2.0 with that archive.

## Paper estimand and reproduction

The paper reports the **observed AOI assignment versus no treatment
anywhere**, over 100 days, conditional on the retained pretreatment history.
The sign is `saved = N(no treatment) - N(observed assignment)`.

| Method | No treatment | Observed assignment | Saved |
|---|---:|---:|---:|
| Naive | 961.576 | 1013.084 | -51.508 |
| SEM | 1105.422 | 930.222 | 175.200 |

From the package root, after unpacking the current results archive:

```bash
export PPDISENTANGLE_OUTPUT_ROOT=/path/to/PPDisentangle-output
Rscript inst/oklahoma/paper/oklahoma_paper_assets.R --contrast observed
```

This generates figures, LaTeX tables, full-precision `expected_counts.csv`
(including all partition sensitivities), bootstrap summaries, and a build
manifest under `oklahoma/paper/generated/`. It does not refit the models.
The script stops if the requested contrast is unavailable.

The saved object retains legacy all-or-nothing top-level fields, including
`config$ATE_CONTRAST`. The paper renderer explicitly selects
`ate$by_contrast$observed` and `ate_total_mean_observed`. The all-or-nothing
results are an additional estimand; pass `--contrast all_or_nothing` to
render them separately. Legacy HTML reports and scenario launchers can
use that additional estimand and are not the paper reproduction entrypoint.

## Archived design and fitting settings

- Intervention date: 18 March 2015; post-treatment data end on 24 June 2015,
  before the next regional directive. The regulatory action concerned well
  depth and conditional disposal-volume reductions; it was not a uniform
  immediate 50% cut across the AOI.
- Magnitude threshold 2.5; 2472 pretreatment events, split into 1236 KDE
  training events and 1236 events retained for estimation/history; 818
  post-treatment events over 98 days. Counterfactual evaluation uses 100 days.
- Partition: 77 counties; 9 treated by centroid inside the OCC AOI. This is
  an analysis assignment rule, not a statement that every well in a county
  faced the same regulatory requirement.
- KDE: Scott's isotropic bandwidth, trained on the first half of pretreatment
  events. Each component's shape has spatial mean one on its observed support.
- Bivariate C/D-only fits (also exposed through E/F aliases): naive and SEM,
  shared kernel shape parameters estimated, spatial magnitude exponent
  `gamma=0`, Gutenberg-Richter `beta_gr=2.92`, temporal truncation 250 days.
  Finite truncation permits temporal `p<1`; spatial `q>1` is enforced.
- SEM: one outer round, 2000 inner iterations, 20 proposals per iteration,
  zero retained labellings, sequential Bernoulli initialization, start from
  the corresponding naive C fit, and a monotone complete-data likelihood
  gate. Sequential MAP proposals are off.
- Counterfactuals: 500 simulations per contrast, retaining the pretreatment
  history. Fitted background density is kept on an observed-support
  reference; immigrant totals scale when intervention support changes
  (`density_reference=TRUE`). This differs from a fixed immigrant budget.
- Bootstrap: 512 attempted replicates per fit, full refitting of the matching
  bivariate law, conditional on retained observed pretreatment history;
  SEM uses 500 inner iterations. There are 507 usable naive and 512 usable
  SEM replicates. After excluding failures/explosive refits, the entire
  distribution is translated to have the fitted contrast as its mean.
  Reported quantiles are therefore recentered quantiles. SDs are 75.18 and
  39.52; recentering leaves SDs unchanged.
- Partition sensitivity: 500 inner SEM iterations on grids at 1, 2 and 5
  times the estimated range, plus AOI-versus-remainder. Naive contrasts
  remain negative and SEM contrasts positive; SEM magnitude ranges from
  23.696 to 350.038. No temporal-cutoff sensitivity results are in this run.

The saved configuration is the authoritative record of this fit. Defaults
in `oklahoma_analysis.R` and the exploratory/cluster launch scripts are
configurable and need not reproduce this exact archived run. Starting a
new fit is different from reproducing the published figures from saved fits.
Do not compare likelihoods across temporal truncation settings.

## Prepared inputs

The tracked directory `oklahoma_induced_seismicity_data_regional20150318/`
contains the USGS event snapshot, OCC AOI geometry, projected coordinates,
window geometry and metadata used by the analysis. Refreshing these inputs
from live services can change the data and is unnecessary for reproduction.

`oklahoma_counties_2022.geojson` freezes the 77 Census cartographic county
boundaries returned by `tigris::counties(state="OK", cb=TRUE, year=2022)`.
It supplies the existing paper maps without a live Census download. Source:
[US Census 2022 cartographic boundaries](https://www.census.gov/geographies/mapping-files/time-series/geo/cartographic-boundary.html).
The maps use EPSG:5070. `Oklahoma_data_and_viz.R` is the data-acquisition
script; `oklahoma_analysis.R` performs fitting; `ate_bivariate.R` evaluates
contrasts; `paper/oklahoma_paper_assets.R` is the publication renderer.

See [the archive guide](../zenodo/README.md) for full figure reproduction
and [the cluster guide](../nesi/README.md) for launching new analyses.
