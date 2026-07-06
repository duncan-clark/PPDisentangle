# PPDisentangle

R package: **causal inference for spatiotemporal point processes** using a stochastic EM approach with ETAS/Hawkes outcomes, process labelling, and treatment-effect summaries (including spillover).

## Repository layout

| Path | Contents |
|------|----------|
| `R/` | R package source: SEM, labelling, Hawkes/ETAS fitting, ATE summaries, plotting, and standard errors. |
| `src/` | Rcpp implementations for simulation and likelihood speedups. |
| `man/` | Generated package documentation from roxygen2. |
| `tests/testthat/` | Unit and consistency tests for simulation, likelihoods, labelling, ATE, and standard errors. |
| `inst/oklahoma/` | Oklahoma induced-seismicity application scripts, report source, and small prepared input snapshot. |
| `inst/sim_study/` | Simulation-study drivers and cluster helpers. |
| `artifacts/` | Versioned publication artifacts separated from the package source tree. |
| `output/`, `cluster_output/` | Local analysis/HPC outputs; ignored by git. |

## Oklahoma application

See [`inst/oklahoma/README.md`](inst/oklahoma/README.md).

## Simulation study

See [`inst/sim_study/README.md`](inst/sim_study/README.md).

## Artifact policy

Source code, tests, package documentation, analysis scripts, and the small prepared Oklahoma input snapshot are kept in the package tree. Generated publication figures/tables live under `artifacts/` so they can be versioned without being bundled into package builds. Local run outputs, rendered reports, compiled objects, and R check directories are ignored.

## Installation (development tree)

From the package directory:

```r
# install.packages("devtools")
devtools::install()
```

## Testing

```r
devtools::test()
```

License: **MIT** — see `DESCRIPTION` and `LICENSE`.
