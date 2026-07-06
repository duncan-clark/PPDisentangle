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
| `docs/paper/` | Overleaf manuscript snapshot for context (not compiled from this repo). |
| `../PPDisentangle-output/` | **Local only** — all analysis RDS, logs, and paper assets; see [`inst/OUTPUT.md`](inst/OUTPUT.md). |

## Oklahoma application

See [`inst/oklahoma/README.md`](inst/oklahoma/README.md).

## Simulation study

See [`inst/sim_study/README.md`](inst/sim_study/README.md).

## Output policy

The repository tracks **source code, tests, analysis scripts, and small prepared inputs** only. Generated results live in a sibling folder `PPDisentangle-output/` (or `PPDISENTANGLE_OUTPUT_ROOT`) so they survive git clones, branch switches, and repo restructures. Archive large result sets on Zenodo separately from the software release; see [`inst/OUTPUT.md`](inst/OUTPUT.md).

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
