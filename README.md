# PPDisentangle

R package: **causal inference for spatiotemporal point processes** using a stochastic EM approach with ETAS/Hawkes outcomes, process labelling, and treatment-effect summaries (including spillover).

## Data and code availability (for papers and reproducibility)

A ready-to-use **Data and Code Availability** paragraph (e.g. for *JRSS-B*), plus file paths inside this package, is in:

**[`inst/DATA_AND_CODE_AVAILABILITY.md`](inst/DATA_AND_CODE_AVAILABILITY.md)**

After installation, open it from R:

```r
file.show(system.file("DATA_AND_CODE_AVAILABILITY.md", package = "PPDisentangle"))
```

Add your public repository URL and, when available, a Zenodo DOI or release tag in the manuscript; the bundled text marks where to do that.

## Oklahoma application

See [`inst/oklahoma/README.md`](inst/oklahoma/README.md).

## Installation (development tree)

From the package directory:

```r
# install.packages("devtools")
devtools::install()
```

License: **MIT** — see `DESCRIPTION` and `LICENSE`.
