# Data and code availability (JRSS-B and general use)

This file is part of the **PPDisentangle** R package. After installation, you can open it in R with:

```r
file.show(system.file("DATA_AND_CODE_AVAILABILITY.md", package = "PPDisentangle"))
```

Or browse it in the package source at `inst/DATA_AND_CODE_AVAILABILITY.md`.

---

## Statement for copy-paste into a manuscript

**Data and code availability.** Earthquake event data for Oklahoma were obtained from the U.S. Geological Survey (USGS) Earthquake Hazards Program event web service (`https://earthquake.usgs.gov/fdsnws/event/1/query`) for the spatio-temporal window, magnitude threshold, and geographic extent described in the manuscript and recorded in the analysis metadata (see `inst/oklahoma/oklahoma_induced_seismicity_data_regional20150318/metadata.json` in the **PPDisentangle** repository). Regulatory Area-of-Interest (AOI) geometry for Oklahoma Corporation Commission directive AOI_20150318 was retrieved from the Oklahoma Corporation Commission public GIS service (`https://gis.occ.ok.gov/server/rest/services/PUBLIC/DIRECTIVE_AOIs/MapServer`). County boundaries used for spatial partitioning were taken from the U.S. Census Bureau (accessed in reproducible form via the R package **tigris**).

The code implementing the bivariate ETAS models, stochastic EM labelling, kernel-density background estimation, treatment-effect summaries, sensitivity analyses, and reporting workflows is provided as the R package **PPDisentangle** (this repository). Application-specific scripts and documentation for the Oklahoma study live under `inst/oklahoma/` (including `oklahoma_analysis.R`, `run_nesi.sh`, and Quarto report sources). A processed event table, AOI geometry, and configuration metadata sufficient to reproduce the Oklahoma analyses are bundled under `inst/oklahoma/oklahoma_induced_seismicity_data_regional20150318/`.

Software is implemented in **R**; key dependencies include **spatstat**, **sf**, **Rcpp**, and other packages declared in the package `DESCRIPTION`. The package is distributed under the **MIT** license (see the `LICENSE` file and `DESCRIPTION`). Third-party data services (USGS, OCC) are subject to their respective terms of use.

**Archival version.** For citation and long-term access, replace this sentence in your paper with a persistent identifier for a tagged release (e.g. GitHub release or Zenodo DOI) once you have published one.

---

## Pointing collaborators and reviewers to the package

- **Package root / clone URL:** use the Git repository where **PPDisentangle** is hosted (add the URL to your paper’s `DESCRIPTION` `URL` field when ready).
- **Oklahoma case study:** `inst/oklahoma/README.md`
- **Paper figures / LaTeX fragments:** `inst/oklahoma/paper/README.md`
