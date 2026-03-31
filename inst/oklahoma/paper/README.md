# Oklahoma paper assets (LaTeX + R)

This folder separates **writing** from **computed outputs**:

| File | Role |
|------|------|
| `oklahoma_application_section.tex` | Main-text LaTeX section (`sec:application`). |
| `oklahoma_application_appendix.tex` | Supplement listing HTML-report material not in the main section. |
| `oklahoma_paper_assets.R` | Builds all Oklahoma figures and `\input`-able table `.tex` fragments. |
| `generated/` | Populated by the R script (tables, CSV summary, manifest). Do not edit by hand. |

## Build

From the **repository root**:

```bash
Rscript inst/oklahoma/paper/oklahoma_paper_assets.R
```

Options:

```text
--input       Path to oklahoma_results*.rds (default: output/oklahoma/oklahoma_results_job5226010.rds)
--plots-dir   PDF output directory (default: plots/oklahoma)
--tex-dir     Generated LaTeX directory (default: inst/oklahoma/paper/generated)
--data-dir    Oklahoma CSV bundle (default: inst/oklahoma/oklahoma_induced_seismicity_data_regional20150318)
```

Or in R: `setwd("<repo>"); source("inst/oklahoma/paper/oklahoma_paper_assets.R")`.

## LaTeX layout

- Figures are written under `plots/oklahoma/` (e.g. `ATE_diff.pdf`, `cumulative_count.pdf`, plus ECDF / running-mean / simulation histogram PDFs).
- `oklahoma_application_section.tex` uses paths relative to the repo root as the LaTeX working directory. If your main document lives elsewhere, set `\graphicspath` or change the `\includegraphics` and `\input` paths.

The preamble of your main document should include `\usepackage{booktabs}` for the generated tables.

## Supersedes

The previous `oklahoma_publication_ef.Rmd` / `oklahoma_publication_ef.R` workflow is replaced by this layout.
