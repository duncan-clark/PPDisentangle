# Oklahoma paper assets (builder + generated outputs)

Prose for the manuscript lives in your main document (e.g. Overleaf). This folder holds the **R driver** and **regenerated artefacts** only:

| Path | Role |
|------|------|
| `oklahoma_paper_assets.R` | Builds Oklahoma figures (PDFs) and `\input`-able table `.tex` fragments. |
| `generated/` | Tables, CSV summary, manifest RDS — **do not edit by hand**. |
| `generated/figures/` | PDFs from the same run (partition, point patterns, SEM traces, bootstrap plots, …). |

## Build

From the **repository root**:

```bash
Rscript inst/oklahoma/paper/oklahoma_paper_assets.R
```

Options:

```text
--input       Path to results `.rds` (default: first existing of
              output/oklahoma/for_paper.rds, inst/oklahoma/paper/for_paper.rds)
--plots-dir   PDF output directory (default: inst/oklahoma/paper/generated/figures)
--tex-dir     Generated LaTeX directory (default: inst/oklahoma/paper/generated)
--data-dir    Oklahoma CSV bundle (default: inst/oklahoma/oklahoma_induced_seismicity_data_regional20150318)
```

Or in R: `setwd("<repo>"); source("inst/oklahoma/paper/oklahoma_paper_assets.R")`.

## LaTeX / Overleaf

- Default figure output: `generated/figures/` (alongside `generated/*.tex`).
- Point `\includegraphics` and `\input` at copies of those paths in your project, or set `\graphicspath` if the layout differs.
- Optional: `--plots-dir` to write PDFs elsewhere (e.g. a local `figures/` mirror).

The preamble of your main document should include `\usepackage{booktabs}` for the generated tables.

## Supersedes

The previous `oklahoma_publication_ef.Rmd` / `oklahoma_publication_ef.R` workflow is replaced by this layout.
