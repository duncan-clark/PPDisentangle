# Oklahoma paper assets builder

Prose for the manuscript lives on **Overleaf**. A context snapshot is in
`docs/paper/paper.tex`. This folder holds the **R driver** only; regenerated
publication outputs go under `../PPDisentangle-output/oklahoma/paper/generated/`.

| Path | Role |
|------|------|
| `oklahoma_paper_assets.R` | Builds Oklahoma figures (PDFs) and `\input`-able table `.tex` fragments. |
| `../../../PPDisentangle-output/oklahoma/paper/generated/` | Tables, CSV summary, manifest RDS — **do not edit by hand**. |
| `../../../PPDisentangle-output/oklahoma/paper/generated/figures/` | PDFs from the same run. |

## Build

From the **repository root**:

```bash
Rscript inst/oklahoma/paper/oklahoma_paper_assets.R
```

Options:

```text
--input       Path to results `.rds` (default: first existing of
              ../PPDisentangle-output/oklahoma/for_paper.rds,
              inst/oklahoma/paper/for_paper.rds)
--plots-dir   PDF output directory (default: PPDisentangle-output/.../figures)
--tex-dir     Generated LaTeX directory (default: PPDisentangle-output/.../generated)
--data-dir    Oklahoma CSV bundle (default: inst/oklahoma/oklahoma_induced_seismicity_data_regional20150318)
```

The preamble of your main document should include `\usepackage{booktabs}` for the generated tables.
