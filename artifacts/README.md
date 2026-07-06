# Artifacts

This directory holds versioned publication artifacts that are generated from
analysis scripts but should not live inside the R package source tree.

| Path | Contents | Regeneration command |
|------|----------|----------------------|
| `oklahoma/paper/generated/` | Oklahoma paper figures, LaTeX table fragments, CSV summaries, and build manifest. | `Rscript inst/oklahoma/paper/oklahoma_paper_assets.R` |
| `sim_study/generated/` | Simulation-study paper figures and LaTeX table fragments. | `Rscript inst/sim_study/plot_time_sweep_publication.R` |
| `paper/paper.tex` | Archived manuscript TeX snapshot with paths adjusted to this artifact layout. | Maintained externally, then copied here when needed. |

Local run outputs belong in `output/` or `cluster_output/`; those directories
are ignored by git.
