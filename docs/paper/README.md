# Manuscript context (Overleaf is canonical)

`paper.tex`, `revision.tex`, `robustness.tex`, and other files in this directory
are **local context snapshots** for developers and AI agents. They are
**gitignored** and must not be committed or pushed. The live manuscript is edited
and compiled on **Overleaf**.

Generated figures and `\input`-able table fragments for Overleaf are produced by:

```bash
Rscript inst/oklahoma/paper/oklahoma_paper_assets.R
Rscript inst/sim_study/plot_time_sweep_publication.R
```

Outputs land under `../PPDisentangle-output/` (see `inst/OUTPUT.md`). Copy those
assets into your Overleaf project as needed.

This snapshot may lag behind Overleaf; update it occasionally when notation or
figure references change materially.
