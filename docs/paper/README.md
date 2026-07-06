# Manuscript context (Overleaf is canonical)

`paper.tex` here is a **read-only context snapshot** for developers and AI
agents working on this repository. The live manuscript is edited and compiled on
**Overleaf** — do not treat this file as something to build from the repo.

Generated figures and `\input`-able table fragments for Overleaf are produced by:

```bash
Rscript inst/oklahoma/paper/oklahoma_paper_assets.R
Rscript inst/sim_study/plot_time_sweep_publication.R
```

Outputs land under `../PPDisentangle-output/` (see `inst/OUTPUT.md`). Copy those
assets into your Overleaf project as needed.

This snapshot may lag behind Overleaf; update it occasionally when notation or
figure references change materially.
