# Manuscript context (Overleaf is canonical)

`paper.tex`, `revision.tex`, and other files in this directory
are **local context snapshots** for developers and AI agents. They are
**gitignored** and must not be committed or pushed. The live manuscript is edited
and compiled on **Overleaf**.

**Canonical robustness appendix prose for the paper is the inlined block in
`revision.tex`** (label `app:simulation_robustness`). Local helpers such as
`robustness.tex` and auto-generated `simulation_robustness_appendix.tex` are
optional drafts / figure-pipeline byproducts — do not treat them as the
manuscript source of truth.

Robustness appendix figures are produced by:

```bash
Rscript inst/sim_study/sim_study_robustness.R --replot robustness_7568933
# or a fresh run:
Rscript inst/sim_study/sim_study_robustness.R --sims 32 --target-points 2500
```

Outputs:
- `PPDisentangle-output/sim_study/generated/robustness/figures/` (PDFs)

Copy the PDFs to `plots/sim_study/robustness/` on Overleaf (paths used by
`\RobustnessIncludeFig{...}` in `revision.tex`).

Other generated assets:

```bash
Rscript inst/oklahoma/paper/oklahoma_paper_assets.R
Rscript inst/sim_study/plot_time_sweep_publication.R
Rscript inst/zenodo/reproduce_paper_figures.R
```

Outputs land under `../PPDisentangle-output/` (see `inst/OUTPUT.md`). Copy those
assets into your Overleaf project as needed.

This snapshot may lag behind Overleaf; update it occasionally when notation or
figure references change materially.
