# Manuscript context (Overleaf is canonical)

`paper.tex`, `revision.tex`, and other files in this directory
are **local context snapshots** for developers and AI agents. They are
**gitignored** and must not be committed or pushed. The live manuscript is edited
and compiled on **Overleaf**.

Robustness appendix figures and the Overleaf-ready fragment are produced by:

```bash
Rscript inst/sim_study/sim_study_robustness.R --sims 32 --target-points 2500
```

Outputs:
- `PPDisentangle-output/sim_study/generated/robustness/figures/` (PDFs)
- `PPDisentangle-output/sim_study/generated/robustness/simulation_robustness_appendix.tex`

Copy the `.tex` file and PDFs to `plots/sim_study/robustness/` on Overleaf, then
`\input{plots/sim_study/robustness/simulation_robustness_appendix.tex}` inside
`\section{Additional simulation results}` in `revision.tex`.

Local PDF preview (from `generated/robustness/`):

```bash
pdflatex -jobname=robustness_standalone '\def\robustnessstandalone{}\input{simulation_robustness_appendix.tex}'
```

Other generated assets:

```bash
Rscript inst/oklahoma/paper/oklahoma_paper_assets.R
Rscript inst/sim_study/plot_time_sweep_publication.R
```

Outputs land under `../PPDisentangle-output/` (see `inst/OUTPUT.md`). Copy those
assets into your Overleaf project as needed.

This snapshot may lag behind Overleaf; update it occasionally when notation or
figure references change materially.
