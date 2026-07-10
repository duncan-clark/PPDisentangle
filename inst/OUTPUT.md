# Generated outputs (outside this repository)

All analysis outputs live in a **sibling directory** next to this git checkout:

```text
/path/to/
  PPDisentangle/           ← this repo (code + scripts only)
  PPDisentangle-output/    ← Oklahoma RDS, sim study RDS, paper figures, logs
```

Default path: `../PPDisentangle-output/` relative to the package root.

Override for a custom location:

```bash
export PPDISENTANGLE_OUTPUT_ROOT=/path/to/my/outputs
```

## Layout

| Path under `PPDisentangle-output/` | Contents | Regeneration |
|------|----------|--------------|
| `oklahoma/` | Full analysis RDS, exploratory plots, rendered reports (`oklahoma_report.html`) | `inst/oklahoma/oklahoma_analysis.R` or `Rscript inst/oklahoma/render_oklahoma_report.R` |
| `oklahoma/plots/` | Exploratory PNGs from the Oklahoma pipeline | same |
| `oklahoma/paper/generated/` | Paper figures (PDF), LaTeX table fragments | `Rscript inst/oklahoma/paper/oklahoma_paper_assets.R` |
| `sim_study/` | Cluster/local run RDS and logs | `bash inst/sim_study/run_nesi.sh` |
| `sim_study/generated/` | Simulation-study paper figures and tables | `Rscript inst/sim_study/plot_time_sweep_publication.R` |

## R helpers

After loading the package:

```r
PPDisentangle::pp_output_root()
PPDisentangle::pp_output_path("oklahoma")
```

## Zenodo

- **Software:** GitHub release of PPDisentangle.
- **Paper outputs:** one companion archive with the canonical simulation and
  Oklahoma artifacts (see `inst/zenodo/`).

Canonical jobs in the archive:

| Component | Identifier |
|-----------|------------|
| Main simulation | `sim_study/paper/main_5228509/` (`*_summary_FOR_PAPER.rds` + raw `tmp*.rds`) |
| Robustness appendix | `sim_study/paper/robustness_7568933/` |
| Oklahoma application | `oklahoma/for_paper.rds` |

Manuscript robustness appendix text: Overleaf / `docs/paper/revision.tex` (local
snapshot is gitignored).

```bash
# pack when ready (excludes local-only refreshes, logs, duplicate summary.rds,
# and robustness_standalone.pdf preview builds)
bash inst/zenodo/pack_zenodo.sh

# reproduce figures after unpacking the Zenodo tarball
export PPDISENTANGLE_OUTPUT_ROOT=/path/to/PPDisentangle-output
Rscript inst/zenodo/reproduce_paper_figures.R
```

Session pin: `inst/zenodo/sessionInfo.txt` (refresh with
`Rscript inst/zenodo/capture_session.R`).

A context snapshot of the manuscript lives in `docs/paper/paper.tex` (Overleaf
is canonical; the repo copy is not compiled here).

## NeSI cluster workflow

Submit from your laptop, compute on Mahuika, fetch results locally:

```bash
bash inst/nesi/submit.sh sim_study --sims 100
bash inst/nesi/fetch.sh sim_study JOB_ID
```

Cluster outputs default to `/nesi/project/uoo04008/PPDisentangle-output`.
See [`inst/nesi/README.md`](nesi/README.md).
