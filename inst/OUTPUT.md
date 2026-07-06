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
| `oklahoma/` | Full analysis RDS, exploratory plots, rendered reports | `inst/oklahoma/oklahoma_analysis.R` |
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
- **Simulation results:** tarball of `PPDisentangle-output/sim_study/`.
- **Oklahoma results:** tarball of `PPDisentangle-output/oklahoma/`.

A context snapshot of the manuscript lives in `docs/paper/paper.tex` (Overleaf
is canonical; the repo copy is not compiled here).
