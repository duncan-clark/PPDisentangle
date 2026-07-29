# Zenodo companion archive helpers

Paper figures are regenerated from frozen simulation/application outputs that
live **outside** the git repo (Zenodo deposit).

## Canonical paper jobs

| Component | Path under `PPDisentangle-output/` |
|-----------|--------------------------------------|
| Main simulation | `sim_study/paper/main_5228509/` |
| Illustrative Hawkes realisation | `.../simulated_hawkes_hawkes_process.pdf` (`fig:pp_realiz`) |
| Robustness appendix | `sim_study/paper/robustness_merged_tcal/` |
| Oklahoma | `oklahoma/for_paper.rds` |

The robustness archive is the production **67-scenario** suite:
time-calibrated K-separation / SNR (`7767397`), K×spatial-range surface
(`7839871`), and the remaining families from `7762623` / earlier resume
sources. Legacy `paper/robustness_7568933/` is superseded and not packed.

Local-only material lives under `sim_study/local/` and is not packed.

Manuscript robustness appendix text is maintained in Overleaf /
`docs/paper/revision.tex` (gitignored local snapshot). Generated
`simulation_robustness_appendix.tex` is optional.

## Reproduce figures from a Zenodo unpack

```bash
export PPDISENTANGLE_OUTPUT_ROOT=/path/to/unpacked/PPDisentangle-output
Rscript inst/zenodo/reproduce_paper_figures.R
```

This regenerates main-sim figures, robustness figure PDFs (+ optional appendix
`.tex` fragment), and Oklahoma paper assets. It does **not** compile
`robustness.pdf`.

A `sessionInfo()` snapshot is written beside the outputs and compared against
`inst/zenodo/sessionInfo.txt` in this repo.

## Session / dependency pin

- Package dependencies: `DESCRIPTION`
- Frozen R session snapshot: [`sessionInfo.txt`](sessionInfo.txt)
- Refresh the snapshot after intentional dependency upgrades:

```bash
Rscript inst/zenodo/capture_session.R
```

## Pack a Zenodo tarball (maintainer; when ready)

```bash
bash inst/zenodo/pack_zenodo.sh
```
