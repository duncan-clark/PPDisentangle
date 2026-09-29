# Reproduce the paper outputs

Use PPDisentangle **v0.2.0** with the updated companion results deposit.
Software version history: https://doi.org/10.5281/zenodo.21221021.
Results v2: https://doi.org/10.5281/zenodo.23040032.
The July 29 v0.1.0 / results v1 pair predates the current Oklahoma fits.

## Install and reproduce

Clone the software at tag `v0.2.0`, install it and the plotting dependencies,
and unpack the results archive outside the source checkout:

```r
# From the cloned package directory:
# install.packages("remotes")
remotes::install_local(".", dependencies = TRUE)
```

```bash
export PPDISENTANGLE_OUTPUT_ROOT=/path/to/unpacked/PPDisentangle-output
Rscript inst/zenodo/reproduce_paper_figures.R
Rscript inst/zenodo/validate_paper_outputs.R "$PPDISENTANGLE_OUTPUT_ROOT"
```

The reproduction script regenerates figures and tables from saved fits;
it does not run the simulations, SEM fits, or bootstrap again. The illustrative
Hawkes realisation is a frozen PDF copied from the archive. The manuscript
and bibliography are maintained in Overleaf and are not distributed or
compiled by this repository. A local `docs/paper/revised.tex` snapshot is
ignored by Git. The optional generated robustness TeX is a reporting aid,
not the manuscript source.

## Frozen inputs and conventions

| Component | Path below the output root |
|---|---|
| Main simulation, job 5228509 | `sim_study/paper/main_5228509/` |
| Robustness, 67 scenarios | `sim_study/paper/robustness_merged_tcal/` |
| Oklahoma, job 8804859 + ATE backfill | `oklahoma/for_paper.rds` |
| Regenerated simulation assets | `sim_study/generated/` |
| Regenerated Oklahoma assets | `oklahoma/paper/generated/` |

The 67 scenarios comprise separation (7), background intensity (5),
separation × spatial range (25), kernel misspecification (16), assignment (5),
effect modification (3), and geometry (6). Excluding the 25-cell spatial
surface leaves the 42 scenarios in the six other families.

Oklahoma reproduction explicitly selects **observed-vs-none** over 100 days,
including the partition table. Legacy top-level ATE fields remain
all-or-nothing for provenance; they are not the paper estimates. Expected
saved counts are -51.508 (naive) and +175.200 (SEM). Of 512 bootstrap attempts,
507 and 512 are retained, respectively. Bootstrap distributions and all
reported quantiles are recentered to the corresponding fitted contrast.
See [the Oklahoma guide](../oklahoma/README.md) for the full archived recipe.

Simulation robustness uses the **all-or-nothing** contrast and the archived
fitting/generation conventions. The known control parameters assist SEM
initialization and label updates as well as contrast evaluation. The count
benchmark in label proposals is a stationary mean-rate approximation.
Component immigrant rates in the simulations are total budgets on each
component's assigned support; changing support changes immigrant density.
Nominal count calibration does not guarantee identical realized sample sizes.

The univariate all-treated simulator iterates over components present in
the post-treatment assignment; hence it omits control-component descendants
from pretreatment events in that world. Its simulation reference ('truth')
uses the same convention. These frozen outputs therefore do not implement
an intervention retaining those control descendants. The Oklahoma bivariate
counterfactual code retains the full supplied pretreatment history. This
release preserves the fitted scientific results and makes that distinction
explicit; it does not rerun the simulation study under a different intervention.

## Provenance and packaging

Each archive contains `release_manifest.json`, a per-file `MD5SUMS` list,
and a reproduction-session snapshot. The session snapshot is diagnostic,
not a dependency lockfile; installing packages today may select other versions.
The repository's `inst/zenodo/sessionInfo.txt` is a reference snapshot.

```bash
# From the package root, using the canonical local saved outputs:
bash inst/zenodo/pack_zenodo.sh --out /path/to/results.tar.gz
```

Packaging copies frozen inputs to a temporary tree, regenerates all derived
assets there, checks the canonical counts/settings and scenario count, then
writes provenance/checksums and compresses the tree. It excludes local refreshes,
logs, obsolete reports, LaTeX build products and macOS metadata. It does not
modify the source result files. Set `RSCRIPT` or `GIT` when these executables
are not on the default path. Publish only a package built from the committed
release source. Keep the software tag and results version together.
