# Performance Review and Test Harness

This folder contains scripts for a broad speedup review workflow without changing core algorithms.

## Hotspot Coverage

- **SEM:** proposal scoring, weight calculation, and outer parameter optimization cost.
- **Simulation/Likelihood:** Hawkes simulation and likelihood precompute effects.
- **Decode:** end-to-end quick-check timing path (optional).

## Scripts

- `benchmark_hotspots.R`
  - Runs component benchmarks and correctness guards.
  - Writes JSON output for one run.
- `compare_benchmarks.R`
  - Compares baseline vs candidate JSON files.
  - Reports speedup and correctness pass/fail.
- `run_ab_perf.sh`
  - Creates temporary worktrees for two refs.
  - Runs benchmark suite on each and compares.

## Typical Usage

Run a single benchmark file on current checkout:

```bash
Rscript inst/perf/benchmark_hotspots.R --quick --out=inst/perf/results/current.json
```

Compare two refs end-to-end:

```bash
bash inst/perf/run_ab_perf.sh main HEAD --quick --min-speedup 1.00
```

Include decode timing in the suite:

```bash
bash inst/perf/run_ab_perf.sh main HEAD --quick --run-decode --min-speedup 1.05
```

## Pass Criteria

- **Correctness:** every component reports `correctness_pass = TRUE`.
- **Speed:** candidate `primary_time_sec` must satisfy `baseline / candidate >= min-speedup`.

The harness is designed so that when optimization changes are introduced later, you can show:
1. measurable speedups, and
2. no functional breakage.
