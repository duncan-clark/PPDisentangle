# A/B Benchmark Report

- Baseline: `/Users/dac6/Library/CloudStorage/GoogleDrive-clark.a.duncan@gmail.com/My Drive/my_laptop/Documents/Academics/computation/r_code/CausalPointProcess/PPDisentangle/inst/perf/results/base_full.json`
- Candidate: `/Users/dac6/Library/CloudStorage/GoogleDrive-clark.a.duncan@gmail.com/My Drive/my_laptop/Documents/Academics/computation/r_code/CausalPointProcess/PPDisentangle/inst/perf/results/cand_full.json`
- Minimum required speedup: `0.000`

## Results

| Component | Baseline (s) | Candidate (s) | Speedup | Base OK | Cand OK |
|---|---:|---:|---:|:---:|:---:|
| likelihood | 0.1130 | 0.1140 | 0.991 | yes | yes |
| simulation | 0.0460 | 0.0440 | 1.045 | yes | yes |
| sem | 2.0230 | 2.0250 | 0.999 | yes | yes |

## Verdict

- Correctness pass: **TRUE**
- Speedup pass: **TRUE**
- Overall pass: **TRUE**
