# A/B Benchmark Report

- Baseline: `/Users/dac6/Library/CloudStorage/GoogleDrive-clark.a.duncan@gmail.com/My Drive/my_laptop/Documents/Academics/computation/r_code/CausalPointProcess/PPDisentangle/inst/perf/results/base_quick.json`
- Candidate: `/Users/dac6/Library/CloudStorage/GoogleDrive-clark.a.duncan@gmail.com/My Drive/my_laptop/Documents/Academics/computation/r_code/CausalPointProcess/PPDisentangle/inst/perf/results/cand_quick.json`
- Minimum required speedup: `1.000`

## Results

| Component | Baseline (s) | Candidate (s) | Speedup | Base OK | Cand OK |
|---|---:|---:|---:|:---:|:---:|
| likelihood | 0.0290 | 0.0290 | 1.000 | yes | yes |
| simulation | 0.0300 | 0.0320 | 0.938 | yes | yes |
| sem | 0.6660 | 0.6980 | 0.954 | yes | yes |

## Verdict

- Correctness pass: **TRUE**
- Speedup pass: **FALSE**
- Overall pass: **FALSE**
