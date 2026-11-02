# Benchmark suite `scale`

Versions:

- sabertooth: `Sabertooth v0.3.0`
- sabertooth_exh: `Sabertooth v0.3.0 (easysearch: exhaustive Smith-Waterman, no prefilter)`

## scale


| run | wall s (min–max, n) | peak RSS MiB | AUC1 | TP@1e-3 | FP@1e-3 | FDP@1e-3 | FP@E10 | FP@rec 0.02 | FP@rec 0.03 | FP@rec 0.04 | FP@rec 0.05 | FP@rec 0.06 | FP@rec 0.07 | FP@rec 0.08 | FP@rec 0.10 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| st_s5.7_cold | 16.62 (16.50–16.74, 3) | 896 | 0.1900 | 4125 | 2 | 0.0005 | 4846 | 0 | 0 | 1 | 1 | 1 | 5 | 186 | n/r |
| st_s5.7_mem512m | 28.90 (28.24–28.90, 3) | 231 | 0.1900 | 4125 | 2 | 0.0005 | 4846 | 0 | 0 | 1 | 1 | 1 | 5 | 186 | n/r |
| st_s5.7_mem256m | 41.15 (40.22–41.95, 3) | 90 | 0.1900 | 4125 | 2 | 0.0005 | 4846 | 0 | 0 | 1 | 1 | 1 | 5 | 186 | n/r |

FP@rec x = false positives at the first E-value threshold reaching recall x (n/r = not reached).
