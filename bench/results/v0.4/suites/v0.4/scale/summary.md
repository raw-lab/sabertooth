# Benchmark suite `scale`

Versions:

- sabertooth: `Sabertooth v0.4.0`
- sabertooth_exh: `Sabertooth v0.4.0 (easysearch: exhaustive Smith-Waterman, no prefilter)`

## scale


| run | wall s (min–max, n) | peak RSS MiB | AUC1 | TP@1e-3 | FP@1e-3 | FDP@1e-3 | FP@E10 | FP@rec 0.02 | FP@rec 0.03 | FP@rec 0.04 | FP@rec 0.05 | FP@rec 0.06 | FP@rec 0.07 | FP@rec 0.08 | FP@rec 0.10 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| st_s5.7_cold | 4.46 (4.45–4.58, 3) | 473 | 0.1930 | 4130 | 2 | 0.0005 | 4695 | 0 | 0 | 1 | 1 | 1 | 5 | 277 | n/r |
| st_s5.7_mem512m | 6.99 (6.83–7.45, 3) | 261 | 0.1930 | 4130 | 2 | 0.0005 | 4695 | 0 | 0 | 1 | 1 | 1 | 5 | 277 | n/r |
| st_s5.7_mem256m | 19.62 (19.53–19.78, 3) | 136 | 0.1930 | 4130 | 2 | 0.0005 | 4695 | 0 | 0 | 1 | 1 | 1 | 5 | 277 | n/r |

FP@rec x = false positives at the first E-value threshold reaching recall x (n/r = not reached).
