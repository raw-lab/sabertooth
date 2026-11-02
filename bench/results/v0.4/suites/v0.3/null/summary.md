# Benchmark suite `null`

Versions:

- sabertooth: `Sabertooth v0.3.0`
- sabertooth_exh: `Sabertooth v0.3.0 (easysearch: exhaustive Smith-Waterman, no prefilter)`

## null_q


| run | wall s | peak RSS MiB | hits E≤1e-3 | E≤1e-2 | E≤0.1 | E≤1 | obs/exp E≤0.1 | obs/exp E≤1 |
|---|---|---|---|---|---|---|---|---|
| st_s5.7 | 2.32 | 59 | 1 | 4 | 37 | 362 | 0.370 | 0.362 |
| st_s7.5 | 6.44 | 67 | 1 | 6 | 55 | 516 | 0.550 | 0.516 |
| st_s5.7_plain | 2.94 | 64 | 2 | 10 | 72 | 569 | 0.720 | 0.569 |

## null_db


| run | wall s | peak RSS MiB | hits E≤1e-3 | E≤1e-2 | E≤0.1 | E≤1 | obs/exp E≤0.1 | obs/exp E≤1 |
|---|---|---|---|---|---|---|---|---|
| st_s5.7 | 8.82 | 62 | 3 | 13 | 169 | 1448 | 0.422 | 0.362 |
| st_s7.5 | 18.03 | 74 | 4 | 15 | 229 | 2098 | 0.573 | 0.524 |
| st_s5.7_plain | 8.78 | 64 | 9 | 28 | 274 | 2174 | 0.685 | 0.543 |

FP@rec x = false positives at the first E-value threshold reaching recall x (n/r = not reached).
