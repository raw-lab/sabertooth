# Benchmark suite `null`

Versions:

- sabertooth: `Sabertooth v0.3.0`
- baseline: `Sabertooth v0.2.0`
- mmseqs: `MMseqs2 8cc5ce367b5638c4306c2d7cfc652dd099a4643f`
- diamond: `diamond version 2.2.8`
- sabertooth_exh: `Sabertooth v0.3.0 (easysearch: exhaustive Smith-Waterman, no prefilter)`

## null_q


| run | wall s | peak RSS MiB | hits E≤1e-3 | E≤1e-2 | E≤0.1 | E≤1 | obs/exp E≤0.1 | obs/exp E≤1 |
|---|---|---|---|---|---|---|---|---|
| v02_default | 9.30 | 294 | 4 | 31 | 208 | 1429 | 2.080 | 1.429 |
| st_s5.7 | 1.32 | 59 | 1 | 4 | 37 | 362 | 0.370 | 0.362 |
| st_s7.5 | 3.71 | 67 | 1 | 6 | 55 | 516 | 0.550 | 0.516 |
| st_s5.7_plain | 1.91 | 64 | 2 | 10 | 72 | 569 | 0.720 | 0.569 |
| mm_s5.7 | 1.44 | 881 | 0 | 1 | 13 | 76 | 0.130 | 0.076 |
| dm_sensitive | 2.03 | 28 | 0 | 2 | 3 | 9 | 0.030 | 0.009 |
| st_exhaustive | 28.75 | 15 | 1 | 9 | 89 | 931 | 0.890 | 0.931 |
| st_exhaustive_plain | 27.81 | 18 | 2 | 21 | 177 | 1515 | 1.770 | 1.515 |

## null_db


| run | wall s | peak RSS MiB | hits E≤1e-3 | E≤1e-2 | E≤0.1 | E≤1 | obs/exp E≤0.1 | obs/exp E≤1 |
|---|---|---|---|---|---|---|---|---|
| v02_default | 29.93 | 296 | 18 | 113 | 796 | 5765 | 1.990 | 1.441 |
| st_s5.7 | 5.68 | 62 | 3 | 13 | 169 | 1448 | 0.422 | 0.362 |
| st_s7.5 | 14.06 | 74 | 4 | 15 | 229 | 2098 | 0.573 | 0.524 |
| st_s5.7_plain | 6.18 | 64 | 9 | 28 | 274 | 2174 | 0.685 | 0.543 |
| mm_s5.7 | 4.52 | 886 | 3 | 5 | 35 | 285 | 0.088 | 0.071 |
| dm_sensitive | 4.23 | 38 | 3 | 5 | 23 | 56 | 0.057 | 0.014 |
| st_exhaustive | 109.90 | 16 | 5 | 27 | 356 | 3762 | 0.890 | 0.941 |
| st_exhaustive_plain | 114.99 | 19 | 12 | 53 | 589 | 5772 | 1.472 | 1.443 |

FP@rec x = false positives at the first E-value threshold reaching recall x (n/r = not reached).
