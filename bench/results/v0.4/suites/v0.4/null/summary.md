# Benchmark suite `null`

Versions:

- sabertooth: `Sabertooth v0.4.0`
- mmseqs: `MMseqs2 8cc5ce367b5638c4306c2d7cfc652dd099a4643f`
- diamond: `diamond version 2.2.8`
- sabertooth_exh: `Sabertooth v0.4.0 (easysearch: exhaustive Smith-Waterman, no prefilter)`

## null_q


| run | wall s | peak RSS MiB | hits E≤1e-3 | E≤1e-2 | E≤0.1 | E≤1 | obs/exp E≤0.1 | obs/exp E≤1 |
|---|---|---|---|---|---|---|---|---|
| st_s5.7 | 2.15 | 48 | 1 | 5 | 44 | 415 | 0.440 | 0.415 |
| st_s7.5 | 6.90 | 53 | 1 | 5 | 61 | 610 | 0.610 | 0.610 |
| st_s5.7_plain | 2.76 | 63 | 2 | 11 | 91 | 666 | 0.910 | 0.666 |
| mm_s5.7 | 1.30 | 884 | 0 | 1 | 13 | 76 | 0.130 | 0.076 |
| dm_sensitive | 2.85 | 34 | 0 | 2 | 3 | 9 | 0.030 | 0.009 |
| st_exhaustive | 31.24 | 15 | 1 | 9 | 89 | 931 | 0.890 | 0.931 |
| st_exhaustive_plain | 33.30 | 18 | 2 | 21 | 177 | 1515 | 1.770 | 1.515 |

## null_db


| run | wall s | peak RSS MiB | hits E≤1e-3 | E≤1e-2 | E≤0.1 | E≤1 | obs/exp E≤0.1 | obs/exp E≤1 |
|---|---|---|---|---|---|---|---|---|
| st_s5.7 | 7.77 | 51 | 4 | 11 | 191 | 1610 | 0.477 | 0.403 |
| st_s7.5 | 25.58 | 57 | 5 | 18 | 273 | 2463 | 0.682 | 0.616 |
| st_s5.7_plain | 8.77 | 63 | 8 | 31 | 304 | 2501 | 0.760 | 0.625 |
| mm_s5.7 | 5.95 | 886 | 3 | 5 | 35 | 285 | 0.088 | 0.071 |
| dm_sensitive | 5.76 | 37 | 3 | 5 | 23 | 56 | 0.057 | 0.014 |
| st_exhaustive | 123.88 | 17 | 5 | 27 | 356 | 3762 | 0.890 | 0.941 |
| st_exhaustive_plain | 126.55 | 19 | 12 | 53 | 589 | 5772 | 1.472 | 1.443 |

FP@rec x = false positives at the first E-value threshold reaching recall x (n/r = not reached).
