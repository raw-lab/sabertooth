<!-- table:before_after_v04 -->
| | v0.3.0 | v0.4.0 | change |
|---|---|---|---|
| *held-out SCOPe40 (2.6 M residues): 4,000 queries, 2 threads* | | | |
| wall time, s | 12.07 | 12.02 | 1.0× faster |
| peak RSS, MiB | 67 | 58 | 1.2× less |
| sensitivity (AUC1) | 0.2017 | 0.2146 | +6.4% |
| true / false positives at E ≤ 1e-3 | 20,238 / 19 | 21,619 / 19 | +6.8% TP |
| false positives at matched recall 0.08 | 16 | 12 | 25% fewer |
| false positives at matched recall 0.10 | 419 | 89 | 4.7× fewer |
| low-complexity-insert queries: TP / FP at E ≤ 1e-3 | 3922 / 7 | 4127 / 7 |  |
| shuffled queries: chance hits per query at E ≤ 1 (≤ 1 = calibrated) | 0.362 | 0.415 |  |
| *68 M residues (374,323 sequences): 1,000 held-out queries* | | | |
| wall time, s | 21.3 | 7.3 | 2.9× faster |
| peak RSS, MiB | 888 | 472 | 1.9× less |
| sensitivity (AUC1) | 0.1837 | 0.1849 | +0.7% |
| true / false positives at E ≤ 1e-3 | 4889 / 4 | 4999 / 4 |  |
| false positives at matched recall 0.10 | 450 | 367 | 18% fewer |
| makedb + createindex (single build): s / peak RSS MiB / on disk MiB | 23.3 / 1450 / 921 | 23.1 / 1430 / 529 | 1.7× smaller |
| on a 1 GiB machine (cgroup): wall time, s | 29.7 | 7.3 | 4.1× faster |
| *265 M residues (1,454,323 sequences): 1,000 held-out queries* | | | |
| wall time, s | 69.2 | 20.7 | 3.4× faster |
| peak RSS, MiB | 2705 | 1419 | 1.9× less |
| sensitivity (AUC1) | 0.1775 | 0.1777 | +0.1% |
| true / false positives at E ≤ 1e-3 | 4651 / 3 | 4737 / 3 |  |
| false positives at matched recall 0.10 | 1750 | 2243 | **28% more** |
| makedb + createindex (single build): s / peak RSS MiB / on disk MiB | 88.9 / 1522 / 2964 | 99.7 / 1524 / 1705 | 1.7× smaller |
<!-- /table -->

<!-- table:head_to_head -->
| database | Sabertooth 0.4 default: wall s / RSS MiB / AUC1 | MMseqs2 `-s 5.7`: wall s / RSS MiB / AUC1 | Sabertooth speed-up | Sabertooth AUC1 | Sabertooth 0.4 `-s 7.5` | MMseqs2 `-s 7.5` | Sabertooth speed-up | Sabertooth AUC1 |
|---|---|---|---|---|---|---|---|---|
| 2.6 M (SCOPe40) | 3.04 / 54 / 0.2095 | 2.38 / 884 / 0.1577 | 0.78× | +33% | 8.14 / 62 / 0.2418 | 11.00 / 1764 / 0.2302 | 1.35× | +5.1% |
| 17 M | 3.28 / 225 / 0.1957 | 3.20 / 984 / 0.1561 | 0.97× | +25% | 10.70 / 252 / 0.2251 | 21.69 / 1977 / 0.2239 | 2.03× | +0.6% |
| 34 M | 5.02 / 308 / 0.1913 | 4.36 / 1097 / 0.1553 | 0.87× | +23% | 17.08 / 360 / 0.2195 | 31.20 / 2205 / 0.2206 | 1.83× | −0.5% |
| 68 M | 7.30 / 472 / 0.1849 | 6.38 / 1329 / 0.1541 | 0.87× | +20% | 26.71 / 531 / 0.2120 | 39.83 / 2673 / 0.2168 | 1.49× | −2.2% |
| 134 M | 10.79 / 790 / 0.1822 | 10.80 / 1772 / 0.1529 | 1.00× | +19% | 48.68 / 871 / 0.2059 | 49.82 / 3529 / 0.2124 | 1.02× | −3.1% |
| 265 M | 20.66 / 1419 / 0.1777 | 18.15 / 2653 / 0.1519 | 0.88× | +17% | 90.92 / 1495 / 0.2004 | 61.85 / 5259 / 0.2078 | 0.68× | −3.5% |
<!-- /table -->

<!-- table:headline_v04 -->
| | v0.3 default | v0.4 `-s 2` | **v0.4 default** | v0.4 `-s 7.5` | MMseqs2 `-s 5.7` | MMseqs2 `-s 7.5` | DIAMOND `--sensitive` | DIAMOND `--ultra-sensitive` | exhaustive SW |
|---|---|---|---|---|---|---|---|---|---|
| *held-out SCOPe40: 4,000 queries × 14,323 domains* |  |  |  |  |  |  |  |  |  |
| wall time, s (median of 3; exhaustive SW: 1 run) | 12.07 | 2.07 | 12.02 | 32.15 | 8.40 | 41.81 | 7.07 | 20.72 | 127.64 |
| peak RSS, MiB | 67 | 49 | 58 | 66 | 886 | 1772 | 38 | 80 | 29 |
| AUC1 (sensitivity to 1st FP) | 0.2017 | 0.0812 | 0.2146 | 0.2434 | 0.1595 | 0.2335 | 0.1660 | 0.2456 | 0.2636 |
| true positives at E ≤ 1e-3 | 20238 | 6970 | 21619 | 24709 | 16661 | 25936 | 20124 | 29504 | 25854 |
| false positives at E ≤ 1e-3 | 19 | 8 | 19 | 20 | 12 | 18 | 14 | 16 | 22 |
| FP at matched recall 0.07 | 8 | — | 7 | 7 | 18 | 6 | 7 | 6 | 7 |
| FP at matched recall 0.08 | 16 | — | 12 | 7 | 2053 | 6 | 11 | 6 | 7 |
| FP at matched recall 0.10 | 419 | — | 89 | 20 | — | 11 | — | 7 | 17 |
| *low-complexity control (1,000 LC-insert queries × SCOPe40 + decoys)* |  |  |  |  |  |  |  |  |  |
| TP / FP at E ≤ 1e-3 | 3922 / 7 | 1443 / 4 | 4127 / 7 | 4600 / 12 | 3471 / 55 | 5101 / 714 | 4004 / 3 | 5520 / 5 | — |
| *null: 1,000 shuffled queries (expected chance hits per query = E)* |  |  |  |  |  |  |  |  |  |
| chance hits per query at E ≤ 1 | 0.362 | — | 0.415 | 0.610 | 0.076 | — | 0.009 | — | 0.931 |
| chance hits per query at E ≤ 0.01 | 0.0040 | — | 0.0050 | 0.0050 | 0.0010 | — | 0.0020 | — | 0.0090 |
<!-- /table -->

<!-- table:null_v04 -->
| run | workload | hits E≤1e-3 | E≤0.01 | E≤0.1 | E≤1 | observed/expected at E≤0.1 | at E≤1 |
|---|---|---|---|---|---|---|---|
| v0.3 st_s5.7 | null_q | 1 | 4 | 37 | 362 | 0.370 | 0.362 |
| v0.3 st_s7.5 | null_q | 1 | 6 | 55 | 516 | 0.550 | 0.516 |
| st_s5.7 | null_q | 1 | 5 | 44 | 415 | 0.440 | 0.415 |
| st_s7.5 | null_q | 1 | 5 | 61 | 610 | 0.610 | 0.610 |
| st_s5.7_plain | null_q | 2 | 11 | 91 | 666 | 0.910 | 0.666 |
| st_exhaustive | null_q | 1 | 9 | 89 | 931 | 0.890 | 0.931 |
| st_exhaustive_plain | null_q | 2 | 21 | 177 | 1515 | 1.770 | 1.515 |
| mm_s5.7 | null_q | 0 | 1 | 13 | 76 | 0.130 | 0.076 |
| dm_sensitive | null_q | 0 | 2 | 3 | 9 | 0.030 | 0.009 |
| v0.3 st_s5.7 | null_db | 3 | 13 | 169 | 1448 | 0.422 | 0.362 |
| v0.3 st_s7.5 | null_db | 4 | 15 | 229 | 2098 | 0.573 | 0.524 |
| st_s5.7 | null_db | 4 | 11 | 191 | 1610 | 0.477 | 0.403 |
| st_s7.5 | null_db | 5 | 18 | 273 | 2463 | 0.682 | 0.616 |
| st_s5.7_plain | null_db | 8 | 31 | 304 | 2501 | 0.760 | 0.625 |
| st_exhaustive | null_db | 5 | 27 | 356 | 3762 | 0.890 | 0.941 |
| st_exhaustive_plain | null_db | 12 | 53 | 589 | 5772 | 1.472 | 1.443 |
| mm_s5.7 | null_db | 3 | 5 | 35 | 285 | 0.088 | 0.071 |
| dm_sensitive | null_db | 3 | 5 | 23 | 56 | 0.057 | 0.014 |
<!-- /table -->

<!-- table:presets_v04 -->
| `-s` | v0.4 k-mer score / cap | SCOPe40: wall s | AUC1 | TP / FP@1e-3 | v0.3 same `-s`: wall s | AUC1 | 68 M: wall s | AUC1 | TP / FP@1e-3 | v0.3 same `-s`: wall s | AUC1 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 32 / 64 | 0.44 | 0.0664 | 1287 / 1 | 0.91 | 0.1434 | 0.37 | 0.0628 | 1183 / 1 | 7.22 | 0.1362 |
| 2 | 30 / 128 | 0.39 | 0.0810 | 1646 / 1 | 1.14 | 0.1631 | 0.44 | 0.0755 | 1492 / 1 | 10.32 | 0.1537 |
| 4 | 26 / 256 | 0.72 | 0.1503 | 3120 / 1 | 1.50 | 0.1798 | 1.16 | 0.1380 | 2740 / 1 | 13.90 | 0.1685 |
| 5.7 (default) | 22 / 512 | 2.10 | 0.2144 | 4719 / 1 | 2.13 | 0.2032 | 5.05 | 0.1930 | 4130 / 2 | 16.59 | 0.1900 |
| 6.5 | 20 / 512 | 2.98 | 0.2306 | 5003 / 1 | 2.84 | 0.2165 | 8.78 | 0.2078 | 4371 / 2 | 22.99 | 0.2012 |
| 7.5 | 18 / 1024 | 5.35 | 0.2433 | 5228 / 1 | 4.78 | 0.2326 | 18.80 | 0.2186 | 4566 / 2 | 40.91 | 0.2140 |
| 8.5 | 18 / 2048 | 7.27 | 0.2471 | 5304 / 1 | 7.83 | 0.2375 | 27.78 | 0.2219 | 4622 / 2 | 61.31 | 0.2190 |
| 9 | 16 / 4096 | 14.78 | 0.2513 | 5371 / 2 | 9.90 | 0.2376 | 58.55 | 0.2241 | 4684 / 3 | 78.42 | 0.2198 |
<!-- /table -->

<!-- table:scale_limits_v04 -->
| run (68 M residues) | machine | shards | wall s | peak RSS MiB | cgroup peak MiB | AUC1 |
|---|---|---|---|---|---|---|
| v0.4 default | cold page cache | 1 | 4.5 (4.5–4.6) | 473 | — | 0.1930 |
| v0.4 default | 512 MiB RAM, `--memory-limit 384M` | 3 | 7.0 (6.8–7.5) | 261 | 512 | 0.1930 |
| v0.4 default | 256 MiB RAM, `--memory-limit 192M` | 7 | 19.6 (19.5–19.8) | 136 | 256 | 0.1930 |
| v0.3 default | cold page cache | 1 | 16.6 (16.5–16.7) | 896 | — | 0.1900 |
| v0.3 default | 512 MiB RAM, `--memory-limit 384M` | 5 | 28.9 (28.2–28.9) | 231 | 512 | 0.1900 |
| v0.3 default | 256 MiB RAM, `--memory-limit 192M` | 19 | 41.2 (40.2–42.0) | 90 | 256 | 0.1900 |
<!-- /table -->

<!-- table:mem1g_v04 -->
| run (68 M residues, 1 GiB RAM) | memory option | shards | wall s | peak RSS MiB | cgroup peak MiB | AUC1 |
|---|---|---|---|---|---|---|
| v0.4 default | `--memory-limit 1G` | 1 | 7.3 (6.5–7.9) | 475 | 472 | 0.1849 |
| v0.3 default | `--memory-limit 1G` | 2 | 29.7 (29.6–29.8) | 486 | 919 | 0.1837 |
| MMseqs2 `-s 5.7` | none (see text) | — | 259.8 (254.1–265.3) | 1022 | 1024 | 0.1541 |
| DIAMOND `--sensitive` | none (see text) | — | 55.7 (55.5–55.7) | 268 | 336 | 0.1635 |
<!-- /table -->

<!-- table:eval_eval_v04 -->
**held-out homology (4,000 queries; median of 3 runs after a warm-up, exhaustive SW one run)**

| run | wall s | RSS MiB | AUC1 | TP@1e-3 | FP@1e-3 | FDP@1e-3 | FP@E10 | FP@rec 0.06 | FP@rec 0.07 | FP@rec 0.08 | FP@rec 0.10 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| v0.3 st_s5.7 | 12.07 (11.92–12.39) | 67 | 0.2017 | 20238 | 19 | 0.0009 | 16577 | 6 | 8 | 16 | 419 |
| v0.3 st_s7.5 | 29.08 (28.36–29.38) | 80 | 0.2313 | 23709 | 20 | 0.0008 | 24992 | 6 | 7 | 8 | 22 |
| st_s2 | 2.07 (2.00–2.07) | 49 | 0.0812 | 6970 | 8 | 0.0011 | 1965 | — | — | — | — |
| st_s5.7 | 12.02 (11.37–12.02) | 58 | 0.2146 | 21619 | 19 | 0.0009 | 18528 | 6 | 7 | 12 | 89 |
| st_s7.5 | 32.15 (31.59–33.41) | 66 | 0.2434 | 24709 | 20 | 0.0008 | 29544 | 6 | 7 | 7 | 20 |
| mm_s5.7 | 8.40 (8.39–9.06) | 886 | 0.1595 | 16661 | 12 | 0.0007 | 2604 | 7 | 18 | 2053 | — |
| mm_s7.5 | 41.81 (40.71–42.92) | 1772 | 0.2335 | 25936 | 18 | 0.0007 | 12157 | 6 | 6 | 6 | 11 |
| dm_sensitive | 7.07 (7.06–7.20) | 38 | 0.1660 | 20124 | 14 | 0.0007 | 244 | 6 | 7 | 11 | — |
| dm_ultra | 20.72 (19.62–21.31) | 80 | 0.2456 | 29504 | 16 | 0.0005 | 2575 | 6 | 6 | 6 | 7 |
| st_exhaustive | 127.64 | 29 | 0.2636 | 25854 | 22 | 0.0009 | 52142 | 6 | 7 | 7 | 17 |
<!-- /table -->

<!-- table:eval_lc_v04 -->
**low-complexity inserts vs SCOPe40 + decoys (single runs)**

| run | wall s | RSS MiB | AUC1 | TP@1e-3 | FP@1e-3 | FDP@1e-3 | FP@E10 | FP@rec 0.04 | FP@rec 0.06 | FP@rec 0.07 | FP@rec 0.08 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| v0.3 st_s5.7 | 4.17 | 111 | 0.1675 | 3922 | 7 | 0.0018 | 2194 | 3 | 4 | 14 | 45 |
| v0.3 st_s7.5 | 10.93 | 132 | 0.1897 | 4422 | 11 | 0.0025 | 3255 | 3 | 3 | 11 | 22 |
| st_s2 | 0.54 | 81 | 0.0685 | 1443 | 4 | 0.0028 | 293 | — | — | — | — |
| st_s5.7 | 3.36 | 87 | 0.1739 | 4127 | 7 | 0.0017 | 2359 | 3 | 3 | 8 | 19 |
| st_s7.5 | 9.84 | 96 | 0.1962 | 4600 | 12 | 0.0026 | 3580 | 3 | 3 | 5 | 17 |
| mm_s5.7 | 2.65 | 906 | 0.1343 | 3471 | 55 | 0.0156 | 632 | 17 | 90 | 336 | — |
| mm_s7.5 | 13.73 | 1811 | 0.1726 | 5101 | 714 | 0.1228 | 4726 | 100 | 231 | 345 | 573 |
| dm_sensitive | 6.28 | 46 | 0.1437 | 4004 | 3 | 0.0007 | 68 | 3 | 3 | 5 | — |
| dm_ultra | 13.12 | 89 | 0.1936 | 5520 | 5 | 0.0009 | 382 | 3 | 3 | 3 | 3 |
<!-- /table -->

<!-- table:eval_frag_v04 -->
**fragments, 30–60% of the domain (single runs)**

| run | wall s | RSS MiB | AUC1 | TP@1e-3 | FP@1e-3 | FDP@1e-3 | FP@E10 | FP@rec 0.03 | FP@rec 0.04 | FP@rec 0.05 | FP@rec 0.06 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| v0.3 st_s5.7 | 1.18 | 61 | 0.1171 | 2272 | 4 | 0.0018 | 3827 | 4 | 5 | 47 | 621 |
| v0.3 st_s7.5 | 3.52 | 69 | 0.1387 | 2710 | 4 | 0.0015 | 5739 | 3 | 4 | 5 | 39 |
| st_s2 | 0.50 | 46 | 0.0413 | 738 | 3 | 0.0040 | 393 | — | — | — | — |
| st_s5.7 | 1.09 | 50 | 0.1237 | 2325 | 4 | 0.0017 | 4095 | 3 | 5 | 43 | 405 |
| st_s7.5 | 3.44 | 55 | 0.1467 | 2790 | 4 | 0.0014 | 6607 | 3 | 4 | 5 | 29 |
| mm_s5.7 | 0.87 | 884 | 0.0796 | 1827 | 3 | 0.0016 | 517 | 3 | 409 | — | — |
| mm_s7.5 | 4.82 | 1761 | 0.1324 | 2920 | 4 | 0.0014 | 2664 | 3 | 3 | 4 | 31 |
| dm_sensitive | 2.78 | 35 | 0.0929 | 2315 | 3 | 0.0013 | 26 | 3 | 3 | — | — |
| dm_ultra | 5.79 | 56 | 0.1320 | 3422 | 3 | 0.0009 | 264 | 3 | 3 | 3 | 4 |
<!-- /table -->

<!-- table:eval_multi_v04 -->
**two-domain fusions (single runs)**

| run | wall s | RSS MiB | AUC1 | TP@1e-3 | FP@1e-3 | FDP@1e-3 | FP@E10 | FP@rec 0.06 | FP@rec 0.07 | FP@rec 0.08 | FP@rec 0.10 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| v0.3 st_s5.7 | 6.17 | 70 | 0.1402 | 9683 | 12 | 0.0012 | 4422 | 3 | 5 | 12 | 469 |
| v0.3 st_s7.5 | 15.39 | 80 | 0.1621 | 11109 | 12 | 0.0011 | 6637 | 3 | 3 | 6 | 16 |
| st_s2 | 1.31 | 50 | 0.0503 | 3456 | 6 | 0.0017 | 583 | — | — | — | — |
| st_s5.7 | 6.30 | 61 | 0.1488 | 10319 | 11 | 0.0011 | 5045 | 3 | 5 | 8 | 108 |
| st_s7.5 | 15.69 | 71 | 0.1704 | 11548 | 12 | 0.0010 | 7729 | 3 | 3 | 5 | 13 |
| mm_s5.7 | 4.58 | 885 | 0.1059 | 7871 | 6 | 0.0008 | 732 | 5 | 12 | — | — |
| mm_s7.5 | 21.22 | 1767 | 0.1656 | 11991 | 7 | 0.0006 | 3635 | 3 | 3 | 4 | 8 |
| dm_sensitive | 4.66 | 35 | 0.1114 | 9399 | 5 | 0.0005 | 100 | 3 | 4 | 9 | — |
| dm_ultra | 11.58 | 59 | 0.1755 | 13384 | 7 | 0.0005 | 969 | 3 | 3 | 3 | 5 |
<!-- /table -->

<!-- table:series_v04 -->
Queries: `q_eval1k.fa` (sha256 e153440d7e45…), 2 threads unless stated, median of 3 measured runs after 1 warm-up run(s); peak RSS: the largest of the measured runs.

| binary | version | sha256 |
|---|---|---|
| DIAMOND | diamond version 2.2.8 | `d3d72578fa97e87c…` |
| MMseqs2 | 8cc5ce367b5638c4306c2d7cfc652dd099a4643f | `b7ef6e0e33df5dd4…` |
| Sabertooth 0.4 | Sabertooth v0.4.0 | `48c9d9d45bd5a526…` |
| Sabertooth 0.3 | Sabertooth v0.3.0 | `549806e86771203e…` |

### 2.6 M residues (SCOPe40 only)

| configuration | wall s (min–max) | CPU s | peak RSS MiB | AUC1 | TP / FP at E ≤ 1e-3 | FP at recall 0.06 / 0.08 / 0.10 |
|---|---|---|---|---|---|---|
| Sabertooth 0.4 `-s 1` | 0.51 (0.46–0.59) | 0.9 | 48 | 0.0672 | 1626 / 1 | – / – / – |
| Sabertooth 0.4 `-s 2` | 0.54 (0.51–0.56) | 1.0 | 48 | 0.0832 | 1869 / 3 | – / – / – |
| Sabertooth 0.4 `-s 4` | 1.04 (1.03–1.12) | 2.0 | 50 | 0.1432 | 3405 / 4 | 11 / – / – |
| Sabertooth 0.4 default (`-s 5.7`) | 3.04 (2.85–3.22) | 5.8 | 54 | 0.2095 | 5813 / 5 | 3 / 3 / 7 |
| Sabertooth 0.4 `-s 7.5` | 8.14 (7.63–8.66) | 15.5 | 62 | 0.2418 | 6629 / 6 | 3 / 3 / 5 |
| Sabertooth 0.3 default | 3.32 (3.20–3.34) | 6.3 | 64 | 0.1989 | 5443 / 6 | 3 / 3 / 22 |
| MMseqs2 `-s 1` | 0.29 (0.28–0.29) | 0.5 | 884 | 0.0213 | 385 / 0 | – / – / – |
| MMseqs2 `-s 4` | 0.66 (0.66–0.73) | 1.3 | 884 | 0.0917 | 2158 / 2 | – / – / – |
| MMseqs2 default (`-s 5.7`) | 2.38 (2.32–2.46) | 4.4 | 884 | 0.1577 | 4483 / 3 | 3 / 58 / – |
| MMseqs2 `-s 7.5` | 11.00 (10.95–11.07) | 20.3 | 1764 | 0.2302 | 7007 / 5 | 3 / 3 / 3 |
| DIAMOND `--sensitive` | 3.42 (3.37–3.54) | 6.3 | 38 | 0.1647 | 5482 / 4 | 3 / 3 / – |

| database build (single run) | wall s | peak RSS MiB | on disk MiB (all files) | seed index MiB | note |
|---|---|---|---|---|---|
| Sabertooth 0.4 | 0.6 | 62 | 35 | 31.4 | makedb 0.0 s + createindex 0.5 s |
| Sabertooth 0.3 | 0.6 | 63 | 45 | 41.1 | makedb 0.0 s + createindex 0.6 s |
| MMseqs2 | 4.4 | 889 | 882 | – | createdb 0.1 s + createindex 4.4 s |
| DIAMOND | 0.1 | 21 | 3 | – | makedb 0.1 s |

| same `-s` | Sabertooth speed-up | Sabertooth AUC1 vs MMseqs2 | MMseqs2 / Sabertooth peak RSS |
|---|---|---|---|
| `-s 1` | 0.56× | +215.7% | 18.41× |
| `-s 4` | 0.64× | +56.2% | 17.73× |
| `-s 5.7` | 0.78× | +32.8% | 16.28× |
| `-s 7.5` | 1.35× | +5.1% | 28.36× |

### 17 M residues (SCOPe40 + filler)

| configuration | wall s (min–max) | CPU s | peak RSS MiB | AUC1 | TP / FP at E ≤ 1e-3 | FP at recall 0.06 / 0.08 / 0.10 |
|---|---|---|---|---|---|---|
| Sabertooth 0.4 `-s 1` | 0.44 (0.42–0.48) | 0.8 | 214 | 0.0658 | 1538 / 1 | – / – / – |
| Sabertooth 0.4 `-s 2` | 0.51 (0.50–0.52) | 0.9 | 214 | 0.0811 | 1770 / 3 | – / – / – |
| Sabertooth 0.4 `-s 4` | 1.00 (0.98–1.24) | 1.9 | 216 | 0.1365 | 3186 / 3 | 143 / – / – |
| Sabertooth 0.4 default (`-s 5.7`) | 3.28 (3.27–3.43) | 6.4 | 225 | 0.1957 | 5325 / 5 | 3 / 3 / 50 |
| Sabertooth 0.4 `-s 7.5` | 10.70 (10.37–10.79) | 20.3 | 252 | 0.2251 | 6051 / 5 | 3 / 3 / 6 |
| Sabertooth 0.3 default | 8.88 (8.67–9.00) | 16.9 | 271 | 0.1905 | 5119 / 4 | 3 / 3 / 117 |
| MMseqs2 `-s 1` | 0.29 (0.25–0.30) | 0.5 | 980 | 0.0212 | 377 / 0 | – / – / – |
| MMseqs2 `-s 4` | 0.87 (0.85–1.08) | 1.6 | 980 | 0.0912 | 2101 / 2 | – / – / – |
| MMseqs2 default (`-s 5.7`) | 3.20 (3.12–3.21) | 6.2 | 984 | 0.1561 | 4325 / 3 | 3 / 242 / – |
| MMseqs2 `-s 7.5` | 21.69 (21.31–22.44) | 41.2 | 1977 | 0.2239 | 6605 / 4 | 3 / 3 / 3 |
| DIAMOND `--sensitive` | 15.56 (15.53–16.03) | 28.9 | 85 | 0.1640 | 5296 / 5 | 3 / 3 / – |

| database build (single run) | wall s | peak RSS MiB | on disk MiB (all files) | seed index MiB | note |
|---|---|---|---|---|---|
| Sabertooth 0.4 | 4.0 | 287 | 227 | 204.8 | makedb 0.2 s + createindex 3.8 s |
| Sabertooth 0.3 | 3.2 | 292 | 257 | 235.5 | makedb 0.2 s + createindex 3.1 s |
| MMseqs2 | 9.2 | 988 | 1013 | – | createdb 0.3 s + createindex 8.9 s |
| DIAMOND | 0.8 | 39 | 19 | – | makedb 0.8 s |

| same `-s` | Sabertooth speed-up | Sabertooth AUC1 vs MMseqs2 | MMseqs2 / Sabertooth peak RSS |
|---|---|---|---|
| `-s 1` | 0.67× | +210.1% | 4.59× |
| `-s 4` | 0.87× | +49.6% | 4.53× |
| `-s 5.7` | 0.97× | +25.4% | 4.37× |
| `-s 7.5` | 2.03× | +0.6% | 7.85× |

### 34 M residues (SCOPe40 + filler)

| configuration | wall s (min–max) | CPU s | peak RSS MiB | AUC1 | TP / FP at E ≤ 1e-3 | FP at recall 0.06 / 0.08 / 0.10 |
|---|---|---|---|---|---|---|
| Sabertooth 0.4 `-s 1` | 0.47 (0.45–0.47) | 0.9 | 291 | 0.0648 | 1506 / 1 | – / – / – |
| Sabertooth 0.4 `-s 2` | 0.54 (0.53–0.54) | 1.0 | 291 | 0.0800 | 1731 / 3 | – / – / – |
| Sabertooth 0.4 `-s 4` | 1.16 (1.12–1.21) | 2.3 | 294 | 0.1344 | 3110 / 3 | 387 / – / – |
| Sabertooth 0.4 default (`-s 5.7`) | 5.02 (4.76–5.07) | 9.2 | 308 | 0.1913 | 5181 / 4 | 3 / 3 / 115 |
| Sabertooth 0.4 `-s 7.5` | 17.08 (15.87–17.24) | 31.1 | 360 | 0.2195 | 5867 / 4 | 3 / 3 / 8 |
| Sabertooth 0.3 default | 14.98 (14.58–15.81) | 28.3 | 475 | 0.1868 | 5007 / 4 | 3 / 3 / 235 |
| MMseqs2 `-s 1` | 0.26 (0.26–0.26) | 0.5 | 1089 | 0.0212 | 376 / 0 | – / – / – |
| MMseqs2 `-s 4` | 1.02 (0.99–1.14) | 1.9 | 1090 | 0.0911 | 2080 / 2 | – / – / – |
| MMseqs2 default (`-s 5.7`) | 4.36 (4.30–4.43) | 8.3 | 1097 | 0.1553 | 4258 / 3 | 3 / 458 / – |
| MMseqs2 `-s 7.5` | 31.20 (30.89–31.26) | 58.9 | 2205 | 0.2206 | 6457 / 3 | 3 / 3 / 3 |
| DIAMOND `--sensitive` | 28.93 (28.76–29.62) | 54.4 | 146 | 0.1638 | 5216 / 3 | 3 / 3 / – |

| database build (single run) | wall s | peak RSS MiB | on disk MiB (all files) | seed index MiB | note |
|---|---|---|---|---|---|
| Sabertooth 0.4 | 6.8 | 543 | 326 | 282.1 | makedb 0.3 s + createindex 6.5 s |
| Sabertooth 0.3 | 6.6 | 554 | 474 | 430.0 | makedb 0.3 s + createindex 6.3 s |
| MMseqs2 | 15.6 | 1101 | 1164 | – | createdb 0.5 s + createindex 15.1 s |
| DIAMOND | 1.4 | 65 | 37 | – | makedb 1.4 s |

| same `-s` | Sabertooth speed-up | Sabertooth AUC1 vs MMseqs2 | MMseqs2 / Sabertooth peak RSS |
|---|---|---|---|
| `-s 1` | 0.56× | +205.5% | 3.74× |
| `-s 4` | 0.88× | +47.6% | 3.71× |
| `-s 5.7` | 0.87× | +23.2% | 3.56× |
| `-s 7.5` | 1.83× | -0.5% | 6.13× |

### 68 M residues (SCOPe40 + filler)

| configuration | wall s (min–max) | CPU s | peak RSS MiB | AUC1 | TP / FP at E ≤ 1e-3 | FP at recall 0.06 / 0.08 / 0.10 |
|---|---|---|---|---|---|---|
| Sabertooth 0.4 `-s 1` | 0.48 (0.47–0.50) | 0.9 | 449 | 0.0641 | 1459 / 1 | – / – / – |
| Sabertooth 0.4 `-s 2` | 0.64 (0.63–0.66) | 1.2 | 450 | 0.0788 | 1678 / 4 | – / – / – |
| Sabertooth 0.4 `-s 4` | 1.68 (1.62–1.73) | 3.2 | 454 | 0.1299 | 3001 / 4 | – / – / – |
| Sabertooth 0.4 default (`-s 5.7`) | 7.30 (6.93–7.30) | 13.9 | 472 | 0.1849 | 4999 / 4 | 3 / 4 / 367 |
| Sabertooth 0.4 `-s 7.5` | 26.71 (26.60–29.90) | 51.8 | 531 | 0.2120 | 5674 / 4 | 3 / 3 / 14 |
| Sabertooth 0.3 default | 21.28 (20.72–21.33) | 41.3 | 888 | 0.1837 | 4889 / 4 | 3 / 4 / 450 |
| MMseqs2 `-s 1` | 0.30 (0.30–0.36) | 0.5 | 1314 | 0.0212 | 374 / 0 | – / – / – |
| MMseqs2 `-s 4` | 1.30 (1.27–1.32) | 2.4 | 1316 | 0.0906 | 2054 / 2 | – / – / – |
| MMseqs2 default (`-s 5.7`) | 6.38 (6.08–6.39) | 12.3 | 1329 | 0.1541 | 4191 / 3 | 3 / 900 / – |
| MMseqs2 `-s 7.5` | 39.83 (39.72–39.97) | 76.9 | 2673 | 0.2168 | 6286 / 4 | 3 / 3 / 3 |
| DIAMOND `--sensitive` | 55.08 (54.90–55.18) | 105.7 | 272 | 0.1635 | 5156 / 3 | 3 / 3 / – |
| Sabertooth 0.4 default, 1 GiB RAM | 7.27 (6.47–7.91) | 13.3 | 475 (cgroup peak 472) | 0.1849 | 4999 / 4 | 3 / 4 / 367 |
| Sabertooth 0.3 default, 1 GiB RAM | 29.71 (29.61–29.84) | 55.2 | 486 (cgroup peak 919) | 0.1837 | 4889 / 4 | 3 / 4 / 450 |
| MMseqs2 default, 1 GiB RAM | 259.80 (254.11–265.34) | 313.8 | 1022 (cgroup peak 1024) | 0.1541 | 4191 / 3 | 3 / 900 / – |
| DIAMOND `--sensitive`, 1 GiB RAM | 55.67 (55.53–55.73) | 106.0 | 268 (cgroup peak 336) | 0.1635 | 5156 / 3 | 3 / 3 / – |

| database build (single run) | wall s | peak RSS MiB | on disk MiB (all files) | seed index MiB | note |
|---|---|---|---|---|---|
| Sabertooth 0.4 | 23.1 | 1430 | 529 | 440.4 | makedb 0.7 s + createindex 22.4 s |
| Sabertooth 0.3 | 23.3 | 1450 | 921 | 832.6 | makedb 0.6 s + createindex 22.7 s |
| MMseqs2 | 25.2 | 1336 | 1473 | – | createdb 1.0 s + createindex 24.3 s |
| DIAMOND | 3.0 | 110 | 75 | – | makedb 3.0 s |

| same `-s` | Sabertooth speed-up | Sabertooth AUC1 vs MMseqs2 | MMseqs2 / Sabertooth peak RSS |
|---|---|---|---|
| `-s 1` | 0.63× | +202.3% | 2.92× |
| `-s 4` | 0.77× | +43.3% | 2.90× |
| `-s 5.7` | 0.87× | +20.0% | 2.82× |
| `-s 7.5` | 1.49× | -2.2% | 5.03× |

### 134 M residues (SCOPe40 + filler)

| configuration | wall s (min–max) | CPU s | peak RSS MiB | AUC1 | TP / FP at E ≤ 1e-3 | FP at recall 0.06 / 0.08 / 0.10 |
|---|---|---|---|---|---|---|
| Sabertooth 0.4 `-s 1` | 0.58 (0.53–0.63) | 1.1 | 752 | 0.0637 | 1424 / 1 | – / – / – |
| Sabertooth 0.4 `-s 2` | 0.64 (0.64–0.65) | 1.2 | 753 | 0.0779 | 1638 / 3 | – / – / – |
| Sabertooth 0.4 `-s 4` | 1.96 (1.87–2.16) | 3.8 | 757 | 0.1277 | 2923 / 3 | – / – / – |
| Sabertooth 0.4 default (`-s 5.7`) | 10.79 (10.48–10.82) | 20.6 | 790 | 0.1822 | 4866 / 3 | 3 / 5 / 881 |
| Sabertooth 0.4 `-s 7.5` | 48.68 (47.36–50.81) | 93.5 | 871 | 0.2059 | 5494 / 3 | 3 / 3 / 21 |
| Sabertooth 0.3 default | 37.54 (37.40–38.27) | 73.1 | 1526 | 0.1806 | 4760 / 3 | 3 / 5 / 870 |
| MMseqs2 `-s 1` | 0.34 (0.34–0.38) | 0.6 | 1742 | 0.0212 | 373 / 0 | – / – / – |
| MMseqs2 `-s 4` | 1.97 (1.95–2.03) | 3.7 | 1746 | 0.0905 | 2041 / 2 | – / – / – |
| MMseqs2 default (`-s 5.7`) | 10.80 (10.51–10.99) | 21.1 | 1772 | 0.1529 | 4144 / 3 | 3 / – / – |
| MMseqs2 `-s 7.5` | 49.82 (49.42–49.82) | 94.4 | 3529 | 0.2124 | 6140 / 3 | 3 / 3 / 3 |
| DIAMOND `--sensitive` | 108.65 (108.52–108.70) | 209.2 | 491 | 0.1626 | 5066 / 3 | 3 / 3 / – |

| database build (single run) | wall s | peak RSS MiB | on disk MiB (all files) | seed index MiB | note |
|---|---|---|---|---|---|
| Sabertooth 0.4 | 40.9 | 1427 | 920 | 743.1 | makedb 1.3 s + createindex 39.6 s |
| Sabertooth 0.3 | 45.1 | 1425 | 1627 | 1450.2 | makedb 1.2 s + createindex 43.8 s |
| MMseqs2 | 51.5 | 1783 | 2066 | – | createdb 2.0 s + createindex 49.5 s |
| DIAMOND | 6.1 | 198 | 147 | – | makedb 6.1 s |

| same `-s` | Sabertooth speed-up | Sabertooth AUC1 vs MMseqs2 | MMseqs2 / Sabertooth peak RSS |
|---|---|---|---|
| `-s 1` | 0.59× | +199.9% | 2.32× |
| `-s 4` | 1.00× | +41.2% | 2.31× |
| `-s 5.7` | 1.00× | +19.1% | 2.24× |
| `-s 7.5` | 1.02× | -3.1% | 4.05× |

### 265 M residues (SCOPe40 + filler)

| configuration | wall s (min–max) | CPU s | peak RSS MiB | AUC1 | TP / FP at E ≤ 1e-3 | FP at recall 0.06 / 0.08 / 0.10 |
|---|---|---|---|---|---|---|
| Sabertooth 0.4 `-s 1` | 0.82 (0.69–0.82) | 1.5 | 1359 | 0.0633 | 1403 / 1 | – / – / – |
| Sabertooth 0.4 `-s 2` | 0.91 (0.87–1.00) | 1.8 | 1361 | 0.0775 | 1614 / 3 | – / – / – |
| Sabertooth 0.4 `-s 4` | 2.99 (2.98–3.08) | 5.8 | 1369 | 0.1256 | 2865 / 3 | – / – / – |
| Sabertooth 0.4 default (`-s 5.7`) | 20.66 (19.24–20.70) | 37.8 | 1419 | 0.1777 | 4737 / 3 | 3 / 5 / 2243 |
| Sabertooth 0.4 `-s 7.5` | 90.92 (89.59–91.06) | 178.5 | 1495 | 0.2004 | 5342 / 3 | 3 / 3 / 32 |
| Sabertooth 0.3 default | 69.25 (69.01–69.42) | 135.2 | 2705 | 0.1775 | 4651 / 3 | 3 / 7 / 1750 |
| MMseqs2 `-s 1` | 0.55 (0.44–0.61) | 0.8 | 2604 | 0.0212 | 370 / 0 | – / – / – |
| MMseqs2 `-s 4` | 3.11 (3.08–3.26) | 6.0 | 2613 | 0.0897 | 2018 / 2 | – / – / – |
| MMseqs2 default (`-s 5.7`) | 18.15 (18.14–18.29) | 35.5 | 2653 | 0.1519 | 4081 / 3 | 3 / – / – |
| MMseqs2 `-s 7.5` | 61.85 (61.61–62.41) | 117.1 | 5259 | 0.2078 | 5977 / 3 | 3 / 3 / 5 |
| DIAMOND `--sensitive` | 215.57 (215.02–216.11) | 414.7 | 932 | 0.1615 | 4990 / 3 | 3 / 3 / – |
| Sabertooth 0.4 default, 1 thread | 35.79 (35.56–35.83) | 35.4 | 1397 | 0.1777 | 4737 / 3 | 3 / 5 / 2243 |
| MMseqs2 default, 1 thread | 38.71 (38.41–38.98) | 38.3 | 2628 | 0.1519 | 4081 / 3 | 3 / – / – |

| database build (single run) | wall s | peak RSS MiB | on disk MiB (all files) | seed index MiB | note |
|---|---|---|---|---|---|
| Sabertooth 0.4 | 99.7 | 1524 | 1705 | 1350.1 | makedb 2.9 s + createindex 96.7 s |
| Sabertooth 0.3 | 88.9 | 1522 | 2964 | 2609.5 | makedb 2.3 s + createindex 86.5 s |
| MMseqs2 | 107.0 | 2660 | 3257 | – | createdb 4.0 s + createindex 103.0 s |
| DIAMOND | 11.8 | 376 | 293 | – | makedb 11.8 s |

| same `-s` | Sabertooth speed-up | Sabertooth AUC1 vs MMseqs2 | MMseqs2 / Sabertooth peak RSS |
|---|---|---|---|
| `-s 1` | 0.67× | +198.7% | 1.92× |
| `-s 4` | 1.04× | +40.0% | 1.91× |
| `-s 5.7` | 0.88× | +17.0% | 1.87× |
| `-s 7.5` | 0.68× | -3.5% | 3.52× |

### Matched sensitivity

Wall time each tool needs to reach the other default's AUC1, interpolated (log wall time linear in AUC1) between that tool's measured `-s` settings; “>” = not reached by its most sensitive measured setting, “≤” = already exceeded by its fastest.

| database | MMseqs2 `-s 5.7`: AUC1 / wall s | Sabertooth 0.4 at that AUC1, s | | Sabertooth 0.4 default: AUC1 / wall s | MMseqs2 at that AUC1, s | |
|---|---|---|---|---|---|---|
| 2.6 M | 0.1577 / 2.38 | 1.31 | Sabertooth 1.8× faster | 0.2095 / 3.04 | 7.11 | MMseqs2 2.3× slower |
| 17 M | 0.1561 / 3.20 | 1.49 | Sabertooth 2.2× faster | 0.1957 / 3.28 | 9.79 | MMseqs2 3.0× slower |
| 34 M | 0.1553 / 4.36 | 1.99 | Sabertooth 2.2× faster | 0.1913 / 5.02 | 12.92 | MMseqs2 2.6× slower |
| 68 M | 0.1541 / 6.38 | 3.21 | Sabertooth 2.0× faster | 0.1849 / 7.30 | 15.70 | MMseqs2 2.2× slower |
| 134 M | 0.1529 / 10.80 | 4.31 | Sabertooth 2.5× faster | 0.1822 / 10.79 | 22.89 | MMseqs2 2.1× slower |
| 265 M | 0.1519 / 18.15 | 7.92 | Sabertooth 2.3× faster | 0.1777 / 20.66 | 32.00 | MMseqs2 1.5× slower |

### Scaling with database size

Wall time (s) per size; exponent b of wall ∝ residues^b (least squares over 17, 34, 68, 134, 265 M).

| configuration | 2.6 M | 17 M | 34 M | 68 M | 134 M | 265 M | b |
|---|---|---|---|---|---|---|---|
| Sabertooth 0.4 `-s 1` | 0.51 | 0.44 | 0.47 | 0.48 | 0.58 | 0.82 | 0.21 |
| Sabertooth 0.4 `-s 2` | 0.54 | 0.51 | 0.54 | 0.64 | 0.64 | 0.91 | 0.20 |
| Sabertooth 0.4 `-s 4` | 1.04 | 1.00 | 1.16 | 1.68 | 1.96 | 2.99 | 0.39 |
| Sabertooth 0.4 default (`-s 5.7`) | 3.04 | 3.28 | 5.02 | 7.30 | 10.79 | 20.66 | 0.65 |
| Sabertooth 0.4 `-s 7.5` | 8.14 | 10.70 | 17.08 | 26.71 | 48.68 | 90.92 | 0.78 |
| Sabertooth 0.3 default | 3.32 | 8.88 | 14.98 | 21.28 | 37.54 | 69.25 | 0.73 |
| MMseqs2 `-s 1` | 0.29 | 0.29 | 0.26 | 0.30 | 0.34 | 0.55 | 0.22 |
| MMseqs2 `-s 4` | 0.66 | 0.87 | 1.02 | 1.30 | 1.97 | 3.11 | 0.47 |
| MMseqs2 default (`-s 5.7`) | 2.38 | 3.20 | 4.36 | 6.38 | 10.80 | 18.15 | 0.64 |
| MMseqs2 `-s 7.5` | 11.00 | 21.69 | 31.20 | 39.83 | 49.82 | 61.85 | 0.37 |
| DIAMOND `--sensitive` | 3.42 | 15.56 | 28.93 | 55.08 | 108.65 | 215.57 | 0.96 |

| configuration | 2.6 M | 17 M | 34 M | 68 M | 134 M | 265 M |
|---|---|---|---|---|---|---|
| Sabertooth 0.4 default (`-s 5.7`) — AUC1 | 0.2095 | 0.1957 | 0.1913 | 0.1849 | 0.1822 | 0.1777 |
| Sabertooth 0.3 default — AUC1 | 0.1989 | 0.1905 | 0.1868 | 0.1837 | 0.1806 | 0.1775 |
| MMseqs2 default (`-s 5.7`) — AUC1 | 0.1577 | 0.1561 | 0.1553 | 0.1541 | 0.1529 | 0.1519 |
| Sabertooth 0.4 `-s 7.5` — AUC1 | 0.2418 | 0.2251 | 0.2195 | 0.2120 | 0.2059 | 0.2004 |
| MMseqs2 `-s 7.5` — AUC1 | 0.2302 | 0.2239 | 0.2206 | 0.2168 | 0.2124 | 0.2078 |
| DIAMOND `--sensitive` — AUC1 | 0.1647 | 0.1640 | 0.1638 | 0.1635 | 0.1626 | 0.1615 |
| Sabertooth 0.4 default (`-s 5.7`) — peak RSS MiB | 54 | 225 | 308 | 472 | 790 | 1419 |
| Sabertooth 0.3 default — peak RSS MiB | 64 | 271 | 475 | 888 | 1526 | 2705 |
| MMseqs2 default (`-s 5.7`) — peak RSS MiB | 884 | 984 | 1097 | 1329 | 1772 | 2653 |
| Sabertooth 0.4 `-s 7.5` — peak RSS MiB | 62 | 252 | 360 | 531 | 871 | 1495 |
| MMseqs2 `-s 7.5` — peak RSS MiB | 1764 | 1977 | 2205 | 2673 | 3529 | 5259 |
| DIAMOND `--sensitive` — peak RSS MiB | 38 | 85 | 146 | 272 | 491 | 932 |
<!-- /table -->

