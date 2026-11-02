# A/B runs (host B)

Written by `bench/run_ab.py` (see `bench/README.md` for the commands): 68 M
residues (`series_68M.fa`, default index), the 1,000 tuning queries, 2
threads; one unmeasured run per variant, then 3 alternating rounds.
`<set>/<variant>.<round>.json` = measure.py record + `--stats-json` of that
run; `<set>/<variant>.eval.json` = evaluation of the variant's hits;
`<set>/summary.json` = medians, ranges, hit hashes.

| set | variants | result (median wall, AUC1) |
|---|---|---|
| `simd` | `avx2` (default) / `avx512` (`SABERTOOTH_SIMD`) | 4.62 / 4.90 s, identical hits |
| `evidence` | `window48` (default) / `whole_diagonal` (`--diag-window 0`) / `xdrop30` (`--ungapped-xdrop 30`, v0.3's evidence) | 4.77 / 5.27 / 7.05 s; AUC1 0.1930 / 0.1934 / 0.1932 |
| `prefetch` | development builds with / without the hit-loop posting prefetch (as in `../dev/ab_*`), on their own full-offset index | 4.92 / 5.12 s, identical hits |

`mmseqs_split/`: MMseqs2's memory options in a 1 GiB memory cgroup on the
68 M database with (`mm68`) and without (`mmnoidx68`) its precomputed index,
held-out queries, one run each: `--split-memory-limit 768M | 1G` exits at once
with "Cannot fit databases into …" (both databases); `--split 0` and
`--split 4` without the index are killed at the limit after about 7 s;
`--split 4` with the index completes (`mm68_split4.json`).

`mem1g_hostb/`: the 1 GiB machine on host B with the series protocol (one
warm-up outside the limit, then 3 runs in the cgroup with the database
evicted before each; held-out queries `q_eval1k.fa`): Sabertooth 0.4
`--memory-limit 1G`, MMseqs2 `-s 5.7` without an option and with `--split 4`
(`*.measure.json`, evaluations in `*.eval.json`): Sabertooth 4.8 s
(4.7–4.9, cgroup peak 469 MiB, AUC1 0.1849); MMseqs2 163.9 s (162.8–167.6)
without an option and 166.1 s (165.6–166.7) with `--split 4`, AUC1 0.1541
both (same hits, different order).
