# How the v0.4 defaults were chosen

Everything here uses the **tuning** queries (`q_tune.fa`, 1,000 queries from
the tuning superfamilies); the held-out `q_eval*.fa` queries were used only for
the final measurements (`bench/results/v0.4/`). Records of every run below
are in `bench/results/v0.4/tuning/` (hit files not kept) and
`bench/results/v0.4/experiments.jsonl` (the development experiment log, one
JSON line per run: command arguments, wall/CPU per stage, counters, AUC1,
TP/FP).

## 1. Exact gate replay (`tune.sh`, `../tune_gate.py`)

A stricter ungapped gate admits a *subset* of the candidates of a looser one,
and a candidate's gapped alignment does not depend on which other candidates
exist. So one search at a low gate with
`SABERTOOTH_TRACE_CANDIDATES=<file>` (one line per candidate: query, target,
ungapped score, DP cells) plus its hit file gives the exact output of every
higher gate: `tune_gate.py` replays gates `G0..46` and reports candidates and
DP cells per query, AUC1, TP/FP at E ≤ 1e-3 and FP at matched recall for each.

```
tune.sh GRID INDEX EXTRA_IDS|- G0 "T:C T:C ..." [--spaced MASK]
```

runs one traced search per (k-mer score T, cap C) and replays the gates.

## 2. Cost model (`frontier.py`, `frontier_cmp.py`, `gate_table.py`, `gate_lambda.py`)

Wall time of a configuration is modelled from its counters and per-unit CPU
costs measured on host A (`../ENVIRONMENT.md`) without concurrent load:
`CPU ≈ k-mers·c_k + postings·c_p + diagonals·c_d + DP cells·c_c + tracebacks·c_t`
(68 M residues: c_k ≈ 70 ns, c_p ≈ 50–55 ns, c_d ≈ 114 ns, c_c ≈ 0.131 ns,
c_t ≈ 90 µs; the gapped term is exact per gate from the replay). Grids whose
wall times were disturbed by other jobs are compared on these counts, never
on their wall times; every timing in the README tables is a separate clean run
(the 265 M stage profile quoted under Scaling is from the traced tuning run
`s265_sp9/t22_c512.stats.json`).

* `frontier.py DIR c_k c_p c_d c_c c_t` — Pareto frontier (modelled CPU vs AUC1)
* `frontier_cmp.py c... DIR1 DIR2` — best AUC1 within CPU budgets, two grids side by side
  (contiguous `scale` vs spaced `scale_sp9`: the spaced seed is ahead by 0.007–0.027
  AUC1 at every budget from 3 to 100 CPU-s per 1,000 queries; SCOPe40: 0.002–0.021).
  All costs in ns, c_t too: `frontier_cmp.py 70 52 114 0.131 90000 ../results/v0.4/tuning/scale ../results/v0.4/tuning/scale_sp9`
* `gate_table.py c... DIR T:C,...` — per-gate cost and AUC1 of chosen configurations
* `gate_lambda.py LAMBDAS T:C,...` — the gate maximising `AUC1 − λ·ln(CPU)` at
  2.6, 17, 68 and 265 M residues (the evidence for `size_adjusted_gate`)

## 3. Decisions

| decision | evidence (tuning queries) |
|---|---|
| spaced seed `111001101` | screen of 55 weight-6 masks of span 7–9 on SCOPe40 at equal settings (`experiments.jsonl`, names `sp2_*`): 111001101 AUC1 0.2189 vs contiguous 0.2078, same time; confirmed on full grids at 2.6 M and 68 M (`*_sp9` vs contiguous) |
| segment evidence within 48 residues of the seed hits | 68 M, v0.3 settings (`w*` records): whole diagonal AUC1 0.1905 / 7.35 s evidence CPU; window 16: 0.1880 / 4.97; 32: 0.1897 / 5.60; **48: 0.1901 / 6.05**; 64: 0.1902 / 6.43 |
| two-hit seeding rejected | development replay of the 68 M candidate trace (script and output not kept; it needs a per-diagonal trace written by a development build): AUC1 0.190 → 0.163 (two hits on a diagonal), 0.163 (two within 40 residues), 0.170 (two hits or one k-mer ≥ 26) |
| presets (`SEQUENCE_PRESETS`) | points of the 68 M spaced-seed frontier (`scale_sp9`), cost roughly doubling per step; T16/C2048 adds nothing over T18/C2048 (cap-bound), T16/C4096 +0.002 AUC1 for 2× the cost (top preset only) |
| gate 42 at 10^8 residues, −2.5 per e-fold below, +1.0 above | `gate_lambda.py 0.02,0.025,0.03 22:512,24:256`: optimal gates 31–34 (2.6 M), 35–38 (17 M), 38–42 (68 M), 41–42 (265 M); the rule gives 33 / 38 / 41 / 43 (at 265 M one above the optimum, slightly toward speed) |
| MMseqs2 reference frontier | `mm_frontier.sh` (`mm_scope40/`, `mm_scale/`, `mm_frontier.log`): `-s 1…8.5` on the same queries |
