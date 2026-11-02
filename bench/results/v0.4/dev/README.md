# Development records (kept as found)

`--stats-json` outputs of development runs that some v0.4 statements cite. They
were written during development, not by a benchmark driver: each file is one
run, its `command` field holds only the subcommand, and the settings are given
here. Host by file time (`../HOSTS.md`). The decisions they supported were
re-measured with shipped drivers where the code still allows it
(`../ab/`, `../layout/`).

| file | host | what | settings |
|---|---|---|---|
| `ab_pf.{1,2,3}.json`, `ab_nopf.{1,2,3}.json`, `ab_pf.sh` | A | hit-loop posting prefetch on / off, alternating runs | two development builds (reporting version 0.3.0) differing only in `PF_POSTINGS` (8 / 0, `src/prefilter.rs`); 68 M database, spaced-seed index with full offsets; tuning queries; `--spaced 111001101 --kmer-score 22 --max-kmers 512 --min-diag-score 41`, 2 threads (`ab_pf.sh`). Re-measured on host B: `../ab/prefetch/` |
| `l_s17.json` | A | 17 M residues, sparse k-mer table | development build (reports version 0.3.0, before the version bump) with that day's gate 37 (the release uses 38 at 17 M: 11,207 vs 11,059 hits), k-mer score 22, tuning queries, 2 threads; `index.source` names the index |
| `l_s17_direct.json` | A | 17 M residues, direct k-mer table with full 4-byte offsets (before the 16-bit grouped offsets) | as above (same hits) |
| `l40.json` | A | SCOPe40 (2.6 M) with a direct table (full offsets) | development build, gate 33 (as the release at 2.6 M); the default layout there is sparse |
| `simd_avx2.json`, `simd_avx512.json` | B | AVX2 vs `SABERTOOTH_SIMD=avx512`, one run each | a v0.4.0 build (reports 0.4.0), 68 M, tuning queries, defaults (gate 41). Superseded by `../ab/simd/` (3 alternating rounds) |

Not kept (cited in `CHANGES.md` / `CHECKLIST-v0.4.md` as development
measurements without a record): the two-hit seeding replay — neither the
script version that applied the cited rules nor its output (it read a
per-diagonal candidate trace written by a development build; the release's
`SABERTOOTH_TRACE_CANDIDATES` trace has prefilter scores only) —, posting
prefetch at seed creation, reading evidence lanes straight from the targets,
and the 17 M / 68 M grouped-vs-full-offset timing ranges other than the one
run above.
