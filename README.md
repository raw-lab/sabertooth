# Sabertooth

```
   ___       _           _              _   _
  / __| __ _| |__  ___ _| |_ ___  ___ _| |_| |_
  \__ \/ _` | '_ \/ -_)  _/ _ \/ _ \  _|  _| ' \
  |___/\__,_|_.__/\___|\__\___/\___/\__|\__|_||_|
       \V/   \V/    fast · light · profile-sensitive
```

Pure-Rust protein homology search in the MMseqs2 mould: sequence and
profile (PSSM) search against a persistent, memory-mapped seed index, with
exact Smith–Waterman–Gotoh alignments and calibrated gapped E-values. One
static binary, three runtime dependencies (`rayon`, `flate2`, `memmap2`).

**v0.4** rebuilds the prefilter for large databases — a spaced seed, SIMD
ungapped evidence, packed seed hits, an ungapped gate that follows the
database size, and a smaller index format — and measures it head-to-head
against MMseqs2 on databases from 2.6 M to 265 M residues. See
[`CHANGES.md`](CHANGES.md) for what changed and why (and what was measured
and rejected), and [`CHECKLIST-v0.4.md`](CHECKLIST-v0.4.md) for the
item-by-item status with evidence. (v0.3, the correctness / false-positive
overhaul of v0.2: [`CHECKLIST.md`](CHECKLIST.md).)

**Before / after** (v0.3.0 vs v0.4.0 defaults; both binaries re-run on the
same 2-vCPU host (host A), 2 threads, held-out queries, median of 3 runs; every number
comes from a run record in `bench/results/v0.4/` — see
[Benchmarks](#benchmarks) for the method and what did *not* improve):

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

**Against MMseqs2 as the database grows** (SCOPe40 + synthetic filler,
1,000 held-out queries, each tool on its own prebuilt index, E ≤ 10 and
≤ 4,000 hits, 2 threads, host A, median of 3 runs; wall s / peak RSS MiB / AUC1):

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

* **Matched sensitivity, less time.** At every size Sabertooth reaches the
  AUC1 of MMseqs2's default 1.8–2.5× faster, and MMseqs2 needs 1.5–3.0×
  longer than Sabertooth's default to reach its AUC1 (interpolated between
  each tool's measured `-s` settings; table "Matched sensitivity" below).
* **Defaults (`-s 5.7` both).** Sabertooth's AUC1 is 17–33% higher and its
  peak memory 1.9–16× lower; MMseqs2's default is up to 28% faster (on the
  2.6 M SCOPe40 alone; 14–15% at 34, 68 and 265 M, within 3% at 17 and
  134 M).
* **Top setting (`-s 7.5`).** Sabertooth is 1.0–2.0× faster up to 134 M
  residues and uses 3.5–28× less memory, but MMseqs2's high-sensitivity
  search scales better with the database (time ∝ residues^0.37 vs ^0.78): at
  265 M MMseqs2 is 1.5× faster and its AUC1 3.7% higher (see limitations).
* **vs v0.3.** 2.7–3.5× faster from 17 M residues up (1.0–1.1× on SCOPe40),
  1.2–1.9× less memory, an index 1.15–1.95× smaller, and at least as sensitive
  at every size (AUC1 +0.1% at 265 M to +5.3% on SCOPe40), with equal or
  fewer false positives at matched recall 0.06–0.10 up to 68 M — but by
  recall 0.10 it admits about as many as v0.3 at 134 M (881 vs 870) and 28%
  more at 265 M (2,243 vs 1,750; see limitations). On a 1 GiB machine at 68 M residues: 7.3 s vs v0.3's 29.7 s
  and MMseqs2's 260 s.

AUC1 = fraction of each query's true homologs ranked above its first false
positive, averaged over queries. "Matched recall 0.08" = the false positives
admitted by the time a tool has found 8% of all true homolog pairs.

---

## Build

```bash
cargo build --release                       # sabertooth binary in target/release/
cargo build --release --features distributed  # + HydraMPP distributed backend
cargo test --release                        # unit + CLI integration tests
cargo test --release --features distributed # + distributed parity tests
./target/release/sabertooth doctor          # self-checks
```

The crate declares `rust-version = "1.89"` (AVX-512 intrinsics); it was
developed and tested with rustc 1.97 (see "Limitations"). No `target-cpu`
flags are needed: AVX2 is detected at run time (SSE2/SSSE3 are the fallbacks;
AVX-512 kernels exist but are opt-in with `SABERTOOTH_SIMD=avx512` because
they were not faster on either benchmark CPU — 6% slower on host B; other
architectures use the exact scalar kernels).

## Quick start

```bash
sabertooth makedb proteins.fasta.gz protDB          # sequence database (MMseqs2 format)
sabertooth createindex protDB --memory-limit 2G     # persistent seed index protDB.stidx
sabertooth search queries.fasta protDB --out hits.m8 --stats-json stats.json
```

## Databases and indexes — two separate steps

| step | command | writes | notes |
|---|---|---|---|
| sequence database | `makedb <in.fasta[.gz]\|db> <db>` | `<db>`, `<db>.index`, `.dbtype`, `<db>_h`, `.lookup`, `.source` | MMseqs2-compatible amino-acid DB (readable by `mmseqs`) |
| seed index | `createindex <db\|fasta> [out.stidx]` | `<db>.stidx` | persistent, versioned, checksummed k-mer index + compact residue arena |

* **Inputs** for every command: FASTA, gzip FASTA (streamed), `-` (stdin),
  an MMseqs2 database prefix, and a `.stidx` seed index (its sequences and
  identifiers; an index keeps no header descriptions). (A database made
  by `mmseqs createdb` gives the same hits as its FASTA, but `createdb`
  shuffles entries by default, so ties may be listed in another order;
  `--shuffle 0` keeps them identical.) Empty inputs,
  malformed FASTA (with the offending line), corrupt gzip, compressed / split /
  profile / nucleotide MMseqs2 databases and inconsistent `.index` files are
  rejected with an explanation instead of being searched as "no hits".
* **Index use.** `search` uses `<target>.stidx` automatically when it exists,
  was built with the search's own seed and masking (its defaults — the spaced
  seed `111001101` for sequence search, k = 5 for profile search, masking on
  — or the explicit `--k` / `--spaced` / `--mask`), and is not stale (size and
  mtime of the source files are recorded — no rescan). An index built
  differently is never adopted silently, so results do not depend on which
  indexes exist (a v0.3 index, contiguous k = 6, is used with `--k 6`).
  Otherwise it builds a temporary index (removed afterwards) and says so on
  stderr (and in `--stats-json` → `index.source`); `--require-index` turns
  that into an error. `indexinfo <file> --verify`
  checks every section checksum.
* **Build memory.** `createindex` streams sequences into the output file and
  builds postings in partitions sized to half of `--memory-limit` (default:
  a 1 GiB partition buffer); the postings, k-mer table and residues are
  identical for any budget (`small_budget_forces_multiple_partitions_with_identical_output`).
  Writes are atomic (temp file + rename): an interrupted build never leaves an
  index that opens. Build time, peak memory and size at every database size
  of the series are in the [scaling tables](#benchmarks); builds under
  1 GiB / 512 MiB / 256 MiB limits were measured for v0.3
  (`bench/results/build.json`, other host) and the partitioning is unchanged.
* **Size.** 5.3 bytes per residue for the 265 M-residue benchmark database
  (1,350 MiB), 6.8 for 68 M (440 MiB; v0.3: 833 MiB), 12.6 for the 2.6 M
  SCOPe40 (31 MiB): residues (1 byte), postings bit-packed into 4 bytes
  (target and position), and either the sparse k-mer table (8 bytes per
  distinct k-mer, small databases) or the direct table (16-bit offsets for
  all 20⁶ codes plus group bases: 132 MB, from ~10 M residues), whichever
  `--index-layout auto` picks.

## Searching

```bash
sabertooth search  <query> <target> [options]   # prefiltered sequence search
sabertooth easysearch <query> <target>          # exhaustive Smith–Waterman (no prefilter)
sabertooth profilesearch <msa> <target>         # PSSM search (per-profile calibrated)
```

### Sensitivity presets

`-s` sets the similar-k-mer threshold and the per-window k-mer cap together
(explicit `--kmer-score`, `--max-kmers` override); the ungapped gate follows
the database size (`--min-diag-score` fixes it). The anchors are points of the
measured speed–sensitivity frontier on the tuning queries against the
68 M-residue database (`bench/tuning/README.md`). Measured on host B (2.1 GHz;
see Benchmarks) with the tuning queries (1,000), single runs, 2 threads,
against SCOPe40 (2.6 M residues, gate 33) and the 68 M-residue database
(gate 41), with v0.3 at the same `-s` for comparison. The presets were chosen
on these queries, so the table is in-sample; held-out numbers for `-s 2`,
`5.7` and `7.5` are in the Benchmarks tables.

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

### Pipeline and what each knob means

1. **Seeds.** For each query window the exact top-`--max-kmers` similar
   k-mers scoring ≥ `--kmer-score` (ties broken deterministically), plus the
   query's own k-mer. The default seed is the spaced pattern `111001101`
   (weight 6 over a span of 9; `--spaced` sets another, `--k N` a contiguous
   one). Low-complexity windows (SEG 12/2.2/2.5) are not seeded (`--mask`),
   and low-complexity target regions are not indexed. Each k-mer is resolved
   by one lookup in the index's direct k-mer table (large indexes) or by
   prefix-batched probes of its sparse table (small ones), prefetched ahead.
2. **Hits** are gathered in bounded target blocks (`--hit-block`) as packed
   64-bit keys (target, diagonal, query position) grouped by an LSD radix
   sort; per target up to `--prefilter-diags` diagonals are kept.
3. **Ungapped evidence.** Each kept diagonal gets the score of its best
   ungapped segment within `--diag-window` (48) residues of its seed hits,
   computed 32 diagonals at a time with AVX2 (exact, rescored when the 8-bit
   lanes saturate); co-linear segments of a target are chained (gap-costed).
   `--ungapped-xdrop x` switches to v0.3's X-drop extension from every seed
   hit instead.
4. **Gate.** Targets whose evidence reaches the ungapped gate are candidates.
   The default gate follows the database size (42 at 10⁸ residues, lower for
   smaller and slightly higher for larger databases, see `--min-diag-score`):
   chance candidates grow with the database while the score a weak alignment
   needs to be significant grows only with its logarithm. `--max-candidates` /
   `--max-seed-occ` are explicit caps (off by default, counted in the stats
   when they bind).
5. **Gapped score.** Exact striped Smith–Waterman–Gotoh: a `u8` pass (32 AVX2
   lanes) and, when its headroom is exceeded, `i16` (16 lanes); scalar only for
   `open ≤ extend` or `i16` saturation.
6. **Traceback** (only for candidates whose score passes `--evalue`): end cell
   from the SIMD pass, start from a SIMD pass over the reversed prefixes, then a
   banded traceback widened until it reproduces the score — exact, with fixed
   tie rules, so the output never depends on kernel, threads, memory limit or
   sharding.
7. **Statistics.** `E = K·(m−ℓ)·(N−n·ℓ)·e^(−λS)` with gapped `λ, K` and the
   length correction `ℓ` fitted on random sequences for the scoring scheme in
   use (BLOSUM62, 11/1, with or without composition-bias correction; table in
   `bench/results/calibration_table.txt`, reproducible with
   `sabertooth calibrate --table --samples 20000 --seed 24301`). `N`, `n` are
   always the whole database. **An E-value is the expected number of chance
   alignments with at least this score in a database of this size — not a
   probability, p-value, FDR or q-value.** Profile searches fit `λ, K` for each
   PSSM by simulation.
8. **False-positive controls** (sequence queries, on by default):
   composition-bias correction (`--comp-bias-corr`, MMseqs2's local
   correction with its own calibrated statistics) and hard masking of
   low-complexity query positions in scoring (`--mask-scoring`, SEG 15/2.3/2.3:
   their score rows become the `X` row). Coverage filters `--min-qcov` /
   `--min-tcov` are applied after traceback.

`--max-hits` (default 300) only limits output per query. Output is BLAST tab
(`query target pident alnlen mismatch gapopen qstart qend tstart tend evalue
bits`, 1-based); `--extra-cols` appends raw score, lengths, coverages and the
CIGAR. Queries are streamed in batches (`--query-batch`) and results are
written in input order; duplicate query identifiers stay separate records.
Errors (unreadable or malformed input, a corrupt index or shard) end the run
with a message and exit status 1; an `--out` file may then be incomplete.

### Memory and threads

`--threads` (default: all cores) parallelises across queries, or across
candidates for batches smaller than the thread count — never both.
`--memory-limit` sizes threads, query batches and per-worker buffers to a
budget. The seed index is memory-mapped: its pages are file-backed and
reclaimable, and peak RSS includes the index pages a search touches.

**Less RAM than the index (out-of-core search).** With `--memory-limit`,
`search` also checks whether the mapped index fits next to the planned
working memory (¾ of the budget minus worker scratch). If not, the index is
cut into residue-balanced shards that each fit (persistent, in
`<index>.shards-<N>/` or `--shard-dir`, built once and reused while the
manifest matches), and every query batch is searched against one shard at a
time; hits are merged per query with whole-database E-values. The cost is
one seed generation per shard and one read of the shards per query batch;
the alternative — mapping an index larger than RAM — pages the index in and
out on every query (v0.3 on a 512 MiB machine: 23× slower than with shards,
`bench/results/scale.json`; MMseqs2 on the 1 GiB machine below: 41× its
unlimited time). `--shards N` sets
the number of shards explicitly (`--shards 1`: never shard). Nothing in the
output depends on any of these settings (tested byte-for-byte; the only
exceptions are the off-by-default caps `--max-candidates` / `--max-seed-occ`,
which then apply per shard).

`--stats` / `--stats-json` report stage CPU times (seed generation, hit
collection, diagonal evidence, gapped scoring, traceback), every prefilter
and alignment counter, the memory plan, peak RSS and page faults.

## Distributed search (`--features distributed`)

```bash
sabertooth search q.fasta protDB --distributed --shards 8     # local cores
sabertooth search q.fasta protDB --distributed --shards 32 --hydra-host   # head
sabertooth search q.fasta protDB --distributed --shards 32 --hydra-client HEAD_IP  # workers
```

The target index is cut into residue-balanced shards with their own
persistent indexes (`--shard-dir`, reused while the manifest matches); query
batches go to every shard with the complete parameter set; workers compute
E-values for the whole database; hits are merged per query in the local
engine's order and re-capped. **Output is byte-identical to a local search**
(tested with 1, 2, 3, 5 and 7 shards). On a cluster the shard directory must be on shared
storage. `--gpus N` reserves GPU slots for scheduling/pinning only — there is
no GPU kernel; all alignment runs on CPUs.

## Benchmarks

Methodology, datasets (with checksums) and commands:
[`bench/README.md`](bench/README.md). Hosts: the tuning, the database-size
series (with its 1 GiB-machine runs) and the SCOPe40 held-out suites ran on
a 2-vCPU Xeon @ 2.8 GHz VM (host A); after a container restart the preset
table, the 512 MiB / 256 MiB / cold-cache table, the index-layout comparison
and the A/B runs (including MMseqs2's memory options on a 1 GiB machine) ran
on a 2-vCPU Xeon @ 2.1 GHz VM (host B). Every timed
row comes from the host its table or section names (`bench/ENVIRONMENT.md`,
`bench/results/v0.4/HOSTS.md`; the untimed null table mixes both); timings
from different hosts are not comparable — the same development binaries and
search ran 1.3–1.4× faster on host B
(`bench/results/v0.4/dev/ab_pf.*` vs `ab/prefetch/`). 8 GiB RAM, no swap,
2 threads unless marked; E ≤ 10 and ≤ 4,000 hits per query for every tool
(full precision–recall range), each tool on its own prebuilt database and
index. Timed runs are medians of 3 measured runs after an unmeasured warm-up
(min–max shown; peak RSS is the largest of the measured runs) unless marked
as single runs (the low-complexity, fragment and two-domain workloads,
exhaustive Smith–Waterman, the preset table and database builds);
memory-limited runs start with the database evicted from the page cache.
v0.3.0, MMseqs2 18-8cc5c and DIAMOND 2.2.8 were re-run on the same host as
the v0.4 numbers they are compared with (the v0.3 release was measured on
another VM; its records stay in `bench/results/*.json` as history).
Machine-readable records: `bench/results/v0.4/` (series, suites per version,
layout and A/B runs, tuning, development records). The tables below are
generated from them by `bench/readme_tables_v04.py` and
`bench/series_tables.py`.

### Scaling: Sabertooth vs MMseqs2 from 2.6 M to 265 M residues

Every database is SCOPe40 (14,323 domains, the only true homologs) plus
order-2 Markov proteins with SCOPe composition and lengths (filler hits are
false positives); 1,000 held-out queries; host A. "Same `-s`" compares settings of the
same name, which the two tools do not calibrate to each other; "matched
sensitivity" interpolates each tool's own settings to the other's AUC1.

<!-- table:series_v04 -->
Queries: `q_eval1k.fa` (sha256 e153440d7e45…), 2 threads unless stated, median of 3 measured runs after 1 warm-up run(s); peak RSS: the largest of the measured runs.

| binary | version | sha256 |
|---|---|---|
| DIAMOND | diamond version 2.2.8 | `d3d72578fa97e87c…` |
| MMseqs2 | 8cc5ce367b5638c4306c2d7cfc652dd099a4643f | `b7ef6e0e33df5dd4…` |
| Sabertooth 0.4 | Sabertooth v0.4.0 | `48c9d9d45bd5a526…` |
| Sabertooth 0.3 | Sabertooth v0.3.0 | `549806e86771203e…` |

#### 2.6 M residues (SCOPe40 only)

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

#### 17 M residues (SCOPe40 + filler)

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

#### 34 M residues (SCOPe40 + filler)

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

#### 68 M residues (SCOPe40 + filler)

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

#### 134 M residues (SCOPe40 + filler)

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

#### 265 M residues (SCOPe40 + filler)

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

#### Matched sensitivity

Wall time each tool needs to reach the other default's AUC1, interpolated (log wall time linear in AUC1) between that tool's measured `-s` settings; “>” = not reached by its most sensitive measured setting, “≤” = already exceeded by its fastest.

| database | MMseqs2 `-s 5.7`: AUC1 / wall s | Sabertooth 0.4 at that AUC1, s | | Sabertooth 0.4 default: AUC1 / wall s | MMseqs2 at that AUC1, s | |
|---|---|---|---|---|---|---|
| 2.6 M | 0.1577 / 2.38 | 1.31 | Sabertooth 1.8× faster | 0.2095 / 3.04 | 7.11 | MMseqs2 2.3× slower |
| 17 M | 0.1561 / 3.20 | 1.49 | Sabertooth 2.2× faster | 0.1957 / 3.28 | 9.79 | MMseqs2 3.0× slower |
| 34 M | 0.1553 / 4.36 | 1.99 | Sabertooth 2.2× faster | 0.1913 / 5.02 | 12.92 | MMseqs2 2.6× slower |
| 68 M | 0.1541 / 6.38 | 3.21 | Sabertooth 2.0× faster | 0.1849 / 7.30 | 15.70 | MMseqs2 2.2× slower |
| 134 M | 0.1529 / 10.80 | 4.31 | Sabertooth 2.5× faster | 0.1822 / 10.79 | 22.89 | MMseqs2 2.1× slower |
| 265 M | 0.1519 / 18.15 | 7.92 | Sabertooth 2.3× faster | 0.1777 / 20.66 | 32.00 | MMseqs2 1.5× slower |

#### Scaling with database size

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

* **Scaling.** Wall time grows as residues^b with b = 0.65 for Sabertooth's
  default and 0.64 for MMseqs2's (v0.3: 0.73, DIAMOND `--sensitive`: 0.96);
  at `-s 7.5` b = 0.78 for Sabertooth and 0.37 for MMseqs2. Sabertooth's
  database-proportional work is the ungapped evidence of every diagonal with
  a chance seed hit (single-hit seeding — two-hit seeding would remove most
  of it but lost 14% AUC1 at 68 M in a development replay and was rejected):
  in the traced tuning run at 265 M residues (tuning queries, default k-mer
  settings, `bench/results/v0.4/tuning/s265_sp9/t22_c512.stats.json`) diagonal
  evidence took 31.6 CPU-s per 1,000 queries — more than
  seed generation, hit collection, gapped alignment and traceback together at
  the default gate.
* **Memory.** Peak RSS is mostly the memory-mapped index: Sabertooth 0.4
  1,419 MiB at 265 M residues (index 1,350 MiB) against MMseqs2 2,653 MiB
  (`-s 5.7`) and 5,259 MiB (`-s 7.5`); DIAMOND is the lightest (932 MiB at
  265 M) but its time grows linearly (216 s at 265 M, 10× Sabertooth's
  default, with a lower AUC1: 0.162 vs 0.178).
* **Index builds** take about as long as v0.3's and less than MMseqs2's
  (99.7 vs 88.9 vs 107.0 s at 265 M); database + index on disk is 58% of
  v0.3's and 52% of MMseqs2's (1,705 vs 2,964 vs 3,257 MiB at 265 M).
* **Threads (265 M).** 1 → 2 threads: Sabertooth 35.8 → 20.7 s (1.73×),
  MMseqs2 38.7 → 18.2 s (2.13×); single-threaded Sabertooth's default is 8%
  faster than MMseqs2's.
* **Sensitivity vs database size.** Every tool's AUC1 falls as filler is
  added (the first false positive comes earlier). Sabertooth's falls more
  (default −15% from 2.6 M to 265 M, MMseqs2 −4%) partly by design: its gate
  rises with the database to keep the time down. `--min-diag-score` sets a
  lower gate when sensitivity matters more than time (tuning queries, 265 M,
  default k-mer settings: gate 34 instead of 43 gives AUC1 +3.7% for +81%
  modelled CPU, `bench/results/v0.4/tuning/s265_sp9`).

### Held-out homology search on SCOPe40 (all tools)

Host A (timings); the null rows are accuracy only.

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

* **v0.3 → v0.4** on SCOPe40 (2.6 M residues, where the gate is the same as
  v0.3's): the same time (12.0 vs 12.1 s), AUC1 +6.4% (0.215 vs 0.202),
  1,381 more true positives at E ≤ 1e-3 with the same 19 false positives,
  and 4.7× fewer false positives at matched recall 0.10 (89 vs 419); `-s 7.5`
  gains 5.2% AUC1 for 11% more time.
* **vs MMseqs2**: the default is slower than MMseqs2 `-s 5.7` on this small
  database (12.0 vs 8.4 s) but far more sensitive (0.215 vs 0.160) with 15×
  less memory; `-s 7.5` is 1.3× faster than MMseqs2 `-s 7.5` (32.2 vs
  41.8 s) and more sensitive (0.243 vs 0.234) with 27× less memory.
* **vs DIAMOND**: `--sensitive` is faster (7.1 s) and lighter (38 MiB) than
  the default but less sensitive (0.166); `--ultra-sensitive` is the most
  sensitive prefiltered search here (0.246) and faster than `-s 7.5` (see
  limitations).
* **Exhaustive Smith–Waterman** with Sabertooth's scoring reaches AUC1 0.264:
  the prefilter costs 19% of that at the default and 8% at `-s 7.5`.
* The first false positives are the same few cross-fold SCOP pairs in every
  run (disputable labels), so differences of a few FPs at E ≤ 1e-3 carry
  little signal; FP at matched recall 0.08–0.10 separates the settings.

### False-positive controls

Host A.

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

Every LC-insert query has a 30–60-residue homopolymer / tandem-repeat /
biased stretch; the database adds 2,000 low-complexity and 2,000
compositionally biased decoys. The controls (seed masking, composition-bias
correction, hard masking of low-complexity query positions in scoring) are
unchanged from v0.3; their per-control effect is in v0.3's ablations
(`bench/results/tune.json`, `CHECKLIST.md` 6.12).

### Fragments and two-domain queries

Host A.

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

### E-value calibration on null data

Accuracy only (no timings): host A, except v0.3's `null_db` rows (host B).

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

`null_q`: 1,000 composition-preserving shuffles of real queries vs SCOPe40
(expected chance hits with E ≤ x: 1,000·x); `null_db`: the 4,000 held-out
queries vs a shuffled SCOPe40 (expected 4,000·x). Exhaustive Smith–Waterman
measures the statistics themselves (observed/expected ≈ 0.9: calibrated,
slightly conservative); the prefiltered search is more conservative because
the prefilter drops most chance alignments. The statistics are unchanged
from v0.3.

### Limited RAM and cold cache (68 M residues)

Host B, tuning queries (as in v0.3); the database files are evicted before
each of 3 measured runs, after a one-query run that builds the out-of-core
shards:

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

The 1 GiB machine, from the database-size series (host A, held-out queries,
database evicted before each of 3 measured runs; the unmeasured warm-up runs
outside the limit and builds v0.3's shards):

<!-- table:mem1g_v04 -->
| run (68 M residues, 1 GiB RAM) | memory option | shards | wall s | peak RSS MiB | cgroup peak MiB | AUC1 |
|---|---|---|---|---|---|---|
| v0.4 default | `--memory-limit 1G` | 1 | 7.3 (6.5–7.9) | 475 | 472 | 0.1849 |
| v0.3 default | `--memory-limit 1G` | 2 | 29.7 (29.6–29.8) | 486 | 919 | 0.1837 |
| MMseqs2 `-s 5.7` | none (see text) | — | 259.8 (254.1–265.3) | 1022 | 1024 | 0.1541 |
| DIAMOND `--sensitive` | none (see text) | — | 55.7 (55.5–55.7) | 268 | 336 | 0.1635 |
<!-- /table -->

* **1 GiB machine**: v0.4 7.3 s without shards — its 440 MiB index fits —
  against v0.3's 29.7 s (2 out-of-core shards), MMseqs2's 260 s (its 1.4 GiB
  index is paged in and out: 41× its unlimited time) and DIAMOND's 55.7 s
  (268 MiB; it needs no option). MMseqs2's memory options do not help
  (host B, held-out queries, `bench/results/v0.4/ab/mmseqs_split/` and
  `ab/mem1g_hostb/`): `--split-memory-limit` 768M or 1G makes it refuse to
  start ("Cannot fit databases into …", with or without its precomputed
  index; likewise in v0.3's runs), `--split 0` or `--split 4` without the
  index are killed at the limit, and `--split 4` with the index takes as long
  as no option (166 vs 164 s, median of 3, same hits) — against 4.8 s for
  Sabertooth 0.4 on that host with the same protocol.
* **512 MiB / 256 MiB machines**: v0.4 7.0 s (3 shards) / 19.6 s (7 shards)
  against v0.3's 28.9 s (5) / 41.2 s (19), with the same hits as the
  unlimited run (same AUC1, TP and FP; sharded output is byte-identical by
  construction and test, `out_of_core_shards_give_byte_identical_output`).
  v0.4 uses more of the budget: peak RSS 261 / 136 MiB against v0.3's
  231 / 90 MiB, within the same limits.
* **Cold page cache**: 4.5 s vs v0.3's 16.6 s.

## Trade-offs and limitations

Measured, not hidden:

* **MMseqs2 leads at the top of the sensitivity range on large databases.**
  At `-s 7.5` MMseqs2's time grows with the database as residues^0.37 against
  ^0.78 for Sabertooth: at 265 M residues it is 1.5× faster (61.9 vs 90.9 s)
  with a 3.7% higher AUC1 (0.208 vs 0.200); at 68 M Sabertooth `-s 7.5` is
  1.5× faster but 2.2% less sensitive. MMseqs2 `-s 7.5` also finds more true
  positives at E ≤ 1e-3 at every size (5,977 vs 5,342 at 265 M) and admits
  fewer false positives by recall 0.10 (3–5 vs 5–32). Sabertooth's
  high-sensitivity cost at scale is the ungapped evidence of every chance seed
  diagonal (work proportional to the database; two-hit seeding, which would
  avoid most of it, lost 14% AUC1 at 68 M in a development replay — see the
  last items). Its most sensitive preset, `-s 9`, reaches AUC1 0.251 on the
  SCOPe40 tuning queries, where MMseqs2 `-s 8.5` reached 0.264 (tuning
  records).
* **More false positives at high recall on the largest database.** By the
  time the default has found 10% of all true pairs it has admitted fewer
  false positives than v0.3 up to 68 M (7 / 50 / 115 / 367 vs 22 / 117 / 235 /
  450 at 2.6 / 17 / 34 / 68 M), about as many at 134 M (881 vs 870) and 28%
  more at 265 M (2,243 vs 1,750): the gate that rises with the database drops
  weak true hits, so that recall needs a looser E-value. At recall 0.06–0.08
  it is equal or better at every size; MMseqs2's default does not reach
  recall 0.10 at any size.
* **DIAMOND `--ultra-sensitive`** is faster and more sensitive than Sabertooth
  `-s 7.5` on held-out SCOPe40 (20.7 vs 32.2 s, AUC1 0.246 vs 0.243) and on
  two-domain queries (11.6 vs 15.7 s, 0.176 vs 0.170; single runs);
  Sabertooth `-s 7.5` leads on fragments (3.4 vs 5.8 s, 0.147 vs 0.132) and
  low-complexity-insert queries (9.8 vs 13.1 s, 0.196 vs 0.194; single
  runs). DIAMOND's time grows linearly with the database (216 s for
  `--sensitive` at 265 M residues).
* **The fastest settings.** MMseqs2 `-s 1` is faster than Sabertooth `-s 1`
  at every size (0.26–0.55 vs 0.44–0.82 s for 1,000 queries) but finds a third
  as much (AUC1 0.021 vs 0.063–0.067), and MMseqs2 `-s 4` is faster than
  Sabertooth `-s 4` up to 68 M (Sabertooth at 0.64–0.88× its speed, level
  from 134 M) while Sabertooth's AUC1 is 40–56% higher. At the same `-s 5.7`
  MMseqs2's default is up to 28% faster than Sabertooth's (14–15% at 34, 68 and 265 M; within 3% at
  17 and 134 M) — Sabertooth's default is set more sensitive.
* **The `-s` levels moved.** v0.4's `-s 1…4` are much faster and less
  sensitive than v0.3's same levels (v0.3 `-s 2` on held-out SCOPe40: AUC1
  0.154 in 5.9 s; v0.4 `-s 2`: 0.081 in 2.1 s); scripts that relied on v0.3's
  `-s 2` sensitivity want v0.4 `-s 4`–`5` (preset table). At the top, v0.4
  `-s 6.5 / 7.5 / 9` are 5–49% slower than v0.3's same levels on SCOPe40 (and
  4.6–6.5% more sensitive; tuning queries, host B).
* **Small databases gain sensitivity, not speed.** On SCOPe40 (2.6 M
  residues) the v0.4 default is as fast as v0.3's (12.0 vs 12.1 s for the
  4,000 held-out queries; MMseqs2's default is faster there, 8.4 s); the
  speed-ups grow with the database because v0.4 removes database-proportional
  work.
* **Sensitivity falls with database size.** By design the gate rises with
  the database to bound the time; at 265 M residues the default's AUC1 is
  only 0.1% above v0.3's. `--min-diag-score` lowers it when sensitivity
  matters more than time (see Scaling).
* **More memory than a sparse index on mid-size databases.** The direct k-mer
  table (132 MB for weight 6) is chosen once the sparse tables would exceed
  64 MiB, i.e. from roughly 10 M residues: at 17 M residues the default index
  is 205 MiB and the search peaks at 223 MiB, against 172 / 198 MiB with
  `createindex --index-layout sparse` — which is 1.9× slower (4.72 vs 2.46 s;
  tuning queries, host B, identical hits, `bench/results/v0.4/layout/`). At
  68 M the direct index is both smaller and faster (440 vs 582 MiB, 4.9 vs
  7.4 s). Under a memory limit v0.4 also uses more of the budget than v0.3
  (Limited RAM above).
* **Index builds** are not faster (single builds: 23.1 vs 23.3 s for v0.3 at
  68 M residues, 99.7 vs 88.9 s at 265 M) and peak at ~1.4× the partition
  buffer (1,430 MiB with the default 1 GiB buffer at 68 M); `--memory-limit`
  bounds them as in v0.3.
* **Synthetic large databases.** The series adds order-2 Markov filler to
  SCOPe40: the number of true homologs stays fixed while the database grows,
  and the filler has no internal homology. Real large databases (UniRef) add
  homologs and low-complexity families too; the gate and presets were tuned
  on this construction and are unverified there. Largest database run:
  265 M residues.
* **Two small machines, 2 threads.** Thread scaling beyond 2 is unverified
  (1 → 2 threads: 1.73× for Sabertooth vs 2.13× for MMseqs2 at 265 M); the
  benchmark host changed during the runs, so each table or section states its
  host.
  The AVX-512 kernels are not faster than AVX2 — within noise on host A
  (18.36 vs 18.52 s, a development run in `experiments.jsonl`), 6% slower on
  host B (4.90 vs 4.62 s, `bench/results/v0.4/ab/simd/`) — so they are opt-in;
  non-x86 CPUs get the exact scalar kernels (no NEON).
* **Out-of-core costs time.** Each shard repeats seed generation and each
  query batch re-reads the shards: on a 256 MiB machine the 68 M search takes
  19.6 s (7 shards) against 4.5 s with enough RAM (host B, cold cache). The
  planner's headroom rule is unchanged from v0.3 (not re-tuned).
* **Evidence from development runs.** Three decisions rest on development
  measurements whose records were not kept: rejecting two-hit seeding (an
  exact replay of a 68 M candidate trace gave AUC1 0.190 → 0.163–0.170;
  neither the replay script nor its output was kept, and it needs a
  per-diagonal trace that this release does not write), and not adopting posting
  prefetch at seed creation or evidence lanes read straight from the targets
  (both slower). The hit-loop prefetch, the evidence window, the index layouts
  and the AVX-512 decision were re-measured with shipped drivers
  (`bench/results/v0.4/ab/`, `layout/`); the surviving development records
  are in `bench/results/v0.4/dev/`.
* **Known issues.** `sabertooth help` lists `SABERTOOTH_SIMD=avx2|sse2|scalar`
  without `avx512` (accepted, see Build), and the error for a negative
  `--ungapped-xdrop` still says "0 = whole-diagonal scoring" (0 selects the
  default windowed evidence; `--diag-window 0` is the whole diagonal). Both
  are left as built so that the shipped source rebuilds the benchmarked
  binary bit for bit.
* **Distributed mode** is tested on one host only (byte-identical output); no
  multi-node timing.
* **Seed-guided windows** (`--window 1`) and **profile search** are not part
  of the benchmarks.
* **MSRV.** `rust-version = "1.89"` is declared; only rustc 1.97 was used.
* **No GPU computation.** `--gpus` only reserves and pins devices.

## Layout

```
src/
  seqio.rs fasta.rs mmdb.rs   input formats, validation, MMseqs2 DB reader/writer
  index.rs                    persistent seed index (build, open, lookup)
  prefilter.rs mask.rs        seeds, top-N k-mers, hit blocks, gate, SEG masking
  ungapped.rs                 batched SIMD ungapped segment evidence (AVX2/SSSE3/AVX-512/scalar)
  simd.rs align.rs            striped SW kernels (SSE2/AVX2/AVX-512, u8/i16), canonical traceback
  compbias.rs evalue.rs       composition-bias correction, calibrated statistics
  search.rs                   staged search engine (sequence + profile)
  distributed.rs              HydraMPP backend (feature `distributed`)
  main.rs metrics.rs          CLI, memory plan, stats JSON
  profile.rs msa.rs matrix.rs alphabet.rs nucl.rs cluster.rs taxonomy.rs
tests/cli_db_equivalence.rs   CLI integration and parity tests
bench/                        dataset generators, benchmark drivers (run_suite.py, run_series.py, run_layouts.py, run_ab.py),
                              evaluation, tables, results (bench/results/v0.4/ for this release)
bench/tuning/                 gate replay, cost model and frontier scripts behind the v0.4 defaults
```

## License / provenance

Reimplements algorithms of MMseqs2 (Steinegger & Söding, *Nat. Biotechnol.*
2017) from their published description and open-source formulas; the bundled
`blosum62.out` comes from the MMseqs2 distribution (GPL-3.0), so treat this
code as GPL-3.0. (`Cargo.toml` declares `license = "MIT"`; the mismatch is
pre-existing and left for the owner to resolve.) Earlier design notes: `COMPARISON.md`, `HARDENING.md`,
`STEPS.md`, `VALIDATION.md` (historical, v0.1–v0.2).
