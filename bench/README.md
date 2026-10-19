# Benchmarks

Everything here is reproducible from the scripts in this directory. Results
used in the top-level README live in `bench/results/` (machine-readable JSON
plus Markdown summaries); raw hit files are not committed (they are large and
regenerable).

## Environment

See `ENVIRONMENT.md`: 2 vCPU Intel Xeon @ 2.1 GHz KVM guest (AVX2 + AVX-512
present), 8 GiB RAM, no swap, rustc 1.97.0, release profile (thin LTO,
`codegen-units = 1`, no `target-cpu=native`). Tools: Sabertooth v0.2.0
baseline (untouched source, same compiler), Sabertooth v0.3.0, MMseqs2
18-8cc5c (`mmseqs version` = 8cc5ce367b5638c4306c2d7cfc652dd099a4643f),
DIAMOND 2.2.8. **Every timing is a 2-thread run on this one machine**; numbers
on other hardware will differ, ratios less so.

## Data

`make_datasets.py` builds every dataset deterministically (seed 20261018) from
an ASTRAL SCOPe 40%-identity domain set: `db/astral.fa` of
github.com/ivujevic/tachyon at commit `fd5ed981e01064b5da01f4bbddb5ad0a3019eb1e`
(saved as `astral_tachyon.fa`; 14,323 domains, 2,611,087 residues, SHA-256
`beafd78dfd81fc251f0bffc38e905febb5b637a490d9459a66499f5a5e2d71bf`; headers
carry the SCOP sccs labels; the exact SCOPe release is not stated in that
repository). `ds/MANIFEST.tsv` lists the SHA-256, record and residue counts
of each generated file. Running it twice gives byte-identical files (checked:
the scaling run regenerated all files and the manifests agree).

| file | what | used for |
|---|---|---|
| `scope40.fa` | the curated database (truth = SCOP sccs labels) | homology workloads |
| `split.tsv` | superfamily → tune / eval (SHA-256 parity of the sf id) | tuning vs held-out separation |
| `q_tune.fa`, `q_eval.fa` | 1,000 / 4,000 queries from tune / eval superfamilies (each with ≥1 true homolog) | parameter choice / final evaluation |
| `q_frag*.fa` | 30–60% fragments of queries | short queries |
| `q_multi*.fa` | two-domain fusions of different folds | multidomain queries |
| `q_lc*.fa` | queries with a 30–60 aa low-complexity insert (homopolymers, short tandem repeats, biased stretches) | low-complexity FP control |
| `q_null*.fa` | composition-preserving shuffles of queries (no homologs) | E-value calibration |
| `db_decoy.fa` | scope40 + windowed shuffles (local composition kept) + 2,000 LC decoys + 2,000 compositionally biased decoys (TM-, acidic-, basic-, Pro-, Ser/Thr-, Cys-rich) | FP controls |
| `db_null.fa` | composition-preserving shuffle of every scope40 domain | E-value calibration |
| `ds_scale/scale_db.fa` | scope40 + 360,000 order-2 Markov proteins of SCOPe composition and length distribution (374,323 sequences, 68,310,060 residues; SHA-256 in `ds_scale/scale_db.sha256`) | large-database scaling |

Truth (MMseqs2 / foldseek SCOP benchmark convention, `eval_scop.py`): a hit
is a **true positive** if the target shares a superfamily with the query
(self-hits excluded), a **false positive** if the target is in a different
fold or is any synthetic decoy/filler, and is **ignored** if same fold but a
different superfamily. Windowed shuffles of the query's own superfamily are
ignored (they keep residual similarity). Composite queries are TPs for any
source domain's superfamily and FPs only if the fold differs from all
sources.

**Tuning vs evaluation.** Every default and preset (k-mer thresholds, caps,
gates, scoring mask, X-drop) was chosen on the `*_tune` files only, and
nothing was changed after seeing `eval` results. The `eval` suite uses the
held-out superfamilies. Its Sabertooth runs were executed twice with
byte-identical outputs: in the first full pass, and again with the delivered
build (which adds the out-of-core mode and lint fixes, outside the default
search path); the MMseqs2, DIAMOND and v0.2 records are from the first pass
(same machine, unchanged binaries). The same holds for the `tune` and `null`
suites, with three additions made after the first pass: the `kadane` and
`diag1` ablations (run once, after the held-out suite), and the `plain`
ablation, re-run with all three false-positive controls off (its first
version left the scoring mask on). The `scale` and `build` suites were run
entirely with the second-pass build. The delivered source differs from that
build only outside the default search path (see `ENVIRONMENT.md`); the
delivered binary reproduced eleven benchmark outputs byte for byte (eval,
lc, frag, multi, tune, lc_tune `nomask`, tune `plain`, null_db, scale, scale
out-of-core, scale `-s 7.5`). All binaries are listed in `ENVIRONMENT.md`.

## Metrics (all with explicit denominators)

Per run `eval_scop.py` reports: the TP/FP pair denominators; at E-value
thresholds 1e-10…10 the TP and FP counts, recall (TP / all TP pairs), recall
split into same-family and superfamily-only ("distant") pairs, precision,
**FDP** = FP / (TP + FP) and **FPR** = FP / all FP pairs; **FP at matched
recall** (the FP count at the first E-value threshold reaching each recall
level — the main "false positives at matched sensitivity" number); a
precision–recall curve (about 200 points); **AUC1** (per query, the fraction of
its TPs ranked above its first FP, averaged — threshold-free and independent
of E-value calibration); and AUC1 by query length. Null workloads report the
observed number of hits with E ≤ x per query against the expected x
(E-value calibration: observed/expected ≤ 1 is conservative).

`measure.py` records wall time, peak RSS (`ru_maxrss` of the child tree),
user/sys CPU and page faults per repetition (min/median/max over `--repeat`);
`--warmup N` first runs the command N times unmeasured (warm page cache,
independent of what ran before); `--evict` drops files from the page cache
first (cold start); `--cgroup-limit` runs inside a fresh cgroup-v1 memory
cgroup whose limit counts page cache, i.e. a machine with that much RAM (peak
usage and failcnt recorded; the database files are evicted first so that
their pages are charged to that cgroup). Sabertooth runs also save
`--stats-json` (stage CPU times, every prefilter/alignment counter, memory
plan, shards). Runs marked `prepare` (the out-of-core runs) first execute the
same command with one query under the same limits, which builds the
persistent shards; that step is recorded separately (`prepare` in the run
record), and the measured repetitions reuse the shards.

All tools report E ≤ 10 with up to 4,000 hits per query so the curves cover
the same range, and each searches its own prebuilt database (build time,
peak RSS and size recorded separately). For MMseqs2 this means
`--max-seqs 4000` (its default of 300 also limits how many prefilter
candidates are aligned, so it would be faster and less sensitive); for
DIAMOND `-k 4000` (default 25). Sabertooth's `--max-hits` limits output only. The v0.2 baseline reads FASTA: its
MMseqs2-database path silently searched zero targets (see `CHANGES.md`).

## Running

```bash
# datasets (deterministic)
python3 -I bench/make_datasets.py --scope astral_tachyon.fa --out DS
python3 -I bench/make_datasets.py --scope astral_tachyon.fa --out DS_SCALE --markov 360000
cat DS_SCALE/scope40.fa DS_SCALE/markov_360000.fa > DS_SCALE/scale_db.fa
grep '^>' DS_SCALE/markov_360000.fa | sed 's/^>//; s/ .*//' > DS_SCALE/scale.extra.ids

# suites (each writes OUT/<suite>/{summary.md,summary.json,<run>.json}); tools optional
T="--sabertooth target/release/sabertooth --baseline OLD/target/release/sabertooth \
   --mmseqs mmseqs --diamond diamond --threads 2"
python3 bench/run_suite.py --data DS --out OUT --suite tune  $T
python3 bench/run_suite.py --data DS --out OUT --suite eval  $T --repeat 3
python3 bench/run_suite.py --data DS --out OUT --suite null  $T
python3 bench/run_suite.py --data DS_SCALE --out OUT --suite scale $T --repeat 3 --warmup 1
python3 bench/run_suite.py --data DS_SCALE --out OUT --suite build $T   # createindex under limits
python3 bench/plot_results.py OUT bench/results          # figures + compact JSON + tables.md
python3 bench/readme_tables.py OUT > bench/results/readme_tables.md
python3 bench/readme_tables.py OUT --fill README.md      # refresh the README tables
# the delivered tree also holds every run record: OUT/<suite>/<run>__<workload>.json
# (not *.stats.json, which is embedded in them) copied to bench/results/runs/<suite>/,
# OUT/build/*.json to bench/results/runs/build/ and OUT/build/summary.json to
# bench/results/build.json
```

Suites: `tune` (sabertooth `-s` presets and ablations on the tuning split),
`eval` (held-out split, every tool at two operating points, exhaustive
Smith–Waterman as the sensitivity ceiling of the scoring scheme), `null`
(E-value calibration), `scale` (68 M-residue database: presets, 1 vs 2
threads, cold page cache, frequent-seed caps, 1 GiB / 512 MiB / 256 MiB
machines via cgroups — Sabertooth with `--memory-limit` = ¾ of the machine,
which selects the out-of-core shards, plus the 512 MiB machine with
`--shards 1`; MMseqs2 with and without `--split-memory-limit` — and every
tool at its default operating point), `build` (`createindex` on the 68 M-residue
database by default and under 1 GiB / 512 MiB / 256 MiB limits).
`--only NAME…` re-runs selected runs and keeps the other records of the suite
in its summary.

Peak RSS (`ru_maxrss`) counts resident file-backed pages as well, some of
which may be charged to another cgroup (shared libraries, pages cached
earlier), so it can slightly exceed the cgroup peak; the cgroup peak is the
binding number for the memory-limited runs.
