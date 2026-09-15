# Sabertooth

```
   ___       _           _              _   _
  / __| __ _| |__  ___ _| |_ ___  ___ _| |_| |_
  \__ \/ _` | '_ \/ -_)  _/ _ \/ _ \  _|  _| ' \
  |___/\__,_|_.__/\___|\__\___/\___/\__|\__|_||_|
       \V/   \V/    fast · light · profile-sensitive
```

A pure-Rust reimplementation of the **profile / PSSM sensitivity core** of
[MMseqs2](https://github.com/soedinglab/MMseqs2), built to serve as the search
backend for a Rust rewrite of profile-scanning tools such as geNomad.

> **Scope, stated honestly.** MMseqs2 is ~75,000 lines of C++ in `src/` alone,
> spanning clustering, taxonomy, iterative search workflows, MSA generation,
> vectorized kernels, and much more. Sabertooth is **~2,700 lines of Rust** that
> faithfully reimplements the specific path that matters for remote-homology
> profile search — and nothing it doesn't. It is **not** a drop-in replacement
> for all of MMseqs2. What it *does* implement, it implements to match: the
> substitution-matrix reconstruction and PSSM construction are **cell-for-cell
> identical** to an independent implementation of MMseqs2's documented formulas
> (see [Validation](#validation)).

---

## What's implemented

The complete sensitivity-critical profile-search pipeline:

| Stage | Module | Notes |
|-------|--------|-------|
| Amino-acid alphabet | `alphabet.rs` | MMseqs2 residue order `ACDEFGHIKLMNPQRSTVWY`, `X` sentinel, ambiguity folding |
| Substitution matrix | `matrix.rs` | Reconstructs joint-prob `P`, pseudocount matrix `R`, and integer scores from `blosum62.out` using MMseqs2's exact `exp(λ·S)·p_a·p_b` derivation |
| Profile / PSSM | `profile.rs` | Henikoff position-based sequence weights, per-column `Neff_M`, substitution pseudocounts, log-odds PSSM (bit-factor 8) — MMseqs2 `PSSMCalculator` defaults `pca=1.0, pcb=1.5` |
| MSA input | `profile.rs` | a3m (lowercase inserts dropped) and aligned-FASTA parsing |
| K-mer prefilter | `prefilter.rs` | Exact + similar-k-mer generation via best-first branch-and-bound, plus an **ungapped diagonal (Kadane) gate** mirroring MMseqs2's ungapped prefilter |
| Local alignment | `align.rs` | Smith–Waterman–Gotoh affine-gap DP with full traceback, generic over sequence and profile scorers |
| Statistics | `evalue.rs` | Karlin–Altschul bit scores and E-values; analytic ungapped λ for the matrix, principled `λ = ln2/8` for log-odds PSSM scores |
| Search driver | `search.rs` | Rayon-parallel search, E-value + coverage filtering, BLAST tab (`.m8`) output |
| CLI | `main.rs` | `createdb`, `msa2profile`, `search`, `profilesearch`, `align`, plus `doctor`/`version`/`info`/`help` |

### What's intentionally **not** implemented

So there's no ambiguity about what you're getting:

- **Clustering / linclust**, **taxonomy assignment**, and the **cascaded/iterative
  search** workflow (`mmseqs search`'s multi-round profile bootstrapping).
- **MSA generation** (`result2msa`) — Sabertooth *consumes* MSAs, it doesn't build them.
- **SIMD-striped Smith–Waterman.** Alignment uses a clean scalar DP. It is correct
  and parallelised across target sequences with rayon, but a single alignment is not
  vectorized the way MMseqs2's kernels are.
- **Profile–profile** and **nucleotide / translated** search.
- **The MMseqs2 on-disk database format.** Sabertooth reads and writes FASTA.
- **ALP-calibrated gapped statistics.** MMseqs2 links the ALP library to fit gapped
  Gumbel parameters. Sabertooth uses an analytically solved ungapped λ for the matrix
  and the exact `ln2/8` λ for PSSM log-odds (both principled), which is sufficient for
  ranking and thresholding but is not the same calibration machinery.

These are noted again, with rationale, in [`COMPARISON.md`](COMPARISON.md).

---

## Building

Requires a Rust toolchain. Developed and tested against **rustc 1.75.0** (the crate
pins `rust-version = "1.75"`).

```bash
cargo build --release
```

The only direct dependency is [`rayon`](https://crates.io/crates/rayon) for
data-parallelism. On rustc 1.75 you may need to hold `rayon-core` at a
1.75-compatible release (newer `rayon-core` raises the MSRV):

```bash
cargo update -p rayon-core --precise 1.12.1
```

Run the test suite (20 unit + integration tests):

```bash
cargo test --release
```

Confirm the build is healthy with the self-check:

```bash
./target/release/sabertooth doctor
```

---

## Usage

### Build a profile (PSSM) from an MSA

```bash
sabertooth msa2profile family.a3m --out family.pssm
```

Emits a per-column PSSM table (20 amino-acid log-odds scores + `Neff` per position).

### Profile search

```bash
sabertooth profilesearch family.a3m targets.fasta --out hits.m8 --evalue 1e-3
```

### Sequence search

```bash
sabertooth search queries.fasta targets.fasta --out hits.m8 --threads 8
```

### Pairwise alignment

```bash
sabertooth align query.fasta targets.fasta
```

Output for the search commands is BLAST tabular (`.m8`):
`query  target  pident  alnlen  mismatch  gapopen  qstart  qend  tstart  tend  evalue  bits`.

### Key options

| Option | Meaning | Default |
|--------|---------|---------|
| `--evalue <f>` | Max E-value reported | `1e-3` |
| `--k <int>` | K-mer length | 6 (seq), 5 (profile) |
| `--kmer-score <int>` | Similar-k-mer score threshold | 25 (seq), 40 (profile) |
| `--min-diag-score <int>` | Ungapped diagonal gate | 15 (seq), 30 (profile) |
| `--min-cov <f>` | Min query coverage 0..1 | 0 |
| `--max-hits <int>` | Max hits per query (0 = all) | 300 |
| `--pca` / `--pcb` / `--bias` | PSSM pseudocount admixture / score bias | 1.0 / 1.5 / 0.0 |
| `--threads <n>` | Worker threads | all cores |

---

## Validation

Because a full MMseqs2 binary could not be built in the development sandbox
(no CMake, and the `simde` / `gzstream` git submodules were absent), fidelity was
established by **differential testing against an independent second
implementation** of MMseqs2's documented formulas, written from scratch in Python
(`compare/mmseqs_reference.py`) as a wholly separate code path.

Running `compare/compare.py` diffs Sabertooth's actual CLI output against that
reference, cell by cell:

| Quantity | Cells compared | Exact match | Max &#124;diff&#124; |
|----------|---------------:|------------:|---------------------:|
| Integer substitution matrix | 400 | **400 / 400 (100%)** | 0 |
| PSSM — ungapped MSA (`family.a3m`) | 2620 | **2620 / 2620 (100%)** | 0 |
| PSSM — gapped MSA (`gapped.fasta`) | 1000 | **1000 / 1000 (100%)** | 0 |

Two independent implementations agreeing to the integer on every cell — across the
matrix reconstruction, sequence weighting, `Neff`, pseudocounts, and log-odds
rounding — is strong evidence the port is faithful. (This validates fidelity *to the
documented algorithm*. A byte-diff against a compiled `mmseqs` binary is the natural
next check and is not claimed here.)

Reproduce:

```bash
python3 compare/compare.py target/release/sabertooth demo/family.a3m
python3 compare/compare.py target/release/sabertooth demo/gapped.fasta
```

### Worked example: remote-homology sensitivity

The `demo/` directory contains a 4-sequence protein-kinase MSA and a target set
mixing diverged kinase homologs with unrelated decoys. Profile search recovers the
remote homolog that a naive threshold would miss:

```
$ sabertooth profilesearch demo/family.a3m demo/targets.fasta --evalue 1e-3
  true_homolog_diverged    E=3.31e-33   bits=124.6   id=26.5%
  another_kinase_homolog   E=9.63e-11   bits=50.0    id=16.2%
  (decoys: E = 0.29 – 0.90, correctly excluded)
```

A true homolog is detected at **26.5% sequence identity** with `E = 3×10⁻³³`, while
unrelated decoys land at `E ≈ 0.3–0.9` and are filtered — the profile sensitivity
that makes PSSM search worthwhile.

---

## Architecture notes

- **Correct profile statistics.** A PSSM value is `8·log2(p/p_back)`. For pure
  log-odds scores the Karlin–Altschul equation `Σ p_back·e^{λM} = 1` is solved
  exactly by `λ = ln2 / bit_factor = ln2/8`. Using the substitution matrix's own λ
  (from a ×2-bit scale) on ×8-scale PSSM scores would massively inflate significance;
  Sabertooth uses the correct per-scale λ, which is what keeps decoys non-significant.
- **Prefilter that doesn't drop homologs.** Similar-k-mer generation is a best-first
  branch-and-bound so that capping the number of k-mers per position keeps the
  *highest-scoring* ones (including the exact seed), and a cheap ungapped Kadane pass
  along the seed diagonal gates candidates before the expensive gapped alignment —
  the same division of labour MMseqs2 uses.
- **Light footprint.** One direct dependency, a ~7.6 MB static binary, no build-time
  code generation.

---

## Layout

```
src/
  alphabet.rs   amino-acid encoding
  matrix.rs     substitution-matrix reconstruction
  profile.rs    PSSM construction from MSAs
  prefilter.rs  k-mer prefilter + ungapped diagonal gate
  align.rs      Smith–Waterman–Gotoh alignment
  evalue.rs     Karlin–Altschul statistics
  search.rs     parallel search driver
  banner.rs     ASCII-art banners
  main.rs       CLI
  data/blosum62.out   MMseqs2's exact BLOSUM62 data file
compare/        independent Python reference + diff harness
demo/           example MSA, targets, and gapped test data
```

## License / provenance

Reimplements the algorithms of MMseqs2 (Steinegger & Söding, *Nat. Biotechnol.*
2017) from their published description and open-source formulas. The bundled
`blosum62.out` originates from the MMseqs2 distribution (GPL-3.0); treat this
reimplementation as GPL-3.0 accordingly.
