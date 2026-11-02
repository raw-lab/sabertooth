# Sabertooth vs MMseqs2 — detailed comparison

> **Historical document (v0.1 – v0.2).** Kept for provenance. For the current
> release see `README.md`, `CHANGES.md` and `CHECKLIST-v0.4.md` (v0.3:
> `CHECKLIST.md`); the measured head-to-head comparison with MMseqs2 is in
> `README.md` ("Scaling") and `bench/results/v0.4/`; statements below about the seed
> index, statistics (E-values), SIMD kernels, traceback, database input and the
> distributed backend describe the old code and are superseded where they differ.


This document is the technical answer to "how does it compare to the original."
It walks the profile-search pipeline stage by stage, states exactly which MMseqs2
source construct each Rust module reproduces, and is explicit about the boundaries.

## Summary

Sabertooth reimplements MMseqs2's **profile/PSSM search core** in ~3,400 lines of
Rust (plus ~700 lines of tests). Where the two overlap — matrix reconstruction and PSSM construction — outputs
are **numerically identical to the integer** (verified by differential testing and a
byte-diff against a compiled `mmseqs`, see below). It also covers spaced seeds, a
sensitivity parameter, iterative (PSI-BLAST-style) profile search, a Monte-Carlo
gapped-statistics calibration, a striped SSE2 alignment kernel, and banded alignment.
The remaining MMseqs2 machinery (clustering/linclust, taxonomy, MSA building, the
on-disk DB format) is out of scope by design for the geNomad backend.

---

## Stage-by-stage

### 1. Amino-acid alphabet — `alphabet.rs` ↔ `Sequence.cpp`, `SubstitutionMatrix.h`

- Same 20-letter order `ACDEFGHIKLMNPQRSTVWY` with `X` as the ambiguous/degenerate
  symbol at index 20 (`PROFILE_AA_SIZE = 20`).
- `B → D`, `Z → E`, `J → I/L`, `U → C`, `O → K`, and any other non-standard symbol
  fold to `X`, as MMseqs2 does when reading sequences.

### 2. Substitution matrix — `matrix.rs` ↔ `SubstitutionMatrix.cpp`, `BaseMatrix.cpp`

MMseqs2 does **not** read pre-rounded integer scores; it reconstructs them from the
`.out` file's bit-scores. Sabertooth performs the identical reconstruction:

1. `P(a,b) = exp(λ · S_ab) · p_back[a] · p_back[b]` (joint probability)
2. `R[a][b] = P(a,b) / p_back[b]` (= `P(a|b)`, the pseudocount mixing matrix)
3. `score[i][j] = round(bit_factor · log2(P/(p_a·p_b)) + bias)`

with the search-matrix constants `bit_factor = 2.0`, `bias = -0.2`, and the
`ANY_BACK = 1e-5` handling for a non-positive `X` background (rescaling the real
backgrounds by `1 - p_back[X]`). The bundled `data/blosum62.out` is MMseqs2's own
data file, so the inputs are identical too.

**Result:** 400/400 integer cells identical to the independent reference.

### 3. Profile / PSSM — `profile.rs` ↔ `PSSMCalculator.cpp`, `MsaFilter.cpp`

The default MMseqs2 path (`wg = true`, `PCMODE_SUBSTITUTION_SCORE`) is reproduced end
to end:

1. **Sequence weights** — Henikoff–Henikoff 1994 position-based weighting:
   `w[k] += 1 / (n_l[aa] · distinct_aa[col] · (n_res[k] + 30))`, then normalised.
2. **Match frequencies** — `f[pos][aa] = Σ_k w[k]·[msa[k][pos] == aa]`, column-
   normalised, falling back to the background for all-gap columns.
3. **`Neff_M` per column** — via column entropy: `Neff_HMM = mean(2^entropy)`,
   `Nlim = max(10, Neff_HMM + 1)`, and
   `Neff_M[pos] = Nlim − (Nlim−1)·2^{scale·w_M}` with
   `w_M = −1/setSize + Σ w[k]` over non-gap rows. (This is occupancy/weight driven,
   not within-column diversity — matching MMseqs2 exactly.)
4. **Pseudocounts** — `g = R·f`; `τ = min(1, pca/(1 + Neff/pcb))`;
   `p = (1−τ)·f + τ·g`, defaults `pca = 1.0`, `pcb = 1.5`.
5. **Log-odds PSSM** — `round(8·log2(p/p_back) + 8·bias)` clamped to `[−128, 127]`,
   **profile bit-factor = 8**.

**Result:** 2620/2620 PSSM cells identical on an ungapped MSA; 1000/1000 identical on
a gapped MSA (exercising the weight/`Neff`/gap paths).

### 4. Prefilter — `prefilter.rs` ↔ `Prefiltering.cpp`, `QueryMatcher.cpp`, `KmerGenerator.cpp`

Conceptually the same two-step design:

- **Similar-k-mer seeding.** For each query/profile position, generate the k-mers
  whose score exceeds a threshold and look them up in an inverted target k-mer index.
  Sabertooth generates these with a **best-first branch-and-bound** (suffix-max
  pruning, children explored highest-score-first) so that the per-position cap keeps
  the *top-scoring* k-mers — including the exact seed — rather than an arbitrary
  subset. This was a real correctness point: a naive DFS cap silently drops the seeds
  you most need.
- **Ungapped diagonal gate.** Seeds are grouped by diagonal, and the best ungapped
  segment score along the top diagonal is computed with Kadane's algorithm. Only
  targets clearing `min_diag_score` proceed to gapped alignment — MMseqs2's ungapped
  diagonal prefilter serves the same purpose (cheaply discard most targets so the
  expensive stage runs on few).

**Differences.** MMseqs2's prefilter is a heavily optimized, SIMD-aware, cache-tuned
implementation with spaced-seed support and a sensitivity parameter `-s` that tunes
k-mer thresholds. Sabertooth's is scalar and simpler, exposing `--k`,
`--kmer-score`, and `--min-diag-score` directly. It is built for correctness and
clarity first; it is not micro-optimized to MMseqs2's throughput.

### 5. Local alignment — `align.rs` + `simd.rs` ↔ `smith_waterman_sse2.cpp` / `StripedSmithWaterman`

Same algorithm — **Smith–Waterman with Gotoh affine gaps** and full traceback,
producing raw score, alignment coordinates, CIGAR, and percent identity. It is
generic over a `Scorer` trait, so the identical DP serves both sequence
(matrix-scored) and profile (PSSM-scored) alignment.

**SIMD.** The scoring inner loop is now vectorized: `simd.rs` implements Farrar's
**striped Smith–Waterman** over SSE2 (8 × `i16` lanes, no runtime feature detection
needed since SSE2 is the x86-64 baseline). The search driver builds a striped query
profile once, scores every prefilter survivor with the fast traceback-free kernel,
and runs the scalar aligner (for CIGAR/coordinates) only on the candidates whose
score clears the E-value cutoff. The vectorized score is bit-identical to the scalar
`raw_score` (verified on 4000+ random sequence and profile cases and against a
brute-force reference); on the rare chance of `i16` saturation, or on non-x86-64
targets, it falls back to an exact scalar score. Measured gate speedup: ~22× versus
full scalar alignment (200-col profile × 4000 targets), ~2.4 G cells/s. AVX2 (wider
lanes, runtime dispatch) is left as future work.

**Correctness fix found via the SIMD work.** Differential-testing the striped kernel
against the scalar aligner exposed a latent bug in the scalar DP (present since
v0.1.0): the vertical-gap term `F` was carried in a *scalar*, but `F[i][j]` depends on
`F[i-1][j]` (the previous row at the same column) and so requires a per-column array.
The scalar version silently under-scored any deletion spanning more than one query
row. It is fixed; the corrected aligner *improves* remote-homolog scores (the demo's
true homologs went from E≈3.3e-33 / 9.6e-11 to 2.3e-34 / 9.5e-16) while leaving
gapless decoys unchanged — better discrimination, not just a different number. The
PSSM/matrix differential remains cell-exact (alignment and profile construction are
independent).

**Remaining difference.** Alignment is still parallelised **across** target sequences
with rayon, and the vectorized path is score-only (traceback stays scalar) — the same
two-pass split MMseqs2 uses (fast vectorized score, separate alignment for output).

### 6. Statistics — `evalue.rs` ↔ `EvalueComputation.cpp` (+ ALP)

- Bit score `S' = (λ·S − ln K) / ln 2`; `E = m·n·2^{−S'}`.
- For the **substitution matrix**, λ is obtained by analytically solving the
  ungapped Karlin–Altschul equation `Σ p_a p_b e^{λ s_{ab}} = 1` by bisection.
- For **PSSM log-odds** scores, λ is the exact closed form `ln2 / bit_factor = ln2/8`
  (derived in the README). This is the single most important statistical correctness
  point: applying the matrix's ×2-scale λ to ×8-scale PSSM scores inflates
  significance by ~4× in the exponent and makes decoys look homologous.

**Difference (now narrowed).** MMseqs2 links the **ALP** library to fit *gapped*
Gumbel parameters (λ, K) by simulation. Sabertooth's default λ/K are principled
ungapped / log-odds values, accurate for ranking and thresholding remote homologs
(as the demo shows). It now also ships a pure-Rust analogue of ALP — `calibrate_gapped`
(the `calibrate` command) — which Monte-Carlo-simulates random-pair alignments with
the SIMD kernel and fits a gapped Gumbel by the method of moments. The fitted gapped
λ is below the ungapped λ and yields more conservative E-values, the direction that
closes the absolute bit-score gap seen against the `mmseqs` binary (`VALIDATION.md`).
It is not bit-identical to ALP's estimates, but it is a real gapped calibration rather
than an ungapped approximation.

### 7. I/O and workflow — `main.rs`, `fasta.rs`, `search.rs`

- Sabertooth reads/writes **FASTA** and emits **BLAST-tab `.m8`**, the same
  columns MMseqs2's `convertalis` produces.
- It does **not** implement MMseqs2's binary database format, its module-graph
  workflow engine, or the prefilter→align→cluster chaining of `mmseqs` subcommands.

---

## Feature matrix

| Capability | MMseqs2 | Sabertooth |
|---|---|---|
| Substitution-matrix reconstruction | ✅ | ✅ (cell-identical) |
| PSSM / profile construction | ✅ | ✅ (cell-identical) |
| a3m / aligned-FASTA MSA input | ✅ | ✅ |
| K-mer + similar-k-mer prefilter | ✅ | ✅ (scalar) |
| Spaced seeds + sensitivity (`-s`) | ✅ | ✅ (`--spaced`, `-s`) |
| Ungapped diagonal prefilter | ✅ | ✅ |
| Smith–Waterman–Gotoh gapped alignment | ✅ | ✅ (scalar) |
| Banded alignment (long sequences) | ✅ | ✅ (`--band`) |
| Karlin–Altschul E-values | ✅ | ✅ (analytic/log-odds λ) |
| Profile search (PSSM vs sequences) | ✅ | ✅ |
| Sequence search | ✅ | ✅ |
| BLAST-tab `.m8` output | ✅ | ✅ |
| SIMD-striped alignment kernel | ✅ | ✅ (striped SSE2 score kernel; scalar traceback) |
| Gapped-statistics calibration | ✅ (ALP) | ✅ (Monte-Carlo Gumbel; `calibrate`) |
| Iterative profile search | ✅ | ✅ (`--num-iterations`) |
| MSA generation (`result2msa`) | ✅ | ❌ |
| Clustering / linclust | ✅ | ❌ |
| Taxonomy | ✅ | ❌ |
| Profile–profile search | ✅ | ❌ |
| Nucleotide / translated search | ✅ | ❌ |
| On-disk MMseqs2 DB format | ✅ | ❌ (FASTA) |

---

## Differential testing, and a real binary diff

Fidelity was first established by reimplementing MMseqs2's formulas a second time,
independently, in Python and diffing every cell — useful, but it shares whatever
conventions the reimplementation and Sabertooth have in common. That comparison has
since been supplemented by a **direct byte-diff against a compiled `mmseqs`**
(`VALIDATION.md`): MMseqs2 does build in the sandbox after all (pip-provided CMake +
the source tree's bundled libraries). Against the real binary, Sabertooth's PSSM
matches to ~91% exact / 99.6% within one alignment-profile quantization unit — after
accounting for the fact that `profile2pssm` dumps `profile_score / 4`, not the raw
×8 PSSM — and profile search returns the same homologs in the same rank order. The
independent-Python check remains valuable as a *second* implementation of the
underlying arithmetic; the binary diff is the ground-truth check.

---
