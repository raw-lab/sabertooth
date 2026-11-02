//! K-mer prefilter: seeds → bounded, target-blocked hit aggregation →
//! multi-diagonal ungapped / chained evidence → ranked candidates.
//!
//! ## Seeds
//!
//! For every query window the generator emits **exactly the top
//! `max_kmers_per_pos`** similar k-mers whose summed score (substitution matrix
//! for sequence queries, PSSM for profiles) reaches `kmer_score` — ties at the
//! boundary broken by a fixed rank order — so a per-window cap keeps the
//! genuinely best seeds. (v0.2 used a depth-first traversal with locally
//! sorted children; its first `cap` k-mers were not the global top `cap`.)
//! The enumeration is depth-first with suffix-maximum pruning; when the cap is
//! exceeded, one counting pass over partial-sum histograms gives the exact
//! boundary score, after which only k-mers above it (plus the first ties) are
//! kept (`top_kmers`; a best-first heap implementation serves as fallback and
//! test oracle). For sequence queries the query's own (exact) k-mer is always
//! seeded when it contains no `X` and no masked position, even if it is below
//! the threshold or the cap is reached.
//!
//! Lookups go through a k-mer presence bitmap first; k-mers sharing their
//! first `w − 1` residues are resolved by one bucket probe, and probes are
//! software-pipelined with prefetches.
//!
//! ## Bounded hit processing
//!
//! v0.2 materialised every `(target, diagonal)` hit of a query and sorted the
//! whole vector. Here posting lists (sorted by target) are consumed in **target
//! blocks** sized so a block holds at most `hit_block` hits: each seed keeps a
//! cursor into its posting list and a block takes the postings of targets in
//! `[T_b, T_{b+1})`. A target's evidence therefore never straddles a block
//! boundary — it is complete when its block is evaluated — and the hit buffer
//! (reused across queries, per worker) stays bounded. A single target with more
//! hits than `hit_block` grows the buffer and is counted as an overrun.
//!
//! ## Evidence per target
//!
//! Hits are grouped by diagonal; on the `max_diagonals` best-supported
//! diagonals the ungapped segment around each seed is found by X-drop
//! extension (`xdrop`, default 30 raw units; `--ungapped-xdrop 0` scores the
//! best segment anywhere on the diagonal by Kadane instead), and pairs of
//! segments that are co-linear in query and target are chained with an affine
//! gap penalty, so an indel-split homolog is not reduced to its single most
//! frequent diagonal.
//! A target passes when its best single-segment or chained score reaches
//! `min_diag_score`. Optional caps (`max_seed_occ`, `max_candidates`) are explicit
//! sensitivity trade-offs and every time they bind it is counted in
//! [`PrefilterStats`].

use crate::align::{GapCosts, ProfileScorer, ScoreTable, SeqScorer};
use crate::alphabet::{ALPHABET_SIZE, PROFILE_AA_SIZE};
use crate::fasta::SeqDb;
use crate::index::SeedIndex;
use crate::mask::MaskParams;
use crate::matrix::SubstitutionMatrix;
use crate::profile::Profile;
use crate::ungapped::{DiagKernel, DiagProfile, Lane};
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::time::Instant;

/// Number of real residues used for k-mer coding.
const A: u64 = PROFILE_AA_SIZE as u64; // 20

/// Largest supported k-mer length (`20^12 < 2^64`; k > ~9 is not useful for
/// protein sensitivity).
pub const MAX_K: usize = 12;
/// Smallest useful k-mer length.
pub const MIN_K: usize = 2;

/// Validate a requested k-mer length, returning a clear error if unsupported.
pub fn check_k(k: usize) -> Result<(), String> {
    if !(MIN_K..=MAX_K).contains(&k) {
        Err(format!(
            "k-mer length {} is out of the supported range {}..={} \
             (codes are packed into a base-20 u64)",
            k, MIN_K, MAX_K
        ))
    } else {
        Ok(())
    }
}

/// A (possibly spaced) seed pattern. `offsets` are the positions within a window
/// of length `span` that contribute residues to the k-mer code; `weight` is the
/// number of such positions (the effective k). The window always begins and ends
/// on a care position so the seed diagonal is defined at the window start.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeedPattern {
    offsets: Vec<usize>,
    span: usize,
}

impl SeedPattern {
    /// Contiguous seed of length `k`.
    pub fn contiguous(k: usize) -> Self {
        SeedPattern {
            offsets: (0..k).collect(),
            span: k,
        }
    }

    /// Parse a `0/1` mask such as `"110101"`. Must start and end with `1`.
    pub fn from_mask(mask: &str) -> Result<Self, String> {
        let bytes = mask.as_bytes();
        if bytes.is_empty() {
            return Err("empty seed mask".into());
        }
        let mut offsets = Vec::new();
        for (i, &b) in bytes.iter().enumerate() {
            match b {
                b'1' => offsets.push(i),
                b'0' => {}
                _ => return Err(format!("seed mask must be 0/1, got '{}'", b as char)),
            }
        }
        if offsets.is_empty() {
            return Err("seed mask has no care ('1') positions".into());
        }
        if offsets[0] != 0 || *offsets.last().unwrap() != bytes.len() - 1 {
            return Err("seed mask must start and end with '1'".into());
        }
        Ok(SeedPattern {
            offsets,
            span: bytes.len(),
        })
    }

    /// Number of coded positions (the effective k).
    #[inline]
    pub fn weight(&self) -> usize {
        self.offsets.len()
    }

    /// Total window length spanned by the pattern.
    #[inline]
    pub fn span(&self) -> usize {
        self.span
    }

    /// Care positions within the window.
    #[inline]
    pub fn offsets(&self) -> &[usize] {
        &self.offsets
    }

    /// Render back to a `0/1` mask.
    pub fn mask_string(&self) -> String {
        let mut m = vec![b'0'; self.span];
        for &o in &self.offsets {
            m[o] = b'1';
        }
        String::from_utf8(m).unwrap()
    }

    /// The default seed of sequence searches and `createindex`
    /// ([`DEFAULT_SEQUENCE_SEED`]).
    pub fn sequence_default() -> Self {
        Self::from_mask(DEFAULT_SEQUENCE_SEED).expect("valid default seed")
    }

    /// The default seed of profile searches (contiguous 5-mer).
    pub fn profile_default() -> Self {
        Self::contiguous(5)
    }
}

/// Default seed pattern of sequence searches and `createindex`: weight 6 over
/// a span of 9. At equal modeled cost it ranks more true homologs above the
/// first false positive than the contiguous 6-mer of v0.3 on both the SCOPe
/// tuning split and the 68 M-residue scale database (`bench/README.md`);
/// `--k 6` selects the contiguous seed.
pub const DEFAULT_SEQUENCE_SEED: &str = "111001101";

/// Map a sensitivity value (higher = more sensitive) to a similar-k-mer score
/// threshold, mirroring MMseqs2's `-s` in spirit (monotone; `s ≈ 5.7`
/// reproduces the default thresholds). `scale` is 2 for sequence scores and 8
/// for PSSM scores.
pub fn kmer_score_for_sensitivity(s: f64, scale: i32) -> i32 {
    let s = s.clamp(1.0, 9.0);
    let base = if scale >= 8 { 40.0 } else { 25.0 };
    let per_unit = 0.75 * scale as f64;
    let thr = base + (5.7 - s) * per_unit;
    thr.round().max(4.0) as i32
}

/// Sequence-search sensitivity presets: `-s` → `(kmer_score, max_kmers,
/// min_diag_score)`, the gate given for a [`GATE_REF_RESIDUES`]-residue
/// database (the CLI moves it with the database size, see
/// [`size_adjusted_gate`]). Anchors are points of the measured
/// speed–sensitivity frontier of the default spaced seed on the SCOPe tuning
/// queries against the 68 M-residue scale database (`bench/README.md`);
/// values between anchors are interpolated (k-mer score linearly, cap
/// geometrically, gate linearly). Lower k-mer thresholds and larger caps
/// generate more seeds (sensitivity, cost).
pub const SEQUENCE_PRESETS: [(f64, i32, usize, i32); 10] = [
    (1.0, 32, 64, 42),
    (2.0, 30, 128, 42),
    (3.0, 28, 128, 42),
    (4.0, 26, 256, 42),
    (5.0, 24, 256, 42),
    (5.7, 22, 512, 42),
    (6.5, 20, 512, 42),
    (7.5, 18, 1024, 42),
    (8.5, 18, 2048, 42),
    (9.0, 16, 4096, 42),
];

/// Interpolated sequence preset for sensitivity `s` (clamped to 1..9).
pub fn sequence_preset(s: f64) -> (i32, usize, i32) {
    let s = s.clamp(1.0, 9.0);
    let t = &SEQUENCE_PRESETS;
    let mut k = 0;
    while k + 2 < t.len() && s > t[k + 1].0 {
        k += 1;
    }
    let (s0, k0, c0, g0) = t[k];
    let (s1, k1, c1, g1) = t[k + 1];
    let f = ((s - s0) / (s1 - s0)).clamp(0.0, 1.0);
    let kmer = (k0 as f64 + f * (k1 - k0) as f64).round() as i32;
    let cap = ((c0 as f64).ln() + f * ((c1 as f64).ln() - (c0 as f64).ln())).exp();
    let cap = ((cap / 32.0).round() as usize).max(1) * 32;
    let gate = (g0 as f64 + f * (g1 - g0) as f64).round() as i32;
    (kmer, cap, gate)
}

/// Database size (residues) at which the [`SEQUENCE_PRESETS`] gates apply.
pub const GATE_REF_RESIDUES: f64 = 1.0e8;
/// Change of the ungapped gate (raw score units) per e-fold of database size
/// below [`GATE_REF_RESIDUES`] …
pub const GATE_SLOPE_BELOW: f64 = 2.5;
/// … and above it.
pub const GATE_SLOPE_ABOVE: f64 = 1.0;
/// Bounds of the size-adjusted gate.
pub const GATE_MIN: i32 = 20;
pub const GATE_MAX: i32 = 60;

/// The ungapped gate for a database of `db_residues` residues, given the gate
/// `gate` of a [`GATE_REF_RESIDUES`]-residue database: `gate + slope ·
/// ln(db_residues / GATE_REF_RESIDUES)` with [`GATE_SLOPE_BELOW`] for smaller
/// and [`GATE_SLOPE_ABOVE`] for larger databases, rounded and clamped to
/// `GATE_MIN..=GATE_MAX`.
///
/// Chance seed matches — and with them the diagonals to score and the
/// candidates sent to gapped alignment — grow linearly with the database,
/// while the per-query work that does not depend on it (similar k-mer
/// generation) stays fixed. On small databases gapped alignment is a small
/// part of the run, so a lower gate buys sensitivity cheaply; from about
/// 10^8 residues the database-proportional work dominates and the best gate
/// rises only slowly (the E-value a weak alignment would need also grows
/// with the database). Fitted on the SCOPe tuning split at 2.6, 17, 68 and
/// 265 M residues (`bench/README.md`). The CLI applies it to the preset gate
/// unless `--min-diag-score` is given; library users call it themselves.
pub fn size_adjusted_gate(gate: i32, db_residues: usize) -> i32 {
    let x = (db_residues.max(1) as f64 / GATE_REF_RESIDUES).ln();
    let slope = if x < 0.0 { GATE_SLOPE_BELOW } else { GATE_SLOPE_ABOVE };
    ((gate as f64 + slope * x).round() as i32).clamp(GATE_MIN, GATE_MAX)
}

/// Prefilter configuration. One definition serves the CLI, the library and
/// the distributed backend, so their defaults cannot drift apart.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "distributed", derive(serde::Serialize, serde::Deserialize))]
pub struct PrefilterParams {
    /// K-mer length (seed weight).
    pub k: usize,
    /// Minimum summed score for a similar k-mer to be generated.
    pub kmer_score: i32,
    /// Maximum similar k-mers generated per query window (best first).
    pub max_kmers_per_pos: usize,
    /// Minimum number of seed hits (any diagonal) for a target.
    pub min_hits: usize,
    /// Minimum ungapped (single segment or chained) score for a target to
    /// survive the prefilter, in raw score units.
    pub min_diag_score: i32,
    /// Minimum seed hits on one diagonal (1 = a single hit suffices; 2 =
    /// MMseqs2-style double hit).
    pub diag_hits: usize,
    /// Diagonals per target scored by ungapped extension.
    pub max_diagonals: usize,
    /// Keep at most this many candidates per query, best first (0 = all).
    pub max_candidates: usize,
    /// Skip k-mers whose posting list is longer than this (0 = never).
    pub max_seed_occ: usize,
    /// Always seed the query's own k-mer (sequence queries).
    pub exact_seeds: bool,
    /// Hits held per target block (bounds seed-hit working memory).
    pub hit_block: usize,
    /// Low-complexity masking of query seed windows.
    pub mask: Option<MaskParams>,
    /// Ungapped evidence per diagonal: `0` = best segment anywhere on the
    /// diagonal (Kadane over its full length); `x > 0` = X-drop extension
    /// (drop-off `x`, raw units) from each seed hit, so the cost is the
    /// segment length rather than the diagonal length.
    pub xdrop: i32,
    /// Compute each candidate's target region (`t_lo..=t_hi`, needed only by
    /// seed-guided windows). Off: diagonal-segment evidence reports
    /// `0..=u32::MAX` (the whole target) and skips the boundary computation.
    pub regions: bool,
    /// Segment evidence (`xdrop = 0`) is the best ungapped segment within
    /// this many residues of the outermost seed hits on each diagonal
    /// (0 = anywhere on the diagonal).
    pub diag_window: usize,
}

impl PrefilterParams {
    /// Sequence-search defaults (the CLI defaults).
    pub fn sequence_default() -> Self {
        // the `-s 5.7` preset; the gate is that of a 10^8-residue database
        // (the CLI adjusts it with `size_adjusted_gate`)
        PrefilterParams {
            k: 6,
            kmer_score: 22,
            max_kmers_per_pos: 512,
            min_hits: 1,
            min_diag_score: 42,
            diag_hits: 1,
            max_diagonals: 4,
            max_candidates: 0,
            max_seed_occ: 0,
            exact_seeds: true,
            hit_block: 1 << 20,
            mask: None,
            // segment evidence: the best ungapped segment within 48 residues
            // of the seed hits on each diagonal, scored in SIMD batches
            // (`xdrop > 0` selects v0.3's scalar X-drop extension instead)
            xdrop: 0,
            regions: false,
            diag_window: 48,
        }
    }

    /// Profile-search defaults (the CLI defaults).
    pub fn profile_default() -> Self {
        PrefilterParams {
            k: 5,
            kmer_score: 40,
            min_diag_score: 30,
            exact_seeds: false,
            ..Self::sequence_default()
        }
    }

    /// Validate parameter combinations.
    pub fn validate(&self) -> Result<(), String> {
        check_k(self.k)?;
        if self.max_kmers_per_pos == 0 {
            return Err("--max-kmers must be at least 1".into());
        }
        if self.max_diagonals == 0 {
            return Err("--prefilter-diags must be at least 1".into());
        }
        if self.diag_hits == 0 {
            return Err("--diag-hits must be at least 1".into());
        }
        if self.hit_block < 1024 {
            return Err("--hit-block must be at least 1024".into());
        }
        Ok(())
    }
}

impl Default for PrefilterParams {
    fn default() -> Self {
        Self::sequence_default()
    }
}

/// A prefilter candidate: a target plus seeding evidence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Candidate {
    pub target_id: u32,
    /// Seed hits on this target.
    pub hits: u32,
    /// Diagonal (`target_pos − query_pos`) of the best ungapped segment.
    pub best_diag: i64,
    /// Prefilter score: best ungapped segment or chain (raw units).
    pub diag_score: i32,
    /// Target region covered by the supporting segments (0-based, inclusive).
    pub t_lo: u32,
    pub t_hi: u32,
}

/// Counters for one or more prefilter runs (summable).
#[derive(Clone, Copy, Debug, Default)]
#[cfg_attr(feature = "distributed", derive(serde::Serialize, serde::Deserialize))]
pub struct PrefilterStats {
    pub queries: u64,
    pub windows: u64,
    pub windows_skipped: u64,
    pub kmers_generated: u64,
    pub windows_at_cap: u64,
    pub exact_seeds_added: u64,
    pub lookups: u64,
    pub prefix_probes: u64,
    pub presence_rejects: u64,
    pub seeds: u64,
    pub postings_visited: u64,
    pub seeds_over_occ_cap: u64,
    pub postings_skipped_occ: u64,
    pub blocks: u64,
    pub hit_buffer_peak: u64,
    pub hit_buffer_overruns: u64,
    pub targets_hit: u64,
    pub targets_evaluated: u64,
    pub diagonals_scored: u64,
    pub xdrop_extensions: u64,
    pub chains_used: u64,
    pub candidates: u64,
    pub candidates_capped: u64,
    pub ns_seeds: u64,
    pub ns_hits: u64,
    pub ns_evidence: u64,
}

impl PrefilterStats {
    pub fn add(&mut self, o: &PrefilterStats) {
        macro_rules! sum { ($($f:ident),*) => { $( self.$f += o.$f; )* } }
        sum!(
            queries, windows, windows_skipped, kmers_generated, windows_at_cap, exact_seeds_added,
            lookups, prefix_probes, presence_rejects, seeds, postings_visited, seeds_over_occ_cap,
            postings_skipped_occ, blocks, hit_buffer_overruns, targets_hit, targets_evaluated,
            diagonals_scored, xdrop_extensions, chains_used, candidates, candidates_capped, ns_seeds, ns_hits, ns_evidence
        );
        self.hit_buffer_peak = self.hit_buffer_peak.max(o.hit_buffer_peak);
    }
}

// ---------------------------------------------------------------------------
// exact similar-k-mer enumeration (top `cap` per window)
// ---------------------------------------------------------------------------

/// Residues of one window position sorted by descending score (ties: residue).
type SortedRow = [(i32, u8); PROFILE_AA_SIZE];

fn sorted_row(row: &[i32; ALPHABET_SIZE]) -> SortedRow {
    let mut v = [(0i32, 0u8); PROFILE_AA_SIZE];
    for (r, slot) in v.iter_mut().enumerate() {
        *slot = (row[r], r as u8);
    }
    v.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
    v
}

/// Reusable memory for [`top_kmers`].
#[derive(Default)]
struct KmerWork {
    dist: Vec<u64>,
    tmp: Vec<u64>,
    heap: BinaryHeap<HeapItem>,
}

/// Score ranges wider than this (summed over the window) use the best-first
/// fallback instead of the counting pass (never reached with `i8` matrices
/// and PSSMs, whose summed range is at most `12 * 255`).
const MAX_DP_RANGE: i64 = 1 << 16;

/// The top-`cap` similar k-mers of one window.
///
/// Canonical order: score descending, then the rank vector (each position's
/// residues ranked by score descending, residue ascending) lexicographically
/// ascending. The returned **set** is exactly the first `cap` k-mers with
/// score ≥ `threshold` in that order; the order of `out` is unspecified.
/// Returns `true` iff more than `cap` k-mers reach the threshold.
///
/// Method: depth-first enumeration in rank-lexicographic order with exact
/// suffix-maximum pruning (rows are sorted, so the first failing rank ends a
/// level). While at most `cap` k-mers have been found nothing else happens.
/// When the cap is exceeded, one counting pass over partial-sum histograms
/// (exact `u64` counts) gives the boundary score `t*` (the cap-th best score)
/// and the number of k-mers strictly above it; the collected list is reduced
/// to those above `t*` plus the first ties in DFS order, and the remaining
/// enumeration only accepts scores above `t*` (ties until the quota is full).
/// Each window therefore visits at most `2·cap + 1` leaves, instead of
/// paying a heap operation per emitted k-mer.
fn top_kmers(rows: &[SortedRow], threshold: i32, cap: usize, work: &mut KmerWork, out: &mut Vec<(u64, i32)>) -> bool {
    out.clear();
    let w = rows.len();
    if w == 0 || cap == 0 {
        return false;
    }
    assert!(w <= MAX_K, "seed weight {} exceeds {}", w, MAX_K);
    let mut smax = [0i32; MAX_K + 1];
    let mut range: i64 = 0;
    for i in (0..w).rev() {
        smax[i] = smax[i + 1] + rows[i][0].0;
        range += (rows[i][0].0 - rows[i][PROFILE_AA_SIZE - 1].0) as i64;
    }
    if smax[0] < threshold {
        return false;
    }
    if range > MAX_DP_RANGE {
        return top_kmers_heap(rows, threshold, cap, &mut work.heap, out);
    }
    let mut thr = threshold;
    let mut capped = false;
    let mut tstar = i32::MIN;
    let mut ties_left: u64 = 0;
    let last = w - 1;
    let mut rank = [0u8; MAX_K];
    let mut part = [0i32; MAX_K];
    let mut code = [0u64; MAX_K];
    let mut d = 0usize;
    loop {
        if d == last {
            let base_s = part[last];
            let base_c = code[last] * A;
            for &(s, aa) in rows[last].iter() {
                let ns = base_s + s;
                if ns < thr {
                    break;
                }
                out.push((base_c + aa as u64, ns));
                if capped {
                    if ns == tstar {
                        // thr == tstar only while ties remain
                        ties_left -= 1;
                        if ties_left == 0 {
                            thr = tstar + 1;
                        }
                    }
                } else if out.len() > cap {
                    capped = true;
                    let (t, above) = cap_boundary(rows, &smax[..=w], thr, cap, work);
                    tstar = t;
                    let need = cap as u64 - above;
                    let mut kept = 0u64;
                    out.retain(|&(_, sc)| {
                        if sc > tstar {
                            true
                        } else if sc == tstar && kept < need {
                            kept += 1;
                            true
                        } else {
                            false
                        }
                    });
                    ties_left = need - kept;
                    thr = if ties_left == 0 { tstar + 1 } else { tstar };
                }
            }
            if last == 0 {
                break;
            }
            d = last - 1;
            rank[d] += 1;
            continue;
        }
        let r = rank[d] as usize;
        if r < PROFILE_AA_SIZE {
            let (s, aa) = rows[d][r];
            let ns = part[d] + s;
            if ns + smax[d + 1] >= thr {
                part[d + 1] = ns;
                code[d + 1] = code[d] * A + aa as u64;
                d += 1;
                rank[d] = 0;
                continue;
            }
        }
        // level exhausted (later ranks score lower)
        if d == 0 {
            break;
        }
        d -= 1;
        rank[d] += 1;
    }
    debug_assert!(!capped || out.len() == cap);
    capped
}

/// Exact boundary of the top-`cap` list among k-mers scoring ≥ `thr`:
/// `(t, above)` where `t` is the largest score with at least `cap` k-mers
/// scoring ≥ `t` and `above` is the number scoring > `t` (< cap). Partial
/// sums that cannot reach `thr` are pruned, so each histogram spans at most
/// `max_score − thr + 1` bins. Requires at least `cap` k-mers ≥ `thr`.
fn cap_boundary(rows: &[SortedRow], smax: &[i32], thr: i32, cap: usize, work: &mut KmerWork) -> (i32, u64) {
    let KmerWork { dist, tmp, .. } = work;
    dist.clear();
    dist.push(1);
    let (mut lo, mut hi) = (0i32, 0i32);
    for (i, row) in rows.iter().enumerate() {
        let nhi = hi + row[0].0;
        let nlo = (lo + row[PROFILE_AA_SIZE - 1].0).max(thr.saturating_sub(smax[i + 1]));
        debug_assert!(nlo <= nhi);
        tmp.clear();
        tmp.resize((nhi - nlo + 1) as usize, 0);
        for (off, &cnt) in dist.iter().enumerate() {
            if cnt == 0 {
                continue;
            }
            let s = lo + off as i32;
            for &(v, _) in row.iter() {
                let t = s + v;
                if t < nlo {
                    break;
                }
                tmp[(t - nlo) as usize] += cnt;
            }
        }
        std::mem::swap(dist, tmp);
        lo = nlo;
        hi = nhi;
    }
    let mut above = 0u64;
    let mut t = hi;
    loop {
        let c = dist[(t - lo) as usize];
        if above + c >= cap as u64 {
            return (t, above);
        }
        above += c;
        assert!(t > lo, "cap_boundary: fewer than cap k-mers reach the threshold");
        t -= 1;
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
struct HeapItem {
    score: i32,
    code: u64,
    ranks: [u8; MAX_K],
    last: u8,
}

impl Ord for HeapItem {
    fn cmp(&self, o: &Self) -> Ordering {
        // max-heap by score; ties: lexicographically smaller rank vector first
        self.score.cmp(&o.score).then(o.ranks.cmp(&self.ranks))
    }
}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

/// Best-first (k-best sums) enumeration in the canonical order of
/// [`top_kmers`]: emits up to `cap` k-mers with score ≥ `threshold`, in order.
/// Every rank vector is generated once (children increment only coordinates
/// ≥ the parent's last incremented one) and scores never increase along an
/// increment, so the max-heap pops k-mers in global order. Used as the
/// fallback for extreme score ranges and as the test oracle.
fn top_kmers_heap(
    rows: &[SortedRow],
    threshold: i32,
    cap: usize,
    heap: &mut BinaryHeap<HeapItem>,
    out: &mut Vec<(u64, i32)>,
) -> bool {
    out.clear();
    heap.clear();
    let w = rows.len();
    let code_of = |ranks: &[u8; MAX_K]| -> u64 {
        let mut c = 0u64;
        for i in 0..w {
            c = c * A + rows[i][ranks[i] as usize].1 as u64;
        }
        c
    };
    let ranks = [0u8; MAX_K];
    let s0: i32 = rows.iter().map(|r| r[0].0).sum();
    if s0 < threshold || cap == 0 {
        return false;
    }
    heap.push(HeapItem {
        score: s0,
        code: code_of(&ranks),
        ranks,
        last: 0,
    });
    while let Some(it) = heap.pop() {
        if out.len() >= cap {
            return true;
        }
        out.push((it.code, it.score));
        for j in it.last as usize..w {
            let r = it.ranks[j] as usize;
            if r + 1 >= PROFILE_AA_SIZE {
                continue;
            }
            let ns = it.score - rows[j][r].0 + rows[j][r + 1].0;
            if ns < threshold {
                continue; // later ranks at j only score lower
            }
            let mut nr = it.ranks;
            nr[j] += 1;
            heap.push(HeapItem {
                score: ns,
                code: code_of(&nr),
                ranks: nr,
                last: j as u8,
            });
        }
    }
    false
}

// ---------------------------------------------------------------------------
// per-worker scratch
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct Seed {
    qpos: u32,
    cur: usize,
    end: usize,
}

/// Reusable per-worker prefilter memory.
#[derive(Default)]
pub struct PrefilterScratch {
    rows: Vec<SortedRow>,
    window_rows: Vec<SortedRow>,
    kw: KmerWork,
    kmers: Vec<(u64, i32)>,
    seeds: Vec<Seed>,
    /// packed seed hits of the current target block (and radix-sort buffer)
    keys: Vec<u64>,
    keys_tmp: Vec<u64>,
    ev: EvidenceScratch,
    dprof: DiagProfile,
    cands: Vec<Candidate>,
    /// per-prefix posting ranges (see [`SeedIndex::prefix_ranges`])
    ranges: [(usize, usize); 20],
    /// k-mer runs of the current window that need an index probe
    runs: Vec<(u32, u32)>,
}

/// Look-ahead distances (in probed runs) of the seed-lookup prefetches.
const PF_LOOKUP: usize = 16;
/// Look-ahead (in seeds) of the posting-list prefetches of the hit loop.
const PF_POSTINGS: usize = 8;
const PF_ENTRY: usize = 4;

impl PrefilterScratch {
    pub fn new() -> Self {
        Self::default()
    }
    /// Bytes reserved (memory accounting).
    pub fn heap_bytes(&self) -> usize {
        self.rows.capacity() * std::mem::size_of::<SortedRow>()
            + self.seeds.capacity() * std::mem::size_of::<Seed>()
            + (self.keys.capacity() + self.keys_tmp.capacity()) * 8
            + self.ev.heap_bytes()
            + self.cands.capacity() * std::mem::size_of::<Candidate>()
            + self.kmers.capacity() * 16
            + (self.kw.dist.capacity() + self.kw.tmp.capacity()) * 8
            + self.kw.heap.capacity() * std::mem::size_of::<HeapItem>()
    }
}

/// Query description for the prefilter.
pub struct PrefilterQuery<'a> {
    /// Seeding scores per query column (`[col][residue]`).
    pub seed_rows: &'a [[i32; ALPHABET_SIZE]],
    /// Query residues for exact seeds (sequence queries).
    pub residues: Option<&'a [u8]>,
    /// Scores for ungapped extension (may include composition-bias correction).
    pub ext_rows: &'a [[i32; ALPHABET_SIZE]],
    /// Query columns excluded from seeding (low-complexity mask).
    pub seed_mask: Option<&'a [bool]>,
}

/// Best ungapped segment on diagonal `d` that contains the seed starting at
/// query position `qs0` (`span` residues), by X-drop extension in both
/// directions: `(score, q_start, q_end)`.
#[inline]
fn xdrop_segment(rows: &[[i32; ALPHABET_SIZE]], t: &[u8], d: i64, qs0: usize, span: usize, x: i32) -> (i32, usize, usize) {
    let m = rows.len() as i64;
    let n = t.len() as i64;
    let q_lo = (-d).max(0) as usize; // first query position on the diagonal
    let q_hi = m.min(n - d) as usize; // one past the last
    let sc = |q: usize| -> i32 {
        let tr = t[(q as i64 + d) as usize] as usize;
        if tr < ALPHABET_SIZE {
            rows[q][tr]
        } else {
            0
        }
    };
    let s_end = (qs0 + span).min(q_hi);
    let mut seed = 0i32;
    for q in qs0..s_end {
        seed += sc(q);
    }
    // right
    let (mut cur, mut best, mut be) = (seed, seed, s_end.saturating_sub(1).max(qs0));
    let mut q = s_end;
    while q < q_hi {
        cur += sc(q);
        if cur > best {
            best = cur;
            be = q;
        } else if cur < best - x {
            break;
        }
        q += 1;
    }
    // left (gain over the right-extended segment)
    let (mut cur_l, mut best_l, mut bs) = (0i32, 0i32, qs0);
    let mut q = qs0;
    while q > q_lo {
        q -= 1;
        cur_l += sc(q);
        if cur_l > best_l {
            best_l = cur_l;
            bs = q;
        } else if cur_l < best_l - x {
            break;
        }
    }
    (best + best_l, bs, be)
}

/// Run the prefilter for one query. Candidates are returned best first
/// (prefilter score descending, then target id ascending).
pub fn prefilter_query(
    q: &PrefilterQuery,
    index: &SeedIndex,
    params: &PrefilterParams,
    gaps: GapCosts,
    sc: &mut PrefilterScratch,
    stats: &mut PrefilterStats,
) -> Vec<Candidate> {
    let pattern = index.pattern().clone();
    assert_eq!(
        pattern.weight(),
        params.k,
        "index seed weight {} does not match prefilter k {} (use PrefilterParams synced to the index)",
        pattern.weight(),
        params.k
    );
    let m = q.seed_rows.len();
    let span = pattern.span();
    stats.queries += 1;
    sc.cands.clear();
    if m < span || index.is_empty() {
        return Vec::new();
    }
    let t_seed = Instant::now();
    // ---- seeds -----------------------------------------------------------
    sc.rows.clear();
    sc.rows.extend(q.seed_rows.iter().map(sorted_row));
    sc.seeds.clear();
    let mut total_postings: usize = 0;
    let occ_cap = params.max_seed_occ;
    let direct = index.is_direct();
    // the presence bitmap serves the sparse layout only (built once per
    // index on first use; the direct table answers absence by itself)
    let presence = if direct { None } else { Some(index.presence()) };
    #[inline]
    fn push_seed(
        seeds: &mut Vec<Seed>,
        stats: &mut PrefilterStats,
        total: &mut usize,
        occ_cap: usize,
        qpos: u32,
        (a, b): (usize, usize),
    ) {
        let len = b - a;
        if occ_cap > 0 && len > occ_cap {
            stats.seeds_over_occ_cap += 1;
            stats.postings_skipped_occ += len as u64;
        } else {
            seeds.push(Seed { qpos, cur: a, end: b });
            *total += len;
        }
    }
    for start in 0..=(m - span) {
        stats.windows += 1;
        // eligibility: no X (sequence) and no masked care position
        let mut skip = false;
        let mut exact: u64 = 0;
        for &off in pattern.offsets() {
            let c = start + off;
            if let Some(mk) = q.seed_mask {
                if mk[c] {
                    skip = true;
                    break;
                }
            }
            if let Some(res) = q.residues {
                let r = res[c];
                if r as usize >= PROFILE_AA_SIZE {
                    skip = true;
                    break;
                }
                exact = exact * A + r as u64;
            }
        }
        if skip {
            stats.windows_skipped += 1;
            continue;
        }
        let qpos = start as u32;
        sc.window_rows.clear();
        for &off in pattern.offsets() {
            sc.window_rows.push(sc.rows[start + off]);
        }
        let at_cap = top_kmers(&sc.window_rows, params.kmer_score, params.max_kmers_per_pos, &mut sc.kw, &mut sc.kmers);
        if at_cap {
            stats.windows_at_cap += 1;
        }
        let nk = sc.kmers.len();
        let mut have_exact = false;
        if direct {
            // direct layout (dense k-mer spaces): one pair of offsets per
            // k-mer, fetched PF_LOOKUP k-mers ahead; a presence test would
            // rarely reject anything
            for r in 0..nk.min(PF_LOOKUP) {
                index.prefetch_code(sc.kmers[r].0);
            }
            for k in 0..nk {
                if k + PF_LOOKUP < nk {
                    index.prefetch_code(sc.kmers[k + PF_LOOKUP].0);
                }
                let c = sc.kmers[k].0;
                have_exact |= c == exact;
                let rg = index.direct_range(c);
                if rg.1 > rg.0 {
                    push_seed(&mut sc.seeds, stats, &mut total_postings, occ_cap, qpos, rg);
                }
            }
            stats.prefix_probes += nk as u64;
        } else {
            // sparse layout: resolve in runs sharing the first w−1 residues
            // (the generator emits k-mers grouped that way) — a presence-bitmap
            // test per k-mer, then one bucket probe per run with a present k-mer
            sc.runs.clear();
            let mut rejects = 0u64;
            let mut i = 0;
            while i < nk {
                let prefix = sc.kmers[i].0 / A;
                let mut j = i;
                let mut any = false;
                while j < nk && sc.kmers[j].0 / A == prefix {
                    let c = sc.kmers[j].0;
                    have_exact |= c == exact;
                    let p = presence.is_some_and(|pr| pr.may_contain(c));
                    any |= p;
                    rejects += !p as u64;
                    j += 1;
                }
                if any {
                    sc.runs.push((i as u32, j as u32));
                }
                i = j;
            }
            stats.presence_rejects += rejects;
            // probes with a software-pipelined look-ahead (bucket entries
            // PF_LOOKUP runs ahead, k-mer table entries PF_ENTRY ahead) — the
            // probes are independent random accesses into the index
            let nr = sc.runs.len();
            let run_prefix = |r: usize, kmers: &[(u64, i32)], runs: &[(u32, u32)]| kmers[runs[r].0 as usize].0 / A;
            for r in 0..nr.min(PF_LOOKUP) {
                index.prefetch_prefix_bucket(run_prefix(r, &sc.kmers, &sc.runs));
            }
            for r in 0..nr {
                if r + PF_LOOKUP < nr {
                    index.prefetch_prefix_bucket(run_prefix(r + PF_LOOKUP, &sc.kmers, &sc.runs));
                }
                if r + PF_ENTRY < nr {
                    index.prefetch_prefix_entries(run_prefix(r + PF_ENTRY, &sc.kmers, &sc.runs));
                }
                let (a, b) = sc.runs[r];
                let prefix = sc.kmers[a as usize].0 / A;
                let mask = index.prefix_ranges(prefix, &mut sc.ranges);
                stats.prefix_probes += 1;
                for k in a as usize..b as usize {
                    let c = sc.kmers[k].0;
                    let rr = (c % A) as usize;
                    if (mask >> rr) & 1 == 1 {
                        push_seed(&mut sc.seeds, stats, &mut total_postings, occ_cap, qpos, sc.ranges[rr]);
                    }
                }
            }
        }
        stats.lookups += nk as u64;
        let add_exact = q.residues.is_some() && params.exact_seeds && !have_exact;
        if add_exact {
            stats.exact_seeds_added += 1;
            stats.lookups += 1;
            stats.prefix_probes += 1;
            let rg = index.posting_range(exact);
            if rg.1 > rg.0 {
                push_seed(&mut sc.seeds, stats, &mut total_postings, occ_cap, qpos, rg);
            }
        }
        stats.kmers_generated += (nk + add_exact as usize) as u64;
    }
    stats.seeds += sc.seeds.len() as u64;
    stats.postings_visited += total_postings as u64;
    stats.ns_seeds += t_seed.elapsed().as_nanos() as u64;
    if sc.seeds.is_empty() {
        return Vec::new();
    }

    // ---- hits, by target block ------------------------------------------
    // A hit is packed into one u64: (target − block start, diagonal + m,
    // query position), most significant first, so one radix sort groups the
    // hits of a target and, within it, of each diagonal (query positions
    // ascending). Blocks are added when a block's target range would not fit
    // the key.
    let n_targets = index.len() as u64;
    let qbits = bits_for(m as u64 - 1);
    let dbits = bits_for(m as u64 + index.max_seq_len() as u64);
    assert!(qbits + dbits <= 62, "query/target too long for packed seed hits");
    let sh_t = qbits + dbits;
    let max_block_targets = 1u64 << (64 - sh_t).min(62);
    let blocks = (total_postings.div_ceil(params.hit_block).max(1) as u64).max(n_targets.div_ceil(max_block_targets));
    let kadane_mode = params.xdrop == 0;
    let dk = if kadane_mode { crate::ungapped::diag_kernel() } else { DiagKernel::Scalar };
    let simd = dk != DiagKernel::Scalar && sc.dprof.rebuild(q.ext_rows);
    let dk = if simd { dk } else { DiagKernel::Scalar };
    let mut out: Vec<Candidate> = Vec::new();
    let mut ns_hits = 0u64;
    let mut ns_ev = 0u64;
    for b in 0..blocks {
        let t_lo = b * n_targets / blocks;
        let t_hi = (b + 1) * n_targets / blocks;
        if t_hi <= t_lo {
            continue;
        }
        let th = Instant::now();
        stats.blocks += 1;
        sc.keys.clear();
        let mb = m as u64;
        let ns = sc.seeds.len();
        for si in 0..ns {
            // the posting lists are random accesses: fetch a few seeds ahead
            if si + PF_POSTINGS < ns {
                index.prefetch_posting(sc.seeds[si + PF_POSTINGS].cur);
            }
            let s = &mut sc.seeds[si];
            let (mut k, end) = (s.cur, s.end);
            let qp = s.qpos as u64;
            while k < end {
                let (tid, tpos) = index.posting(k);
                if tid >= t_hi {
                    break;
                }
                sc.keys.push(((tid - t_lo) << sh_t) | ((tpos + mb - qp) << qbits) | qp);
                k += 1;
            }
            s.cur = k;
        }
        let nh = sc.keys.len() as u64;
        stats.hit_buffer_peak = stats.hit_buffer_peak.max(nh);
        if nh as usize > params.hit_block {
            stats.hit_buffer_overruns += 1;
        }
        // whole-diagonal evidence does not need the query positions in order
        let lo_bit = if kadane_mode { qbits } else { 0 };
        radix_sort_u64(&mut sc.keys, &mut sc.keys_tmp, lo_bit, sh_t + bits_for(t_hi - t_lo - 1));
        ns_hits += th.elapsed().as_nanos() as u64;
        let te = Instant::now();
        let lay = KeyLayout {
            qbits,
            dmask: (1u64 << dbits) - 1,
            sh_t,
            t_lo,
            m: m as i64,
        };
        let PrefilterScratch { keys, ev, dprof, .. } = sc;
        if kadane_mode {
            evaluate_block_kadane(keys, &lay, q.ext_rows, dk, dprof, index, params, gaps, span, stats, ev, &mut out);
        } else {
            evaluate_block_xdrop(keys, &lay, q.ext_rows, index, params, gaps, span, stats, &mut out);
        }
        ns_ev += te.elapsed().as_nanos() as u64;
    }
    stats.ns_hits += ns_hits;
    stats.ns_evidence += ns_ev;

    // ---- rank + optional candidate budget ---------------------------------
    let order = |a: &Candidate, b: &Candidate| b.diag_score.cmp(&a.diag_score).then(a.target_id.cmp(&b.target_id));
    if params.max_candidates > 0 && out.len() > params.max_candidates {
        out.select_nth_unstable_by(params.max_candidates - 1, order);
        stats.candidates_capped += (out.len() - params.max_candidates) as u64;
        out.truncate(params.max_candidates);
    }
    out.sort_unstable_by(order);
    stats.candidates += out.len() as u64;
    out
}

/// Bits needed to represent every value in `0..=x` (0 for `x = 0`).
#[inline]
fn bits_for(x: u64) -> u32 {
    64 - x.leading_zeros()
}

/// LSD radix sort of `v` by bits `[lo_bit, hi_bit)` (11-bit digits; digits
/// that are equal in every key are skipped). Stable, so keys equal in the
/// sorted bits keep their input order.
fn radix_sort_u64(v: &mut Vec<u64>, tmp: &mut Vec<u64>, lo_bit: u32, hi_bit: u32) {
    const D: u32 = 11;
    const NB: usize = 1 << D;
    let n = v.len();
    if n < 2 || hi_bit <= lo_bit {
        return;
    }
    if n <= 48 {
        // tiny: insertion sort on the masked key (stable)
        let mask = if hi_bit >= 64 { u64::MAX } else { (1u64 << hi_bit) - 1 } & !((1u64 << lo_bit) - 1);
        for i in 1..n {
            let x = v[i];
            let kx = x & mask;
            let mut j = i;
            while j > 0 && v[j - 1] & mask > kx {
                v[j] = v[j - 1];
                j -= 1;
            }
            v[j] = x;
        }
        return;
    }
    tmp.clear();
    tmp.resize(n, 0);
    let mut cnt = [0u32; NB];
    let mut shift = lo_bit;
    while shift < hi_bit {
        let width = (hi_bit - shift).min(D);
        let dmask = (1u64 << width) - 1;
        cnt[..1 << width].fill(0);
        for &x in v.iter() {
            cnt[((x >> shift) & dmask) as usize] += 1;
        }
        if cnt[..1 << width].iter().any(|&c| c as usize == n) {
            shift += width;
            continue; // one bucket: this digit is constant
        }
        let mut sum = 0u32;
        for c in cnt[..1 << width].iter_mut() {
            let t = *c;
            *c = sum;
            sum += t;
        }
        for &x in v.iter() {
            let d = ((x >> shift) & dmask) as usize;
            // SAFETY-free indexing: cnt[d] < n by construction
            tmp[cnt[d] as usize] = x;
            cnt[d] += 1;
        }
        std::mem::swap(v, tmp);
        shift += width;
    }
}

/// Field layout of a packed seed hit (see `prefilter_query`).
struct KeyLayout {
    qbits: u32,
    dmask: u64,
    sh_t: u32,
    t_lo: u64,
    m: i64,
}

impl KeyLayout {
    #[inline]
    fn tid(&self, k: u64) -> u64 {
        (k >> self.sh_t) + self.t_lo
    }
    #[inline]
    fn diag(&self, k: u64) -> i64 {
        ((k >> self.qbits) & self.dmask) as i64 - self.m
    }
    #[inline]
    fn qpos(&self, k: u64) -> usize {
        (k & ((1u64 << self.qbits) - 1)) as usize
    }
    /// (target, diagonal) identity of a hit
    #[inline]
    fn td(&self, k: u64) -> u64 {
        k >> self.qbits
    }
}

#[derive(Clone, Copy, Default)]
struct DiagRun {
    diag: i64,
    count: u32,
    /// index of the run's first hit in the block's key array
    start: u32,
    /// smallest / largest query position of the run's seed hits
    qmin: u32,
    qmax: u32,
}

/// The `kmax` best-supported diagonals of one target's (sorted) hits: count
/// descending, then diagonal ascending.
#[inline]
fn top_diagonals(group: &[u64], base: usize, lay: &KeyLayout, kmax: usize, top: &mut Vec<DiagRun>) {
    top.clear();
    let mut a = 0;
    while a < group.len() {
        let key = lay.td(group[a]);
        let q0 = lay.qpos(group[a]) as u32;
        let (mut qmin, mut qmax) = (q0, q0);
        let mut b = a + 1;
        while b < group.len() && lay.td(group[b]) == key {
            let q = lay.qpos(group[b]) as u32;
            qmin = qmin.min(q);
            qmax = qmax.max(q);
            b += 1;
        }
        let run = DiagRun {
            diag: lay.diag(group[a]),
            count: (b - a) as u32,
            start: (base + a) as u32,
            qmin,
            qmax,
        };
        a = b;
        // insert keeping `top` sorted, size ≤ kmax
        let pos = top.iter().position(|t| run.count > t.count).unwrap_or(top.len());
        if pos < kmax {
            top.insert(pos, run);
            if top.len() > kmax {
                top.pop();
            }
        }
    }
}

/// Best single segment / co-linear chain of two segments, v0.3 rules:
/// `segs[k] = (score, q_start, q_end, diag)`. Returns `(score, used)`.
fn best_or_chain(segs: &[(i32, usize, usize, i64)], gaps: GapCosts, stats: &mut PrefilterStats) -> (i32, [usize; 2], usize) {
    let ns = segs.len();
    let mut best = (0i32, 0usize);
    for (k, sg) in segs.iter().enumerate() {
        if sg.0 > best.0 {
            best = (sg.0, k);
        }
    }
    let mut chain = (i32::MIN, 0usize, 0usize);
    for x in 0..ns {
        for y in 0..ns {
            if x == y {
                continue;
            }
            let (sx, _qsx, qex, dx) = segs[x];
            let (sy, qsy, _qey, dy) = segs[y];
            if sx <= 0 || sy <= 0 || dx == dy {
                continue;
            }
            let tex = qex as i64 + dx;
            let tsy = qsy as i64 + dy;
            if qex < qsy && tex < tsy {
                let g = gaps.cost((dy - dx).unsigned_abs() as usize);
                let c = sx + sy - g;
                if c > chain.0 {
                    chain = (c, x, y);
                }
            }
        }
    }
    if chain.0 > best.0 {
        stats.chains_used += 1;
        (chain.0, [chain.1, chain.2], best.1)
    } else {
        (best.0, [best.1, best.1], best.1)
    }
}

/// Target region covered by the used segments.
#[inline]
fn region(segs: &[(i32, usize, usize, i64)], used: [usize; 2]) -> (u32, u32) {
    let mut lo = u32::MAX;
    let mut hi = 0u32;
    for &u in &used {
        let (_, qs, qe, d) = segs[u];
        lo = lo.min((qs as i64 + d).max(0) as u32);
        hi = hi.max((qe as i64 + d).max(0) as u32);
    }
    (lo, hi)
}

/// X-drop evidence (`--ungapped-xdrop > 0`): one extension per seed hit not
/// already covered, scalar.
#[allow(clippy::too_many_arguments)]
fn evaluate_block_xdrop(
    keys: &[u64],
    lay: &KeyLayout,
    ext_rows: &[[i32; ALPHABET_SIZE]],
    index: &SeedIndex,
    params: &PrefilterParams,
    gaps: GapCosts,
    span: usize,
    stats: &mut PrefilterStats,
    out: &mut Vec<Candidate>,
) {
    let kmax = params.max_diagonals;
    let mut top: Vec<DiagRun> = Vec::with_capacity(kmax + 1);
    let mut i = 0;
    while i < keys.len() {
        let tid = lay.tid(keys[i]);
        let mut j = i + 1;
        while j < keys.len() && lay.tid(keys[j]) == tid {
            j += 1;
        }
        stats.targets_hit += 1;
        let group = &keys[i..j];
        let base = i;
        i = j;
        // look-ahead: fetch offsets 16 hits ahead and residues 8 hits ahead
        if let Some(&h) = keys.get(j + 16) {
            index.prefetch_offset(lay.tid(h) as usize);
        }
        if let Some(&h) = keys.get(j + 8) {
            let tpos = (lay.qpos(h) as i64 + lay.diag(h)).max(0) as usize;
            index.prefetch_residues(lay.tid(h) as usize, tpos);
        }
        if group.len() < params.min_hits {
            continue;
        }
        top_diagonals(group, base, lay, kmax, &mut top);
        if (top[0].count as usize) < params.diag_hits {
            continue;
        }
        stats.targets_evaluated += 1;
        let t = index.seq(tid as usize);
        let ns = top.len().min(16);
        let mut segs: [(i32, usize, usize, i64); 16] = [(0, 0, 0, 0); 16];
        for (k, run) in top.iter().take(ns).enumerate() {
            let mut seg = (i32::MIN, 0usize, 0usize);
            let mut covered_to: i64 = -1;
            for &h in &keys[run.start as usize..(run.start + run.count) as usize] {
                let qp = lay.qpos(h);
                if (qp as i64) <= covered_to {
                    continue;
                }
                let r = xdrop_segment(ext_rows, t, run.diag, qp, span, params.xdrop);
                stats.xdrop_extensions += 1;
                covered_to = r.2 as i64;
                if r.0 > seg.0 {
                    seg = r;
                }
            }
            stats.diagonals_scored += 1;
            segs[k] = (seg.0, seg.1, seg.2, run.diag);
        }
        let (score, used, bi) = best_or_chain(&segs[..ns], gaps, stats);
        if score < params.min_diag_score {
            continue;
        }
        let (lo, hi) = region(&segs[..ns], used);
        out.push(Candidate {
            target_id: tid as u32,
            hits: group.len() as u32,
            best_diag: segs[bi].3,
            diag_score: score,
            t_lo: lo,
            t_hi: hi,
        });
    }
}

/// One diagonal awaiting its batched score.
#[derive(Clone, Copy)]
struct Job {
    tid: u32,
    diag: i32,
    /// query positions to score
    lo: u32,
    hi: u32,
}

/// One target awaiting batched diagonal scores.
#[derive(Clone, Copy)]
struct TargetRec {
    tid: u32,
    hits: u32,
    j0: u32,
    j1: u32,
}

/// Diagonals scored per evidence chunk: bounds the evidence buffers
/// independently of the hit-block size.
const EV_CHUNK: usize = 1 << 14;

/// Reusable buffers of the segment evidence.
#[derive(Default)]
pub(crate) struct EvidenceScratch {
    top: Vec<DiagRun>,
    targets: Vec<TargetRec>,
    jobs: Vec<Job>,
    /// jobs ordered by window start (lanes of a batch share query positions)
    order: Vec<u32>,
    counts: Vec<u32>,
    scores: Vec<i32>,
    lane_out: Vec<i32>,
    diag: crate::ungapped::DiagScratch,
}

impl EvidenceScratch {
    fn heap_bytes(&self) -> usize {
        self.targets.capacity() * std::mem::size_of::<TargetRec>()
            + self.jobs.capacity() * std::mem::size_of::<Job>()
            + (self.order.capacity() + self.counts.capacity() + self.scores.capacity() + self.lane_out.capacity()) * 4
    }
}

/// Segment evidence (`--ungapped-xdrop 0`, the default): the best ungapped
/// segment on each kept diagonal (within `diag_window` of its seed hits, or
/// anywhere on it), scored in SIMD batches (see [`crate::ungapped`]);
/// segment boundaries are computed only where a chain could change the score
/// and, with `regions`, for passing targets. With `diag_window = 0` the
/// results equal v0.3's scalar whole-diagonal mode.
#[allow(clippy::too_many_arguments)]
fn evaluate_block_kadane(
    keys: &[u64],
    lay: &KeyLayout,
    ext_rows: &[[i32; ALPHABET_SIZE]],
    dk: DiagKernel,
    dprof: &DiagProfile,
    index: &SeedIndex,
    params: &PrefilterParams,
    gaps: GapCosts,
    span: usize,
    stats: &mut PrefilterStats,
    ev: &mut EvidenceScratch,
    out: &mut Vec<Candidate>,
) {
    let kmax = params.max_diagonals;
    let mq = ext_rows.len();
    ev.targets.clear();
    ev.jobs.clear();
    // group hits by target, choose diagonals; score in chunks
    let mut i = 0;
    while i < keys.len() {
        let tid = lay.tid(keys[i]);
        let mut j = i + 1;
        while j < keys.len() && lay.tid(keys[j]) == tid {
            j += 1;
        }
        stats.targets_hit += 1;
        let group = &keys[i..j];
        let base = i;
        i = j;
        if let Some(&h) = keys.get(j + 16) {
            index.prefetch_offset(lay.tid(h) as usize);
        }
        if group.len() < params.min_hits {
            continue;
        }
        if group.len() == 1 {
            // the common case: a single seed hit
            let h = group[0];
            let q = lay.qpos(h) as u32;
            ev.top.clear();
            ev.top.push(DiagRun {
                diag: lay.diag(h),
                count: 1,
                start: base as u32,
                qmin: q,
                qmax: q,
            });
        } else {
            top_diagonals(group, base, lay, kmax, &mut ev.top);
        }
        if (ev.top[0].count as usize) < params.diag_hits {
            continue;
        }
        stats.targets_evaluated += 1;
        let j0 = ev.jobs.len() as u32;
        let len = index.seq_len(tid as usize);
        for run in ev.top.iter().take(16) {
            let (mut lo, mut hi) = crate::ungapped::diag_range(mq, len, run.diag);
            if params.diag_window > 0 {
                lo = lo.max((run.qmin as usize).saturating_sub(params.diag_window));
                hi = hi.min(run.qmax as usize + span + params.diag_window);
            }
            ev.jobs.push(Job {
                tid: tid as u32,
                diag: run.diag as i32,
                lo: lo as u32,
                hi: hi.max(lo) as u32,
            });
        }
        ev.targets.push(TargetRec {
            tid: tid as u32,
            hits: group.len() as u32,
            j0,
            j1: ev.jobs.len() as u32,
        });
        if ev.jobs.len() >= EV_CHUNK {
            score_and_gate(ext_rows, dk, dprof, index, params, gaps, stats, ev, out);
        }
    }
    score_and_gate(ext_rows, dk, dprof, index, params, gaps, stats, ev, out);
}

/// Score the pending diagonals of `ev` in lane batches, then chain, gate and
/// emit candidates for the pending targets; clears both lists.
#[allow(clippy::too_many_arguments)]
fn score_and_gate(
    ext_rows: &[[i32; ALPHABET_SIZE]],
    dk: DiagKernel,
    dprof: &DiagProfile,
    index: &SeedIndex,
    params: &PrefilterParams,
    gaps: GapCosts,
    stats: &mut PrefilterStats,
    ev: &mut EvidenceScratch,
    out: &mut Vec<Candidate>,
) {
    let njobs = ev.jobs.len();
    if njobs == 0 {
        ev.targets.clear();
        return;
    }
    stats.diagonals_scored += njobs as u64;
    let mq = ext_rows.len();
    let nl = dk.lanes();
    ev.scores.clear();
    ev.scores.resize(njobs, 0);
    ev.lane_out.clear();
    ev.lane_out.resize(nl, 0);
    // ---- order jobs by window start (counting sort) so the lanes of a batch
    // cover the same query positions
    ev.counts.clear();
    ev.counts.resize(mq + 1, 0);
    for jb in &ev.jobs {
        ev.counts[jb.lo as usize] += 1;
    }
    let mut acc = 0u32;
    for c in ev.counts.iter_mut() {
        let t = *c;
        *c = acc;
        acc += t;
    }
    ev.order.clear();
    ev.order.resize(njobs, 0);
    for (j, jb) in ev.jobs.iter().enumerate() {
        let slot = &mut ev.counts[jb.lo as usize];
        ev.order[*slot as usize] = j as u32;
        *slot += 1;
    }
    // ---- score in lane batches, prefetching the next batch's residues
    let prefetch_batch = |order: &[u32], jobs: &[Job]| {
        for &j in order {
            let jb = &jobs[j as usize];
            if jb.lo < jb.hi {
                let (start, _) = index.seq_span(jb.tid as usize);
                let a = (jb.lo as i64 + jb.diag as i64) as usize;
                let b = (jb.hi as i64 + jb.diag as i64) as usize;
                index.prefetch_residue_range(start + a, start + b);
            }
        }
    };
    let mut lanes: Vec<Lane> = Vec::with_capacity(nl);
    prefetch_batch(&ev.order[..nl.min(njobs)], &ev.jobs);
    let mut c0 = 0;
    while c0 < njobs {
        let c1 = (c0 + nl).min(njobs);
        let p1 = (c1 + nl).min(njobs);
        prefetch_batch(&ev.order[c1..p1], &ev.jobs);
        lanes.clear();
        for &j in &ev.order[c0..c1] {
            let jb = &ev.jobs[j as usize];
            lanes.push(Lane {
                target: index.seq(jb.tid as usize),
                diag: jb.diag as i64,
                lo: jb.lo as usize,
                hi: jb.hi as usize,
            });
        }
        crate::ungapped::score_batch(dk, dprof, ext_rows, &lanes, &mut ev.lane_out, &mut ev.diag);
        for (i, &j) in ev.order[c0..c1].iter().enumerate() {
            ev.scores[j as usize] = ev.lane_out[i];
        }
        c0 = c1;
    }
    // ---- chains, gate, candidates
    let mut segs: [(i32, usize, usize, i64); 16] = [(0, 0, 0, 0); 16];
    let mut have_pos = [false; 16];
    for tr in &ev.targets {
        let js = tr.j0 as usize..tr.j1 as usize;
        let ns = js.len();
        if ns == 1 && ev.scores[tr.j0 as usize] < params.min_diag_score {
            continue; // the common case: one diagonal, below the gate
        }
        let mut best = (0i32, 0usize);
        for (k, jj) in js.clone().enumerate() {
            let sc = ev.scores[jj];
            segs[k] = (sc, 0, 0, ev.jobs[jj].diag as i64);
            have_pos[k] = false;
            if sc > best.0 {
                best = (sc, k);
            }
        }
        let t = index.seq(tr.tid as usize);
        let mut pos_of = |k: usize, segs: &mut [(i32, usize, usize, i64); 16]| {
            if !have_pos[k] {
                let jb = &ev.jobs[tr.j0 as usize + k];
                let (s2, qs, qe) = crate::ungapped::kadane_segment(ext_rows, t, segs[k].3, jb.lo as usize, jb.hi as usize);
                debug_assert_eq!(s2, segs[k].0);
                segs[k].1 = qs;
                segs[k].2 = qe;
                have_pos[k] = true;
            }
        };
        // chains of two co-linear segments: only pairs whose chain value
        // could exceed the best single segment need boundaries
        let mut chain = (i32::MIN, 0usize, 0usize);
        if ns > 1 {
            for x in 0..ns {
                for y in 0..ns {
                    if x == y {
                        continue;
                    }
                    let (sx, dx) = (segs[x].0, segs[x].3);
                    let (sy, dy) = (segs[y].0, segs[y].3);
                    if sx <= 0 || sy <= 0 || dx == dy {
                        continue;
                    }
                    let c = sx + sy - gaps.cost((dy - dx).unsigned_abs() as usize);
                    if c <= best.0 || c <= chain.0 {
                        continue; // cannot change the result
                    }
                    pos_of(x, &mut segs);
                    pos_of(y, &mut segs);
                    let (qex, qsy) = (segs[x].2, segs[y].1);
                    if qex < qsy && (qex as i64 + dx) < (qsy as i64 + dy) {
                        chain = (c, x, y);
                    }
                }
            }
        }
        let (score, used) = if chain.0 > best.0 {
            stats.chains_used += 1;
            (chain.0, [chain.1, chain.2])
        } else {
            (best.0, [best.1, best.1])
        };
        if score < params.min_diag_score {
            continue;
        }
        let (lo, hi) = if params.regions {
            for &u in &used {
                pos_of(u, &mut segs);
            }
            region(&segs[..ns], used)
        } else {
            (0, u32::MAX)
        };
        out.push(Candidate {
            target_id: tr.tid,
            hits: tr.hits,
            best_diag: segs[best.1].3,
            diag_score: score,
            t_lo: lo,
            t_hi: hi,
        });
    }
    ev.jobs.clear();
    ev.targets.clear();
}

// ---------------------------------------------------------------------------
// convenience wrappers (library API compatibility)
// ---------------------------------------------------------------------------

/// The seed index type (kept under its v0.2 name).
pub type KmerIndex = SeedIndex;

impl SeedIndex {
    /// Build a transient contiguous-seed index over an in-memory database.
    pub fn build(db: &SeqDb, k: usize) -> SeedIndex {
        Self::build_spaced(db, SeedPattern::contiguous(k))
    }

    /// Build a transient index with an arbitrary seed pattern (no masking).
    pub fn build_spaced(db: &SeqDb, pattern: SeedPattern) -> SeedIndex {
        let opts = crate::index::BuildOptions::new(pattern);
        SeedIndex::build_transient(db, &opts)
            .expect("building a transient seed index failed")
            .0
    }
}

/// Seeding rows of a scorer.
fn rows_of<S: crate::align::Scorer>(s: &S) -> Vec<[i32; ALPHABET_SIZE]> {
    ScoreTable::from_scorer(s).rows().to_vec()
}

/// Prefilter a plain sequence query (no masking, no composition correction).
pub fn prefilter_sequence(
    query: &[u8],
    _db: &SeqDb,
    index: &SeedIndex,
    mat: &SubstitutionMatrix,
    params: PrefilterParams,
) -> Vec<Candidate> {
    let rows = rows_of(&SeqScorer { query, mat });
    let q = PrefilterQuery {
        seed_rows: &rows,
        residues: Some(query),
        ext_rows: &rows,
        seed_mask: None,
    };
    let mut sc = PrefilterScratch::new();
    let mut st = PrefilterStats::default();
    prefilter_query(&q, index, &params, GapCosts::sequence_default(), &mut sc, &mut st)
}

/// Prefilter a profile (PSSM) query.
pub fn prefilter_profile(profile: &Profile, _db: &SeqDb, index: &SeedIndex, params: PrefilterParams) -> Vec<Candidate> {
    let rows = rows_of(&ProfileScorer { profile });
    let q = PrefilterQuery {
        seed_rows: &rows,
        residues: None,
        ext_rows: &rows,
        seed_mask: None,
    };
    let mut sc = PrefilterScratch::new();
    let mut st = PrefilterStats::default();
    prefilter_query(&q, index, &params, GapCosts::profile_default(), &mut sc, &mut st)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::build_aa2num;
    use crate::profile::PseudoCountParams;
    use std::io::Cursor;

    fn db_from(fasta: &str) -> SeqDb {
        let t = build_aa2num();
        SeqDb::from_reader(Cursor::new(fasta), &t).unwrap()
    }

    fn enc(s: &str) -> Vec<u8> {
        let t = build_aa2num();
        s.bytes().map(|b| t[b as usize]).collect()
    }

    /// Rank vector of `code` (canonical tie order of `top_kmers`).
    fn rank_vector(sorted: &[SortedRow], code: u64) -> Vec<u8> {
        let w = sorted.len();
        let mut digits = vec![0u8; w];
        let mut x = code;
        for i in (0..w).rev() {
            digits[i] = (x % 20) as u8;
            x /= 20;
        }
        (0..w)
            .map(|i| sorted[i].iter().position(|&(_, r)| r == digits[i]).unwrap() as u8)
            .collect()
    }

    #[test]
    fn top_kmers_is_exact_global_order() {
        // brute force over all 20^w k-mers for small w and random score rows,
        // including narrow score ranges that produce many ties at the boundary
        let mut rng = crate::evalue::SplitMix64::new(11);
        let mut kw = KmerWork::default();
        let mut capped_cases = 0;
        for w in 2..=4usize {
            for case in 0..40 {
                let spread = if case % 2 == 0 { 15 } else { 4 };
                let rows: Vec<[i32; ALPHABET_SIZE]> = (0..w)
                    .map(|_| {
                        let mut r = [0i32; ALPHABET_SIZE];
                        for x in r.iter_mut() {
                            *x = rng.below(spread) as i32 - (spread as i32) / 2;
                        }
                        r
                    })
                    .collect();
                let sorted: Vec<SortedRow> = rows.iter().map(sorted_row).collect();
                let thr = rng.below(10) as i32 - 3;
                let mut all: Vec<(u64, i32)> = Vec::new();
                for c in 0..20u64.pow(w as u32) {
                    let mut x = c;
                    let mut sc = 0;
                    for i in (0..w).rev() {
                        sc += rows[i][(x % 20) as usize];
                        x /= 20;
                    }
                    if sc >= thr {
                        all.push((c, sc));
                    }
                }
                let canon = |v: &mut Vec<(u64, i32)>| {
                    let mut keyed: Vec<(i32, Vec<u8>, u64, i32)> =
                        v.iter().map(|&(c, s)| (-s, rank_vector(&sorted, c), c, s)).collect();
                    keyed.sort_unstable();
                    v.clear();
                    v.extend(keyed.into_iter().map(|(_, _, c, s)| (c, s)));
                };
                canon(&mut all);
                for cap in [1usize, 2, 5, 50, 333, 100000] {
                    let want = &all[..all.len().min(cap)];
                    // best-first oracle: exact canonical order
                    let mut heap = BinaryHeap::new();
                    let mut out_h = Vec::new();
                    let at_cap_h = top_kmers_heap(&sorted, thr, cap, &mut heap, &mut out_h);
                    assert_eq!(&out_h[..], want, "heap w={} thr={} cap={}", w, thr, cap);
                    assert_eq!(at_cap_h, all.len() > cap);
                    // depth-first + exact boundary: same set
                    let mut out = Vec::new();
                    let at_cap = top_kmers(&sorted, thr, cap, &mut kw, &mut out);
                    canon(&mut out);
                    assert_eq!(&out[..], want, "dfs w={} thr={} cap={}", w, thr, cap);
                    assert_eq!(at_cap, all.len() > cap);
                    capped_cases += at_cap as usize;
                }
            }
        }
        assert!(capped_cases > 100, "too few capped cases exercised: {}", capped_cases);
    }

    #[test]
    fn top_kmers_matches_oracle_on_blosum_windows() {
        // real matrix rows, k = 6, thresholds/caps used in practice
        let m = SubstitutionMatrix::blosum62();
        let mut rng = crate::evalue::SplitMix64::new(5);
        let mut kw = KmerWork::default();
        let mut heap = BinaryHeap::new();
        let (mut a, mut b) = (Vec::new(), Vec::new());
        for _ in 0..300 {
            let q: Vec<u8> = (0..6).map(|_| rng.below(20) as u8).collect();
            let rows = rows_of(&SeqScorer { query: &q, mat: &m });
            let sorted: Vec<SortedRow> = rows.iter().map(sorted_row).collect();
            for (thr, cap) in [(25, 256), (18, 1024), (12, 64), (30, 8)] {
                let ca = top_kmers(&sorted, thr, cap, &mut kw, &mut a);
                let cb = top_kmers_heap(&sorted, thr, cap, &mut heap, &mut b);
                assert_eq!(ca, cb);
                a.sort_unstable();
                b.sort_unstable();
                assert_eq!(a, b, "query {:?} thr {} cap {}", q, thr, cap);
            }
        }
    }

    #[test]
    fn old_dfs_generator_was_not_globally_ordered() {
        // Documents the v0.2 defect: with a cap, depth-first emission (best child
        // first) is not the global top list. Rows where position 0's runner-up
        // combined with position 1's best beats position 0's best with position
        // 1's runner-up.
        let mut r0 = [-10i32; ALPHABET_SIZE];
        r0[0] = 10;
        r0[1] = 9;
        let mut r1 = [-10i32; ALPHABET_SIZE];
        r1[0] = 10;
        r1[1] = 0;
        let sorted = vec![sorted_row(&r0), sorted_row(&r1)];
        let mut out = Vec::new();
        top_kmers(&sorted, 0, 2, &mut KmerWork::default(), &mut out);
        out.sort_unstable();
        // global top-2: (0,0)=20, (1,0)=19 — DFS would give (0,0)=20, (0,1)=10
        assert_eq!(out, vec![(0, 20), (20, 19)]); // codes 0·20+0 and 1·20+0
    }

    #[test]
    fn exact_kmer_is_found() {
        let m = SubstitutionMatrix::blosum62();
        let db = db_from(">t1\nGGGGGMKVLLACDEFGHIKGGGG\n>t2\nPPPPPPPPPPPPPP\n");
        let idx = KmerIndex::build(&db, 6);
        let q = enc("MKVLLACDEFGHIK");
        let cands = prefilter_sequence(&q, &db, &idx, &m, PrefilterParams::default());
        assert!(!cands.is_empty());
        assert_eq!(cands[0].target_id, 0);
    }

    #[test]
    fn exact_seed_below_threshold_is_kept() {
        // "AAAAAA" self-scores 24 < 25: only the exact-seed rule finds it
        let m = SubstitutionMatrix::blosum62();
        let db = db_from(">t\nWWWWAAAAAAWWWW\n");
        let idx = KmerIndex::build(&db, 6);
        let q = enc("AAAAAA");
        let mut p = PrefilterParams::default();
        p.kmer_score = 25;
        p.min_diag_score = 1;
        p.exact_seeds = false;
        assert!(prefilter_sequence(&q, &db, &idx, &m, p).is_empty());
        p.exact_seeds = true;
        assert_eq!(prefilter_sequence(&q, &db, &idx, &m, p).len(), 1);
    }

    #[test]
    fn similar_kmer_reaches_diverged_target() {
        let m = SubstitutionMatrix::blosum62();
        let db = db_from(">t\nGGGGMKILLSCDEFGHLKGGGG\n");
        let idx = KmerIndex::build(&db, 6);
        let q = enc("MKVLLACDEFGHIK");
        let mut p = PrefilterParams::default();
        p.kmer_score = 15;
        let cands = prefilter_sequence(&q, &db, &idx, &m, p);
        assert!(!cands.is_empty(), "similar-kmer seeding failed to find diverged target");
    }

    #[test]
    fn profile_prefilter_finds_target() {
        let m = SubstitutionMatrix::blosum62();
        let msa = ">q\nMKVLLACDEFGHIK\n>h\nMKILLSCDEFGHLK\n";
        let prof = Profile::from_msa_str(msa, &m, PseudoCountParams::default()).unwrap();
        let db = db_from(">t\nGGGGMKVLLACDEFGHIKGGGG\n");
        let idx = KmerIndex::build(&db, 6);
        let mut p = PrefilterParams::profile_default();
        p.k = 6;
        let cands = prefilter_profile(&prof, &db, &idx, p);
        assert!(!cands.is_empty());
    }

    #[test]
    fn k_out_of_range_is_rejected() {
        assert!(check_k(1).is_err());
        assert!(check_k(MAX_K + 1).is_err());
        assert!(check_k(6).is_ok());
        assert!(check_k(MAX_K).is_ok());
    }

    #[test]
    fn long_kmers_do_not_collide_under_u64() {
        let db = db_from(">a\nACDEFGHIKLMNPQRSTVWY\n");
        let idx = KmerIndex::build(&db, 8);
        assert_eq!(idx.distinct_kmers(), 13);
    }

    #[test]
    fn diagonal_mode_is_selected() {
        let m = SubstitutionMatrix::blosum62();
        let db = db_from(">t\nMKVLLACDEFGHIKMKVLLAC\n");
        let idx = KmerIndex::build(&db, 6);
        let q = enc("MKVLLACDEFGHIK");
        let cands = prefilter_sequence(&q, &db, &idx, &m, PrefilterParams::default());
        assert!(!cands.is_empty());
        assert!(cands[0].diag_score > 0);
        assert_eq!(cands[0].best_diag, 0);
    }

    #[test]
    fn second_diagonal_and_chain_rescue_indel_split_homolog() {
        // Target: a repeat of a short query motif on a wrong diagonal (most
        // frequent), plus the real homology split by a 4-residue insertion.
        // The most-frequent-diagonal rule (v0.2) scores only the repeat.
        let m = SubstitutionMatrix::blosum62();
        let q = "MKVLLACDEFGHIKLMNPQRSTVWYHHHHWWCCMM";
        let t = format!("{}{}{}", "QRSTQRSTQRSTQRSTQRSTQRSTQRST", "MKVLLACDEFGHIKLMNP", "GGGGQRSTVWYHHHHWWCCMM");
        let db = db_from(&format!(">t\n{}\n", t));
        let idx = KmerIndex::build(&db, 4);
        let mut p = PrefilterParams::default();
        p.k = 4;
        p.kmer_score = 18;
        let qn = enc(q);
        let mut p1 = p;
        p1.max_diagonals = 1;
        let one = prefilter_sequence(&qn, &db, &idx, &m, p1);
        let four = prefilter_sequence(&qn, &db, &idx, &m, p);
        assert_eq!(four.len(), 1);
        assert!(four[0].diag_score > one.first().map(|c| c.diag_score).unwrap_or(0));
    }

    #[test]
    fn target_blocks_give_identical_candidates() {
        let m = SubstitutionMatrix::blosum62();
        let mut rng = crate::evalue::SplitMix64::new(4);
        let aa = b"ACDEFGHIKLMNPQRSTVWY";
        let core: String = (0..80).map(|_| aa[rng.below(20)] as char).collect();
        let mut fasta = String::new();
        for i in 0..400 {
            let mut s: Vec<u8> = core.bytes().collect();
            for _ in 0..(i % 30) {
                let p = rng.below(s.len());
                s[p] = aa[rng.below(20)];
            }
            let pre: String = (0..rng.below(50)).map(|_| aa[rng.below(20)] as char).collect();
            fasta.push_str(&format!(">s{}\n{}{}\n", i, pre, String::from_utf8(s).unwrap()));
        }
        let db = db_from(&fasta);
        let idx = KmerIndex::build(&db, 5);
        let rows = rows_of(&SeqScorer { query: &enc(&core), mat: &m });
        let qn = enc(&core);
        let q = PrefilterQuery { seed_rows: &rows, residues: Some(&qn), ext_rows: &rows, seed_mask: None };
        let mut p = PrefilterParams::default();
        p.k = 5;
        p.kmer_score = 18;
        let mut big = PrefilterStats::default();
        let mut small = PrefilterStats::default();
        let a = prefilter_query(&q, &idx, &p, GapCosts::sequence_default(), &mut PrefilterScratch::new(), &mut big);
        p.hit_block = 1024;
        let b = prefilter_query(&q, &idx, &p, GapCosts::sequence_default(), &mut PrefilterScratch::new(), &mut small);
        assert!(small.blocks > 4, "expected several target blocks, got {}", small.blocks);
        assert_eq!(big.blocks, 1);
        assert!(small.hit_buffer_peak < big.hit_buffer_peak);
        assert_eq!(a, b);
        // caps are counted
        p.max_candidates = 10;
        let mut st = PrefilterStats::default();
        let c = prefilter_query(&q, &idx, &p, GapCosts::sequence_default(), &mut PrefilterScratch::new(), &mut st);
        assert_eq!(c.len(), 10);
        assert_eq!(st.candidates_capped as usize, a.len() - 10);
        assert_eq!(&c[..], &a[..10]);
        p.max_candidates = 0;
        p.max_seed_occ = 3;
        let mut st = PrefilterStats::default();
        prefilter_query(&q, &idx, &p, GapCosts::sequence_default(), &mut PrefilterScratch::new(), &mut st);
        assert!(st.seeds_over_occ_cap > 0 && st.postings_skipped_occ > 0);
    }

    #[test]
    fn query_mask_removes_low_complexity_seeds() {
        let m = SubstitutionMatrix::blosum62();
        let db = db_from(&format!(">lc\nMMMM{}MMMM\n", "Q".repeat(40)));
        let idx = KmerIndex::build(&db, 6);
        let qn = enc(&format!("WWWW{}WWWW", "Q".repeat(40)));
        let rows = rows_of(&SeqScorer { query: &qn, mat: &m });
        let mask = crate::mask::low_complexity_mask(&qn, MaskParams::default());
        let unmasked = PrefilterQuery { seed_rows: &rows, residues: Some(&qn), ext_rows: &rows, seed_mask: None };
        let masked = PrefilterQuery { seed_rows: &rows, residues: Some(&qn), ext_rows: &rows, seed_mask: Some(&mask) };
        let p = PrefilterParams::default();
        let g = GapCosts::sequence_default();
        let mut s1 = PrefilterStats::default();
        let mut s2 = PrefilterStats::default();
        assert!(!prefilter_query(&unmasked, &idx, &p, g, &mut PrefilterScratch::new(), &mut s1).is_empty());
        assert!(prefilter_query(&masked, &idx, &p, g, &mut PrefilterScratch::new(), &mut s2).is_empty());
        assert!(s2.windows_skipped > 30);
    }

    #[test]
    fn seed_pattern_parsing() {
        let p = SeedPattern::from_mask("1101011").unwrap();
        assert_eq!(p.weight(), 5);
        assert_eq!(p.span(), 7);
        assert_eq!(p.mask_string(), "1101011");
        assert!(SeedPattern::from_mask("0110").is_err());
        assert!(SeedPattern::from_mask("1100").is_err());
        assert!(SeedPattern::from_mask("11x1").is_err());
        assert!(SeedPattern::from_mask("").is_err());
        let c = SeedPattern::contiguous(6);
        assert_eq!(c.weight(), 6);
        assert_eq!(c.span(), 6);
    }

    #[test]
    fn spaced_seed_finds_what_contiguous_misses() {
        let m = SubstitutionMatrix::blosum62();
        let q = enc("MKVLLACDEFGHIK");
        let db = db_from(">tgt\nLLACDPF\n");
        let mut strict = PrefilterParams::default();
        strict.kmer_score = 30;
        strict.min_diag_score = 1;
        strict.exact_seeds = false;
        strict.k = 6;
        let idx_c = KmerIndex::build(&db, 6);
        let c_cont = prefilter_sequence(&q, &db, &idx_c, &m, strict);
        assert!(!c_cont.iter().any(|c| c.diag_score > 0), "contiguous seed should miss");
        let pattern = SeedPattern::from_mask("1111101").unwrap();
        strict.k = pattern.weight();
        let idx_s = KmerIndex::build_spaced(&db, pattern);
        let c_spaced = prefilter_sequence(&q, &db, &idx_s, &m, strict);
        assert!(c_spaced.iter().any(|c| c.diag_score > 0), "spaced seed should find the diverged target");
    }

    #[test]
    fn sensitivity_lowers_threshold_monotonically() {
        let mut prev = i32::MAX;
        for i in 1..=9 {
            let thr = kmer_score_for_sensitivity(i as f64, 2);
            assert!(thr <= prev);
            prev = thr;
        }
        assert!(kmer_score_for_sensitivity(5.7, 8) > kmer_score_for_sensitivity(5.7, 2));
    }

    #[test]
    fn presets_are_monotone_and_default_is_s57() {
        let mut prev = (i32::MAX, 0usize);
        for i in 0..=80 {
            let s = 1.0 + i as f64 * 0.1;
            let (k, c, g) = sequence_preset(s);
            assert!(k <= prev.0 && c >= prev.1, "s={} -> {:?}", s, (k, c, g));
            assert!(c % 32 == 0 && (GATE_MIN..=GATE_MAX).contains(&g));
            prev = (k, c);
        }
        for &(s, k, c, g) in SEQUENCE_PRESETS.iter() {
            assert_eq!(sequence_preset(s), (k, c, g));
        }
        let d = PrefilterParams::sequence_default();
        assert_eq!(sequence_preset(5.7), (d.kmer_score, d.max_kmers_per_pos, d.min_diag_score));
    }

    #[test]
    fn size_adjusted_gate_follows_database_size() {
        let r = GATE_REF_RESIDUES as usize;
        assert_eq!(size_adjusted_gate(42, r), 42);
        let mut prev = i32::MIN;
        for e in -14..=8 {
            let n = (GATE_REF_RESIDUES * (e as f64 / 2.0).exp()) as usize;
            let g = size_adjusted_gate(42, n);
            assert!(g >= prev, "gate must not fall as the database grows");
            assert!((GATE_MIN..=GATE_MAX).contains(&g));
            prev = g;
        }
        // two e-folds below / above the reference
        let e2 = std::f64::consts::E.powi(2);
        assert_eq!(size_adjusted_gate(42, (GATE_REF_RESIDUES / e2) as usize), (42.0 - 2.0 * GATE_SLOPE_BELOW).round() as i32);
        assert_eq!(size_adjusted_gate(42, (GATE_REF_RESIDUES * e2) as usize), (42.0 + 2.0 * GATE_SLOPE_ABOVE).round() as i32);
        // tiny and empty databases are clamped
        assert_eq!(size_adjusted_gate(42, 0), GATE_MIN);
        assert_eq!(size_adjusted_gate(42, 10), GATE_MIN);
        // the measured operating points (bench/README.md)
        assert_eq!(size_adjusted_gate(42, 2_611_087), 33);
        assert_eq!(size_adjusted_gate(42, 68_310_060), 41);
        assert_eq!(size_adjusted_gate(42, 265_123_580), 43);
    }

    #[test]
    fn defaults_are_shared_and_valid() {
        assert_eq!(PrefilterParams::default(), PrefilterParams::sequence_default());
        PrefilterParams::sequence_default().validate().unwrap();
        PrefilterParams::profile_default().validate().unwrap();
    }
}

#[cfg(test)]
mod microbench {
    use super::*;
    use crate::alphabet::build_aa2num;

    /// `cargo test --release microbench -- --ignored --nocapture`
    #[test]
    #[ignore]
    fn microbench_enumeration_vs_lookup() {
        let m = SubstitutionMatrix::blosum62();
        let path = std::env::var("ST_BENCH_INDEX").unwrap_or_default();
        if path.is_empty() {
            eprintln!("set ST_BENCH_INDEX=<file.stidx>");
            return;
        }
        let idx = SeedIndex::open(std::path::Path::new(&path)).unwrap();
        let t = build_aa2num();
        let q: Vec<u8> = "MADENKLKLGSGSFGEVFLVKHKESGNIYAMKKIRLEVEEEGVPSTAIREISLLKELKHENIVRLYDVVHSDKKLTLVFEYLDKDLKKYMDSCPEFAKDPRLIKMFLYQILRGLAYCHRRKVLHRDLKPQNLLIDRRTNSLKLADFGLARAFGVPLRTYTHEVVTLWYRAPEILLG".bytes().map(|b| t[b as usize]).collect();
        let rows = rows_of(&SeqScorer { query: &q, mat: &m });
        let sorted: Vec<SortedRow> = rows.iter().map(sorted_row).collect();
        let mut heap = BinaryHeap::new();
        let mut kw = KmerWork::default();
        let mut out = Vec::new();
        for (thr, cap) in [(25, 256), (18, 1024), (15, 4096)] {
            let mut codes = Vec::new();
            let mut n = 0usize;
            let t0 = std::time::Instant::now();
            for _ in 0..20 {
                for st in 0..q.len() - 6 {
                    top_kmers_heap(&sorted[st..st + 6], thr, cap, &mut heap, &mut out);
                    n += out.len();
                }
            }
            let th = t0.elapsed();
            let mut n2 = 0usize;
            let t0 = std::time::Instant::now();
            for _ in 0..20 {
                for st in 0..q.len() - 6 {
                    top_kmers(&sorted[st..st + 6], thr, cap, &mut kw, &mut out);
                    n2 += out.len();
                    codes.extend(out.iter().map(|x| x.0));
                }
            }
            let te = t0.elapsed();
            assert_eq!(n, n2);
            let t1 = std::time::Instant::now();
            let mut hits = 0usize;
            for &c in &codes {
                hits += idx.postings(c).len();
            }
            let tl = t1.elapsed();
            eprintln!(
                "thr {} cap {}: {} kmers, best-first {:.1} ns/kmer, depth-first+boundary {:.1} ns/kmer, lookup {:.1} ns/kmer ({} postings)",
                thr,
                cap,
                n,
                th.as_nanos() as f64 / n as f64,
                te.as_nanos() as f64 / n as f64,
                tl.as_nanos() as f64 / codes.len() as f64,
                hits
            );
        }
    }
}
