//! Local alignment (Smith–Waterman with Gotoh affine gaps) with an **exact
//! three-state traceback**.
//!
//! A single dynamic-programming core is driven by a pair-scoring closure so the
//! same logic serves sequence↔sequence (substitution matrix), profile↔sequence
//! (PSSM) and profile↔profile (column co-emission) alignment.
//!
//! ## Why three states
//!
//! With affine gaps the optimal path is a walk over three matrices: `H` (any
//! state), `E` (gap consuming target residues, CIGAR `I`) and `F` (gap consuming
//! query residues, CIGAR `D`). Storing a single "best direction" per `H` cell —
//! what v0.2 did — cannot tell a gap *extension* from a gap *opening*, so the
//! traceback could leave a gap early (`24M2I1M2I9M` for a true `24M4I...`) and
//! report coordinates, identities and gap counts of a path that does not have
//! the reported score. Here every cell records the `H` source **and** whether
//! `E`/`F` were extended or opened, and traceback walks the state machine, so
//! the returned CIGAR always rescores to exactly `raw_score`
//! ([`rescore_alignment`] checks this in tests).
//!
//! ## The canonical method ([`align_canonical`], used by every sequence and PSSM search)
//!
//! 1. vectorised score + end-cell pass ([`crate::simd`]); the end cell is the
//!    first optimum in column-major order (smallest target position, then
//!    smallest query position);
//! 2. vectorised pass over the reversed prefixes, whose own column-major-first
//!    optimum is the start of an optimal alignment ending at that cell;
//! 3. anchored traceback of the start/end rectangle inside a diagonal band
//!    (radius [`DEFAULT_BAND_RADIUS`]) that is doubled until it reproduces the
//!    optimal score — exact, `O(n + band·L)` memory.
//!
//! Inside the band `H` prefers diagonal, then `E`, then `F`, and `E`/`F` prefer
//! opening over extending when equal. Every step is deterministic and identical
//! for the SIMD and scalar kernels, so the reported alignment never depends on
//! the CPU, thread count, memory limit or sharding.
//!
//! Closure-driven alignment ([`align_pairs`], used for profile–profile
//! scoring, which has no residue table) keeps the older scheme: one full
//! traceback pass (1 byte per cell) up to [`FULL_TRACEBACK_MAX_CELLS`], else a
//! linear-memory score/end pass, an anchored reverse start search and the same
//! widening banded traceback ([`align_exact_banded`]); its end cell is the
//! first maximum in row-major order.

use crate::alphabet::{ALPHABET_SIZE, PROFILE_AA_SIZE};
use crate::matrix::SubstitutionMatrix;
use crate::profile::Profile;
use crate::simd::{QueryProfile, SimdScratch};

/// Abstraction over "the score of placing this target residue against query
/// column `qpos`". Implemented by a plain sequence, a PSSM, and [`ScoreTable`].
pub trait Scorer {
    /// Number of query columns.
    fn query_len(&self) -> usize;
    /// Score of aligning `target_residue` (internal index) at query column `qpos`.
    fn score(&self, qpos: usize, target_residue: u8) -> i32;
    /// Query residue at a column (for identity accounting), if meaningful.
    fn query_residue(&self, qpos: usize) -> u8;
}

/// Sequence-vs-sequence scorer backed by the integer substitution matrix.
pub struct SeqScorer<'a> {
    pub query: &'a [u8],
    pub mat: &'a SubstitutionMatrix,
}

impl<'a> Scorer for SeqScorer<'a> {
    #[inline]
    fn query_len(&self) -> usize {
        self.query.len()
    }
    #[inline]
    fn score(&self, qpos: usize, target_residue: u8) -> i32 {
        self.mat.score(self.query[qpos], target_residue)
    }
    #[inline]
    fn query_residue(&self, qpos: usize) -> u8 {
        self.query[qpos]
    }
}

/// Profile-vs-sequence scorer backed by a PSSM.
pub struct ProfileScorer<'a> {
    pub profile: &'a Profile,
}

impl<'a> Scorer for ProfileScorer<'a> {
    #[inline]
    fn query_len(&self) -> usize {
        self.profile.query_len
    }
    #[inline]
    fn score(&self, qpos: usize, target_residue: u8) -> i32 {
        self.profile.score(qpos, target_residue)
    }
    #[inline]
    fn query_residue(&self, qpos: usize) -> u8 {
        self.profile.query_num[qpos]
    }
}

/// Precomputed per-column score rows (`rows[qpos][residue]`, residues 0..=20),
/// optionally with a per-column additive bias (composition-bias correction).
/// Target residues outside the alphabet score 0.
#[derive(Clone, Debug)]
pub struct ScoreTable {
    rows: Vec<[i32; ALPHABET_SIZE]>,
    qres: Vec<u8>,
}

impl ScoreTable {
    /// Tabulate any scorer.
    pub fn from_scorer<S: Scorer>(s: &S) -> ScoreTable {
        Self::with_bias(s, None)
    }

    /// Tabulate a scorer, adding `bias[qpos]` to every residue score of column
    /// `qpos`.
    pub fn with_bias<S: Scorer>(s: &S, bias: Option<&[i32]>) -> ScoreTable {
        let m = s.query_len();
        let mut rows = vec![[0i32; ALPHABET_SIZE]; m];
        let mut qres = Vec::with_capacity(m);
        for (c, row) in rows.iter_mut().enumerate() {
            let b = bias.map(|b| b[c]).unwrap_or(0);
            for (r, slot) in row.iter_mut().enumerate() {
                *slot = s.score(c, r as u8) + b;
            }
            qres.push(s.query_residue(c));
        }
        ScoreTable { rows, qres }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.rows.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    /// Score row of query column `qpos`.
    #[inline]
    pub fn row(&self, qpos: usize) -> &[i32; ALPHABET_SIZE] {
        &self.rows[qpos]
    }

    /// All rows.
    pub fn rows(&self) -> &[[i32; ALPHABET_SIZE]] {
        &self.rows
    }

    /// Replace the rows of the given query columns (hard masking).
    pub fn set_rows(&mut self, cols: &[usize], row: &[i32; ALPHABET_SIZE]) {
        for &c in cols {
            self.rows[c] = *row;
        }
    }

    /// Largest and smallest score in the table.
    pub fn score_range(&self) -> (i32, i32) {
        let mut lo = i32::MAX;
        let mut hi = i32::MIN;
        for r in &self.rows {
            for &v in r.iter() {
                lo = lo.min(v);
                hi = hi.max(v);
            }
        }
        (lo, hi)
    }
}

impl Scorer for ScoreTable {
    #[inline]
    fn query_len(&self) -> usize {
        self.rows.len()
    }
    #[inline]
    fn score(&self, qpos: usize, target_residue: u8) -> i32 {
        let r = target_residue as usize;
        if r < ALPHABET_SIZE {
            self.rows[qpos][r]
        } else {
            0
        }
    }
    #[inline]
    fn query_residue(&self, qpos: usize) -> u8 {
        self.qres[qpos]
    }
}

/// Affine gap penalties (expressed in the same score units as the scorer).
/// A gap of length `L` costs `open + (L-1)·extend`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "distributed", derive(serde::Serialize, serde::Deserialize))]
pub struct GapCosts {
    /// Cost to open a gap (includes the first gap position).
    pub open: i32,
    /// Cost to extend a gap by one position.
    pub extend: i32,
}

impl GapCosts {
    /// Sequence-mode defaults (half-bit matrix scale): 11 / 1, blastp-like.
    pub fn sequence_default() -> Self {
        GapCosts { open: 11, extend: 1 }
    }
    /// Profile-mode defaults, scaled to the ×8 PSSM units (≈ 4× sequence).
    pub fn profile_default() -> Self {
        GapCosts { open: 44, extend: 4 }
    }
    /// Cost of a gap run of `len` positions.
    #[inline]
    pub fn cost(&self, len: usize) -> i32 {
        if len == 0 {
            0
        } else {
            self.open + (len as i32 - 1) * self.extend
        }
    }
}

/// The result of a local alignment. All coordinates are 0-based and inclusive.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Alignment {
    pub raw_score: i32,
    pub query_start: usize,
    pub query_end: usize,
    pub target_start: usize,
    pub target_end: usize,
    /// Alignment columns (aligned pairs + gap columns).
    pub aln_len: usize,
    /// Aligned pairs with identical (real) residues.
    pub identities: usize,
    /// Aligned pairs with different residues.
    pub mismatches: usize,
    /// Gap columns.
    pub gaps: usize,
    /// Gap runs (BLAST `gapopen`).
    pub gap_opens: usize,
    /// Run-length CIGAR over `M` (aligned pair), `I` (target residue against a
    /// query gap) and `D` (query residue against a target gap).
    pub cigar: String,
}

impl Alignment {
    /// Percent identity over alignment columns (BLAST `pident`).
    pub fn pct_identity(&self) -> f64 {
        if self.aln_len == 0 {
            0.0
        } else {
            100.0 * self.identities as f64 / self.aln_len as f64
        }
    }

    /// Fraction of the query covered by the alignment.
    pub fn query_cov(&self, qlen: usize) -> f64 {
        (self.query_end + 1 - self.query_start) as f64 / qlen.max(1) as f64
    }

    /// Fraction of the target covered by the alignment.
    pub fn target_cov(&self, tlen: usize) -> f64 {
        (self.target_end + 1 - self.target_start) as f64 / tlen.max(1) as f64
    }

    /// Shift target coordinates by `offset` (alignment computed on a window
    /// starting at target position `offset`).
    pub fn shift_target(&mut self, offset: usize) {
        self.target_start += offset;
        self.target_end += offset;
    }
}

// ---------------------------------------------------------------------------
// traceback encoding
// ---------------------------------------------------------------------------

const SRC_STOP: u8 = 0;
const SRC_DIAG: u8 = 1;
const SRC_E: u8 = 2; // gap consuming target (CIGAR I)
const SRC_F: u8 = 3; // gap consuming query  (CIGAR D)
const E_EXT: u8 = 4; // E[i][j] extended from E[i][j-1]
const F_EXT: u8 = 8; // F[i][j] extended from F[i-1][j]

const NEG: i32 = i32::MIN / 4;

/// Reusable working memory for alignment (one per worker thread).
pub struct AlignScratch {
    tb: Vec<u8>,
    h_prev: Vec<i32>,
    h_cur: Vec<i32>,
    f_col: Vec<i32>,
    ops: Vec<u8>,
    /// striped-kernel working memory for the score/end passes
    pub simd: SimdScratch,
    /// reversed target prefix (start pass)
    rev_t: Vec<u8>,
    /// Problems with more cells than this use the banded exact path instead of
    /// one full traceback matrix (bounds traceback memory per worker).
    pub full_tb_limit: usize,
}

impl Default for AlignScratch {
    fn default() -> Self {
        AlignScratch {
            tb: Vec::new(),
            h_prev: Vec::new(),
            h_cur: Vec::new(),
            f_col: Vec::new(),
            ops: Vec::new(),
            simd: SimdScratch::new(),
            rev_t: Vec::new(),
            full_tb_limit: FULL_TRACEBACK_MAX_CELLS,
        }
    }
}

impl AlignScratch {
    pub fn new() -> Self {
        Self::default()
    }
    /// Scratch with a custom full-traceback cell limit.
    pub fn with_limit(cells: usize) -> Self {
        AlignScratch {
            full_tb_limit: cells.max(1024),
            ..Self::default()
        }
    }
    /// Bytes currently reserved (memory accounting).
    pub fn heap_bytes(&self) -> usize {
        self.tb.capacity() + (self.h_prev.capacity() + self.h_cur.capacity() + self.f_col.capacity()) * 4
            + self.ops.capacity()
            + self.rev_t.capacity()
    }
}

/// Above this many cells the one-pass full traceback is replaced by the
/// score → anchored start → banded traceback scheme (4 Mi cells = 4 MiB of
/// traceback per worker). The choice depends only on the problem size, never
/// on memory settings, so every execution mode reports the same alignment
/// (the two methods can pick different co-optimal paths).
pub const FULL_TRACEBACK_MAX_CELLS: usize = 1 << 22;

/// Build an [`Alignment`] from reversed ops collected during traceback.
fn finish_alignment<I: Fn(usize, usize) -> bool>(
    ops_rev: &mut [u8],
    score: i32,
    q_start: usize,
    t_start: usize,
    identical: I,
) -> Alignment {
    ops_rev.reverse();
    let ops = &ops_rev[..];
    let mut i = q_start;
    let mut j = t_start;
    let mut identities = 0;
    let mut mismatches = 0;
    let mut gaps = 0;
    let mut gap_opens = 0;
    let mut prev = 0u8;
    for &op in ops {
        match op {
            b'M' => {
                if identical(i, j) {
                    identities += 1;
                } else {
                    mismatches += 1;
                }
                i += 1;
                j += 1;
            }
            b'I' => {
                gaps += 1;
                if prev != b'I' {
                    gap_opens += 1;
                }
                j += 1;
            }
            _ => {
                gaps += 1;
                if prev != b'D' {
                    gap_opens += 1;
                }
                i += 1;
            }
        }
        prev = op;
    }
    Alignment {
        raw_score: score,
        query_start: q_start,
        query_end: i - 1,
        target_start: t_start,
        target_end: j - 1,
        aln_len: ops.len(),
        identities,
        mismatches,
        gaps,
        gap_opens,
        cigar: run_length_encode(ops),
    }
}

/// One-pass local DP with full three-state traceback over an `m × n` problem.
/// Returns `None` if no positive-scoring alignment exists.
fn local_full_tb<P, I>(
    m: usize,
    n: usize,
    pair: P,
    identical: I,
    gaps: GapCosts,
    band: Option<(i64, i64)>,
    s: &mut AlignScratch,
) -> Option<Alignment>
where
    P: Fn(usize, usize) -> i32,
    I: Fn(usize, usize) -> bool,
{
    if m == 0 || n == 0 {
        return None;
    }
    let w = n + 1;
    s.tb.clear();
    s.tb.resize((m + 1) * w, SRC_STOP);
    s.h_prev.clear();
    s.h_prev.resize(w, 0);
    s.h_cur.clear();
    s.h_cur.resize(w, 0);
    s.f_col.clear();
    s.f_col.resize(w, NEG);
    let (mut best, mut bi, mut bj) = (0i32, 0usize, 0usize);
    for i in 1..=m {
        let mut e = NEG;
        s.h_cur[0] = 0;
        let row = i * w;
        for j in 1..=n {
            if let Some((center, radius)) = band {
                if (j as i64 - i as i64 - center).abs() > radius {
                    // outside the band: unreachable for every state
                    s.h_cur[j] = 0;
                    s.f_col[j] = NEG;
                    e = NEG;
                    s.tb[row + j] = SRC_STOP;
                    continue;
                }
            }
            let diag = s.h_prev[j - 1] + pair(i - 1, j - 1);
            // E: gap consuming target (horizontal)
            let e_open = s.h_cur[j - 1] - gaps.open;
            let e_ext = e - gaps.extend;
            let mut bits = 0u8;
            if e_ext > e_open {
                e = e_ext;
                bits |= E_EXT;
            } else {
                e = e_open;
            }
            // F: gap consuming query (vertical)
            let f_open = s.h_prev[j] - gaps.open;
            let f_ext = s.f_col[j] - gaps.extend;
            let f = if f_ext > f_open {
                bits |= F_EXT;
                f_ext
            } else {
                f_open
            };
            s.f_col[j] = f;
            let mut v = diag;
            let mut src = SRC_DIAG;
            if e > v {
                v = e;
                src = SRC_E;
            }
            if f > v {
                v = f;
                src = SRC_F;
            }
            if v <= 0 {
                v = 0;
                src = SRC_STOP;
            }
            s.h_cur[j] = v;
            s.tb[row + j] = bits | src;
            if v > best {
                best = v;
                bi = i;
                bj = j;
            }
        }
        std::mem::swap(&mut s.h_prev, &mut s.h_cur);
    }
    if best <= 0 {
        return None;
    }
    // walk the state machine
    s.ops.clear();
    let (mut i, mut j) = (bi, bj);
    let mut state = 0u8; // 0 = H, 1 = E, 2 = F
    loop {
        let cell = s.tb[i * w + j];
        match state {
            0 => match cell & 3 {
                SRC_DIAG => {
                    s.ops.push(b'M');
                    i -= 1;
                    j -= 1;
                }
                SRC_E => state = 1,
                SRC_F => state = 2,
                _ => break,
            },
            1 => {
                s.ops.push(b'I');
                if cell & E_EXT == 0 {
                    state = 0;
                }
                j -= 1;
            }
            _ => {
                s.ops.push(b'D');
                if cell & F_EXT == 0 {
                    state = 0;
                }
                i -= 1;
            }
        }
        if i == 0 || j == 0 {
            // a path can only reach the border through H = 0 (stop)
            break;
        }
    }
    let mut ops = std::mem::take(&mut s.ops);
    let aln = finish_alignment(&mut ops, best, i, j, identical);
    s.ops = ops;
    Some(aln)
}

/// Linear-memory local score with the (first, row-major) end cell.
/// Returns `(score, i_end, j_end)` with 0-based inclusive end coordinates;
/// score 0 means no positive alignment.
pub fn local_score_end<P>(m: usize, n: usize, pair: P, gaps: GapCosts, s: &mut AlignScratch) -> (i32, usize, usize)
where
    P: Fn(usize, usize) -> i32,
{
    if m == 0 || n == 0 {
        return (0, 0, 0);
    }
    let w = n + 1;
    s.h_prev.clear();
    s.h_prev.resize(w, 0);
    s.h_cur.clear();
    s.h_cur.resize(w, 0);
    s.f_col.clear();
    s.f_col.resize(w, NEG);
    let (mut best, mut bi, mut bj) = (0i32, 0usize, 0usize);
    for i in 1..=m {
        let mut e = NEG;
        s.h_cur[0] = 0;
        for j in 1..=n {
            let diag = s.h_prev[j - 1] + pair(i - 1, j - 1);
            e = (s.h_cur[j - 1] - gaps.open).max(e - gaps.extend);
            let f = (s.h_prev[j] - gaps.open).max(s.f_col[j] - gaps.extend);
            s.f_col[j] = f;
            let v = diag.max(e).max(f).max(0);
            s.h_cur[j] = v;
            if v > best {
                best = v;
                bi = i;
                bj = j;
            }
        }
        std::mem::swap(&mut s.h_prev, &mut s.h_cur);
    }
    if best <= 0 {
        (0, 0, 0)
    } else {
        (best, bi - 1, bj - 1)
    }
}

/// Find the start of the optimal local alignment that ends at `(ie, je)` with
/// score `score` (the first row-major maximum), by an anchored DP over the
/// reversed prefixes: the reversed path must begin at `(ie, je)` and is never
/// restarted. Returns 0-based `(is, js)`. Linear memory.
///
/// Pruning is exact: for the first row-major optimum every proper suffix of the
/// forward alignment scores > 0 (otherwise a prefix would score ≥ `score` and end
/// earlier in row-major order), so reversed cells with value ≤ 0 cannot lie on it.
fn anchored_start<P>(ie: usize, je: usize, pair: P, gaps: GapCosts, score: i32, s: &mut AlignScratch) -> Option<(usize, usize)>
where
    P: Fn(usize, usize) -> i32,
{
    let m = ie + 1;
    let n = je + 1;
    let w = n + 1;
    s.h_prev.clear();
    s.h_prev.resize(w, NEG);
    s.h_cur.clear();
    s.h_cur.resize(w, NEG);
    s.f_col.clear();
    s.f_col.resize(w, NEG);
    // reversed cell (a, b), 1-based ↔ query ie+1-a, target je+1-b (0-based)
    s.h_prev[0] = 0; // the anchor corner is the only start
    for a in 1..=m {
        let mut e = NEG;
        s.h_cur[0] = NEG;
        let qi = ie + 1 - a;
        let mut alive = false;
        for b in 1..=n {
            let tj = je + 1 - b;
            let d = s.h_prev[b - 1];
            let diag = if d > NEG / 2 { d + pair(qi, tj) } else { NEG };
            e = (s.h_cur[b - 1] - gaps.open).max(e - gaps.extend);
            let f = (s.h_prev[b] - gaps.open).max(s.f_col[b] - gaps.extend);
            // the forward alignment's first column is an aligned pair
            if diag == score {
                return Some((qi, tj));
            }
            let mut v = diag.max(e).max(f);
            if v <= 0 {
                v = NEG;
            } else {
                alive = true;
            }
            s.h_cur[b] = v;
            s.f_col[b] = if f <= 0 { NEG } else { f };
            if e <= 0 {
                e = NEG;
            }
        }
        if !alive {
            break;
        }
        std::mem::swap(&mut s.h_prev, &mut s.h_cur);
    }
    None
}

/// Anchored (global within the rectangle) traceback restricted to a diagonal
/// band. The rectangle is query `qs..=qe` × target `ts..=te`; the path starts
/// with the aligned pair `(qs, ts)` and ends at `(qe, te)`. Returns the
/// alignment if the banded optimum equals `score` (which proves it optimal),
/// otherwise `None` (the caller widens the band).
#[allow(clippy::too_many_arguments)]
fn anchored_banded_tb<P, I>(
    qs: usize,
    qe: usize,
    ts: usize,
    te: usize,
    radius: usize,
    pair: &P,
    identical: &I,
    gaps: GapCosts,
    score: i32,
    s: &mut AlignScratch,
) -> Option<Alignment>
where
    P: Fn(usize, usize) -> i32,
    I: Fn(usize, usize) -> bool,
{
    let m = qe + 1 - qs; // rows a = 1..=m
    let n = te + 1 - ts; // cols b = 1..=n
    // band over local diagonals d = b - a covering the segment (0,0)-(m,n)
    let dmin = (0i64).min(n as i64 - m as i64) - radius as i64;
    let dmax = (0i64).max(n as i64 - m as i64) + radius as i64;
    let bw = (dmax - dmin + 1) as usize;
    s.tb.clear();
    s.tb.resize((m + 1) * bw, SRC_STOP);
    let w = n + 1;
    s.h_prev.clear();
    s.h_prev.resize(w, NEG);
    s.h_cur.clear();
    s.h_cur.resize(w, NEG);
    s.f_col.clear();
    s.f_col.resize(w, NEG);
    s.h_prev[0] = 0; // H(0,0): the only start; every other border cell is NEG
    let band_of = |a: usize| -> (usize, usize) {
        let lo = (a as i64 + dmin).max(1) as usize;
        let hi = (a as i64 + dmax).min(n as i64) as usize;
        (lo, hi)
    };
    let mut prev_hi = 0usize; // row 0 holds only column 0
    for a in 1..=m {
        let (lo, hi) = band_of(a);
        // columns entering the band were not computed in row a-1: unreachable
        for b in (prev_hi + 1)..=hi {
            s.h_prev[b] = NEG;
            s.f_col[b] = NEG;
        }
        s.h_cur[lo - 1] = NEG;
        let mut e = NEG;
        let qi = qs + a - 1;
        let row = a * bw;
        for b in lo..=hi {
            let d = s.h_prev[b - 1];
            let diag = if d > NEG / 2 { d + pair(qi, ts + b - 1) } else { NEG };
            let e_open = s.h_cur[b - 1] - gaps.open;
            let e_ext = e - gaps.extend;
            let mut bits = 0u8;
            if e_ext > e_open {
                e = e_ext;
                bits |= E_EXT;
            } else {
                e = e_open;
            }
            let f_open = s.h_prev[b] - gaps.open;
            let f_ext = s.f_col[b] - gaps.extend;
            let f = if f_ext > f_open {
                bits |= F_EXT;
                f_ext
            } else {
                f_open
            };
            s.f_col[b] = f;
            let mut v = diag;
            let mut src = SRC_DIAG;
            if e > v {
                v = e;
                src = SRC_E;
            }
            if f > v {
                v = f;
                src = SRC_F;
            }
            if v <= NEG / 2 {
                v = NEG;
                src = SRC_STOP;
            }
            s.h_cur[b] = v;
            s.tb[row + (b as i64 - a as i64 - dmin) as usize] = bits | src;
        }
        prev_hi = hi;
        std::mem::swap(&mut s.h_prev, &mut s.h_cur);
    }
    if s.h_prev[n] != score {
        return None;
    }
    s.ops.clear();
    let (mut a, mut b) = (m, n);
    let mut state = 0u8;
    while !(a == 0 && b == 0) {
        if a == 0 || b == 0 {
            return None;
        }
        let off = b as i64 - a as i64 - dmin;
        if off < 0 || off as usize >= bw {
            return None;
        }
        let cell = s.tb[a * bw + off as usize];
        match state {
            0 => match cell & 3 {
                SRC_DIAG => {
                    s.ops.push(b'M');
                    a -= 1;
                    b -= 1;
                }
                SRC_E => state = 1,
                SRC_F => state = 2,
                _ => return None,
            },
            1 => {
                s.ops.push(b'I');
                if cell & E_EXT == 0 {
                    state = 0;
                }
                b -= 1;
            }
            _ => {
                s.ops.push(b'D');
                if cell & F_EXT == 0 {
                    state = 0;
                }
                a -= 1;
            }
        }
    }
    if state != 0 {
        return None;
    }
    let mut ops = std::mem::take(&mut s.ops);
    let aln = finish_alignment(&mut ops, score, qs, ts, identical);
    s.ops = ops;
    Some(aln)
}

/// Exact optimal local alignment using linear-memory scoring plus a banded
/// traceback that starts at `band_radius` and doubles until the banded optimum
/// equals the true optimum (it always terminates: a band covering the
/// rectangle is the full problem).
pub fn align_exact_banded<P, I>(m: usize, n: usize, pair: P, identical: I, gaps: GapCosts, band_radius: usize, s: &mut AlignScratch) -> Option<Alignment>
where
    P: Fn(usize, usize) -> i32,
    I: Fn(usize, usize) -> bool,
{
    let (score, ie, je) = local_score_end(m, n, &pair, gaps, s);
    if score <= 0 {
        return None;
    }
    let (is, js) = anchored_start(ie, je, &pair, gaps, score, s)?;
    let rows = ie + 1 - is;
    let cols = je + 1 - js;
    let mut r = band_radius.max(1);
    loop {
        if let Some(a) = anchored_banded_tb(is, ie, js, je, r, &pair, &identical, gaps, score, s) {
            return Some(a);
        }
        if r >= rows + cols {
            // a band this wide covers the rectangle: cannot fail
            debug_assert!(false, "full-width banded traceback failed");
            return None;
        }
        r = (r * 2).min(rows + cols);
    }
}

/// Initial band radius of the canonical traceback (doubled until exact).
pub const DEFAULT_BAND_RADIUS: usize = 16;

/// **Canonical exact local alignment** of a tabulated query against a numeric
/// target — what every sequence/PSSM search reports.
///
/// 1. Forward striped pass: optimal score `S` and the end cell `(ie, je)`, the
///    first optimum in column-major order ([`QueryProfile::sw_score_end_with`]).
/// 2. Reverse striped pass over the reversed prefixes `q[..=ie]`, `t[..=je]`.
///    Every alignment of score `S` inside that rectangle ends at `(ie, je)`
///    (any other end cell would be dominated, i.e. an earlier optimum), so the
///    reversed optimum is `S` and its column-major-first end is a start
///    `(is, js)` of an optimal alignment ending at `(ie, je)` — the one with the
///    latest target start, then the latest query start.
/// 3. Anchored traceback of the rectangle `[is, ie] × [js, je]` inside a
///    diagonal band, widened from `band_radius` until it reproduces `S`
///    (which proves the path optimal).
///
/// Time is two vectorised score passes plus `O(L · band)` scalar traceback
/// cells; memory is `O(n + L · band)`. `qp` must be built from `table`'s rows.
/// Ties are broken deterministically, independent of kernel (SIMD or scalar
/// fallback) and of memory settings.
pub fn align_canonical(
    qp: &QueryProfile,
    table: &ScoreTable,
    target: &[u8],
    gaps: GapCosts,
    band_radius: usize,
    s: &mut AlignScratch,
) -> Option<Alignment> {
    let (score, ie, je) = qp.sw_score_end_with(target, gaps, &mut s.simd)?;
    // reverse pass for the start
    let rows = table.rows();
    let rev_rows: Vec<[i32; ALPHABET_SIZE]> = rows[..=ie].iter().rev().copied().collect();
    // the end-tracking pass runs on the i16 kernels only: skip the u8 profile
    let qpr = QueryProfile::from_rows_opts(rev_rows, crate::simd::kernel_choice(), false);
    s.rev_t.clear();
    s.rev_t.extend(target[..=je].iter().rev());
    let rev_t = std::mem::take(&mut s.rev_t);
    let rev = qpr.sw_score_end_with(&rev_t, gaps, &mut s.simd);
    s.rev_t = rev_t;
    let (is, js) = match rev {
        Some((rs, a, b)) if rs == score => (ie - a, je - b),
        _ => {
            // unreachable by construction; stay exact if it ever happens
            debug_assert!(false, "reverse pass score {:?} != forward {}", rev, score);
            anchored_start(ie, je, |i, j| table.score(i, target[j]), gaps, score, s)?
        }
    };
    let pair = |i: usize, j: usize| table.score(i, target[j]);
    let identical = |i: usize, j: usize| is_identity(table.query_residue(i), target[j]);
    let rows_n = ie + 1 - is;
    let cols_n = je + 1 - js;
    let mut r = band_radius.max(1);
    loop {
        if let Some(a) = anchored_banded_tb(is, ie, js, je, r, &pair, &identical, gaps, score, s) {
            return Some(a);
        }
        if r >= rows_n + cols_n {
            debug_assert!(false, "full-width banded traceback failed");
            return None;
        }
        r = (r * 2).min(rows_n + cols_n);
    }
}

/// Generic exact local alignment over a pair-scoring closure. Uses one full
/// traceback pass when `m·n` is small, otherwise [`align_exact_banded`].
pub fn align_pairs<P, I>(m: usize, n: usize, pair: P, identical: I, gaps: GapCosts, s: &mut AlignScratch) -> Option<Alignment>
where
    P: Fn(usize, usize) -> i32,
    I: Fn(usize, usize) -> bool,
{
    if m == 0 || n == 0 {
        return None;
    }
    if (m + 1).saturating_mul(n + 1) <= s.full_tb_limit {
        local_full_tb(m, n, pair, identical, gaps, None, s)
    } else {
        align_exact_banded(m, n, pair, identical, gaps, 32, s)
    }
}

#[inline]
fn is_identity(q: u8, t: u8) -> bool {
    q == t && (q as usize) < PROFILE_AA_SIZE
}

/// Optimal Smith–Waterman–Gotoh local alignment of a scorer's query against a
/// numeric target sequence. Returns `None` if the best score is non-positive.
pub fn local_align<S: Scorer>(scorer: &S, target: &[u8], gaps: GapCosts) -> Option<Alignment> {
    let mut s = AlignScratch::new();
    local_align_with(scorer, target, gaps, &mut s)
}

/// [`local_align`] with caller-provided scratch (no per-call allocation).
pub fn local_align_with<S: Scorer>(scorer: &S, target: &[u8], gaps: GapCosts, s: &mut AlignScratch) -> Option<Alignment> {
    if scorer.query_len() == 0 || target.is_empty() {
        return None;
    }
    let table = ScoreTable::from_scorer(scorer);
    let qp = QueryProfile::from_rows(table.rows().to_vec());
    align_canonical(&qp, &table, target, gaps, DEFAULT_BAND_RADIUS, s)
}

/// **Exact** optimal local alignment with `O(n + band·L)` traceback memory: the
/// band (initial radius `band_radius`) is widened automatically until the
/// traceback reproduces the optimal score, so the result always equals
/// [`local_align`]'s score. `_center_diag` (the seed diagonal) is accepted for
/// API compatibility; the band is placed on the alignment's own start/end
/// rectangle, which is tighter.
pub fn local_align_banded<S: Scorer>(scorer: &S, target: &[u8], gaps: GapCosts, _center_diag: i64, band_radius: usize) -> Option<Alignment> {
    let mut s = AlignScratch::new();
    local_align_banded_with(scorer, target, gaps, band_radius, &mut s)
}

/// [`local_align_banded`] with caller-provided scratch.
pub fn local_align_banded_with<S: Scorer>(scorer: &S, target: &[u8], gaps: GapCosts, band_radius: usize, s: &mut AlignScratch) -> Option<Alignment> {
    align_exact_banded(
        scorer.query_len(),
        target.len(),
        |i, j| scorer.score(i, target[j]),
        |i, j| is_identity(scorer.query_residue(i), target[j]),
        gaps,
        band_radius,
        s,
    )
}

/// **Heuristic** band-restricted local alignment (the v0.2 `--band` behaviour):
/// only cells within `band_radius` of `center_diag` (`target_pos − query_pos`)
/// are evaluated, so time and memory are `O(m·w)`. The result is the best
/// alignment *inside the band*, which can score below the unrestricted optimum
/// when indels drift the path out of the band. Callers must derive E-values
/// from the returned `raw_score` (never from an unrestricted score). Search
/// does not use this function.
pub fn local_align_in_band<S: Scorer>(scorer: &S, target: &[u8], gaps: GapCosts, center_diag: i64, band_radius: usize) -> Option<Alignment> {
    let m = scorer.query_len();
    let n = target.len();
    if m == 0 || n == 0 {
        return None;
    }
    let mut s = AlignScratch::new();
    local_full_tb(
        m,
        n,
        |i, j| scorer.score(i, target[j]),
        |i, j| is_identity(scorer.query_residue(i), target[j]),
        gaps,
        Some((center_diag, band_radius as i64)),
        &mut s,
    )
}

/// Column co-emission score between column `i` of profile `a` and column `j` of
/// profile `b`: `round(8 · log2(Σ_x P_a[i][x]·P_b[j][x] / f_x))`. It is ~0 when
/// both columns match the background, large-positive when both concentrate on the
/// same residue(s), and negative when they disagree — the HHsearch-style
/// profile column score, on the ×8 scale so the profile gap costs (44/4) apply.
#[inline]
fn pp_col_score(a: &Profile, i: usize, b: &Profile, j: usize, pback: &[f64]) -> i32 {
    let pa = &a.prob[i * PROFILE_AA_SIZE..(i + 1) * PROFILE_AA_SIZE];
    let pb = &b.prob[j * PROFILE_AA_SIZE..(j + 1) * PROFILE_AA_SIZE];
    let mut sum = 0f64;
    for x in 0..PROFILE_AA_SIZE {
        let fx = pback[x].max(1e-6);
        sum += (pa[x] as f64) * (pb[x] as f64) / fx;
    }
    (8.0 * sum.max(1e-6).log2()).round().clamp(-128.0, 127.0) as i32
}

/// Local Smith–Waterman–Gotoh alignment of two **profiles**, scoring column pairs
/// by co-emission ([`pp_col_score`]). Reported identity is the fraction of
/// aligned columns whose consensus residues agree.
pub fn align_profile_profile(a: &Profile, b: &Profile, gaps: GapCosts, pback: &[f64]) -> Option<Alignment> {
    let m = a.query_len;
    let n = b.query_len;
    // precompute the column-pair score matrix once (m·n small for profiles)
    let mut tab = vec![0i32; m * n];
    for i in 0..m {
        for j in 0..n {
            tab[i * n + j] = pp_col_score(a, i, b, j, pback);
        }
    }
    let mut s = AlignScratch::new();
    align_pairs(
        m,
        n,
        |i, j| tab[i * n + j],
        |i, j| a.consensus[i] == b.consensus[j],
        gaps,
        &mut s,
    )
}

/// Compress a per-column op string (`MMMIID...`) into CIGAR (`3M2I1D...`).
fn run_length_encode(ops: &[u8]) -> String {
    let mut out = String::new();
    let mut i = 0;
    while i < ops.len() {
        let c = ops[i];
        let mut run = 1;
        while i + run < ops.len() && ops[i + run] == c {
            run += 1;
        }
        out.push_str(&run.to_string());
        out.push(c as char);
        i += run;
    }
    out
}

/// Parse a CIGAR into `(op, count)` runs.
pub fn cigar_runs(cigar: &str) -> Result<Vec<(u8, usize)>, String> {
    let mut out = Vec::new();
    let mut num = 0usize;
    let mut have = false;
    for ch in cigar.bytes() {
        if ch.is_ascii_digit() {
            num = num
                .checked_mul(10)
                .and_then(|x| x.checked_add((ch - b'0') as usize))
                .ok_or("CIGAR count overflow")?;
            have = true;
        } else {
            if !have || num == 0 {
                return Err(format!("CIGAR op '{}' without a positive count", ch as char));
            }
            if !matches!(ch, b'M' | b'I' | b'D') {
                return Err(format!("unknown CIGAR op '{}'", ch as char));
            }
            out.push((ch, num));
            num = 0;
            have = false;
        }
    }
    if have {
        return Err("trailing count without op".into());
    }
    Ok(out)
}

/// Independently verify an alignment against its sequences: walk the CIGAR from
/// the start coordinates, recompute the score with `scorer` and `gaps`, and
/// check end coordinates, column counts, identities, mismatches, gap columns
/// and gap opens. Returns the recomputed score on success.
pub fn rescore_alignment<S: Scorer>(scorer: &S, target: &[u8], gaps: GapCosts, aln: &Alignment) -> Result<i32, String> {
    let runs = cigar_runs(&aln.cigar)?;
    let (mut i, mut j) = (aln.query_start, aln.target_start);
    let (mut score, mut cols, mut ident, mut mism, mut gapc, mut opens) = (0i32, 0, 0, 0, 0, 0);
    let mut prev = 0u8;
    for &(op, n) in &runs {
        if op == prev {
            return Err("adjacent CIGAR runs with the same op".into());
        }
        match op {
            b'M' => {
                for _ in 0..n {
                    if i >= scorer.query_len() || j >= target.len() {
                        return Err(format!("M out of bounds at q{} t{}", i, j));
                    }
                    score += scorer.score(i, target[j]);
                    if is_identity(scorer.query_residue(i), target[j]) {
                        ident += 1;
                    } else {
                        mism += 1;
                    }
                    i += 1;
                    j += 1;
                }
            }
            b'I' => {
                if j + n > target.len() {
                    return Err("I out of bounds".into());
                }
                score -= gaps.cost(n);
                j += n;
                gapc += n;
                opens += 1;
            }
            _ => {
                if i + n > scorer.query_len() {
                    return Err("D out of bounds".into());
                }
                score -= gaps.cost(n);
                i += n;
                gapc += n;
                opens += 1;
            }
        }
        cols += n;
        prev = op;
    }
    if runs.first().map(|r| r.0) != Some(b'M') || runs.last().map(|r| r.0) != Some(b'M') {
        return Err("local alignment must begin and end with an aligned pair".into());
    }
    let checks = [
        (i, aln.query_end + 1, "query_end"),
        (j, aln.target_end + 1, "target_end"),
        (cols, aln.aln_len, "aln_len"),
        (ident, aln.identities, "identities"),
        (mism, aln.mismatches, "mismatches"),
        (gapc, aln.gaps, "gaps"),
        (opens, aln.gap_opens, "gap_opens"),
    ];
    for (got, want, what) in checks {
        if got != want {
            return Err(format!("{}: recomputed {} but alignment reports {}", what, got, want));
        }
    }
    if score != aln.raw_score {
        return Err(format!(
            "CIGAR rescores to {} but raw_score is {} (cigar {})",
            score, aln.raw_score, aln.cigar
        ));
    }
    Ok(score)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::build_aa2num;
    use crate::profile::{Profile, PseudoCountParams};

    fn encode(seq: &str) -> Vec<u8> {
        let t = build_aa2num();
        seq.bytes().map(|b| t[b as usize]).collect()
    }

    /// Independent reference: textbook Gotoh with three explicit full matrices;
    /// returns only the optimal score (no pointers, no shared code).
    pub(crate) fn reference_score<F: Fn(usize, usize) -> i32>(m: usize, n: usize, s: F, open: i32, ext: i32) -> i32 {
        let ni = i32::MIN / 4;
        let mut h = vec![vec![0i32; n + 1]; m + 1];
        let mut e = vec![vec![ni; n + 1]; m + 1];
        let mut f = vec![vec![ni; n + 1]; m + 1];
        let mut best = 0;
        for i in 1..=m {
            for j in 1..=n {
                e[i][j] = (h[i][j - 1] - open).max(e[i][j - 1] - ext);
                f[i][j] = (h[i - 1][j] - open).max(f[i - 1][j] - ext);
                let d = h[i - 1][j - 1] + s(i - 1, j - 1);
                h[i][j] = d.max(e[i][j]).max(f[i][j]).max(0);
                best = best.max(h[i][j]);
            }
        }
        best
    }

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }
        fn below(&mut self, n: usize) -> usize {
            (self.next() % n as u64) as usize
        }
    }

    fn mutate(rng: &mut Rng, s: &[u8], indel_len: usize) -> Vec<u8> {
        let mut out = Vec::new();
        let mut i = 0;
        while i < s.len() {
            let r = rng.below(100);
            if r < 4 {
                i += 1 + rng.below(indel_len);
                continue;
            }
            if r < 8 {
                for _ in 0..1 + rng.below(indel_len) {
                    out.push(rng.below(20) as u8);
                }
            }
            out.push(if r < 30 { rng.below(20) as u8 } else { s[i] });
            i += 1;
        }
        out
    }

    fn check(sc: &SeqScorer, t: &[u8], gaps: GapCosts) {
        let m = sc.query_len();
        let reference = reference_score(m, t.len(), |i, j| sc.score(i, t[j]), gaps.open, gaps.extend);
        let full = local_align(sc, t, gaps);
        assert_eq!(full.as_ref().map(|a| a.raw_score).unwrap_or(0), reference, "full score");
        if let Some(a) = &full {
            rescore_alignment(sc, t, gaps, a).unwrap_or_else(|e| panic!("full: {}", e));
        }
        for r in [1usize, 3, 16] {
            let banded = local_align_banded(sc, t, gaps, 0, r);
            assert_eq!(banded.as_ref().map(|a| a.raw_score).unwrap_or(0), reference, "banded r={}", r);
            if let Some(a) = &banded {
                rescore_alignment(sc, t, gaps, a).unwrap_or_else(|e| panic!("banded r={}: {}", r, e));
            }
        }
        let simd = crate::simd::QueryProfile::build(sc).sw_score(t, gaps);
        assert_eq!(simd, reference, "simd");
    }

    #[test]
    fn traceback_rescores_exactly_on_random_gapped_pairs() {
        let m = SubstitutionMatrix::blosum62();
        let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
        for gaps in [
            GapCosts { open: 11, extend: 1 },
            GapCosts { open: 5, extend: 1 },
            GapCosts { open: 3, extend: 1 },
            GapCosts { open: 8, extend: 2 },
            GapCosts { open: 2, extend: 2 },
            GapCosts { open: 1, extend: 1 },
        ] {
            for _ in 0..600 {
                let ql = 1 + rng.below(70);
                let q: Vec<u8> = (0..ql).map(|_| rng.below(20) as u8).collect();
                let t = if rng.below(3) > 0 {
                    mutate(&mut rng, &q, 6)
                } else {
                    (0..1 + rng.below(70)).map(|_| rng.below(20) as u8).collect()
                };
                if t.is_empty() {
                    continue;
                }
                check(&SeqScorer { query: &q, mat: &m }, &t, gaps);
            }
        }
    }

    #[test]
    fn long_gaps_repeats_ties_and_boundaries() {
        let m = SubstitutionMatrix::blosum62();
        let gaps = GapCosts::sequence_default();
        let cases: Vec<(&str, String)> = vec![
            // a 30-residue insertion in the target
            ("MKVLLACDEFGHIKLMNPQRSTVWYACDEFGHIKMKV", format!("MKVLLACDEFGHIKLMNPQR{}STVWYACDEFGHIKMKV", "G".repeat(30))),
            // a 25-residue deletion
            ("MKVLLACDEFGHIKLMNPQRWWWWWWWWWWWWWWWWWWWWWWWWWSTVWYACDEFGHIK", "MKVLLACDEFGHIKLMNPQRSTVWYACDEFGHIK".into()),
            // tandem repeats: many equally good alignments
            ("ACDEFACDEFACDEF", "ACDEFACDEFACDEFACDEFACDEF".into()),
            ("WWWW", "WWWWWWWWWW".into()),
            // single residues and boundary-touching matches
            ("W", "W".into()),
            ("W", "AAAAW".into()),
            ("A", "W".into()),
            ("MKVL", "MKVLXXXXMKVL".into()),
            // ambiguous residues on both sides
            ("MKXXLLAXCDE", "MKVVLLAXCDE".into()),
            ("XXXXXX", "XXXXXX".into()),
        ]
        .into_iter()
        .map(|(a, b): (&str, String)| (a, b))
        .collect();
        for (q, t) in &cases {
            let qn = encode(q);
            let tn = encode(t);
            check(&SeqScorer { query: &qn, mat: &m }, &tn, gaps);
        }
        // empty inputs
        let q = encode("MKV");
        assert!(local_align(&SeqScorer { query: &q, mat: &m }, &[], gaps).is_none());
        assert!(local_align(&SeqScorer { query: &[], mat: &m }, &q, gaps).is_none());
    }

    #[test]
    fn v02_traceback_regression_case() {
        // Minimal case found against the v0.2 aligner: it reported raw 28 with
        // CIGAR 2M1D1M2D2M, which rescores to 18 (a 3-residue deletion split
        // into two runs with a spurious aligned pair in between).
        let m = SubstitutionMatrix::blosum62();
        let q = encode("WFMMHGWW");
        let t = encode("WFMWW");
        let sc = SeqScorer { query: &q, mat: &m };
        let gaps = GapCosts::sequence_default();
        let a = local_align(&sc, &t, gaps).unwrap();
        assert_eq!(a.raw_score, reference_score(q.len(), t.len(), |i, j| sc.score(i, t[j]), 11, 1));
        assert_eq!(rescore_alignment(&sc, &t, gaps, &a).unwrap(), a.raw_score, "cigar {}", a.cigar);
        assert_eq!(a.gap_opens, 1, "cigar {}", a.cigar);
        assert_eq!(a.gaps, 3);
    }

    #[test]
    fn long_gap_is_a_single_run() {
        // the v0.2 bug: a long insertion was split into several runs
        let m = SubstitutionMatrix::blosum62();
        let q = encode("MKVLLACDEFGHIKLMNPQRSTVWYMKVLLACDEFGHIK");
        let t = encode("MKVLLACDEFGHIKLMNPQRGGGGSTVWYMKVLLACDEFGHIK");
        let a = local_align(&SeqScorer { query: &q, mat: &m }, &t, GapCosts::sequence_default()).unwrap();
        assert_eq!(a.gap_opens, 1, "cigar {}", a.cigar);
        assert_eq!(a.cigar, "20M4I19M");
        assert_eq!(a.gaps, 4);
        assert_eq!(a.mismatches, 0);
        assert_eq!(a.identities, 39);
    }

    #[test]
    fn identical_sequences_align_fully() {
        let m = SubstitutionMatrix::blosum62();
        let q = encode("MKVLLACDE");
        let t = encode("MKVLLACDE");
        let sc = SeqScorer { query: &q, mat: &m };
        let aln = local_align(&sc, &t, GapCosts::sequence_default()).unwrap();
        assert_eq!(aln.query_start, 0);
        assert_eq!(aln.query_end, q.len() - 1);
        assert_eq!(aln.identities, q.len());
        assert!((aln.pct_identity() - 100.0).abs() < 1e-9);
        assert_eq!(aln.cigar, format!("{}M", q.len()));
        assert!((aln.query_cov(q.len()) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn gap_is_recovered() {
        let m = SubstitutionMatrix::blosum62();
        let q = encode("MKVLLACDEFGHIK");
        let t = encode("MKVLLADEFGHIK");
        let sc = SeqScorer { query: &q, mat: &m };
        let aln = local_align(&sc, &t, GapCosts::sequence_default()).unwrap();
        assert!(aln.raw_score > 0);
        assert!(aln.cigar.contains('D') || aln.cigar.contains('I'));
        rescore_alignment(&sc, &t, GapCosts::sequence_default(), &aln).unwrap();
    }

    #[test]
    fn local_finds_embedded_match() {
        let m = SubstitutionMatrix::blosum62();
        let q = encode("WWWMKVLLACDEWWW");
        let t = encode("GGGGGMKVLLACDEGGGGG");
        let sc = SeqScorer { query: &q, mat: &m };
        let aln = local_align(&sc, &t, GapCosts::sequence_default()).unwrap();
        assert!(aln.identities >= 8);
    }

    #[test]
    fn profile_scorer_aligns() {
        let m = SubstitutionMatrix::blosum62();
        let msa = ">q\nMKVLLACDE\n>h\nMKILLACDE\n>i\nMRVLLSCDE\n";
        let p = Profile::from_msa_str(msa, &m, PseudoCountParams::default()).unwrap();
        let sc = ProfileScorer { profile: &p };
        let t = encode("MKVLLACDE");
        let aln = local_align(&sc, &t, GapCosts::profile_default()).unwrap();
        assert!(aln.raw_score > 0);
        assert_eq!(aln.query_end - aln.query_start + 1, 9);
        rescore_alignment(&sc, &t, GapCosts::profile_default(), &aln).unwrap();
    }

    #[test]
    fn large_problem_uses_banded_path_and_matches_reference() {
        // m·n above FULL_TRACEBACK_MAX_CELLS exercises score → start → band
        let m = SubstitutionMatrix::blosum62();
        let mut rng = Rng(42);
        let q: Vec<u8> = (0..4500).map(|_| rng.below(20) as u8).collect();
        let mut t = mutate(&mut rng, &q[1000..3500], 8);
        let pre: Vec<u8> = (0..2000).map(|_| rng.below(20) as u8).collect();
        t.splice(0..0, pre);
        assert!((q.len() + 1) * (t.len() + 1) > FULL_TRACEBACK_MAX_CELLS);
        let sc = SeqScorer { query: &q, mat: &m };
        let gaps = GapCosts::sequence_default();
        let a = local_align(&sc, &t, gaps).unwrap();
        rescore_alignment(&sc, &t, gaps, &a).unwrap();
        let simd = crate::simd::QueryProfile::build(&sc).sw_score(&t, gaps);
        assert_eq!(a.raw_score, simd);
    }

    #[test]
    fn in_band_heuristic_reports_its_own_score() {
        let m = SubstitutionMatrix::blosum62();
        let q = encode("MKVLLACDEFGHIKLMNPQRSTVWYACDEFGHIK");
        let t = encode("MKVLLACDEFGHIKGGGGGGGGGGLMNPQRSTVWYACDEFGHIK");
        let sc = SeqScorer { query: &q, mat: &m };
        let gaps = GapCosts::sequence_default();
        let full = local_align(&sc, &t, gaps).unwrap();
        let narrow = local_align_in_band(&sc, &t, gaps, 0, 2).unwrap();
        assert!(narrow.raw_score < full.raw_score);
        // the reported alignment is consistent with its (lower) score
        assert_eq!(rescore_alignment(&sc, &t, gaps, &narrow).unwrap(), narrow.raw_score);
    }

    #[test]
    fn profile_profile_aligns_related_families() {
        let m = SubstitutionMatrix::blosum62();
        let msa_a = ">a\nMKVLLACDEFGHIKLMNPQR\n>a2\nMKILLSCDEFGHLKLMNPQR\n";
        let msa_b = ">b\nMKVLLACDEFGHIKLMNPQR\n>b2\nMRVLLACEEFGHIKLLNPQR\n";
        let pa = Profile::from_msa_str(msa_a, &m, PseudoCountParams::default()).unwrap();
        let pb = Profile::from_msa_str(msa_b, &m, PseudoCountParams::default()).unwrap();
        let pback: Vec<f64> = (0..PROFILE_AA_SIZE).map(|i| m.p_back[i]).collect();
        let aln = align_profile_profile(&pa, &pb, GapCosts::profile_default(), &pback).unwrap();
        assert!(aln.raw_score > 0);
        assert!(aln.pct_identity() > 50.0, "{:.1}%", aln.pct_identity());
    }

    #[test]
    fn profile_profile_separates_unrelated() {
        let m = SubstitutionMatrix::blosum62();
        let msa_a = ">a\nMKVLLACDEFGHIKLMNPQR\n>a2\nMKILLSCDEFGHLKLMNPQR\n";
        let msa_related = ">b\nMKVLLACDEFGHIKLMNPQR\n>b2\nMRVLLACEEFGHIKLLNPQR\n";
        let msa_unrel = ">c\nWYWYWYWYWYWYWYWYWYWY\n>c2\nYWYWYWYWYWYWYWYWYWYW\n";
        let pa = Profile::from_msa_str(msa_a, &m, PseudoCountParams::default()).unwrap();
        let pr = Profile::from_msa_str(msa_related, &m, PseudoCountParams::default()).unwrap();
        let pu = Profile::from_msa_str(msa_unrel, &m, PseudoCountParams::default()).unwrap();
        let pback: Vec<f64> = (0..PROFILE_AA_SIZE).map(|i| m.p_back[i]).collect();
        let gaps = GapCosts::profile_default();
        let related = align_profile_profile(&pa, &pr, gaps, &pback).unwrap().raw_score;
        let unrelated = align_profile_profile(&pa, &pu, gaps, &pback).map(|a| a.raw_score).unwrap_or(0);
        assert!(related > unrelated, "related ({}) must outscore unrelated ({})", related, unrelated);
    }

    #[test]
    fn cigar_parser_rejects_garbage() {
        assert!(cigar_runs("3M2I").is_ok());
        assert!(cigar_runs("M").is_err());
        assert!(cigar_runs("3X").is_err());
        assert!(cigar_runs("3M2").is_err());
        assert!(cigar_runs("0M").is_err());
    }
}

#[cfg(test)]
mod microbench {
    use super::*;
    use crate::matrix::SubstitutionMatrix;

    /// `cargo test --release --lib align::microbench -- --ignored --nocapture`
    #[test]
    #[ignore]
    fn microbench_traceback_paths() {
        let m = SubstitutionMatrix::blosum62();
        let mut rng = crate::evalue::SplitMix64::new(7);
        let mut pairs = Vec::new();
        for _ in 0..200 {
            let q: Vec<u8> = (0..180).map(|_| rng.below(20) as u8).collect();
            // homolog: 60% identity with a few indels, embedded in random flanks
            let mut t: Vec<u8> = (0..20).map(|_| rng.below(20) as u8).collect();
            for &c in &q {
                if rng.below(40) == 0 {
                    continue;
                }
                t.push(if rng.below(10) < 6 { c } else { rng.below(20) as u8 });
                if rng.below(50) == 0 {
                    t.push(rng.below(20) as u8);
                }
            }
            t.extend((0..20).map(|_| rng.below(20) as u8));
            pairs.push((q, t));
        }
        let g = GapCosts::sequence_default();
        let mut s = AlignScratch::new();
        let tables: Vec<ScoreTable> = pairs
            .iter()
            .map(|(q, _)| ScoreTable::from_scorer(&crate::align::SeqScorer { query: q, mat: &m }))
            .collect();
        let cells: usize = pairs.iter().map(|(q, t)| q.len() * t.len()).sum();
        let mut best = [f64::MAX; 2];
        let mut sums = [0i64; 2];
        for _rep in 0..7 {
            for (v, slot) in best.iter_mut().enumerate() {
                let t0 = std::time::Instant::now();
                let mut tot = 0i64;
                for (k, (_, t)) in pairs.iter().enumerate() {
                    let tab = &tables[k];
                    let a = match v {
                        0 => local_align_with(tab, t, g, &mut s),
                        _ => local_align_banded_with(tab, t, g, 32, &mut s),
                    };
                    tot += a.map(|a| a.raw_score as i64).unwrap_or(0);
                }
                *slot = slot.min(t0.elapsed().as_nanos() as f64 / cells as f64);
                sums[v] = tot;
            }
        }
        for (v, name) in ["full traceback", "score/start/banded"].iter().enumerate() {
            eprintln!("{}: {:.2} ns/cell (best of 7; checksum {})", name, best[v], sums[v]);
        }
    }
}
