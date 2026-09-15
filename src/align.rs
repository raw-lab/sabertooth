//! Local alignment (Smith–Waterman with Gotoh affine gaps).
//!
//! A single dynamic-programming core is driven by a [`Scorer`] so the exact same
//! optimal local-alignment logic serves both **sequence↔sequence** search (via
//! the integer substitution matrix) and **profile↔sequence** search (via a
//! PSSM). Traceback yields the alignment coordinates, a CIGAR string, percent
//! identity and the raw score used for E-value computation.

use crate::alphabet::PROFILE_AA_SIZE;
use crate::matrix::SubstitutionMatrix;
use crate::profile::Profile;

/// Abstraction over "the score of placing this target residue against query
/// column `qpos`". Implemented by both a plain sequence and a PSSM profile.
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

/// Affine gap penalties (expressed in the same score units as the scorer).
#[derive(Clone, Copy)]
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
}

/// The result of a local alignment.
#[derive(Clone)]
pub struct Alignment {
    pub raw_score: i32,
    pub query_start: usize,
    pub query_end: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub aln_len: usize,
    pub identities: usize,
    pub gaps: usize,
    pub cigar: String,
}

impl Alignment {
    pub fn pct_identity(&self) -> f64 {
        if self.aln_len == 0 {
            0.0
        } else {
            100.0 * self.identities as f64 / self.aln_len as f64
        }
    }
}

/// Banded Smith–Waterman–Gotoh local alignment. Only cells whose diagonal
/// `(j-1) − (i-1)` lies within `band_radius` of `center_diag` (the seed diagonal
/// `target_pos − query_pos` from the prefilter) are evaluated, giving `O(m·w)`
/// time and traceback memory (`w = 2·band_radius+1`) instead of `O(m·n)`. This
/// is the long-sequence path: the full-matrix aligner stores an `(m+1)(n+1)`
/// traceback that is quadratic in memory, which the hardening notes flagged.
///
/// With a band wide enough to contain the optimal path, the result is identical
/// to [`local_align`]. Cells outside the band are treated as the local baseline
/// (0), so the alignment is constrained to stay near the seed diagonal — exactly
/// the assumption the prefilter's diagonal has already established.
pub fn local_align_banded<S: Scorer>(
    scorer: &S,
    target: &[u8],
    gaps: GapCosts,
    center_diag: i64,
    band_radius: usize,
) -> Option<Alignment> {
    let m = scorer.query_len();
    let n = target.len();
    if m == 0 || n == 0 {
        return None;
    }
    let neg_inf = i32::MIN / 4;
    let br = band_radius as i64;
    let bw = 2 * band_radius + 1; // band width in columns

    // Full-width rolling H rows and F column (O(n), not the quadratic term); the
    // quadratic cost was the traceback matrix, which is what we band below.
    let mut h_prev = vec![0i32; n + 2];
    let mut h_cur = vec![0i32; n + 2];
    let mut f_col = vec![neg_inf; n + 2];

    // Banded traceback: m rows × bw columns. Cell (i,j) maps to offset
    // `j − (i + center_diag) + band_radius ∈ [0, bw)`.
    let mut tb = vec![TB_STOP; m * bw];

    let mut best = 0i32;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    // Column band for a given row i (1-based), clamped to [1, n]. Empty if the
    // band lies entirely off the target.
    let col_band = |i: usize| -> (usize, usize) {
        let center = i as i64 + center_diag; // j at the center diagonal
        let lo = (center - br).max(1);
        let hi = (center + br).min(n as i64);
        if lo > hi {
            (1, 0) // empty
        } else {
            (lo as usize, hi as usize)
        }
    };

    let mut prev_hi = 0usize; // right edge filled by the previous row
    for i in 1..=m {
        let (j_lo, j_hi) = col_band(i);
        if j_lo > j_hi {
            // band off the target for this row; nothing to do, advance rows
            h_prev.iter_mut().for_each(|_| {}); // no-op; keep buffers
            std::mem::swap(&mut h_prev, &mut h_cur);
            for v in h_cur.iter_mut() {
                *v = 0;
            }
            prev_hi = 0;
            continue;
        }
        // Reset the left neighbour (out of this row's band) so it reads as 0.
        if j_lo >= 1 {
            h_cur[j_lo - 1] = 0;
        }
        // Newly entered right columns (band shifted right) must not read stale
        // values from an earlier row: clear them to the local/neg baseline.
        let clear_from = prev_hi + 1;
        for j in clear_from.max(j_lo)..=j_hi {
            h_prev[j] = 0;
            f_col[j] = neg_inf;
        }

        let row_base = (i - 1) * bw;
        let center = i as i64 + center_diag;
        let mut e = neg_inf; // horizontal-gap running value
        // E has no carry-in at the left band edge.
        for j in j_lo..=j_hi {
            let s = scorer.score(i - 1, target[j - 1]);
            let diag = h_prev[j - 1] + s;

            e = if j == j_lo {
                h_cur[j - 1] - gaps.open // left neighbour is baseline 0
            } else {
                (h_cur[j - 1] - gaps.open).max(e - gaps.extend)
            };
            f_col[j] = (h_prev[j] - gaps.open).max(f_col[j] - gaps.extend);

            let mut val = diag;
            let mut dir = TB_DIAG;
            if e > val {
                val = e;
                dir = TB_LEFT;
            }
            if f_col[j] > val {
                val = f_col[j];
                dir = TB_UP;
            }
            if val <= 0 {
                val = 0;
                dir = TB_STOP;
            }
            h_cur[j] = val;
            let off = (j as i64 - center + br) as usize; // in [0, bw)
            tb[row_base + off] = dir;

            if val > best {
                best = val;
                best_i = i;
                best_j = j;
            }
        }
        prev_hi = j_hi;
        std::mem::swap(&mut h_prev, &mut h_cur);
    }

    if best <= 0 {
        return None;
    }

    // --- banded traceback ------------------------------------------------
    let mut i = best_i;
    let mut j = best_j;
    let mut ops: Vec<u8> = Vec::new();
    let mut identities = 0usize;
    let mut gap_cols = 0usize;

    while i > 0 && j > 0 {
        let center = i as i64 + center_diag;
        let off = j as i64 - center + br;
        if off < 0 || off as usize >= bw {
            break; // stepped out of the band (shouldn't happen on an optimal path)
        }
        let dir = tb[(i - 1) * bw + off as usize];
        match dir {
            TB_DIAG => {
                ops.push(b'M');
                if scorer.query_residue(i - 1) == target[j - 1]
                    && (scorer.query_residue(i - 1) as usize) < PROFILE_AA_SIZE
                {
                    identities += 1;
                }
                i -= 1;
                j -= 1;
            }
            TB_LEFT => {
                ops.push(b'I');
                gap_cols += 1;
                j -= 1;
            }
            TB_UP => {
                ops.push(b'D');
                gap_cols += 1;
                i -= 1;
            }
            _ => break,
        }
    }

    let query_start = i;
    let target_start = j;
    ops.reverse();
    let aln_len = ops.len();
    let cigar = run_length_encode(&ops);

    Some(Alignment {
        raw_score: best,
        query_start,
        query_end: best_i - 1,
        target_start,
        target_end: best_j - 1,
        aln_len,
        identities,
        gaps: gap_cols,
        cigar,
    })
}


// Traceback direction bits for each cell of H.
const TB_STOP: u8 = 0;
const TB_DIAG: u8 = 1;
const TB_LEFT: u8 = 2; // gap in query (consume target)
const TB_UP: u8 = 3; // gap in target (consume query)

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
/// by co-emission ([`pp_col_score`]). This is Sabertooth's profile–profile search
/// core (cf. HHsearch / MMseqs2 profile-vs-profile). Reported identity is the
/// fraction of aligned columns whose consensus residues agree.
pub fn align_profile_profile(
    a: &Profile,
    b: &Profile,
    gaps: GapCosts,
    pback: &[f64],
) -> Option<Alignment> {
    let m = a.query_len;
    let n = b.query_len;
    if m == 0 || n == 0 {
        return None;
    }
    let neg_inf = i32::MIN / 4;
    let mut h_prev = vec![0i32; n + 1];
    let mut h_cur = vec![0i32; n + 1];
    let mut f_col = vec![neg_inf; n + 1];
    let mut tb = vec![TB_STOP; (m + 1) * (n + 1)];
    let (mut best, mut bi, mut bj) = (0i32, 0usize, 0usize);

    for i in 1..=m {
        let mut e = neg_inf;
        h_cur[0] = 0;
        for j in 1..=n {
            let s = pp_col_score(a, i - 1, b, j - 1, pback);
            let diag = h_prev[j - 1] + s;
            e = (h_cur[j - 1] - gaps.open).max(e - gaps.extend);
            f_col[j] = (h_prev[j] - gaps.open).max(f_col[j] - gaps.extend);
            let mut val = diag;
            let mut dir = TB_DIAG;
            if e > val {
                val = e;
                dir = TB_LEFT;
            }
            if f_col[j] > val {
                val = f_col[j];
                dir = TB_UP;
            }
            if val <= 0 {
                val = 0;
                dir = TB_STOP;
            }
            h_cur[j] = val;
            tb[i * (n + 1) + j] = dir;
            if val > best {
                best = val;
                bi = i;
                bj = j;
            }
        }
        std::mem::swap(&mut h_prev, &mut h_cur);
    }
    if best <= 0 {
        return None;
    }

    let (mut i, mut j) = (bi, bj);
    let mut ops: Vec<u8> = Vec::new();
    let mut identities = 0usize;
    let mut gap_cols = 0usize;
    while i > 0 && j > 0 {
        match tb[i * (n + 1) + j] {
            TB_DIAG => {
                ops.push(b'M');
                if a.consensus[i - 1] == b.consensus[j - 1] {
                    identities += 1;
                }
                i -= 1;
                j -= 1;
            }
            TB_LEFT => {
                ops.push(b'I');
                gap_cols += 1;
                j -= 1;
            }
            TB_UP => {
                ops.push(b'D');
                gap_cols += 1;
                i -= 1;
            }
            _ => break,
        }
    }
    let query_start = i;
    let target_start = j;
    ops.reverse();
    let aln_len = ops.len();
    let cigar = run_length_encode(&ops);
    Some(Alignment {
        raw_score: best,
        query_start,
        query_end: bi - 1,
        target_start,
        target_end: bj - 1,
        aln_len,
        identities,
        gaps: gap_cols,
        cigar,
    })
}

/// Optimal Smith–Waterman–Gotoh local alignment of a scorer's query against a
/// numeric target sequence. Returns `None` if the best score is non-positive.
pub fn local_align<S: Scorer>(
    scorer: &S,
    target: &[u8],
    gaps: GapCosts,
) -> Option<Alignment> {
    let m = scorer.query_len(); // rows (query)
    let n = target.len(); // cols (target)
    if m == 0 || n == 0 {
        return None;
    }

    let neg_inf = i32::MIN / 4;
    // Rolling H rows. The horizontal-gap term E[i][j] depends only on cells to
    // its left in the same row, so it is a scalar carried left→right. The
    // vertical-gap term F[i][j], however, depends on F[i-1][j] — the previous
    // *row* at the same column — so it must be a per-column array carried across
    // rows. (An earlier version used a scalar for F too, which silently
    // under-scored any deletion spanning more than one query row; the striped
    // SIMD kernel and a brute-force reference both caught it.)
    let mut h_prev = vec![0i32; n + 1];
    let mut h_cur = vec![0i32; n + 1];
    let mut f_col = vec![neg_inf; n + 1];

    // Full traceback matrix (m+1) x (n+1). For protein lengths this is fine;
    // very long targets would want a linear-space (Hirschberg) or banded pass.
    let mut tb = vec![TB_STOP; (m + 1) * (n + 1)];

    let mut best = 0i32;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    for i in 1..=m {
        let mut e = neg_inf; // horizontal-gap running value E[i][j-1]
        h_cur[0] = 0;
        for j in 1..=n {
            let s = scorer.score(i - 1, target[j - 1]);
            let diag = h_prev[j - 1] + s;

            // horizontal gap (gap in query, consume target): E[i][j]
            e = (h_cur[j - 1] - gaps.open).max(e - gaps.extend);

            // vertical gap (gap in target, consume query): F[i][j] from row i-1
            f_col[j] = (h_prev[j] - gaps.open).max(f_col[j] - gaps.extend);

            // best of the three, floored at 0 (local)
            let mut val = diag;
            let mut dir = TB_DIAG;
            if e > val {
                val = e;
                dir = TB_LEFT;
            }
            if f_col[j] > val {
                val = f_col[j];
                dir = TB_UP;
            }
            if val <= 0 {
                val = 0;
                dir = TB_STOP;
            }
            h_cur[j] = val;
            tb[i * (n + 1) + j] = dir;

            if val > best {
                best = val;
                best_i = i;
                best_j = j;
            }
        }
        std::mem::swap(&mut h_prev, &mut h_cur);
    }

    if best <= 0 {
        return None;
    }

    // --- traceback -------------------------------------------------------
    let mut i = best_i;
    let mut j = best_j;
    let mut ops: Vec<u8> = Vec::new(); // 'M','I','D'
    let mut identities = 0usize;
    let mut gap_cols = 0usize;

    while i > 0 && j > 0 {
        let dir = tb[i * (n + 1) + j];
        match dir {
            TB_DIAG => {
                ops.push(b'M');
                if scorer.query_residue(i - 1) == target[j - 1]
                    && (scorer.query_residue(i - 1) as usize) < PROFILE_AA_SIZE
                {
                    identities += 1;
                }
                i -= 1;
                j -= 1;
            }
            TB_LEFT => {
                ops.push(b'I'); // insertion in target relative to query
                gap_cols += 1;
                j -= 1;
            }
            TB_UP => {
                ops.push(b'D'); // deletion in target relative to query
                gap_cols += 1;
                i -= 1;
            }
            _ => break, // TB_STOP
        }
    }

    let query_start = i; // 0-based start (i now points before first aligned col)
    let target_start = j;
    ops.reverse();

    let aln_len = ops.len();
    let cigar = run_length_encode(&ops);

    Some(Alignment {
        raw_score: best,
        query_start,
        query_end: best_i - 1,
        target_start,
        target_end: best_j - 1,
        aln_len,
        identities,
        gaps: gap_cols,
        cigar,
    })
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
        out.push_str(&format!("{}{}", run, c as char));
        i += run;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::build_aa2num;
    use crate::fasta::SeqDb;
    use std::io::Cursor;
    use crate::profile::{Profile, PseudoCountParams};

    fn encode(seq: &str) -> Vec<u8> {
        let t = build_aa2num();
        seq.bytes().map(|b| t[b as usize]).collect()
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
    }

    #[test]
    fn gap_is_recovered() {
        let m = SubstitutionMatrix::blosum62();
        // target has a 3-residue deletion in the middle relative to query
        let q = encode("MKVLLACDEFGHIK");
        let t = encode("MKVLLADEFGHIK"); // dropped one residue
        let sc = SeqScorer { query: &q, mat: &m };
        let aln = local_align(&sc, &t, GapCosts::sequence_default()).unwrap();
        assert!(aln.raw_score > 0);
        assert!(aln.cigar.contains('D') || aln.cigar.contains('I'));
    }

    #[test]
    fn local_finds_embedded_match() {
        let m = SubstitutionMatrix::blosum62();
        let q = encode("WWWMKVLLACDEWWW");
        let t = encode("GGGGGMKVLLACDEGGGGG");
        let sc = SeqScorer { query: &q, mat: &m };
        let aln = local_align(&sc, &t, GapCosts::sequence_default()).unwrap();
        // the conserved MKVLLACDE core should be found
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
        // sanity: db reader unaffected
        let db = SeqDb::from_reader(Cursor::new(">x\nMKV\n"), &m.aa2num).unwrap();
        assert_eq!(db.len(), 1);
    }

    // Reference: brute-force affine-gap Smith–Waterman with a full matrix, used
    // to cross-check the banded aligner independently of local_align.
    fn brute_local_score(q: &[u8], t: &[u8], mat: &SubstitutionMatrix, open: i32, ext: i32) -> i32 {
        let (m, n) = (q.len(), t.len());
        let ni = i32::MIN / 4;
        let mut h = vec![vec![0i32; n + 1]; m + 1];
        let mut e = vec![vec![ni; n + 1]; m + 1];
        let mut f = vec![vec![ni; n + 1]; m + 1];
        let mut best = 0;
        for i in 1..=m {
            for j in 1..=n {
                e[i][j] = (h[i][j - 1] - open).max(e[i][j - 1] - ext);
                f[i][j] = (h[i - 1][j] - open).max(f[i - 1][j] - ext);
                let diag = h[i - 1][j - 1] + mat.score(q[i - 1], t[j - 1]) as i32;
                h[i][j] = diag.max(e[i][j]).max(f[i][j]).max(0);
                best = best.max(h[i][j]);
            }
        }
        best
    }

    #[test]
    fn banded_wide_equals_full() {
        let m = SubstitutionMatrix::blosum62();
        // an embedded, gapped homology so the optimal path is non-trivial
        let q = encode("GGGMKVLLACDEFGHIKLMNPQRGGG");
        let t = encode("AAAAAMKVLLACDEFGHKLMNPQRAAAAA"); // one deletion (I dropped)
        let sc = SeqScorer { query: &q, mat: &m };
        let gaps = GapCosts::sequence_default();
        let full = local_align(&sc, &t, gaps).unwrap();
        // wide band (radius >= n) must reproduce the full-matrix result exactly
        let banded =
            local_align_banded(&sc, &t, gaps, full.target_start as i64 - full.query_start as i64, t.len())
                .unwrap();
        assert_eq!(banded.raw_score, full.raw_score);
        assert_eq!(banded.query_start, full.query_start);
        assert_eq!(banded.query_end, full.query_end);
        assert_eq!(banded.target_start, full.target_start);
        assert_eq!(banded.target_end, full.target_end);
        assert_eq!(banded.cigar, full.cigar);
        // and it agrees with an independent brute-force score
        assert_eq!(full.raw_score, brute_local_score(&q, &t, &m, gaps.open, gaps.extend));
    }

    #[test]
    fn banded_narrow_around_diagonal_recovers_alignment() {
        let m = SubstitutionMatrix::blosum62();
        let q = encode("MKVLLACDEFGHIKLMNPQRST");
        let t = encode("MKVLLACDEFGHIKLMNPQRST");
        let sc = SeqScorer { query: &q, mat: &m };
        let gaps = GapCosts::sequence_default();
        // diagonal 0 (self-alignment); a tiny band still captures the full match
        let banded = local_align_banded(&sc, &t, gaps, 0, 2).unwrap();
        let full = local_align(&sc, &t, gaps).unwrap();
        assert_eq!(banded.raw_score, full.raw_score);
        assert_eq!(banded.identities, q.len());
    }

    #[test]
    fn banded_matches_full_across_random_offsets() {
        // For a range of embed offsets (hence diagonals), a band comfortably
        // wider than the indel content reproduces the full-matrix score.
        let m = SubstitutionMatrix::blosum62();
        let core = "MKVLLACDEFGHIKLMNPQR";
        let gaps = GapCosts::sequence_default();
        for pad in [0usize, 3, 7, 11] {
            let q = encode(core);
            let t = encode(&format!("{}{}{}", "A".repeat(pad), core, "A".repeat(pad)));
            let sc = SeqScorer { query: &q, mat: &m };
            let full = local_align(&sc, &t, gaps).unwrap();
            let diag = full.target_start as i64 - full.query_start as i64;
            let banded = local_align_banded(&sc, &t, gaps, diag, 8).unwrap();
            assert_eq!(
                banded.raw_score, full.raw_score,
                "pad {} diag {}",
                pad, diag
            );
        }
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
        assert!(aln.raw_score > 0, "related profiles should score > 0");
        assert!(
            aln.pct_identity() > 50.0,
            "related families should share consensus (got {:.1}%)",
            aln.pct_identity()
        );
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
        let unrelated = align_profile_profile(&pa, &pu, gaps, &pback)
            .map(|a| a.raw_score)
            .unwrap_or(0);
        assert!(
            related > unrelated,
            "related ({}) must outscore unrelated ({})",
            related,
            unrelated
        );
    }
}
