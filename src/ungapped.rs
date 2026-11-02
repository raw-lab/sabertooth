//! Batched ungapped diagonal scoring (prefilter evidence).
//!
//! For every `(target, diagonal)` the prefilter keeps, the evidence is the
//! best ungapped segment anywhere on the diagonal (Kadane's maximum-subarray
//! score, exactly [`kadane_score`]). v0.3 computed it one diagonal at a time
//! (scalar X-drop by default); at 68 M residues that stage was the largest
//! prefilter cost — one random target access and a long scalar loop per
//! diagonal, for ~40,000 diagonals per query.
//!
//! Here diagonals are scored in **batches of SIMD lanes** (16 SSSE3, 32 AVX2,
//! 64 AVX-512BW), the layout MMseqs2 uses for its ungapped stage: every lane
//! holds a different `(target, diagonal)` and all lanes advance along the same
//! query position, so the score row of that position is shared and looked up
//! per lane with byte shuffles. Target residues are first transposed into a
//! position-major buffer (lanes issue their memory accesses together instead
//! of one dependent miss per diagonal).
//!
//! Scores use saturating `u8` arithmetic on biased scores (`score + bias ≥ 0`):
//! `cur = max(0, cur + s)`, `best = max(best, cur)`. A lane is exact while no
//! cell saturated, which is guaranteed when `best ≤ 255 − max(score + bias)`;
//! lanes above that limit are recomputed exactly in scalar code, so results
//! never depend on the kernel (tests compare every kernel with the scalar
//! reference).

use crate::alphabet::ALPHABET_SIZE;
use std::sync::OnceLock;

/// Residue code written outside a lane's valid diagonal range. It scores 0
/// (like every code ≥ 21, see [`kadane_score`]), so before the range the
/// running score stays 0 and after it the score cannot grow: positions outside
/// the diagonal never change the result. Transposition clamps every byte to
/// at most this value.
const SENTINEL: u8 = 21;

/// Query positions `[lo, hi)` at which diagonal `d` (`target_pos −
/// query_pos`) lies inside a query of length `m` and a target of length `n`.
#[inline]
pub fn diag_range(m: usize, n: usize, d: i64) -> (usize, usize) {
    let lo = (-d).max(0);
    let hi = (m as i64).min(n as i64 - d);
    if hi <= lo {
        (0, 0)
    } else {
        (lo as usize, hi as usize)
    }
}

/// Best ungapped segment score on diagonal `d` within query positions
/// `[lo, hi)` (a sub-range of [`diag_range`]), 0 if every segment is
/// non-positive. Residues outside the scoring alphabet score 0 (as in the
/// gapped kernels).
#[inline]
pub fn kadane_score(rows: &[[i32; ALPHABET_SIZE]], t: &[u8], d: i64, lo: usize, hi: usize) -> i32 {
    let (mut best, mut cur) = (0i32, 0i32);
    for q in lo..hi {
        let tr = t[(q as i64 + d) as usize] as usize;
        let s = if tr < ALPHABET_SIZE { rows[q][tr] } else { 0 };
        cur = cur.max(0) + s;
        best = best.max(cur);
    }
    best
}

/// Best ungapped segment within `[lo, hi)` with its query range: `(score,
/// q_start, q_end)` (inclusive; the first maximum in query order; `(0, lo,
/// lo)`-style result when nothing scores above 0). Over the whole diagonal
/// this is v0.3's whole-diagonal evidence, including tie-breaking.
#[inline]
pub fn kadane_segment(rows: &[[i32; ALPHABET_SIZE]], t: &[u8], d: i64, lo: usize, hi: usize) -> (i32, usize, usize) {
    let (mut best, mut bs, mut be) = (0i32, 0usize, 0usize);
    let (mut cur, mut cs) = (0i32, lo);
    for q in lo..hi {
        let tr = t[(q as i64 + d) as usize] as usize;
        let s = if tr < ALPHABET_SIZE { rows[q][tr] } else { 0 };
        if cur <= 0 {
            cur = s;
            cs = q;
        } else {
            cur += s;
        }
        if cur > best {
            best = cur;
            bs = cs;
            be = q;
        }
    }
    (best, bs, be)
}

/// Kernel family for diagonal batches.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DiagKernel {
    Scalar,
    Ssse3,
    Avx2,
    Avx512,
}

impl DiagKernel {
    /// Lanes per batch.
    pub fn lanes(self) -> usize {
        match self {
            DiagKernel::Scalar => 1,
            DiagKernel::Ssse3 => 16,
            DiagKernel::Avx2 => 32,
            DiagKernel::Avx512 => 64,
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            DiagKernel::Scalar => "scalar",
            DiagKernel::Ssse3 => "ssse3-u8x16",
            DiagKernel::Avx2 => "avx2-u8x32",
            DiagKernel::Avx512 => "avx512bw-u8x64",
        }
    }
}

/// Kernel for this machine: AVX2, else SSSE3, else scalar; `SABERTOOTH_SIMD`
/// overrides (`avx512` opts in to the 64-lane kernel, `avx2`, `sse2` → SSSE3
/// when present, `scalar`).
pub fn diag_kernel() -> DiagKernel {
    static K: OnceLock<DiagKernel> = OnceLock::new();
    *K.get_or_init(|| {
        let want = std::env::var("SABERTOOTH_SIMD").unwrap_or_default();
        detect(&want)
    })
}

fn detect(want: &str) -> DiagKernel {
    #[cfg(target_arch = "x86_64")]
    {
        let ssse3 = std::is_x86_feature_detected!("ssse3");
        let avx2 = std::is_x86_feature_detected!("avx2");
        let avx512 = std::is_x86_feature_detected!("avx512bw");
        match want {
            "scalar" => DiagKernel::Scalar,
            "sse2" => {
                if ssse3 {
                    DiagKernel::Ssse3
                } else {
                    DiagKernel::Scalar
                }
            }
            "avx2" => {
                if avx2 {
                    DiagKernel::Avx2
                } else if ssse3 {
                    DiagKernel::Ssse3
                } else {
                    DiagKernel::Scalar
                }
            }
            "avx512" => {
                if avx512 {
                    DiagKernel::Avx512
                } else if avx2 {
                    DiagKernel::Avx2
                } else if ssse3 {
                    DiagKernel::Ssse3
                } else {
                    DiagKernel::Scalar
                }
            }
            // AVX2 by default: the 64-lane AVX-512 kernel measured slower on the
            // benchmark machine (Cascade Lake); `SABERTOOTH_SIMD=avx512` opts in
            _ => {
                if avx2 {
                    DiagKernel::Avx2
                } else if ssse3 {
                    DiagKernel::Ssse3
                } else {
                    DiagKernel::Scalar
                }
            }
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let _ = want;
        DiagKernel::Scalar
    }
}

/// Per-query biased `u8` score tables for the SIMD kernels: for each query
/// position, 16 bytes for residues 0..15 and 16 bytes for 16..31 (residue 20 =
/// `X`; 21..31 score 0, 21 being the sentinel).
pub struct DiagProfile {
    m: usize,
    bias: u8,
    /// lanes whose u8 best exceeds this are recomputed exactly
    limit: u8,
    tables: Vec<[u8; 32]>,
}

impl Default for DiagProfile {
    /// An empty profile (call [`DiagProfile::rebuild`] before use).
    fn default() -> Self {
        DiagProfile {
            m: 0,
            bias: 0,
            limit: 0,
            tables: Vec::new(),
        }
    }
}

/// Minimum `u8` headroom for the SIMD path to be worthwhile.
const MIN_LIMIT: i32 = 64;

impl DiagProfile {
    /// `None` when the score range leaves too little `u8` headroom (e.g.
    /// PSSM queries); callers then use [`kadane_score`].
    pub fn build(rows: &[[i32; ALPHABET_SIZE]]) -> Option<DiagProfile> {
        let mut p = DiagProfile {
            m: 0,
            bias: 0,
            limit: 0,
            tables: Vec::new(),
        };
        if p.rebuild(rows) {
            Some(p)
        } else {
            None
        }
    }

    /// Reuse the allocation for another query; `false` if unsuitable.
    pub fn rebuild(&mut self, rows: &[[i32; ALPHABET_SIZE]]) -> bool {
        let (mut lo, mut hi) = (0i32, 0i32); // residues 21..30 score 0
        for r in rows {
            for &v in r.iter() {
                lo = lo.min(v);
                hi = hi.max(v);
            }
        }
        let bias = -lo;
        let limit = 255 - (hi + bias);
        if rows.is_empty() || bias > 255 || limit < MIN_LIMIT {
            return false;
        }
        self.m = rows.len();
        self.bias = bias as u8;
        self.limit = limit as u8;
        self.tables.clear();
        self.tables.reserve(rows.len());
        for r in rows {
            let mut t = [bias as u8; 32];
            for (res, &v) in r.iter().enumerate() {
                t[res] = (v + bias) as u8;
            }
            self.tables.push(t);
        }
        true
    }

    pub fn query_len(&self) -> usize {
        self.m
    }
}

/// One diagonal to score: the target's residues, the diagonal and the query
/// positions `[lo, hi)` to score (inside [`diag_range`]).
#[derive(Clone, Copy)]
pub struct Lane<'a> {
    pub target: &'a [u8],
    pub diag: i64,
    pub lo: usize,
    pub hi: usize,
}

/// Reusable transposition buffers.
#[derive(Default)]
pub struct DiagScratch {
    /// lane-major residues: `lm[lane * w + i]`
    lm: Vec<u8>,
    /// position-major residues: `pm[i * lanes + lane]`
    pm: Vec<u8>,
    best: Vec<u8>,
}

/// Score a batch of at most `kernel.lanes()` diagonals: `out[i]` = exact
/// [`kadane_score`] of `lanes[i]`. `rows` must be the rows `prof` was built
/// from.
pub fn score_batch(
    kernel: DiagKernel,
    prof: &DiagProfile,
    rows: &[[i32; ALPHABET_SIZE]],
    lanes: &[Lane],
    out: &mut [i32],
    s: &mut DiagScratch,
) {
    let nl = kernel.lanes();
    debug_assert!(lanes.len() <= nl && out.len() >= lanes.len());
    if kernel == DiagKernel::Scalar || lanes.is_empty() {
        for (o, l) in out.iter_mut().zip(lanes) {
            *o = kadane_score(rows, l.target, l.diag, l.lo, l.hi);
        }
        return;
    }
    let m = prof.m;
    // union of the lanes' valid query ranges (the kernel runs only over it)
    let mut q_lo = m;
    let mut q_hi = 0usize;
    let range = |lane: &Lane| -> (usize, usize) {
        debug_assert!({
            let (a, b) = diag_range(m, lane.target.len(), lane.diag);
            lane.lo >= lane.hi || (a <= lane.lo && lane.hi <= b)
        });
        (lane.lo, lane.hi.max(lane.lo))
    };
    for lane in lanes {
        let (i0, i1) = range(lane);
        if i0 < i1 {
            q_lo = q_lo.min(i0);
            q_hi = q_hi.max(i1);
        }
    }
    s.best.clear();
    s.best.resize(nl, 0);
    if q_lo < q_hi {
        // 16-aligned window of positions to transpose
        let a = q_lo & !15;
        let b = (q_hi + 15) & !15;
        let w = b - a;
        // 1) lane-major copy of each lane's diagonal (sentinel outside)
        if s.lm.len() < nl * w {
            s.lm.resize(nl * w, SENTINEL);
        }
        for (l, lane) in lanes.iter().enumerate() {
            let row = &mut s.lm[l * w..(l + 1) * w];
            let (i0, i1) = range(lane);
            if i0 >= i1 {
                row.fill(SENTINEL);
                continue;
            }
            let (r0, r1) = (i0 - a, i1 - a);
            row[..r0].fill(SENTINEL);
            let src = &lane.target[(i0 as i64 + lane.diag) as usize..(i1 as i64 + lane.diag) as usize];
            row[r0..r1].copy_from_slice(src);
            row[r1..].fill(SENTINEL);
        }
        for l in lanes.len()..nl {
            s.lm[l * w..(l + 1) * w].fill(SENTINEL);
        }
        // 2) transpose 16 x 16 blocks to position-major (clamping codes);
        // every byte of pm[..w * nl] is written
        if s.pm.len() < w * nl {
            s.pm.resize(w * nl, SENTINEL);
        }
        // SAFETY: SSE2 is part of the x86-64 baseline; both buffers hold
        // at least w * nl bytes and every block lies inside them.
        #[cfg(target_arch = "x86_64")]
        unsafe {
            let mut src = [std::ptr::null::<u8>(); 16];
            for g in 0..nl / 16 {
                for pb in 0..w / 16 {
                    for (r, sp) in src.iter_mut().enumerate() {
                        *sp = s.lm.as_ptr().add((g * 16 + r) * w + pb * 16);
                    }
                    x86::transpose16(&src, s.pm.as_mut_ptr().add(pb * 16 * nl + g * 16), nl);
                }
            }
        }
        // SAFETY: kernels were selected by runtime feature detection; `pm`
        // holds rows a..b (w rows of nl bytes) and the profile m tables.
        unsafe {
            #[cfg(target_arch = "x86_64")]
            match kernel {
                DiagKernel::Avx512 => x86::kadane_avx512(prof, &s.pm, a, q_lo, q_hi, &mut s.best),
                DiagKernel::Avx2 => x86::kadane_avx2(prof, &s.pm, a, q_lo, q_hi, &mut s.best),
                DiagKernel::Ssse3 => x86::kadane_ssse3(prof, &s.pm, a, q_lo, q_hi, &mut s.best),
                DiagKernel::Scalar => unreachable!(),
            }
        }
    }
    for (i, lane) in lanes.iter().enumerate() {
        let b = s.best[i];
        out[i] = if b > prof.limit {
            // possibly saturated: exact recomputation
            kadane_score(rows, lane.target, lane.diag, lane.lo, lane.hi)
        } else {
            b as i32
        };
    }
}

#[cfg(target_arch = "x86_64")]
mod x86 {
    use super::DiagProfile;
    use std::arch::x86_64::*;

    /// Transpose a 16 x 16 byte block: row `r` (16 bytes at `src[r]`)
    /// becomes column `r` of the output (16 rows at `dst + c·ds`), clamping
    /// every byte to at most the sentinel code.
    #[inline]
    #[target_feature(enable = "sse2")]
    pub unsafe fn transpose16(src: &[*const u8; 16], dst: *mut u8, ds: usize) {
        let cap = _mm_set1_epi8(super::SENTINEL as i8);
        let mut a = [_mm_setzero_si128(); 16];
        for (r, v) in a.iter_mut().enumerate() {
            *v = _mm_min_epu8(_mm_loadu_si128(src[r] as *const __m128i), cap);
        }
        // bytes → 16-bit pairs (rows 2k, 2k+1)
        let mut b = [_mm_setzero_si128(); 16];
        for k in 0..8 {
            b[2 * k] = _mm_unpacklo_epi8(a[2 * k], a[2 * k + 1]); // cols 0..7
            b[2 * k + 1] = _mm_unpackhi_epi8(a[2 * k], a[2 * k + 1]); // cols 8..15
        }
        // 16-bit pairs → 32-bit quads (rows 4j..4j+3)
        let mut c = [_mm_setzero_si128(); 16];
        for j in 0..4 {
            let (lo0, hi0, lo1, hi1) = (b[4 * j], b[4 * j + 1], b[4 * j + 2], b[4 * j + 3]);
            c[4 * j] = _mm_unpacklo_epi16(lo0, lo1); // cols 0..3
            c[4 * j + 1] = _mm_unpackhi_epi16(lo0, lo1); // cols 4..7
            c[4 * j + 2] = _mm_unpacklo_epi16(hi0, hi1); // cols 8..11
            c[4 * j + 3] = _mm_unpackhi_epi16(hi0, hi1); // cols 12..15
        }
        // quads → 64-bit octets (rows 8i..8i+7)
        let mut e = [_mm_setzero_si128(); 16];
        for i in 0..2 {
            for q in 0..4 {
                let x = c[8 * i + q]; // rows 8i..8i+3, col group q
                let y = c[8 * i + 4 + q]; // rows 8i+4..8i+7, col group q
                e[8 * i + 2 * q] = _mm_unpacklo_epi32(x, y); // cols 4q, 4q+1
                e[8 * i + 2 * q + 1] = _mm_unpackhi_epi32(x, y); // cols 4q+2, 4q+3
            }
        }
        // octets → full columns
        for k in 0..8 {
            let x = e[k]; // rows 0..7, cols 2k, 2k+1
            let y = e[8 + k]; // rows 8..15, same cols
            _mm_storeu_si128(dst.add((2 * k) * ds) as *mut __m128i, _mm_unpacklo_epi64(x, y));
            _mm_storeu_si128(dst.add((2 * k + 1) * ds) as *mut __m128i, _mm_unpackhi_epi64(x, y));
        }
    }

    #[target_feature(enable = "avx2")]
    pub unsafe fn kadane_avx2(p: &DiagProfile, buf: &[u8], a: usize, q_lo: usize, q_hi: usize, best: &mut [u8]) {
        const L: usize = 32;
        let vbias = _mm256_set1_epi8(p.bias as i8);
        let v15 = _mm256_set1_epi8(15);
        let mut cur = _mm256_setzero_si256();
        let mut bst = _mm256_setzero_si256();
        let tp = p.tables.as_ptr() as *const u8;
        let bp = buf.as_ptr();
        for i in q_lo..q_hi {
            let t = _mm256_loadu_si256(bp.add((i - a) * L) as *const __m256i);
            let lo = _mm256_broadcastsi128_si256(_mm_loadu_si128(tp.add(i * 32) as *const __m128i));
            let hi = _mm256_broadcastsi128_si256(_mm_loadu_si128(tp.add(i * 32 + 16) as *const __m128i));
            let s_lo = _mm256_shuffle_epi8(lo, t);
            let s_hi = _mm256_shuffle_epi8(hi, t);
            let sc = _mm256_blendv_epi8(s_lo, s_hi, _mm256_cmpgt_epi8(t, v15));
            cur = _mm256_subs_epu8(_mm256_adds_epu8(cur, sc), vbias);
            bst = _mm256_max_epu8(bst, cur);
        }
        _mm256_storeu_si256(best.as_mut_ptr() as *mut __m256i, bst);
    }

    #[target_feature(enable = "avx512bw")]
    pub unsafe fn kadane_avx512(p: &DiagProfile, buf: &[u8], a: usize, q_lo: usize, q_hi: usize, best: &mut [u8]) {
        const L: usize = 64;
        let vbias = _mm512_set1_epi8(p.bias as i8);
        let v15 = _mm512_set1_epi8(15);
        let mut cur = _mm512_setzero_si512();
        let mut bst = _mm512_setzero_si512();
        let tp = p.tables.as_ptr() as *const u8;
        let bp = buf.as_ptr();
        for i in q_lo..q_hi {
            let t = _mm512_loadu_si512(bp.add((i - a) * L) as *const _);
            let lo = _mm512_broadcast_i32x4(_mm_loadu_si128(tp.add(i * 32) as *const __m128i));
            let hi = _mm512_broadcast_i32x4(_mm_loadu_si128(tp.add(i * 32 + 16) as *const __m128i));
            let s_lo = _mm512_shuffle_epi8(lo, t);
            let s_hi = _mm512_shuffle_epi8(hi, t);
            let sc = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(t, v15), s_lo, s_hi);
            cur = _mm512_subs_epu8(_mm512_adds_epu8(cur, sc), vbias);
            bst = _mm512_max_epu8(bst, cur);
        }
        _mm512_storeu_si512(best.as_mut_ptr() as *mut _, bst);
    }

    #[target_feature(enable = "ssse3")]
    pub unsafe fn kadane_ssse3(p: &DiagProfile, buf: &[u8], a: usize, q_lo: usize, q_hi: usize, best: &mut [u8]) {
        const L: usize = 16;
        let vbias = _mm_set1_epi8(p.bias as i8);
        let v15 = _mm_set1_epi8(15);
        let mut cur = _mm_setzero_si128();
        let mut bst = _mm_setzero_si128();
        let tp = p.tables.as_ptr() as *const u8;
        let bp = buf.as_ptr();
        for i in q_lo..q_hi {
            let t = _mm_loadu_si128(bp.add((i - a) * L) as *const __m128i);
            let lo = _mm_loadu_si128(tp.add(i * 32) as *const __m128i);
            let hi = _mm_loadu_si128(tp.add(i * 32 + 16) as *const __m128i);
            let s_lo = _mm_shuffle_epi8(lo, t);
            let s_hi = _mm_shuffle_epi8(hi, t);
            // SSE2 blend: (hi & mask) | (lo & !mask)
            let mk = _mm_cmpgt_epi8(t, v15);
            let sc = _mm_or_si128(_mm_and_si128(mk, s_hi), _mm_andnot_si128(mk, s_lo));
            cur = _mm_subs_epu8(_mm_adds_epu8(cur, sc), vbias);
            bst = _mm_max_epu8(bst, cur);
        }
        _mm_storeu_si128(best.as_mut_ptr() as *mut __m128i, bst);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{ScoreTable, SeqScorer};
    use crate::matrix::SubstitutionMatrix;

    fn kernels() -> Vec<DiagKernel> {
        let mut v = vec![DiagKernel::Scalar];
        #[cfg(target_arch = "x86_64")]
        {
            if std::is_x86_feature_detected!("ssse3") {
                v.push(DiagKernel::Ssse3);
            }
            if std::is_x86_feature_detected!("avx2") {
                v.push(DiagKernel::Avx2);
            }
            if std::is_x86_feature_detected!("avx512bw") {
                v.push(DiagKernel::Avx512);
            }
        }
        v
    }

    #[test]
    fn segment_and_score_agree() {
        let mut rng = crate::evalue::SplitMix64::new(3);
        let m = SubstitutionMatrix::blosum62();
        for _ in 0..300 {
            let q: Vec<u8> = (0..1 + rng.below(80)).map(|_| rng.below(21) as u8).collect();
            let t: Vec<u8> = (0..1 + rng.below(120)).map(|_| rng.below(21) as u8).collect();
            let rows = ScoreTable::from_scorer(&SeqScorer { query: &q, mat: &m }).rows().to_vec();
            for d in -(q.len() as i64) - 2..t.len() as i64 + 2 {
                let (lo, hi) = diag_range(q.len(), t.len(), d);
                assert_eq!(kadane_score(&rows, &t, d, lo, hi), kadane_segment(&rows, &t, d, lo, hi).0);
                // sub-windows
                if hi > lo + 2 {
                    let (a, b) = (lo + 1, hi - 1);
                    assert_eq!(kadane_score(&rows, &t, d, a, b), kadane_segment(&rows, &t, d, a, b).0);
                }
            }
        }
    }

    #[test]
    fn simd_batches_match_scalar_reference() {
        // random and homologous diagonals, every lane count, lanes with empty
        // ranges, strong hits that exceed the u8 limit (exact recomputation)
        let mut rng = crate::evalue::SplitMix64::new(17);
        let mat = SubstitutionMatrix::blosum62();
        let mut sc = DiagScratch::default();
        for case in 0..60 {
            let qlen = 1 + rng.below(if case % 3 == 0 { 400 } else { 90 });
            let q: Vec<u8> = (0..qlen).map(|_| rng.below(21) as u8).collect();
            let rows = ScoreTable::from_scorer(&SeqScorer { query: &q, mat: &mat }).rows().to_vec();
            let prof = DiagProfile::build(&rows).expect("blosum rows fit u8");
            let targets: Vec<Vec<u8>> = (0..70)
                .map(|k| {
                    if k % 5 == 0 {
                        // a copy of the query with a few substitutions: long, high-scoring
                        let mut t = q.clone();
                        for _ in 0..qlen / 10 {
                            let p = rng.below(qlen);
                            t[p] = rng.below(20) as u8;
                        }
                        t
                    } else {
                        (0..1 + rng.below(300)).map(|_| rng.below(21) as u8).collect()
                    }
                })
                .collect();
            let lanes_all: Vec<Lane> = targets
                .iter()
                .enumerate()
                .map(|(k, t)| {
                    let diag = if k % 5 == 0 { 0 } else { rng.below(t.len() + qlen + 4) as i64 - qlen as i64 - 2 };
                    let (mut lo, mut hi) = diag_range(qlen, t.len(), diag);
                    if k % 3 == 1 && hi > lo + 4 {
                        // a window inside the diagonal
                        let a = lo + rng.below(hi - lo);
                        let b = a + rng.below(hi - a) + 1;
                        lo = a;
                        hi = b;
                    }
                    Lane { target: t, diag, lo, hi }
                })
                .collect();
            let want: Vec<i32> = lanes_all.iter().map(|l| kadane_score(&rows, l.target, l.diag, l.lo, l.hi)).collect();
            for k in kernels() {
                let nl = k.lanes();
                let mut out = vec![0i32; nl];
                for (ci, chunk) in lanes_all.chunks(nl).enumerate() {
                    score_batch(k, &prof, &rows, chunk, &mut out, &mut sc);
                    for (j, &o) in out[..chunk.len()].iter().enumerate() {
                        assert_eq!(o, want[ci * nl + j], "kernel {:?} case {} lane {}", k, case, ci * nl + j);
                    }
                }
            }
        }
    }

    #[test]
    fn profiles_with_wide_ranges_are_rejected() {
        let rows = vec![[0i32; ALPHABET_SIZE], {
            let mut r = [0i32; ALPHABET_SIZE];
            r[0] = 200;
            r[1] = -100;
            r
        }];
        assert!(DiagProfile::build(&rows).is_none());
    }
}
