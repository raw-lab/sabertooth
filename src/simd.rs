//! SIMD-vectorized local-alignment scoring (striped Smith–Waterman, Farrar 2007).
//!
//! Most candidates are rejected on score alone, so search first computes a
//! fast, traceback-free **score** with Farrar's striped layout over `i16`
//! lanes: AVX2 (16 lanes) when the CPU has it (runtime detection), otherwise
//! SSE2 (8 lanes, part of the x86-64 baseline). The same kernels, with end-cell
//! tracking, drive the canonical traceback in [`crate::align::align_canonical`].
//! `SABERTOOTH_SIMD=avx2|sse2|scalar` forces a kernel (testing, portability
//! checks); results are identical for every choice.
//!
//! ## Exactness contract
//!
//! [`QueryProfile::sw_score`] returns exactly the optimal Smith–Waterman–Gotoh
//! score (equal to `local_align(...).raw_score`). Differential tests compare it
//! with an independent three-matrix reference. The kernel falls back to an
//! exact scalar computation when the striped algorithm's assumptions fail:
//!
//! * `i16` saturation (score reaches 32767),
//! * a score table entry outside the `i16` range,
//! * `open <= extend`: Farrar's lazy-F termination test (`F − e > H − o`) is
//!   only sufficient when opening costs strictly more than extending,
//! * non-x86-64 targets.
//!
//! v0.2 used the striped path for `open == extend`, where the lazy-F loop stops
//! early (1,333 / 4,000 random pairs scored wrong at open = extend = 1); that
//! case now uses the exact scalar kernel. The lazy-F pass also refreshes `E`
//! from any `H` it raises, so the vector recurrence is literally Gotoh's. (v0.2
//! skipped that refresh; it never changed a score, because a deletion run
//! followed by an insertion run can always be replaced by the equal-cost
//! insertion-then-deletion corner, which the main loop does propagate — a
//! mutation test confirmed no score depends on it.)

use crate::align::{GapCosts, Scorer};
use std::sync::OnceLock;

const NPROF: usize = crate::alphabet::ALPHABET_SIZE; // 21 columns incl. X

/// Which kernel produced a score (for counters / tests).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kernel {
    Avx2I16,
    Sse2I16,
    /// `u8` first pass (exact below its headroom; reported when it decided).
    Avx2U8,
    Sse2U8,
    Scalar,
}

/// `SABERTOOTH_SIMD_U8=0` disables the `u8` first pass (results are identical).
fn u8_enabled() -> bool {
    static U8: OnceLock<bool> = OnceLock::new();
    *U8.get_or_init(|| std::env::var("SABERTOOTH_SIMD_U8").map(|v| v != "0").unwrap_or(true))
}

/// Kernel family to build query profiles for.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum KernelChoice {
    /// AVX2 if the CPU supports it, else SSE2 (default).
    Auto,
    Avx2,
    Sse2,
    Scalar,
}

/// The process-wide kernel choice: `SABERTOOTH_SIMD` (`avx2`, `sse2`,
/// `scalar`, `auto`) or automatic detection. An unavailable request falls
/// back to the best available kernel.
pub fn kernel_choice() -> KernelChoice {
    static CHOICE: OnceLock<KernelChoice> = OnceLock::new();
    *CHOICE.get_or_init(|| match std::env::var("SABERTOOTH_SIMD").as_deref() {
        Ok("avx2") => KernelChoice::Avx2,
        Ok("sse2") => KernelChoice::Sse2,
        Ok("scalar") => KernelChoice::Scalar,
        _ => KernelChoice::Auto,
    })
}

/// Name of the kernel family actually used for profiles built with `choice`.
pub fn kernel_name(choice: KernelChoice) -> &'static str {
    let u8s = if u8_enabled() { " (+u8 first pass)" } else { "" };
    match resolve(choice) {
        Kernel::Avx2I16 | Kernel::Avx2U8 => if u8_enabled() { "avx2-i16x16 (+u8x32 first pass)" } else { "avx2-i16x16" },
        Kernel::Sse2I16 | Kernel::Sse2U8 => if u8s.is_empty() { "sse2-i16x8" } else { "sse2-i16x8 (+u8x16 first pass)" },
        Kernel::Scalar => "scalar",
    }
}

fn resolve(choice: KernelChoice) -> Kernel {
    #[cfg(target_arch = "x86_64")]
    {
        let avx2 = std::is_x86_feature_detected!("avx2");
        match choice {
            KernelChoice::Scalar => Kernel::Scalar,
            KernelChoice::Sse2 => Kernel::Sse2I16,
            KernelChoice::Avx2 | KernelChoice::Auto => {
                if avx2 {
                    Kernel::Avx2I16
                } else {
                    Kernel::Sse2I16
                }
            }
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let _ = choice;
        Kernel::Scalar
    }
}

/// A query (sequence or profile) prepared for fast repeated scoring against many
/// targets. Build once, score many times.
pub struct QueryProfile {
    m: usize,
    /// Per-column score of each residue: `col_scores[col][res]`.
    col_scores: Vec<[i32; NPROF]>,
    kernel: Kernel,
    #[cfg(target_arch = "x86_64")]
    striped: Option<x86::Striped>,
    /// Biased `u8` profile for the first scoring pass (twice the lanes);
    /// `None` when scores are too wide for useful `u8` headroom.
    #[cfg(target_arch = "x86_64")]
    striped8: Option<x86::StripedU8>,
}

/// Reusable SIMD working memory (one per worker thread).
#[derive(Default)]
pub struct SimdScratch {
    #[cfg(target_arch = "x86_64")]
    x86: x86::Scratch,
    h_prev: Vec<i32>,
    h_cur: Vec<i32>,
    f_col: Vec<i32>,
}

impl SimdScratch {
    pub fn new() -> Self {
        Self::default()
    }
}

impl QueryProfile {
    /// Build a query profile from any [`Scorer`] (sequence-, PSSM- or
    /// table-backed).
    pub fn build<S: Scorer>(scorer: &S) -> QueryProfile {
        let m = scorer.query_len();
        let mut col_scores = vec![[0i32; NPROF]; m];
        for (c, row) in col_scores.iter_mut().enumerate() {
            for (r, slot) in row.iter_mut().enumerate() {
                *slot = scorer.score(c, r as u8);
            }
        }
        Self::from_rows(col_scores)
    }

    /// Build from precomputed score rows (process-wide kernel choice).
    pub fn from_rows(col_scores: Vec<[i32; NPROF]>) -> QueryProfile {
        Self::from_rows_with(col_scores, kernel_choice())
    }

    /// Build from precomputed score rows for an explicit kernel family (with
    /// the `u8` first pass unless `SABERTOOTH_SIMD_U8=0`).
    pub fn from_rows_with(col_scores: Vec<[i32; NPROF]>, choice: KernelChoice) -> QueryProfile {
        Self::from_rows_opts(col_scores, choice, u8_enabled())
    }

    /// Build with explicit control of the `u8` first pass (tests, benchmarks).
    pub fn from_rows_opts(col_scores: Vec<[i32; NPROF]>, choice: KernelChoice, use_u8: bool) -> QueryProfile {
        let m = col_scores.len();
        let kernel = resolve(choice);
        #[cfg(target_arch = "x86_64")]
        {
            let (striped, striped8) = match kernel {
                Kernel::Avx2I16 | Kernel::Avx2U8 => (
                    x86::Striped::build(&col_scores, m, 16),
                    if use_u8 { x86::StripedU8::build(&col_scores, m, 32) } else { None },
                ),
                Kernel::Sse2I16 | Kernel::Sse2U8 => (
                    x86::Striped::build(&col_scores, m, 8),
                    if use_u8 { x86::StripedU8::build(&col_scores, m, 16) } else { None },
                ),
                Kernel::Scalar => (None, None),
            };
            QueryProfile { m, col_scores, kernel, striped, striped8 }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            let _ = use_u8;
            QueryProfile { m, col_scores, kernel }
        }
    }

    /// Whether the `u8` first pass is available for this profile.
    pub fn has_u8(&self) -> bool {
        #[cfg(target_arch = "x86_64")]
        {
            self.striped8.is_some()
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            false
        }
    }

    pub fn query_len(&self) -> usize {
        self.m
    }

    /// Kernel family this profile was built for.
    pub fn kernel(&self) -> Kernel {
        self.kernel
    }

    /// Score rows this profile was built from.
    pub fn rows(&self) -> &[[i32; NPROF]] {
        &self.col_scores
    }

    #[cfg(target_arch = "x86_64")]
    #[inline]
    fn simd_ok(&self, gaps: GapCosts) -> bool {
        gaps.open > gaps.extend && gaps.open <= i16::MAX as i32 && gaps.extend >= 0
    }

    /// Optimal Smith–Waterman–Gotoh local-alignment **score** (no traceback).
    pub fn sw_score(&self, target: &[u8], gaps: GapCosts) -> i32 {
        let mut s = SimdScratch::new();
        self.sw_score_with(target, gaps, &mut s).0
    }

    /// [`sw_score`](Self::sw_score) with reusable scratch; also reports which
    /// kernel produced the score.
    pub fn sw_score_with(&self, target: &[u8], gaps: GapCosts, s: &mut SimdScratch) -> (i32, Kernel) {
        if self.m == 0 || target.is_empty() {
            return (0, Kernel::Scalar);
        }
        #[cfg(target_arch = "x86_64")]
        {
            if self.simd_ok(gaps) {
                // SAFETY (all kernel calls): the family was chosen by runtime
                // feature detection; SSE2 is part of the x86-64 baseline.
                if let Some(s8) = &self.striped8 {
                    if gaps.open <= u8::MAX as i32 {
                        let avx2 = self.kernel == Kernel::Avx2I16;
                        let r = unsafe {
                            if avx2 {
                                s8.score_avx2(target, gaps, &mut s.x86)
                            } else {
                                s8.score_sse2(target, gaps, &mut s.x86)
                            }
                        };
                        if let Some(score) = r {
                            return (score, if avx2 { Kernel::Avx2U8 } else { Kernel::Sse2U8 });
                        }
                        // u8 headroom exceeded → i16
                    }
                }
                if let Some(st) = &self.striped {
                    let r = unsafe {
                        match self.kernel {
                            Kernel::Avx2I16 => st.score_avx2(target, gaps, &mut s.x86),
                            _ => st.score_sse2(target, gaps, &mut s.x86),
                        }
                    };
                    if let Some(score) = r {
                        return (score, self.kernel);
                    }
                    // saturated → exact scalar below
                }
            }
        }
        (scalar_gotoh_score(&self.col_scores, self.m, target, gaps, s), Kernel::Scalar)
    }

    /// Optimal score **and end cell** `(score, query_end, target_end)` (0-based,
    /// inclusive), or `None` if no positive-scoring alignment exists.
    ///
    /// The end cell is the first optimum in **column-major order** (smallest
    /// target position, then smallest query position) — the order the striped
    /// kernel visits columns in. The scalar fallback applies the same rule, so
    /// the result never depends on which kernel ran.
    pub fn sw_score_end_with(&self, target: &[u8], gaps: GapCosts, s: &mut SimdScratch) -> Option<(i32, usize, usize)> {
        if self.m == 0 || target.is_empty() {
            return None;
        }
        #[cfg(target_arch = "x86_64")]
        {
            if self.simd_ok(gaps) {
                if let Some(st) = &self.striped {
                    // SAFETY: as in `sw_score_with`.
                    let r = unsafe {
                        match self.kernel {
                            Kernel::Avx2I16 => st.score_end_avx2(target, gaps, self.m, &mut s.x86),
                            _ => st.score_end_sse2(target, gaps, self.m, &mut s.x86),
                        }
                    };
                    match r {
                        Some(Some(x)) => return Some(x),
                        Some(None) => return None, // best score 0
                        None => {}                 // saturated → exact scalar
                    }
                }
            }
        }
        scalar_gotoh_score_end(&self.col_scores, self.m, target, gaps, s)
    }

    /// Exact scalar score + end (the fallback path), exposed for tests.
    pub fn scalar_score_end(&self, target: &[u8], gaps: GapCosts) -> Option<(i32, usize, usize)> {
        let mut s = SimdScratch::new();
        scalar_gotoh_score_end(&self.col_scores, self.m, target, gaps, &mut s)
    }

    /// Exact scalar score (the fallback path), exposed for differential tests.
    pub fn scalar_score(&self, target: &[u8], gaps: GapCosts) -> i32 {
        let mut s = SimdScratch::new();
        scalar_gotoh_score(&self.col_scores, self.m, target, gaps, &mut s)
    }
}

/// Traceback-free scalar Gotoh local-alignment score; O(m·n) time, O(n) memory.
fn scalar_gotoh_score(col_scores: &[[i32; NPROF]], m: usize, target: &[u8], gaps: GapCosts, s: &mut SimdScratch) -> i32 {
    scalar_gotoh_score_end(col_scores, m, target, gaps, s).map(|x| x.0).unwrap_or(0)
}

/// Scalar Gotoh score with the column-major-first optimal end cell.
fn scalar_gotoh_score_end(
    col_scores: &[[i32; NPROF]],
    m: usize,
    target: &[u8],
    gaps: GapCosts,
    s: &mut SimdScratch,
) -> Option<(i32, usize, usize)> {
    let n = target.len();
    let neg = i32::MIN / 4;
    s.h_prev.clear();
    s.h_prev.resize(n + 1, 0);
    s.h_cur.clear();
    s.h_cur.resize(n + 1, 0);
    s.f_col.clear();
    s.f_col.resize(n + 1, neg);
    let (mut best, mut bi, mut bj) = (0i32, 0usize, usize::MAX);
    for (i, row) in col_scores.iter().take(m).enumerate() {
        let mut e = neg;
        s.h_cur[0] = 0;
        for j in 1..=n {
            let tr = target[j - 1] as usize;
            let sc = if tr < NPROF { row[tr] } else { 0 };
            let diag = s.h_prev[j - 1] + sc;
            e = (s.h_cur[j - 1] - gaps.open).max(e - gaps.extend);
            let f = (s.h_prev[j] - gaps.open).max(s.f_col[j] - gaps.extend);
            s.f_col[j] = f;
            let v = diag.max(e).max(f).max(0);
            s.h_cur[j] = v;
            // rows are visited in increasing i, so an equal score in an
            // earlier column is the only tie that moves the end
            if v > best || (v == best && v > 0 && j - 1 < bj) {
                best = v;
                bi = i;
                bj = j - 1;
            }
        }
        std::mem::swap(&mut s.h_prev, &mut s.h_cur);
    }
    if best > 0 {
        Some((best, bi, bj))
    } else {
        None
    }
}

#[cfg(target_arch = "x86_64")]
mod x86 {
    use super::NPROF;
    use crate::align::GapCosts;
    use std::arch::x86_64::*;

    /// Reusable striped DP vectors (one set per register width).
    #[derive(Default)]
    pub struct Scratch {
        h_store: Vec<__m128i>,
        h_load: Vec<__m128i>,
        e_arr: Vec<__m128i>,
        best_col: Vec<__m128i>,
        h_store_y: Vec<__m256i>,
        h_load_y: Vec<__m256i>,
        e_arr_y: Vec<__m256i>,
        best_col_y: Vec<__m256i>,
    }

    /// Biased `u8` striped profile (`score + bias`, all ≥ 0) for the first
    /// scoring pass, with `lanes` = 16 (SSE2) or 32 (AVX2). Cells are exact
    /// while `H + max_profile ≤ 255`; the kernel reports `None` otherwise.
    pub struct StripedU8 {
        lanes: usize,
        seg_len: usize,
        bias: u8,
        /// largest observed `H` for which every cell was exact
        limit: u8,
        data: Vec<Block8>,
    }

    #[repr(C, align(32))]
    #[derive(Clone, Copy)]
    struct Block8([u8; 32]);

    /// Minimum `u8` headroom (`255 − max(score + bias)`) for the pass to pay off.
    const MIN_U8_LIMIT: i32 = 96;

    impl StripedU8 {
        pub fn build(col_scores: &[[i32; NPROF]], m: usize, lanes: usize) -> Option<StripedU8> {
            if m == 0 {
                return None;
            }
            let mut lo = 0i32; // residues outside the alphabet score 0
            let mut hi = 0i32;
            for row in col_scores {
                for &v in row.iter() {
                    lo = lo.min(v);
                    hi = hi.max(v);
                }
            }
            let bias = -lo;
            let limit = 255 - (hi + bias);
            if bias > 255 || limit < MIN_U8_LIMIT {
                return None;
            }
            let seg_len = m.div_ceil(lanes);
            let total = (NPROF + 1) * seg_len * lanes;
            let mut data = vec![Block8([0; 32]); total.div_ceil(32)];
            for res in 0..=NPROF {
                for i in 0..seg_len {
                    for l in 0..lanes {
                        let col = i + l * seg_len;
                        let v = if col < m && res < NPROF { col_scores[col][res] } else { 0 };
                        let k = (res * seg_len + i) * lanes + l;
                        data[k / 32].0[k % 32] = (v + bias) as u8;
                    }
                }
            }
            Some(StripedU8 {
                lanes,
                seg_len,
                bias: bias as u8,
                limit: limit as u8,
                data,
            })
        }

        #[inline]
        fn row(&self, res: u8) -> *const u8 {
            let r = (res as usize).min(NPROF);
            // SAFETY: r ≤ NPROF and data holds NPROF + 1 rows
            unsafe { (self.data.as_ptr() as *const u8).add(r * self.seg_len * self.lanes) }
        }
    }

    #[inline]
    #[target_feature(enable = "sse2")]
    unsafe fn shl1_u8_sse2(v: __m128i) -> __m128i {
        _mm_slli_si128::<1>(v)
    }

    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn shl1_u8_avx2(v: __m256i) -> __m256i {
        _mm256_alignr_epi8::<15>(v, _mm256_permute2x128_si256::<0x08>(v, v))
    }

    #[inline]
    #[target_feature(enable = "sse2")]
    unsafe fn hmax_u8_sse2(v: __m128i) -> u8 {
        let mut m = _mm_max_epu8(v, _mm_srli_si128::<8>(v));
        m = _mm_max_epu8(m, _mm_srli_si128::<4>(m));
        m = _mm_max_epu8(m, _mm_srli_si128::<2>(m));
        m = _mm_max_epu8(m, _mm_srli_si128::<1>(m));
        (_mm_cvtsi128_si32(m) & 0xFF) as u8
    }

    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn hmax_u8_avx2(v: __m256i) -> u8 {
        let lo = _mm256_castsi256_si128(v);
        let hi = _mm256_extracti128_si256::<1>(v);
        hmax_u8_sse2(_mm_max_epu8(lo, hi))
    }

    /// The `u8` score kernel for one register width (saturating unsigned
    /// arithmetic: `H = max(H_diag + p − bias, E, F, 0)`; values below 0 clamp
    /// to 0, which cannot change any `H ≥ 0`, so cells are exact until the
    /// headroom check fails).
    macro_rules! striped_u8_kernel {
        (
            $feat:literal, $V:ty, $L:expr, $full:expr, $score:ident, [$hs:ident, $hl:ident, $ea:ident],
            $setzero:ident, $set1:ident, $adds:ident, $subs:ident, $max:ident, $cmpeq:ident,
            $movemask:ident, $load:ident, $shl1:ident, $hmax:ident
        ) => {
            impl StripedU8 {
                /// `None` when the headroom was exceeded (rescore with `i16`).
                /// Requires `open > extend` and `open ≤ 255`.
                ///
                /// SAFETY: the CPU must support the target feature.
                #[target_feature(enable = $feat)]
                pub unsafe fn $score(&self, target: &[u8], gaps: GapCosts, s: &mut Scratch) -> Option<i32> {
                    debug_assert_eq!(self.lanes, $L);
                    let seg_len = self.seg_len;
                    let v_open = $set1(gaps.open as u8 as i8);
                    let v_extend = $set1(gaps.extend as u8 as i8);
                    let v_bias = $set1(self.bias as i8);
                    let v_zero = $setzero();
                    s.$hs.clear();
                    s.$hs.resize(seg_len, v_zero);
                    s.$hl.clear();
                    s.$hl.resize(seg_len, v_zero);
                    s.$ea.clear();
                    s.$ea.resize(seg_len, v_zero);
                    let h_store = &mut s.$hs;
                    let h_load = &mut s.$hl;
                    let e_arr = &mut s.$ea;
                    let mut v_max = v_zero;
                    let v_limit = $set1(self.limit as i8);
                    for &tb in target {
                        let prof = self.row(tb);
                        // u8: 0 is the floor, so the shifted-in 0 is a valid "−∞" for F
                        let mut v_h = $shl1(h_store[seg_len - 1]);
                        let mut v_f = v_zero;
                        std::mem::swap(h_store, h_load);
                        let hsp = h_store.as_mut_ptr();
                        let hlp = h_load.as_ptr();
                        let eap = e_arr.as_mut_ptr();
                        for i in 0..seg_len {
                            let v_p = $load(prof.add(i * $L) as *const $V);
                            let e = *eap.add(i);
                            v_h = $subs($adds(v_h, v_p), v_bias);
                            v_h = $max(v_h, e);
                            v_h = $max(v_h, v_f);
                            v_max = $max(v_max, v_h);
                            *hsp.add(i) = v_h;
                            let v_h_open = $subs(v_h, v_open);
                            *eap.add(i) = $max($subs(e, v_extend), v_h_open);
                            v_f = $max($subs(v_f, v_extend), v_h_open);
                            v_h = *hlp.add(i);
                        }
                        'lazy: for _ in 0..$L {
                            v_f = $shl1(v_f);
                            for i in 0..seg_len {
                                let v_h = $max(*hsp.add(i), v_f);
                                v_max = $max(v_max, v_h);
                                *hsp.add(i) = v_h;
                                let v_h_open = $subs(v_h, v_open);
                                *eap.add(i) = $max(*eap.add(i), v_h_open);
                                v_f = $subs(v_f, v_extend);
                                // stop when F ≤ H − open in every lane (unsigned)
                                if $movemask($cmpeq($subs(v_f, v_h_open), v_zero)) as u32 == $full {
                                    break 'lazy;
                                }
                            }
                        }
                        // headroom exceeded: later cells may have saturated
                        if $movemask($cmpeq($subs(v_max, v_limit), v_zero)) as u32 != $full {
                            return None;
                        }
                    }
                    Some($hmax(v_max) as i32)
                }
            }
        };
    }

    striped_u8_kernel!(
        "sse2", __m128i, 16, 0xFFFFu32, score_sse2, [h_store, h_load, e_arr],
        _mm_setzero_si128, _mm_set1_epi8, _mm_adds_epu8, _mm_subs_epu8, _mm_max_epu8, _mm_cmpeq_epi8,
        _mm_movemask_epi8, _mm_load_si128, shl1_u8_sse2, hmax_u8_sse2
    );

    striped_u8_kernel!(
        "avx2", __m256i, 32, 0xFFFF_FFFFu32, score_avx2, [h_store_y, h_load_y, e_arr_y],
        _mm256_setzero_si256, _mm256_set1_epi8, _mm256_adds_epu8, _mm256_subs_epu8, _mm256_max_epu8, _mm256_cmpeq_epi8,
        _mm256_movemask_epi8, _mm256_load_si256, shl1_u8_avx2, hmax_u8_avx2
    );

    /// Striped `i16` query profile (Farrar layout) for `lanes` = 8 (SSE2) or
    /// 16 (AVX2). For segment `i` and lane `l` the represented query column is
    /// `i + l·seg_len` (0 past the end). Row `NPROF` is all zeros (target
    /// residues outside the alphabet score 0).
    pub struct Striped {
        lanes: usize,
        seg_len: usize,
        /// flat `i16` profile stored in 32-byte-aligned blocks so every vector
        /// load is aligned to its own width (no cache-line-splitting loads)
        data: Vec<Block>,
    }

    #[repr(C, align(32))]
    #[derive(Clone, Copy)]
    struct Block([i16; 16]);

    impl Striped {
        /// `None` when the query is empty or a score does not fit `i16`
        /// (then the exact scalar path is used).
        pub fn build(col_scores: &[[i32; NPROF]], m: usize, lanes: usize) -> Option<Striped> {
            if m == 0 {
                return None;
            }
            let seg_len = m.div_ceil(lanes);
            let total = (NPROF + 1) * seg_len * lanes;
            let mut flat = vec![0i16; total];
            for res in 0..NPROF {
                for i in 0..seg_len {
                    for l in 0..lanes {
                        let col = i + l * seg_len;
                        let s = if col < m { col_scores[col][res] } else { 0 };
                        if s > i16::MAX as i32 || s < i16::MIN as i32 {
                            return None;
                        }
                        flat[(res * seg_len + i) * lanes + l] = s as i16;
                    }
                }
            }
            let mut data = vec![Block([0; 16]); total.div_ceil(16)];
            for (k, &v) in flat.iter().enumerate() {
                data[k / 16].0[k % 16] = v;
            }
            Some(Striped { lanes, seg_len, data })
        }

        #[inline]
        fn row(&self, res: u8) -> *const i16 {
            let r = (res as usize).min(NPROF);
            // SAFETY: r ≤ NPROF and data holds NPROF + 1 rows of
            // seg_len·lanes values; blocks are contiguous `[i16; 16]`
            unsafe { (self.data.as_ptr() as *const i16).add(r * self.seg_len * self.lanes) }
        }
    }

    #[inline]
    #[target_feature(enable = "sse2")]
    unsafe fn shl1_sse2(v: __m128i) -> __m128i {
        _mm_slli_si128::<2>(v)
    }

    /// Shift all 16 `i16` lanes up by one (lane 0 ← 0), across the 128-bit halves.
    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn shl1_avx2(v: __m256i) -> __m256i {
        _mm256_alignr_epi8::<14>(v, _mm256_permute2x128_si256::<0x08>(v, v))
    }

    /// `i16::MIN` in lane 0, 0 elsewhere (OR-ed into a lane-shifted vector).
    #[inline]
    #[target_feature(enable = "sse2")]
    unsafe fn lane0_min_sse2() -> __m128i {
        _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, i16::MIN)
    }

    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn lane0_min_avx2() -> __m256i {
        _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, i16::MIN)
    }

    /// Horizontal max of 8 packed i16.
    #[inline]
    #[target_feature(enable = "sse2")]
    unsafe fn hmax_sse2(v: __m128i) -> i16 {
        let mut m = v;
        m = _mm_max_epi16(m, _mm_shuffle_epi32::<0b01_00_11_10>(m));
        m = _mm_max_epi16(m, _mm_shuffle_epi32::<0b00_00_00_01>(m));
        m = _mm_max_epi16(m, _mm_shufflelo_epi16::<0b00_00_00_01>(m));
        _mm_extract_epi16::<0>(m) as i16
    }

    /// Horizontal max of 16 packed i16.
    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn hmax_avx2(v: __m256i) -> i16 {
        let lo = _mm256_castsi256_si128(v);
        let hi = _mm256_extracti128_si256::<1>(v);
        hmax_sse2(_mm_max_epi16(lo, hi))
    }

    /// Generates the score and score+end kernels for one register width from a
    /// single body, so SSE2 and AVX2 cannot diverge.
    macro_rules! striped_kernels {
        (
            $feat:literal, $V:ty, $L:expr, $score:ident, $score_end:ident,
            [$hs:ident, $hl:ident, $ea:ident, $bc:ident],
            $setzero:ident, $set1:ident, $adds:ident, $subs:ident, $max:ident,
            $cmpgt:ident, $cmpeq:ident, $movemask:ident, $load:ident, $storeu:ident,
            $shl1:ident, $hmax:ident, $or:ident, $lane0_min:ident
        ) => {
            impl Striped {
                /// Striped local SW score. `None` if the score saturates `i16`
                /// (caller falls back to scalar). Requires `open > extend`.
                ///
                /// SAFETY: the CPU must support the target feature.
                #[target_feature(enable = $feat)]
                pub unsafe fn $score(&self, target: &[u8], gaps: GapCosts, s: &mut Scratch) -> Option<i32> {
                    debug_assert_eq!(self.lanes, $L);
                    let seg_len = self.seg_len;
                    let v_open = $set1(gaps.open as i16);
                    let v_extend = $set1(gaps.extend as i16);
                    let v_zero = $setzero();
                    let v_neg = $set1(i16::MIN);
                    let v_lane0 = $lane0_min();
                    s.$hs.clear();
                    s.$hs.resize(seg_len, v_zero);
                    s.$hl.clear();
                    s.$hl.resize(seg_len, v_zero);
                    s.$ea.clear();
                    s.$ea.resize(seg_len, v_zero);
                    let h_store = &mut s.$hs;
                    let h_load = &mut s.$hl;
                    let e_arr = &mut s.$ea;
                    let mut v_max = v_zero;
                    let sat = $set1(i16::MAX);
                    for &tb in target {
                        let prof = self.row(tb);
                        let mut v_h = $shl1(h_store[seg_len - 1]);
                        // F entering the first segment is not known yet: −∞
                        // (the lazy pass supplies the real values)
                        let mut v_f = v_neg;
                        std::mem::swap(h_store, h_load);
                        // raw pointers: all three arrays hold exactly seg_len vectors
                        let hs = h_store.as_mut_ptr();
                        let hl = h_load.as_ptr();
                        let ea = e_arr.as_mut_ptr();
                        for i in 0..seg_len {
                            let v_p = $load(prof.add(i * $L) as *const $V);
                            let e = *ea.add(i);
                            v_h = $adds(v_h, v_p);
                            v_h = $max(v_h, e);
                            v_h = $max(v_h, v_f);
                            v_h = $max(v_h, v_zero);
                            v_max = $max(v_max, v_h);
                            *hs.add(i) = v_h;
                            let v_h_open = $subs(v_h, v_open);
                            *ea.add(i) = $max($subs(e, v_extend), v_h_open);
                            v_f = $max($subs(v_f, v_extend), v_h_open);
                            v_h = *hl.add(i);
                        }
                        // Lazy-F: propagate F across stripe boundaries. When it
                        // raises H, E for the next column must see the new H.
                        'lazy: for _ in 0..$L {
                            // lane l receives lane l−1; lane 0 (query start) has
                            // no F above it: −∞, not 0 — a 0 there would keep
                            // the termination test below true on most columns
                            v_f = $or($shl1(v_f), v_lane0);
                            for i in 0..seg_len {
                                let v_h = $max(*hs.add(i), v_f);
                                v_max = $max(v_max, v_h);
                                *hs.add(i) = v_h;
                                let v_h_open = $subs(v_h, v_open);
                                *ea.add(i) = $max(*ea.add(i), v_h_open);
                                v_f = $subs(v_f, v_extend);
                                // stop when F cannot improve any H via extension
                                if $movemask($cmpgt(v_f, v_h_open)) == 0 {
                                    break 'lazy;
                                }
                            }
                        }
                    }
                    if $movemask($cmpeq(v_max, sat)) != 0 {
                        None
                    } else {
                        Some($hmax(v_max) as i32)
                    }
                }

                /// Score plus the column-major-first optimal end cell. `None` =
                /// saturated (use the scalar path); `Some(None)` = best score 0.
                ///
                /// SAFETY: the CPU must support the target feature.
                #[target_feature(enable = $feat)]
                pub unsafe fn $score_end(
                    &self,
                    target: &[u8],
                    gaps: GapCosts,
                    m: usize,
                    s: &mut Scratch,
                ) -> Option<Option<(i32, usize, usize)>> {
                    debug_assert_eq!(self.lanes, $L);
                    let seg_len = self.seg_len;
                    let v_open = $set1(gaps.open as i16);
                    let v_extend = $set1(gaps.extend as i16);
                    let v_zero = $setzero();
                    let v_neg = $set1(i16::MIN);
                    let v_lane0 = $lane0_min();
                    s.$hs.clear();
                    s.$hs.resize(seg_len, v_zero);
                    s.$hl.clear();
                    s.$hl.resize(seg_len, v_zero);
                    s.$ea.clear();
                    s.$ea.resize(seg_len, v_zero);
                    s.$bc.clear();
                    s.$bc.resize(seg_len, v_zero);
                    let h_store = &mut s.$hs;
                    let h_load = &mut s.$hl;
                    let e_arr = &mut s.$ea;
                    let best_col = &mut s.$bc;
                    let sat = $set1(i16::MAX);
                    let mut v_max = v_zero;
                    let mut best: i16 = 0;
                    let mut best_j = usize::MAX;
                    for (j, &tb) in target.iter().enumerate() {
                        let prof = self.row(tb);
                        let mut v_h = $shl1(h_store[seg_len - 1]);
                        // F entering the first segment is not known yet: −∞
                        // (the lazy pass supplies the real values)
                        let mut v_f = v_neg;
                        let mut v_col = v_zero;
                        std::mem::swap(h_store, h_load);
                        let hs = h_store.as_mut_ptr();
                        let hl = h_load.as_ptr();
                        let ea = e_arr.as_mut_ptr();
                        for i in 0..seg_len {
                            let v_p = $load(prof.add(i * $L) as *const $V);
                            let e = *ea.add(i);
                            v_h = $adds(v_h, v_p);
                            v_h = $max(v_h, e);
                            v_h = $max(v_h, v_f);
                            v_h = $max(v_h, v_zero);
                            v_col = $max(v_col, v_h);
                            *hs.add(i) = v_h;
                            let v_h_open = $subs(v_h, v_open);
                            *ea.add(i) = $max($subs(e, v_extend), v_h_open);
                            v_f = $max($subs(v_f, v_extend), v_h_open);
                            v_h = *hl.add(i);
                        }
                        'lazy: for _ in 0..$L {
                            // lane l receives lane l−1; lane 0 (query start) has
                            // no F above it: −∞, not 0 — a 0 there would keep
                            // the termination test below true on most columns
                            v_f = $or($shl1(v_f), v_lane0);
                            for i in 0..seg_len {
                                let v_h = $max(*hs.add(i), v_f);
                                v_col = $max(v_col, v_h);
                                *hs.add(i) = v_h;
                                let v_h_open = $subs(v_h, v_open);
                                *ea.add(i) = $max(*ea.add(i), v_h_open);
                                v_f = $subs(v_f, v_extend);
                                if $movemask($cmpgt(v_f, v_h_open)) == 0 {
                                    break 'lazy;
                                }
                            }
                        }
                        v_max = $max(v_max, v_col);
                        // a strictly better column: remember it (first wins ties)
                        if $movemask($cmpgt(v_col, $set1(best))) != 0 {
                            best = $hmax(v_col);
                            best_j = j;
                            best_col.copy_from_slice(h_store);
                        }
                    }
                    if $movemask($cmpeq(v_max, sat)) != 0 {
                        return None;
                    }
                    if best <= 0 {
                        return Some(None);
                    }
                    // smallest query position in the best column holding the optimum
                    let mut buf = [0i16; $L];
                    for q in 0..m {
                        let (seg, lane) = (q % seg_len, q / seg_len);
                        $storeu(buf.as_mut_ptr() as *mut $V, best_col[seg]);
                        if buf[lane] == best {
                            return Some(Some((best as i32, q, best_j)));
                        }
                    }
                    unreachable!("optimum not found in its column")
                }
            }
        };
    }

    striped_kernels!(
        "sse2", __m128i, 8, score_sse2, score_end_sse2,
        [h_store, h_load, e_arr, best_col],
        _mm_setzero_si128, _mm_set1_epi16, _mm_adds_epi16, _mm_subs_epi16, _mm_max_epi16,
        _mm_cmpgt_epi16, _mm_cmpeq_epi16, _mm_movemask_epi8, _mm_load_si128, _mm_storeu_si128,
        shl1_sse2, hmax_sse2, _mm_or_si128, lane0_min_sse2
    );

    striped_kernels!(
        "avx2", __m256i, 16, score_avx2, score_end_avx2,
        [h_store_y, h_load_y, e_arr_y, best_col_y],
        _mm256_setzero_si256, _mm256_set1_epi16, _mm256_adds_epi16, _mm256_subs_epi16, _mm256_max_epi16,
        _mm256_cmpgt_epi16, _mm256_cmpeq_epi16, _mm256_movemask_epi8, _mm256_load_si256, _mm256_storeu_si256,
        shl1_avx2, hmax_avx2, _mm256_or_si256, lane0_min_avx2
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{local_align, GapCosts, ProfileScorer, ScoreTable, SeqScorer};
    use crate::alphabet::build_aa2num;
    use crate::matrix::SubstitutionMatrix;
    use crate::profile::{Profile, PseudoCountParams};

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
        fn range(&mut self, n: usize) -> usize {
            (self.next() % n as u64) as usize
        }
    }

    fn rand_seq(rng: &mut Rng, len: usize, alpha: usize) -> Vec<u8> {
        (0..len).map(|_| rng.range(alpha) as u8).collect()
    }

    const CHOICES: [KernelChoice; 3] = [KernelChoice::Avx2, KernelChoice::Sse2, KernelChoice::Scalar];

    /// the same rows prepared for every kernel family
    fn all_kernels(rows: &[[i32; NPROF]]) -> Vec<QueryProfile> {
        CHOICES.iter().map(|&c| QueryProfile::from_rows_with(rows.to_vec(), c)).collect()
    }

    #[test]
    fn u8_first_pass_is_exact_and_falls_back_on_headroom() {
        let m = SubstitutionMatrix::blosum62();
        let mut rng = Rng(99);
        let mut used_u8 = 0;
        let mut fell_back = 0;
        for round in 0..3000 {
            let ql = 1 + rng.range(140);
            let tl = 1 + rng.range(200);
            let alpha = 2 + rng.range(19);
            let q = rand_seq(&mut rng, ql, alpha);
            // every 5th target embeds the query: high scores exceed the u8 headroom
            let mut t = rand_seq(&mut rng, tl, alpha);
            if round % 5 == 0 {
                let p = rng.range(t.len() + 1);
                t.splice(p..p, q.iter().copied());
            }
            let base = QueryProfile::build(&SeqScorer { query: &q, mat: &m });
            let gaps = [GapCosts::sequence_default(), GapCosts { open: 6, extend: 2 }][round % 2];
            let want = reference(&base.col_scores, &t, gaps.open, gaps.extend);
            for choice in [KernelChoice::Avx2, KernelChoice::Sse2] {
                let with = QueryProfile::from_rows_opts(base.col_scores.clone(), choice, true);
                let without = QueryProfile::from_rows_opts(base.col_scores.clone(), choice, false);
                assert!(with.has_u8() && !without.has_u8());
                let mut sc = SimdScratch::new();
                let (a, ka) = with.sw_score_with(&t, gaps, &mut sc);
                let (b, kb) = without.sw_score_with(&t, gaps, &mut sc);
                assert_eq!(a, want, "u8 path {:?} q={:?} t={:?}", ka, q, t);
                assert_eq!(b, want);
                assert!(!matches!(kb, Kernel::Avx2U8 | Kernel::Sse2U8));
                match ka {
                    Kernel::Avx2U8 | Kernel::Sse2U8 => used_u8 += 1,
                    _ => fell_back += 1,
                }
            }
        }
        assert!(used_u8 > 3000 && fell_back > 300, "u8 {} fallback {}", used_u8, fell_back);
        // wide score tables (PSSM-like) skip the u8 pass entirely
        let wide = vec![[-120i32, 120, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; 30];
        assert!(!QueryProfile::from_rows_opts(wide, KernelChoice::Sse2, true).has_u8());
    }

    #[test]
    fn kernel_families_are_distinct_and_detected() {
        let rows = vec![[1i32; NPROF]; 40];
        let k: Vec<Kernel> = all_kernels(&rows).iter().map(|q| q.kernel()).collect();
        #[cfg(target_arch = "x86_64")]
        {
            if std::is_x86_feature_detected!("avx2") {
                assert_eq!(k, vec![Kernel::Avx2I16, Kernel::Sse2I16, Kernel::Scalar]);
            } else {
                assert_eq!(k, vec![Kernel::Sse2I16, Kernel::Sse2I16, Kernel::Scalar]);
            }
        }
        let _ = k;
    }

    /// independent three-matrix reference (no shared code with the kernels)
    fn reference(rows: &[[i32; NPROF]], t: &[u8], open: i32, ext: i32) -> i32 {
        let (m, n) = (rows.len(), t.len());
        let ni = i32::MIN / 4;
        let mut h = vec![vec![0i32; n + 1]; m + 1];
        let mut e = vec![vec![ni; n + 1]; m + 1];
        let mut f = vec![vec![ni; n + 1]; m + 1];
        let mut best = 0;
        for i in 1..=m {
            for j in 1..=n {
                e[i][j] = (h[i][j - 1] - open).max(e[i][j - 1] - ext);
                f[i][j] = (h[i - 1][j] - open).max(f[i - 1][j] - ext);
                let r = t[j - 1] as usize;
                let sc = if r < NPROF { rows[i - 1][r] } else { 0 };
                h[i][j] = (h[i - 1][j - 1] + sc).max(e[i][j]).max(f[i][j]).max(0);
                best = best.max(h[i][j]);
            }
        }
        best
    }

    #[test]
    fn simd_matches_reference_across_gap_costs_and_compositions() {
        let m = SubstitutionMatrix::blosum62();
        let mut rng = Rng(12345);
        for (open, ext) in [(1, 1), (2, 1), (3, 1), (2, 0), (4, 1), (11, 1), (6, 2), (5, 5)] {
            for _ in 0..800 {
                let ql = 1 + rng.range(90);
                let tl = 1 + rng.range(90);
                let alpha = 2 + rng.range(19); // biased, repeat-rich compositions
                let q = rand_seq(&mut rng, ql, alpha);
                let t = rand_seq(&mut rng, tl, alpha);
                let sc = SeqScorer { query: &q, mat: &m };
                let gaps = GapCosts { open, extend: ext };
                let base = QueryProfile::build(&sc);
                let want = reference(&base.col_scores, &t, open, ext);
                for qp in all_kernels(&base.col_scores) {
                    assert_eq!(qp.sw_score(&t, gaps), want, "{:?} open={} ext={} q={:?} t={:?}", qp.kernel(), open, ext, q, t);
                }
            }
        }
    }

    /// reference end cell: full H matrix, first optimum in column-major order
    fn reference_end(rows: &[[i32; NPROF]], t: &[u8], open: i32, ext: i32) -> Option<(i32, usize, usize)> {
        let (m, n) = (rows.len(), t.len());
        let ni = i32::MIN / 4;
        let mut h = vec![vec![0i32; n + 1]; m + 1];
        let mut e = vec![vec![ni; n + 1]; m + 1];
        let mut f = vec![vec![ni; n + 1]; m + 1];
        for i in 1..=m {
            for j in 1..=n {
                e[i][j] = (h[i][j - 1] - open).max(e[i][j - 1] - ext);
                f[i][j] = (h[i - 1][j] - open).max(f[i - 1][j] - ext);
                let r = t[j - 1] as usize;
                let sc = if r < NPROF { rows[i - 1][r] } else { 0 };
                h[i][j] = (h[i - 1][j - 1] + sc).max(e[i][j]).max(f[i][j]).max(0);
            }
        }
        let mut best = (0, 0, 0);
        for j in 1..=n {
            for i in 1..=m {
                if h[i][j] > best.0 {
                    best = (h[i][j], i - 1, j - 1);
                }
            }
        }
        if best.0 > 0 {
            Some(best)
        } else {
            None
        }
    }

    #[test]
    fn score_end_matches_reference_rule_simd_and_scalar() {
        let m = SubstitutionMatrix::blosum62();
        let mut rng = Rng(424242);
        let mut s = SimdScratch::new();
        let mut simd_used = 0;
        for (open, ext) in [(11, 1), (2, 1), (6, 2), (3, 3), (1, 1)] {
            for round in 0..600 {
                let ql = 1 + rng.range(70);
                let tl = 1 + rng.range(70);
                // tiny alphabets create many tied optima
                let alpha = if round % 3 == 0 { 2 + rng.range(3) } else { 2 + rng.range(19) };
                let q = rand_seq(&mut rng, ql, alpha);
                let mut t = rand_seq(&mut rng, tl, alpha);
                if round % 7 == 0 {
                    t[0] = 20;
                }
                let base = QueryProfile::build(&SeqScorer { query: &q, mat: &m });
                let gaps = GapCosts { open, extend: ext };
                let want = reference_end(&base.col_scores, &t, open, ext);
                for qp in all_kernels(&base.col_scores) {
                    assert_eq!(qp.sw_score_end_with(&t, gaps, &mut s), want, "{:?} open={} ext={} q={:?} t={:?}", qp.kernel(), open, ext, q, t);
                    assert_eq!(want.map(|w| w.0).unwrap_or(0), qp.sw_score(&t, gaps));
                }
                assert_eq!(base.scalar_score_end(&t, gaps), want);
                if open > ext {
                    simd_used += 1;
                }
            }
        }
        assert!(simd_used > 1000);
        // saturating input: scalar fallback, same rule
        let w = build_aa2num()[b'W' as usize];
        let q = vec![w; 3000];
        let qp = QueryProfile::build(&SeqScorer { query: &q, mat: &m });
        let r = qp.sw_score_end_with(&q, GapCosts::sequence_default(), &mut s).unwrap();
        assert_eq!(r, (3000 * m.score(w, w), 2999, 2999));
    }

    #[test]
    fn adjacent_gap_runs_across_stripes_are_scored() {
        // Strongly negative mismatches make adjacent gap runs (D then I, or the
        // equal-cost I then D) optimal across stripe boundaries.
        let mut rows = vec![[-40i32; NPROF]; 40];
        for (c, row) in rows.iter_mut().enumerate() {
            row[c % 20] = 30;
        }
        let tbl_t: Vec<u8> = (0..20u8)
            .chain([5u8, 5, 5])
            .chain(26..40u8)
            .map(|x| x % 20)
            .collect();
        for qp in all_kernels(&rows) {
            let mut rng = Rng(77);
            for (open, ext) in [(44, 4), (20, 1), (12, 2)] {
                let gaps = GapCosts { open, extend: ext };
                assert_eq!(qp.sw_score(&tbl_t, gaps), reference(&rows, &tbl_t, open, ext));
                for _ in 0..400 {
                    let tl = 10 + rng.range(60);
                    let t = rand_seq(&mut rng, tl, 20);
                    assert_eq!(qp.sw_score(&t, gaps), reference(&rows, &t, open, ext), "{:?}", qp.kernel());
                }
            }
        }
    }

    #[test]
    fn simd_matches_scalar_aligner_random() {
        let m = SubstitutionMatrix::blosum62();
        let mut rng = Rng(0x1234_5678_9abc_def0);
        let gaps = GapCosts::sequence_default();
        for _ in 0..2000 {
            let ql = 1 + rng.range(60);
            let q = rand_seq(&mut rng, ql, 20);
            let tl = 1 + rng.range(80);
            let tgt = rand_seq(&mut rng, tl, 20);
            let sc = SeqScorer { query: &q, mat: &m };
            let scalar = local_align(&sc, &tgt, gaps).map(|a| a.raw_score).unwrap_or(0);
            for qp in all_kernels(&QueryProfile::build(&sc).col_scores) {
                assert_eq!(qp.sw_score(&tgt, gaps), scalar);
            }
        }
    }

    #[test]
    fn simd_matches_scalar_profile_random() {
        let m = SubstitutionMatrix::blosum62();
        let t = build_aa2num();
        let mut rng = Rng(0xdead_beef_cafe_babe);
        let gaps = GapCosts::profile_default();
        const AA: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
        for _ in 0..40 {
            let ql = 6 + rng.range(30);
            let rows: Vec<String> = (0..4)
                .map(|_| (0..ql).map(|_| AA[rng.range(20)] as char).collect())
                .collect();
            let msa = format!(">a\n{}\n>b\n{}\n>c\n{}\n>d\n{}\n", rows[0], rows[1], rows[2], rows[3]);
            let prof = Profile::from_msa_str(&msa, &m, PseudoCountParams::default()).unwrap();
            let sc = ProfileScorer { profile: &prof };
            let qp = QueryProfile::build(&sc);
            for _ in 0..50 {
                let tgt: Vec<u8> = (0..1 + rng.range(90)).map(|_| t[AA[rng.range(20)] as usize]).collect();
                let scalar = local_align(&sc, &tgt, gaps).map(|a| a.raw_score).unwrap_or(0);
                assert_eq!(qp.sw_score(&tgt, gaps), scalar);
            }
        }
    }

    #[test]
    fn simd_matches_scalar_with_x_and_out_of_alphabet_residues() {
        let m = SubstitutionMatrix::blosum62();
        let mut rng = Rng(7);
        let q = rand_seq(&mut rng, 30, 20);
        let sc = SeqScorer { query: &q, mat: &m };
        let qp = QueryProfile::build(&sc);
        let gaps = GapCosts::sequence_default();
        let mut tgt = rand_seq(&mut Rng(9), 40, 20);
        tgt[5] = 20; // X
        tgt[15] = 20;
        tgt[20] = 20;
        let scalar = local_align(&sc, &tgt, gaps).map(|a| a.raw_score).unwrap_or(0);
        assert_eq!(qp.sw_score(&tgt, gaps), scalar);
        // a table scorer maps residues >= 21 to 0, as the kernel does
        let tbl = ScoreTable::from_scorer(&sc);
        let mut odd = tgt.clone();
        odd[7] = 200;
        let qp2 = QueryProfile::build(&tbl);
        assert_eq!(qp2.sw_score(&odd, gaps), qp2.scalar_score(&odd, gaps));
    }

    #[test]
    fn saturation_falls_back_to_exact_scalar() {
        // identical 4000-residue W-rich sequences exceed i16::MAX
        let m = SubstitutionMatrix::blosum62();
        let w = build_aa2num()[b'W' as usize];
        let q = vec![w; 4000];
        let sc = SeqScorer { query: &q, mat: &m };
        let qp = QueryProfile::build(&sc);
        let gaps = GapCosts::sequence_default();
        let mut s = SimdScratch::new();
        for qp in all_kernels(&qp.col_scores) {
            let (score, kernel) = qp.sw_score_with(&q, gaps, &mut s);
            assert_eq!(kernel, Kernel::Scalar, "saturating input must use the scalar path");
            assert_eq!(score, 4000 * m.score(w, w));
            assert!(score > i16::MAX as i32);
            // just below the ceiling stays on the profile's kernel and is exact
            let n = (i16::MAX as i32 / m.score(w, w) - 2) as usize;
            let (s2, k2) = qp.sw_score_with(&q[..n], gaps, &mut s);
            assert_eq!(k2, qp.kernel());
            assert_eq!(s2, n as i32 * m.score(w, w));
        }
    }

    #[test]
    fn huge_table_scores_use_scalar() {
        let rows = vec![[40_000i32; NPROF]; 3];
        let qp = QueryProfile::from_rows(rows);
        assert_eq!(qp.sw_score(&[0, 1, 2], GapCosts::sequence_default()), 120_000);
    }
}

#[cfg(test)]
mod microbench {
    use super::*;
    use crate::align::SeqScorer;
    use crate::matrix::SubstitutionMatrix;

    /// `cargo test --release --lib simd::microbench -- --ignored --nocapture`
    #[test]
    #[ignore]
    fn microbench_kernels() {
        let m = SubstitutionMatrix::blosum62();
        let mut x = 0x9E37_79B9_7F4A_7C15u64;
        let mut rnd = |n: u64| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x % n
        };
        for ql in [60usize, 180, 400, 1000] {
            let q: Vec<u8> = (0..ql).map(|_| rnd(20) as u8).collect();
            let targets: Vec<Vec<u8>> = (0..300).map(|_| (0..50 + rnd(600) as usize).map(|_| rnd(20) as u8).collect()).collect();
            let cells: usize = targets.iter().map(|t| t.len() * ql).sum();
            let rows = QueryProfile::build(&SeqScorer { query: &q, mat: &m }).rows().to_vec();
            let mut line = format!("m={:5}:", ql);
            for (choice, u8p) in [(KernelChoice::Sse2, false), (KernelChoice::Avx2, false), (KernelChoice::Sse2, true), (KernelChoice::Avx2, true)] {
                let qp = QueryProfile::from_rows_opts(rows.clone(), choice, u8p);
                let mut s = SimdScratch::new();
                let mut best = f64::MAX;
                let mut chk = 0i64;
                for _ in 0..5 {
                    let t0 = std::time::Instant::now();
                    for t in &targets {
                        chk += qp.sw_score_with(t, GapCosts::sequence_default(), &mut s).0 as i64;
                    }
                    best = best.min(t0.elapsed().as_nanos() as f64 / cells as f64);
                }
                let t0 = std::time::Instant::now();
                for t in &targets {
                    chk += qp.sw_score_end_with(t, GapCosts::sequence_default(), &mut s).map(|r| r.0).unwrap_or(0) as i64;
                }
                let end = t0.elapsed().as_nanos() as f64 / cells as f64;
                line += &format!("  {:?}{} {:.3} ns/cell (+end {:.3})", qp.kernel(), if u8p { "+u8" } else { "" }, best, end);
                let _ = chk;
            }
            eprintln!("{}", line);
        }
    }
}
