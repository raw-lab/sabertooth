//! SIMD-vectorized local-alignment scoring (striped Smith–Waterman, Farrar 2007).
//!
//! The scalar [`crate::align::local_align`] computes the full alignment with
//! traceback — that is what produces the CIGAR and coordinates for a reported
//! hit. But most prefilter survivors are *rejected* by the E-value cutoff, and
//! paying for a full `(m+1)(n+1)` traceback matrix on each of them is wasteful.
//! This module provides a fast, traceback-free **score** using Farrar's striped
//! layout over SSE2 (8 × `i16` lanes), so a search can gate on the score cheaply
//! and only run the scalar aligner on the few candidates that pass.
//!
//! Correctness contract: [`QueryProfile::sw_score`] returns exactly the same raw
//! score as `local_align(...).raw_score` for every input (verified by differential
//! tests against the scalar aligner). On the rare chance of `i16` saturation, or
//! on non-x86-64 targets, it transparently falls back to an equivalent scalar
//! score, so the number is always trustworthy.
//!
//! Only SSE2 is used — it is part of the x86-64 baseline, so no runtime feature
//! detection is needed and the kernel is safe to call on any x86-64 CPU. (AVX2
//! would add lanes but requires runtime dispatch and is left as future work.)

use crate::align::{GapCosts, Scorer};

const NPROF: usize = crate::alphabet::ALPHABET_SIZE; // 21 columns incl. X

/// A query (sequence or profile) prepared for fast repeated scoring against many
/// targets. Build once, score many times.
pub struct QueryProfile {
    m: usize,
    /// Per-column score of each real residue: `col_scores[col][res]`.
    col_scores: Vec<[i32; NPROF]>,
    #[cfg(target_arch = "x86_64")]
    striped: Option<x86::StripedI16>,
}

impl QueryProfile {
    /// Build a query profile from any [`Scorer`] (sequence- or PSSM-backed).
    pub fn build<S: Scorer>(scorer: &S) -> QueryProfile {
        let m = scorer.query_len();
        let mut col_scores = vec![[0i32; NPROF]; m];
        for (c, row) in col_scores.iter_mut().enumerate() {
            for (r, slot) in row.iter_mut().enumerate() {
                *slot = scorer.score(c, r as u8);
            }
        }
        #[cfg(target_arch = "x86_64")]
        {
            let striped = x86::StripedI16::build(&col_scores, m);
            QueryProfile { m, col_scores, striped }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            QueryProfile { m, col_scores }
        }
    }

    /// Optimal Smith–Waterman–Gotoh local-alignment **score** (no traceback).
    /// Always equals `local_align(...).raw_score`.
    pub fn sw_score(&self, target: &[u8], gaps: GapCosts) -> i32 {
        if self.m == 0 || target.is_empty() {
            return 0;
        }
        #[cfg(target_arch = "x86_64")]
        {
            if let Some(striped) = &self.striped {
                // SAFETY: SSE2 is guaranteed on all x86-64 targets.
                if let Some(score) = unsafe { striped.sw_score_sse2(target, gaps) } {
                    return score;
                }
                // i16 saturated → fall through to the exact scalar score.
            }
        }
        scalar_gotoh_score(&self.col_scores, self.m, target, gaps)
    }
}

/// Traceback-free scalar Gotoh local-alignment score over prepared column
/// scores. Portable path and saturation fallback; O(m·n) time, O(n) memory.
/// Byte-for-byte matches `local_align(...).raw_score`.
fn scalar_gotoh_score(col_scores: &[[i32; NPROF]], m: usize, target: &[u8], gaps: GapCosts) -> i32 {
    let n = target.len();
    let neg = i32::MIN / 4;
    let mut h_prev = vec![0i32; n + 1];
    let mut h_cur = vec![0i32; n + 1];
    let mut f_col = vec![neg; n + 1]; // vertical-gap F[i][j], carried across rows
    let mut best = 0i32;
    for i in 0..m {
        let row = &col_scores[i];
        let mut e = neg; // horizontal-gap running value E[i][j-1]
        h_cur[0] = 0;
        for j in 1..=n {
            let tr = target[j - 1] as usize;
            let s = if tr < NPROF { row[tr] } else { 0 };
            let diag = h_prev[j - 1] + s;
            e = (h_cur[j - 1] - gaps.open).max(e - gaps.extend);
            f_col[j] = (h_prev[j] - gaps.open).max(f_col[j] - gaps.extend);
            let mut v = diag;
            if e > v {
                v = e;
            }
            if f_col[j] > v {
                v = f_col[j];
            }
            if v < 0 {
                v = 0;
            }
            h_cur[j] = v;
            if v > best {
                best = v;
            }
        }
        std::mem::swap(&mut h_prev, &mut h_cur);
    }
    best
}

#[cfg(target_arch = "x86_64")]
mod x86 {
    use super::NPROF;
    use crate::align::GapCosts;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    const LANES: usize = 8; // i16 lanes in a 128-bit register

    /// Striped `i16` query profile (Farrar layout). For segment `i` and lane `l`,
    /// the represented query column is `i + l*seg_len` (0 padded past the end).
    pub struct StripedI16 {
        seg_len: usize,
        /// `vprofile[res * seg_len + i]` packed as 8×i16, one lane per stripe.
        vprofile: Vec<[i16; LANES]>,
        /// Zero vector index sentinel for target residue X (>= NPROF): score 0.
        zero_row: Vec<[i16; LANES]>,
        /// True if any prepared score is large enough that overflow is a risk;
        /// then callers should trust the saturation check.
        _big: bool,
    }

    impl StripedI16 {
        pub fn build(col_scores: &[[i32; NPROF]], m: usize) -> Option<StripedI16> {
            if m == 0 {
                return None;
            }
            let seg_len = m.div_ceil(LANES);
            let mut vprofile = vec![[0i16; LANES]; NPROF * seg_len];
            let mut big = false;
            for res in 0..NPROF {
                for i in 0..seg_len {
                    let mut lane = [0i16; LANES];
                    for (l, slot) in lane.iter_mut().enumerate() {
                        let col = i + l * seg_len;
                        let s = if col < m { col_scores[col][res] } else { 0 };
                        if s > i16::MAX as i32 || s < i16::MIN as i32 {
                            big = true;
                        }
                        *slot = s.clamp(i16::MIN as i32, i16::MAX as i32) as i16;
                    }
                    vprofile[res * seg_len + i] = lane;
                }
            }
            let zero_row = vec![[0i16; LANES]; seg_len];
            Some(StripedI16 {
                seg_len,
                vprofile,
                zero_row,
                _big: big,
            })
        }

        /// Striped local SW score over SSE2. Returns `None` if the score
        /// saturates `i16` (caller falls back to scalar).
        ///
        /// SAFETY: caller guarantees SSE2 (always true on x86-64).
        #[target_feature(enable = "sse2")]
        pub unsafe fn sw_score_sse2(&self, target: &[u8], gaps: GapCosts) -> Option<i32> {
            let seg_len = self.seg_len;
            let gap_open = gaps.open.clamp(0, i16::MAX as i32) as i16;
            let gap_extend = gaps.extend.clamp(0, i16::MAX as i32) as i16;

            let v_open = _mm_set1_epi16(gap_open);
            let v_extend = _mm_set1_epi16(gap_extend);
            let v_zero = _mm_setzero_si128();

            // H store/load and E arrays, one 128-bit vector per segment.
            let mut h_store = vec![v_zero; seg_len];
            let mut h_load = vec![v_zero; seg_len];
            let mut e_arr = vec![v_zero; seg_len];

            let mut v_max = v_zero;
            // Track the i16 saturation ceiling to detect overflow.
            let sat = _mm_set1_epi16(i16::MAX);

            for &tb in target {
                let tr = tb as usize;
                let prof: &[[i16; LANES]] = if tr < NPROF {
                    &self.vprofile[tr * seg_len..tr * seg_len + seg_len]
                } else {
                    &self.zero_row[..]
                };

                // vH = last H vector of previous column, shifted left by one lane.
                let mut v_h = _mm_slli_si128(h_store[seg_len - 1], 2);
                let mut v_f = v_zero;

                std::mem::swap(&mut h_store, &mut h_load);

                for i in 0..seg_len {
                    let v_p = load_lane(&prof[i]);
                    // diagonal + match
                    v_h = _mm_adds_epi16(v_h, v_p);
                    // max with E and F, floor at 0 (local)
                    v_h = _mm_max_epi16(v_h, e_arr[i]);
                    v_h = _mm_max_epi16(v_h, v_f);
                    v_h = _mm_max_epi16(v_h, v_zero);
                    v_max = _mm_max_epi16(v_max, v_h);
                    h_store[i] = v_h;

                    // E[i] = max(E[i] - extend, H - open)
                    let v_h_open = _mm_subs_epi16(v_h, v_open);
                    e_arr[i] = _mm_subs_epi16(e_arr[i], v_extend);
                    e_arr[i] = _mm_max_epi16(e_arr[i], v_h_open);

                    // F = max(F - extend, H - open)
                    v_f = _mm_subs_epi16(v_f, v_extend);
                    v_f = _mm_max_epi16(v_f, v_h_open);

                    // load previous column's H for next iteration's diagonal
                    v_h = h_load[i];
                }

                // Lazy-F loop: propagate F across the stripe boundary. A value can
                // cross at most LANES stripes, bounding the outer iterations.
                'lazy: for _ in 0..LANES {
                    v_f = _mm_slli_si128(v_f, 2);
                    for i in 0..seg_len {
                        let mut v_h = h_store[i];
                        v_h = _mm_max_epi16(v_h, v_f);
                        v_max = _mm_max_epi16(v_max, v_h);
                        h_store[i] = v_h;
                        let v_h_open = _mm_subs_epi16(v_h, v_open);
                        v_f = _mm_subs_epi16(v_f, v_extend);
                        // stop when F can no longer improve any H via a new gap
                        let cmp = _mm_cmpgt_epi16(v_f, v_h_open);
                        if _mm_movemask_epi8(cmp) == 0 {
                            break 'lazy;
                        }
                    }
                }
            }

            let best = hmax_epi16(v_max);
            // saturation guard: if we hit the i16 ceiling, the score is unreliable
            let hit_ceiling = _mm_movemask_epi8(_mm_cmpeq_epi16(v_max, sat)) != 0;
            if hit_ceiling {
                None
            } else {
                Some(best as i32)
            }
        }
    }

    #[inline]
    #[target_feature(enable = "sse2")]
    unsafe fn load_lane(lane: &[i16; 8]) -> __m128i {
        _mm_loadu_si128(lane.as_ptr() as *const __m128i)
    }

    /// Horizontal max of 8 packed i16.
    #[inline]
    #[target_feature(enable = "sse2")]
    unsafe fn hmax_epi16(v: __m128i) -> i16 {
        // reduce via shuffles
        let mut m = v;
        m = _mm_max_epi16(m, _mm_shuffle_epi32(m, 0b01_00_11_10)); // hi 64 vs lo 64
        m = _mm_max_epi16(m, _mm_shuffle_epi32(m, 0b00_00_00_01)); // within 64
        m = _mm_max_epi16(m, _mm_shufflelo_epi16(m, 0b00_00_00_01)); // within 32
        _mm_extract_epi16(m, 0) as i16
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{local_align, GapCosts, ProfileScorer, SeqScorer};
    use crate::alphabet::build_aa2num;
    use crate::matrix::SubstitutionMatrix;
    use crate::profile::{Profile, PseudoCountParams};

    // simple deterministic PRNG (xorshift) so tests need no dependencies
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

    fn rand_seq(rng: &mut Rng, len: usize, t: &[u8; 256]) -> Vec<u8> {
        const AA: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
        (0..len).map(|_| t[AA[rng.range(20)] as usize]).collect()
    }

    #[test]
    fn simd_matches_scalar_sequence_random() {
        let m = SubstitutionMatrix::blosum62();
        let t = build_aa2num();
        let mut rng = Rng(0x1234_5678_9abc_def0);
        let gaps = GapCosts::sequence_default();
        let mut checked = 0;
        for _ in 0..2000 {
            let ql = 1 + rng.range(60);
            let tl = 1 + rng.range(80);
            let q = rand_seq(&mut rng, ql, &t);
            let tgt = rand_seq(&mut rng, tl, &t);
            let sc = SeqScorer { query: &q, mat: &m };
            let scalar = local_align(&sc, &tgt, gaps).map(|a| a.raw_score).unwrap_or(0);
            let qp = QueryProfile::build(&sc);
            let simd = qp.sw_score(&tgt, gaps);
            assert_eq!(simd, scalar, "seq mismatch q={:?} t={:?}", q, tgt);
            checked += 1;
        }
        assert_eq!(checked, 2000);
    }

    #[test]
    fn simd_matches_scalar_profile_random() {
        let m = SubstitutionMatrix::blosum62();
        let t = build_aa2num();
        let mut rng = Rng(0xdead_beef_cafe_babe);
        let gaps = GapCosts::profile_default();
        // build one profile from a small random MSA, reuse across targets
        for _ in 0..40 {
            let ql = 6 + rng.range(30);
            let rows: Vec<String> = (0..4)
                .map(|_| {
                    let s = rand_seq(&mut rng, ql, &t);
                    // decode back to letters for the MSA string
                    const AA: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
                    s.iter().map(|&x| AA[x as usize] as char).collect()
                })
                .collect();
            let msa = format!(
                ">a\n{}\n>b\n{}\n>c\n{}\n>d\n{}\n",
                rows[0], rows[1], rows[2], rows[3]
            );
            let prof = Profile::from_msa_str(&msa, &m, PseudoCountParams::default()).unwrap();
            let sc = ProfileScorer { profile: &prof };
            let qp = QueryProfile::build(&sc);
            for _ in 0..50 {
                let tl = 1 + rng.range(90);
                let tgt = rand_seq(&mut rng, tl, &t);
                let scalar = local_align(&sc, &tgt, gaps).map(|a| a.raw_score).unwrap_or(0);
                let simd = qp.sw_score(&tgt, gaps);
                assert_eq!(simd, scalar, "profile mismatch");
            }
        }
    }

    #[test]
    fn simd_matches_scalar_with_x_residues() {
        // targets containing X (index >= 20) must score identically (X → 0).
        let m = SubstitutionMatrix::blosum62();
        let t = build_aa2num();
        let q = rand_seq(&mut Rng(7), 30, &t);
        let sc = SeqScorer { query: &q, mat: &m };
        let qp = QueryProfile::build(&sc);
        let gaps = GapCosts::sequence_default();
        let mut tgt = rand_seq(&mut Rng(9), 40, &t);
        tgt[5] = 20; // X
        tgt[15] = 20; // X
        tgt[20] = 20; // X
        let scalar = local_align(&sc, &tgt, gaps).map(|a| a.raw_score).unwrap_or(0);
        assert_eq!(qp.sw_score(&tgt, gaps), scalar);
    }
}
