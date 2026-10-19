//! Local composition-bias correction for sequence queries.
//!
//! Low-complexity and compositionally biased regions (poly-Q tracts,
//! collagen-like repeats, transmembrane helices, acidic stretches) produce high
//! substitution scores against unrelated proteins with the same bias. A large
//! share of BLAST-style false positives come from such regions.
//!
//! This is the correction MMseqs2 applies by default (`--comp-bias-corr 1`,
//! `SubstitutionMatrix::calcLocalAaBiasCorrection`): for query position `i` with
//! residue `a = q_i`,
//!
//! ```text
//! δ_i = Σ_b p_back(b)·s(a, b)  −  (1/|W_i|) Σ_{j∈W_i, j≠i} s(a, q_j)
//! ```
//!
//! over a window `W_i` of 40 positions centred on `i`, and every score of query
//! column `i` becomes `s(a, ·) + round(δ_i)`. Inside a biased region the local
//! mean self-similarity is high, `δ_i` is strongly negative, and matches that
//! merely reflect shared composition are discounted; in ordinary regions
//! `δ_i ≈ 0`.
//!
//! The correction changes the score distribution, so search uses gapped
//! statistics calibrated **with** the correction (see [`crate::evalue`]).

use crate::alphabet::{ALPHABET_SIZE, PROFILE_AA_SIZE};
use crate::matrix::SubstitutionMatrix;

/// Window length used by MMseqs2.
pub const WINDOW: usize = 40;

/// Per-position additive corrections (already rounded to integer score units).
pub fn local_bias_correction(query: &[u8], mat: &SubstitutionMatrix) -> Vec<i32> {
    local_bias_correction_mode(query, mat, false)
}

/// As [`local_bias_correction`]; with `penalize_only` positive corrections are
/// clipped to 0 (biased regions are discounted, others never boosted).
pub fn local_bias_correction_mode(query: &[u8], mat: &SubstitutionMatrix, penalize_only: bool) -> Vec<i32> {
    let n = query.len();
    let mut out = vec![0i32; n];
    if n == 0 {
        return out;
    }
    // expected score of residue a against the background, per residue
    let mut expected = [0f64; ALPHABET_SIZE];
    for (a, slot) in expected.iter_mut().enumerate() {
        let mut e = 0.0;
        for b in 0..PROFILE_AA_SIZE {
            e += mat.p_back[b] * mat.score(a as u8, b as u8) as f64;
        }
        *slot = e;
    }
    let half = WINDOW / 2;
    for i in 0..n {
        let lo = i.saturating_sub(half);
        let hi = (i + half).min(n);
        let a = query[i];
        let mut sum = 0i64;
        for (j, &b) in query.iter().enumerate().take(hi).skip(lo) {
            if j != i {
                sum += mat.score(a, b) as i64;
            }
        }
        let wl = (hi - lo) as f64;
        let delta = expected[(a as usize).min(ALPHABET_SIZE - 1)] - sum as f64 / wl;
        let d = round_half_away(delta);
        out[i] = if penalize_only { d.min(0) } else { d };
    }
    out
}

#[inline]
fn round_half_away(x: f64) -> i32 {
    if x >= 0.0 {
        (x + 0.5).floor() as i32
    } else {
        -((-x + 0.5).floor() as i32)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::build_aa2num;

    fn enc(s: &str) -> Vec<u8> {
        let t = build_aa2num();
        s.bytes().map(|b| t[b as usize]).collect()
    }

    #[test]
    fn biased_regions_are_discounted_ordinary_ones_are_not() {
        let m = SubstitutionMatrix::blosum62();
        let polyq = enc(&"Q".repeat(60));
        let c = local_bias_correction(&polyq, &m);
        // Q-Q = +5 half-bits; background expectation ≈ -1 → δ ≈ -6
        assert!(c.iter().all(|&d| d <= -5), "{:?}", &c[..5]);
        // a diverse sequence has small corrections on average
        let diverse = enc("MKVLLACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYMKVLLAC");
        let c2 = local_bias_correction(&diverse, &m);
        let mean: f64 = c2.iter().map(|&x| x as f64).sum::<f64>() / c2.len() as f64;
        assert!(mean.abs() < 1.5, "mean correction {}", mean);
    }

    #[test]
    fn short_and_empty_queries() {
        let m = SubstitutionMatrix::blosum62();
        assert!(local_bias_correction(&[], &m).is_empty());
        assert_eq!(local_bias_correction(&enc("W"), &m).len(), 1);
    }
}
