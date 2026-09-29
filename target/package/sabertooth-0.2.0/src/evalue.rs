//! Karlin–Altschul statistics for local alignment significance.
//!
//! MMseqs2 obtains Gumbel parameters (λ, K) from the ALP simulation library for
//! the exact integer matrix + gap costs in use. Sabertooth instead solves for
//! the *ungapped* λ analytically from the score distribution of the matrix it
//! actually uses — so λ automatically tracks the `bit_factor`/`bias` scaling —
//! and uses the standard BLOSUM62 reference `K`. Bit scores and E-values are
//! therefore on a genuine bit scale and rank targets consistently, which is what
//! the prefilter → align → sort pipeline needs.
//!
//! ```text
//! bit_score = (lambda * raw_score - ln K) / ln 2
//! E_value   = K * m * n * exp(-lambda * raw_score) = search_space * 2^(-bit_score)
//! ```

use crate::alphabet::PROFILE_AA_SIZE;
use crate::matrix::SubstitutionMatrix;

/// A small, dependency-free SplitMix64 PRNG for the Monte-Carlo calibration.
struct SplitMix64(u64);
impl SplitMix64 {
    fn new(seed: u64) -> Self {
        SplitMix64(seed)
    }
    #[inline]
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    #[inline]
    fn next_f64(&mut self) -> f64 {
        // 53-bit mantissa in [0, 1)
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Gapped Gumbel (Type-I extreme value) parameters fitted by simulation.
#[derive(Clone, Copy, Debug)]
pub struct GumbelParams {
    /// Scale parameter λ (decay rate of the score distribution's upper tail).
    pub lambda: f64,
    /// Karlin–Altschul K.
    pub k: f64,
    /// Fitted location μ for the simulated sequence length.
    pub mu: f64,
    /// Number of random pairs aligned.
    pub samples: usize,
    /// Mean / std of the simulated optimal local scores (diagnostics).
    pub mean_score: f64,
    pub std_score: f64,
}

/// Calibrate **gapped** local-alignment statistics by Monte-Carlo simulation:
/// align many random sequence pairs drawn from the matrix background
/// distribution and fit a Gumbel distribution to the optimal local scores by
/// the method of moments. This is the pure-Rust analogue of the ALP library
/// MMseqs2 links: the analytic Karlin–Altschul λ is exact only for *ungapped*
/// alignment, and gaps raise the score variance, lowering the effective λ. The
/// fast striped aligner makes the thousands of random alignments cheap.
///
/// Method of moments for a Gumbel(μ, β=1/λ): `mean = μ + γ/λ`,
/// `var = π²/(6λ²)`, so `λ = π/(σ√6)` and `μ = mean − γ/λ`. The local-alignment
/// EVD places the mode at `μ = ln(K·m·n)/λ`, giving `K = e^{λμ}/(m·n)`.
pub fn calibrate_gapped(
    mat: &SubstitutionMatrix,
    gaps: crate::align::GapCosts,
    seq_len: usize,
    num_pairs: usize,
    seed: u64,
) -> GumbelParams {
    use crate::align::SeqScorer;
    use crate::simd::QueryProfile;

    // Cumulative background over the 20 real residues (renormalized).
    let total: f64 = (0..PROFILE_AA_SIZE).map(|i| mat.p_back[i]).sum();
    let mut cum = [0f64; PROFILE_AA_SIZE];
    let mut acc = 0.0;
    for i in 0..PROFILE_AA_SIZE {
        acc += mat.p_back[i] / total;
        cum[i] = acc;
    }
    let sample = |u: f64| -> u8 {
        for (i, &c) in cum.iter().enumerate() {
            if u <= c {
                return i as u8;
            }
        }
        (PROFILE_AA_SIZE - 1) as u8
    };

    let mut rng = SplitMix64::new(seed);
    let mut scores: Vec<f64> = Vec::with_capacity(num_pairs);
    for _ in 0..num_pairs {
        let q: Vec<u8> = (0..seq_len).map(|_| sample(rng.next_f64())).collect();
        let t: Vec<u8> = (0..seq_len).map(|_| sample(rng.next_f64())).collect();
        let sc = SeqScorer { query: &q, mat };
        let qp = QueryProfile::build(&sc);
        scores.push(qp.sw_score(&t, gaps) as f64);
    }

    let n = scores.len().max(1) as f64;
    let mean = scores.iter().sum::<f64>() / n;
    let var = if scores.len() > 1 {
        scores.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / (n - 1.0)
    } else {
        1.0
    };
    let std = var.sqrt().max(1e-9);
    const EULER: f64 = 0.577_215_664_901_532_9;
    let lambda = std::f64::consts::PI / (std * 6f64.sqrt());
    let mu = mean - EULER / lambda;
    let k = (lambda * mu).exp() / (seq_len as f64 * seq_len as f64);
    GumbelParams {
        lambda,
        k,
        mu,
        samples: scores.len(),
        mean_score: mean,
        std_score: std,
    }
}

/// Statistical parameters for significance computation.
pub struct EValueParams {
    /// Ungapped Karlin–Altschul lambda for the matrix's integer scores.
    pub lambda: f64,
    /// Karlin–Altschul K (BLOSUM62 reference default; scale-robust).
    pub k: f64,
    /// Total residues in the target database (n in `m*n`).
    pub db_residues: f64,
    /// ln(K), cached.
    pub log_k: f64,
}

impl EValueParams {
    /// Build parameters for a matrix and target database size.
    pub fn new(mat: &SubstitutionMatrix, db_residues: usize) -> Self {
        let lambda = solve_lambda(mat);
        let k = 0.041; // BLOSUM62 reference K
        EValueParams {
            lambda,
            k,
            db_residues: db_residues.max(1) as f64,
            log_k: k.ln(),
        }
    }

    /// Build parameters for **profile (PSSM) scores**.
    ///
    /// A PSSM value is a pure log-odds quantity scaled by `bit_factor`:
    /// `M = bit_factor · log2(p/pBack)`. For such scores the Karlin–Altschul
    /// equation `Σ pBack·exp(λ·M) = 1` is solved exactly by `λ = ln2 / bit_factor`
    /// (substituting makes `exp(λM) = p/pBack`, and `Σ pBack·(p/pBack) = Σ p = 1`).
    /// With MMseqs2's profile `bit_factor = 8`, `λ = ln2/8 ≈ 0.0866`. This keeps
    /// profile bit scores on a genuine bit scale and — crucially — keeps decoys
    /// non-significant instead of the wildly inflated significance you get from
    /// reusing the ×2 substitution-matrix λ.
    pub fn for_profile(db_residues: usize) -> Self {
        let lambda = std::f64::consts::LN_2 / crate::profile::PROFILE_BIT_FACTOR;
        let k = 0.041;
        EValueParams {
            lambda,
            k,
            db_residues: db_residues.max(1) as f64,
            log_k: k.ln(),
        }
    }

    /// Reconstruct from raw Karlin–Altschul parameters (used by the distributed
    /// backend, where a worker rebuilds the statistics from the shard input).
    pub fn from_raw(lambda: f64, k: f64, db_residues: usize) -> Self {
        EValueParams {
            lambda,
            k,
            db_residues: db_residues.max(1) as f64,
            log_k: k.ln(),
        }
    }

    /// Bit score from a raw integer alignment score.
    #[inline]
    pub fn bit_score(&self, raw_score: i32) -> f64 {
        (self.lambda * raw_score as f64 - self.log_k) / std::f64::consts::LN_2
    }

    /// E-value for a raw score given the query length.
    #[inline]
    pub fn evalue(&self, raw_score: i32, query_len: usize) -> f64 {
        let search_space = query_len.max(1) as f64 * self.db_residues;
        let bits = self.bit_score(raw_score);
        search_space * 2f64.powf(-bits)
    }
}

/// Solve `sum_ij p_i p_j exp(lambda * s_ij) = 1` for the unique positive lambda,
/// using bisection. Requires the expected score to be negative and some
/// positive score to exist (both hold for any real substitution matrix).
fn solve_lambda(mat: &SubstitutionMatrix) -> f64 {
    let n = PROFILE_AA_SIZE;
    let p: Vec<f64> = (0..n).map(|i| mat.p_back[i]).collect();

    // f(lambda) = sum p_i p_j exp(lambda s_ij) - 1 ; monotonically increasing.
    let f = |lambda: f64| -> f64 {
        let mut acc = 0.0;
        for i in 0..n {
            for j in 0..n {
                acc += p[i] * p[j] * (lambda * mat.score[i][j] as f64).exp();
            }
        }
        acc - 1.0
    };

    // Bracket: f(0) = 0 exactly (sum of probs = 1). We want the positive root.
    // f is convex; grows to +inf. Find an upper bound where f > 0.
    let mut lo = 1e-6;
    let mut hi = 2.0;
    // Ensure f(lo) < 0 (just above 0 the derivative is the mean score < 0).
    if f(lo) > 0.0 {
        // Degenerate matrix; fall back to a sane default.
        return 0.3;
    }
    let mut iters = 0;
    while f(hi) < 0.0 && iters < 60 {
        hi *= 2.0;
        iters += 1;
    }
    // Bisection.
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if f(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lambda_is_in_expected_range() {
        // With bit_factor = 2 (half-bit scores), ungapped lambda for BLOSUM62 is
        // ~0.34 on the natural-log scale (i.e. ln2/2 ≈ 0.347).
        let m = SubstitutionMatrix::blosum62();
        let p = EValueParams::new(&m, 1_000_000);
        assert!(p.lambda > 0.25 && p.lambda < 0.45, "lambda = {}", p.lambda);
    }

    #[test]
    fn higher_score_is_more_significant() {
        let m = SubstitutionMatrix::blosum62();
        let p = EValueParams::new(&m, 1_000_000);
        let e_low = p.evalue(40, 300);
        let e_high = p.evalue(200, 300);
        assert!(e_high < e_low);
        assert!(p.bit_score(200) > p.bit_score(40));
    }

    #[test]
    fn gapped_calibration_is_sane_and_below_ungapped() {
        let m = SubstitutionMatrix::blosum62();
        let ungapped = EValueParams::new(&m, 1).lambda;
        let g = calibrate_gapped(
            &m,
            crate::align::GapCosts::sequence_default(),
            200,
            1500,
            0xC0FFEE,
        );
        eprintln!(
            "calibrated: lambda={:.4} k={:.4} mu={:.2} mean={:.2} std={:.2} (ungapped lambda={:.4})",
            g.lambda, g.k, g.mu, g.mean_score, g.std_score, ungapped
        );
        assert!(g.lambda > 0.0 && g.lambda.is_finite(), "lambda={}", g.lambda);
        assert!(g.k > 0.0 && g.k.is_finite(), "k={}", g.k);
        // gaps raise score variance, so the fitted λ is below the ungapped one
        assert!(
            g.lambda < ungapped,
            "gapped λ {} should be < ungapped λ {}",
            g.lambda,
            ungapped
        );
        // reproducible for a fixed seed
        let g2 = calibrate_gapped(
            &m,
            crate::align::GapCosts::sequence_default(),
            200,
            1500,
            0xC0FFEE,
        );
        assert_eq!(g.mean_score, g2.mean_score);
        assert_eq!(g.lambda, g2.lambda);
    }
}
