//! Statistical significance of local alignment scores.
//!
//! ## What changed from v0.2 (and why it matters for false positives)
//!
//! v0.2 reported E-values for **gapped** Smith–Waterman scores using the
//! **ungapped** Karlin–Altschul λ solved analytically for the matrix
//! (λ ≈ 0.32 at the half-bit scale) together with BLAST's *gapped* K = 0.041,
//! and no finite-length correction. Gaps widen the score distribution (the
//! gapped λ is ≈ 0.27), so every E-value was too small — by ~10× at raw score
//! 40 and ~45× at raw score 60 — making unrelated sequences look significant.
//! The `calibrate` command already existed but search never used its result.
//!
//! Now every reported E-value comes from parameters **calibrated for the exact
//! scoring configuration in use**:
//!
//! * Sequence search: λ, K and the finite-size length-adjustment parameters
//!   (α/λ, β; Altschul–Gish / BLAST `BLAST_ComputeLengthAdjustment`) are fitted
//!   by Monte-Carlo simulation of random sequence pairs with the integer matrix
//!   *as reconstructed here* (BLOSUM62, bit factor 2, bias −0.2), gaps 11/1, with
//!   and without composition-bias correction. The values are stored in a
//!   verified table ([`SEQUENCE_TABLE`]) produced by `sabertooth calibrate
//!   --table` (seeded, reproducible) and re-checked by a unit test.
//! * Profile search: each PSSM gets its own fit (HMMER-style): random
//!   background targets with the target database's length distribution are
//!   scored against that PSSM and a Gumbel is fitted by maximum likelihood
//!   ([`calibrate_profile`]). One sequence-derived calibration is never assumed
//!   to hold for every PSSM.
//!
//! ## Definitions (no posterior probabilities)
//!
//! ```text
//! bit_score = (λ·S − ln K) / ln 2
//! E         = K · (m − ℓ) · (N − n·ℓ) · e^(−λ·S)
//! ```
//!
//! `m` = query length, `N` = database residues, `n` = database sequences, `ℓ` =
//! length adjustment (0 when no length model). The E-value is the expected
//! number of alignments with score ≥ S between the query and an unrelated
//! database of this size. It is **not** a posterior probability, a p-value of a
//! single comparison, or a false-discovery rate / q-value. `N` and `n` are always
//! the *global* database size, also when the search is sharded.

use crate::align::{GapCosts, ScoreTable, SeqScorer};
use crate::alphabet::PROFILE_AA_SIZE;
use crate::matrix::SubstitutionMatrix;
use crate::simd::{QueryProfile, SimdScratch};
use rayon::prelude::*;

/// Euler–Mascheroni constant.
pub const EULER: f64 = 0.577_215_664_901_532_9;

/// A small, dependency-free SplitMix64 PRNG for reproducible simulation.
#[derive(Clone)]
pub struct SplitMix64(u64);
impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        SplitMix64(seed)
    }
    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    #[inline]
    pub fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

/// Derive an independent stream for sample `i` of cell `c`.
fn stream_seed(seed: u64, c: u64, i: u64) -> u64 {
    let mut s = SplitMix64::new(seed ^ c.wrapping_mul(0xA24B_AED4_963E_E407) ^ i.wrapping_mul(0x9FB2_1C65_1E98_DF25));
    s.next_u64()
}

/// Sampler for i.i.d. residues from a background distribution.
#[derive(Clone)]
pub struct BackgroundSampler {
    cum: [f64; PROFILE_AA_SIZE],
}

impl BackgroundSampler {
    pub fn new(mat: &SubstitutionMatrix) -> Self {
        let total: f64 = (0..PROFILE_AA_SIZE).map(|i| mat.p_back[i]).sum();
        let mut cum = [0f64; PROFILE_AA_SIZE];
        let mut acc = 0.0;
        for (i, c) in cum.iter_mut().enumerate() {
            acc += mat.p_back[i] / total;
            *c = acc;
        }
        cum[PROFILE_AA_SIZE - 1] = 1.0;
        BackgroundSampler { cum }
    }
    #[inline]
    pub fn sample(&self, rng: &mut SplitMix64) -> u8 {
        let u = rng.next_f64();
        for (i, &c) in self.cum.iter().enumerate() {
            if u < c {
                return i as u8;
            }
        }
        (PROFILE_AA_SIZE - 1) as u8
    }
    pub fn sequence(&self, rng: &mut SplitMix64, len: usize) -> Vec<u8> {
        (0..len).map(|_| self.sample(rng)).collect()
    }
}

// ---------------------------------------------------------------------------
// Gumbel maximum-likelihood fitting
// ---------------------------------------------------------------------------

/// Maximum-likelihood Gumbel (type-I extreme value, maximum) fit:
/// `P(X ≤ x) = exp(−e^{−λ(x−μ)})`.
#[derive(Clone, Copy, Debug)]
pub struct GumbelFit {
    pub lambda: f64,
    pub mu: f64,
    pub n: usize,
}

/// ML location for a fixed λ: `μ = −(1/λ)·ln(mean(e^{−λx}))` (stable form).
pub fn gumbel_mu_given_lambda(xs: &[f64], lambda: f64) -> f64 {
    let xmin = xs.iter().cloned().fold(f64::INFINITY, f64::min);
    let s: f64 = xs.iter().map(|&x| (-lambda * (x - xmin)).exp()).sum();
    xmin - (s / xs.len() as f64).ln() / lambda
}

/// Lawless' ML equations, solved for λ by bisection on
/// `g(λ) = 1/λ − mean + Σx·e^{−λx}/Σe^{−λx}`.
pub fn fit_gumbel_ml(xs: &[f64]) -> Option<GumbelFit> {
    let n = xs.len();
    if n < 10 {
        return None;
    }
    let mean = xs.iter().sum::<f64>() / n as f64;
    let var = xs.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / (n as f64 - 1.0);
    if var.is_nan() || var <= 0.0 {
        return None;
    }
    let xmin = xs.iter().cloned().fold(f64::INFINITY, f64::min);
    let g = |lam: f64| -> f64 {
        let mut s0 = 0.0;
        let mut s1 = 0.0;
        for &x in xs {
            let w = (-lam * (x - xmin)).exp();
            s0 += w;
            s1 += x * w;
        }
        1.0 / lam - mean + s1 / s0
    };
    // moments estimate brackets the root
    let lam0 = std::f64::consts::PI / (var.sqrt() * 6f64.sqrt());
    let (mut lo, mut hi) = (lam0 / 10.0, lam0 * 10.0);
    if g(lo) < 0.0 || g(hi) > 0.0 {
        return None;
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if g(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-12 * hi {
            break;
        }
    }
    let lambda = 0.5 * (lo + hi);
    Some(GumbelFit {
        lambda,
        mu: gumbel_mu_given_lambda(xs, lambda),
        n,
    })
}

// ---------------------------------------------------------------------------
// BLAST length adjustment
// ---------------------------------------------------------------------------

/// Length adjustment ℓ for query length `m` against a database of `n_res`
/// residues in `n_seqs` sequences: the integer part of the fixed point of
/// `ℓ = (α/λ)·(ln K + ln((m − ℓ)(n_res − n_seqs·ℓ))) + β`, bounded so that
/// `K·(m−ℓ)·(n_res−n_seqs·ℓ) ≥ max(m, n_res)` (Altschul & Gish 1996; same
/// procedure as NCBI BLAST's `BLAST_ComputeLengthAdjustment`).
pub fn length_adjustment(k: f64, alpha_d_lambda: f64, beta: f64, m: f64, n_res: f64, n_seqs: f64) -> f64 {
    let log_k = k.ln();
    let a = n_seqs;
    let mb = m * n_seqs + n_res;
    let c = n_res * m - m.max(n_res) / k;
    if c < 0.0 {
        return 0.0;
    }
    let ell_max = 2.0 * c / (mb + (mb * mb - 4.0 * a * c).max(0.0).sqrt());
    let (mut ell_min, mut ell_hi) = (0.0f64, ell_max);
    let mut ell = 0.0f64;
    let mut converged = false;
    for i in 1..=20 {
        let ss = (m - ell) * (n_res - n_seqs * ell);
        let ell_bar = alpha_d_lambda * (log_k + ss.max(1.0).ln()) + beta;
        if ell_bar >= ell {
            ell_min = ell;
            if ell_bar - ell_min <= 1.0 {
                converged = true;
                break;
            }
            if ell_min == ell_hi {
                break;
            }
        } else {
            ell_hi = ell;
        }
        ell = if ell_min <= ell_bar && ell_bar <= ell_hi {
            ell_bar
        } else if i == 1 {
            ell_hi
        } else {
            0.5 * (ell_min + ell_hi)
        };
    }
    let mut out = ell_min.floor();
    if converged {
        let up = ell_min.ceil();
        if up <= ell_hi {
            let ss = (m - up) * (n_res - n_seqs * up);
            if alpha_d_lambda * (log_k + ss.max(1.0).ln()) + beta >= up {
                out = up;
            }
        }
    }
    out.max(0.0)
}

// ---------------------------------------------------------------------------
// E-value parameters
// ---------------------------------------------------------------------------

/// Finite-size (length adjustment) model: `ℓ = (α/λ)·ln(K·m'·N') + β`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LengthModel {
    pub alpha_d_lambda: f64,
    pub beta: f64,
}

/// Calibrated statistics for one scoring configuration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GappedParams {
    pub lambda: f64,
    pub k: f64,
    pub length: Option<LengthModel>,
}

/// One row of the verified sequence-statistics table.
#[derive(Clone, Copy, Debug)]
pub struct TableEntry {
    pub gaps: GapCosts,
    pub comp_bias: bool,
    pub params: GappedParams,
    /// How the entry was produced (command + seed + samples).
    pub provenance: &'static str,
}

/// Verified gapped parameters for the supported sequence-scoring
/// configurations (BLOSUM62 reconstructed at bit factor 2 / bias −0.2).
/// Regenerate with `sabertooth calibrate --table` (values below were produced by
/// that command; see `bench/results/calibration_table.json`). Checked against a
/// fresh, independent-seed calibration by `tests::table_matches_fresh_calibration`.
pub const SEQUENCE_TABLE: &[TableEntry] = &[
    TableEntry {
        gaps: GapCosts { open: 11, extend: 1 },
        comp_bias: false,
        params: GappedParams {
            lambda: 0.3666,
            k: 0.1177,
            length: Some(LengthModel {
                alpha_d_lambda: 2.7,
                beta: -13.0,
            }),
        },
        provenance: "sabertooth 0.3.0 calibrate --table --samples 20000 --seed 24301 (lengths [60, 100, 150, 220, 320, 480, 700, 1000])",
    },
    TableEntry {
        gaps: GapCosts { open: 11, extend: 1 },
        comp_bias: true,
        params: GappedParams {
            lambda: 0.3649,
            k: 0.1076,
            length: Some(LengthModel {
                alpha_d_lambda: 2.4,
                beta: -12.0,
            }),
        },
        provenance: "sabertooth 0.3.0 calibrate --table --samples 20000 --seed 24301 (lengths [60, 100, 150, 220, 320, 480, 700, 1000])",
    },
];

/// Look up the verified table.
pub fn table_params(gaps: GapCosts, comp_bias: bool) -> Option<&'static TableEntry> {
    SEQUENCE_TABLE
        .iter()
        .find(|e| e.gaps == gaps && e.comp_bias == comp_bias)
}

/// Statistical parameters bound to a database size.
#[derive(Clone, Debug)]
pub struct EValueParams {
    /// Gumbel scale λ (natural log units per raw score unit).
    pub lambda: f64,
    /// Karlin–Altschul K.
    pub k: f64,
    /// ln(K), cached.
    pub log_k: f64,
    /// Total residues in the (global) target database.
    pub db_residues: f64,
    /// Number of sequences in the (global) target database.
    pub db_seqs: f64,
    /// Finite-size length model (None: raw `m·N` search space).
    pub length: Option<LengthModel>,
    /// Human-readable provenance of the parameters.
    pub model: String,
}

impl EValueParams {
    /// Build from explicit parameters for a database of `db_residues` residues
    /// in `db_seqs` sequences.
    pub fn from_params(p: GappedParams, db_residues: usize, db_seqs: usize, model: &str) -> Self {
        EValueParams {
            lambda: p.lambda,
            k: p.k,
            log_k: p.k.ln(),
            db_residues: db_residues.max(1) as f64,
            db_seqs: db_seqs.max(1) as f64,
            length: p.length,
            model: model.to_string(),
        }
    }

    /// Calibrated gapped statistics for sequence search with the given gap
    /// costs and composition-bias setting. Errors if the configuration has no
    /// verified table entry (run `sabertooth calibrate` to add one).
    pub fn sequence(gaps: GapCosts, comp_bias: bool, db_residues: usize, db_seqs: usize) -> Result<Self, String> {
        let e = table_params(gaps, comp_bias).ok_or_else(|| {
            format!(
                "no calibrated statistics for gap costs {}/{} with comp-bias-corr {}; \
                 run `sabertooth calibrate --gap-open {} --gap-extend {}` and add the result to SEQUENCE_TABLE",
                gaps.open, gaps.extend, comp_bias as u8, gaps.open, gaps.extend
            )
        })?;
        Ok(Self::from_params(
            e.params,
            db_residues,
            db_seqs,
            &format!(
                "calibrated gapped (BLOSUM62 {}/{}, comp-bias {})",
                gaps.open,
                gaps.extend,
                if comp_bias { "on" } else { "off" }
            ),
        ))
    }

    /// **Uncalibrated, ungapped** parameters (analytic λ of the matrix with
    /// BLAST's K). Kept only for diagnostics (`info`, `doctor`) and API
    /// compatibility; search does not use it for gapped scores.
    pub fn new(mat: &SubstitutionMatrix, db_residues: usize) -> Self {
        let lambda = solve_lambda(mat);
        let k = 0.041;
        EValueParams {
            lambda,
            k,
            log_k: k.ln(),
            db_residues: db_residues.max(1) as f64,
            db_seqs: 1.0,
            length: None,
            model: "ungapped analytic lambda (diagnostic only)".into(),
        }
    }

    /// **Ungapped** profile parameters `λ = ln2/8` (exact for an ideal log-odds
    /// PSSM without gaps). Diagnostic only; profile search calibrates each PSSM
    /// with [`calibrate_profile`].
    pub fn for_profile(db_residues: usize) -> Self {
        let lambda = std::f64::consts::LN_2 / crate::profile::PROFILE_BIT_FACTOR;
        let k = 0.041;
        EValueParams {
            lambda,
            k,
            log_k: k.ln(),
            db_residues: db_residues.max(1) as f64,
            db_seqs: 1.0,
            length: None,
            model: "ungapped ln2/8 (diagnostic only)".into(),
        }
    }

    /// Reconstruct from raw Karlin–Altschul parameters.
    pub fn from_raw(lambda: f64, k: f64, db_residues: usize) -> Self {
        EValueParams {
            lambda,
            k,
            log_k: k.ln(),
            db_residues: db_residues.max(1) as f64,
            db_seqs: 1.0,
            length: None,
            model: "raw".into(),
        }
    }

    /// The gapped parameters this was built from.
    pub fn params(&self) -> GappedParams {
        GappedParams {
            lambda: self.lambda,
            k: self.k,
            length: self.length,
        }
    }

    /// Bit score from a raw integer alignment score.
    #[inline]
    pub fn bit_score(&self, raw_score: i32) -> f64 {
        (self.lambda * raw_score as f64 - self.log_k) / std::f64::consts::LN_2
    }

    /// Effective search space for a query of length `query_len`.
    pub fn search_space(&self, query_len: usize) -> f64 {
        let m = query_len.max(1) as f64;
        match self.length {
            Some(lm) => {
                let ell = length_adjustment(self.k, lm.alpha_d_lambda, lm.beta, m, self.db_residues, self.db_seqs);
                let m_eff = (m - ell).max(1.0 / self.k).max(1.0);
                let n_eff = (self.db_residues - self.db_seqs * ell).max(self.db_residues * 0.01).max(1.0);
                m_eff * n_eff
            }
            None => m * self.db_residues,
        }
    }

    /// E-value for a raw score given the query length.
    #[inline]
    pub fn evalue(&self, raw_score: i32, query_len: usize) -> f64 {
        self.evalue_in_space(raw_score, self.search_space(query_len))
    }

    /// E-value for a raw score in a precomputed search space.
    #[inline]
    pub fn evalue_in_space(&self, raw_score: i32, space: f64) -> f64 {
        self.k * space * (-self.lambda * raw_score as f64).exp()
    }

    /// Smallest raw score whose E-value is ≤ `max_evalue` in `space`.
    pub fn min_score_for(&self, max_evalue: f64, space: f64) -> i32 {
        let s = ((self.k * space / max_evalue).ln() / self.lambda).ceil();
        if s.is_finite() {
            s.max(1.0).min(i32::MAX as f64 / 2.0) as i32
        } else {
            1
        }
    }
}

/// Solve `Σ p_i p_j exp(λ s_ij) = 1` for the ungapped λ (bisection).
fn solve_lambda(mat: &SubstitutionMatrix) -> f64 {
    let n = PROFILE_AA_SIZE;
    let p: Vec<f64> = (0..n).map(|i| mat.p_back[i]).collect();
    let f = |lambda: f64| -> f64 {
        let mut acc = 0.0;
        for i in 0..n {
            for j in 0..n {
                acc += p[i] * p[j] * (lambda * mat.score[i][j] as f64).exp();
            }
        }
        acc - 1.0
    };
    let mut lo = 1e-6;
    let mut hi = 2.0;
    if f(lo) > 0.0 {
        return 0.3;
    }
    let mut iters = 0;
    while f(hi) < 0.0 && iters < 60 {
        hi *= 2.0;
        iters += 1;
    }
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

// ---------------------------------------------------------------------------
// Simulation-based calibration (sequence scoring)
// ---------------------------------------------------------------------------

/// Fit of one simulated (query length, target length) cell.
#[derive(Clone, Debug)]
pub struct CalibrationCell {
    pub m: usize,
    pub n: usize,
    pub samples: usize,
    pub mean: f64,
    pub sd: f64,
    /// Free ML fit of this cell.
    pub lambda_ml: f64,
    pub mu_ml: f64,
    /// ML location with the pooled λ.
    pub mu: f64,
}

/// Result of [`calibrate_sequence`].
#[derive(Clone, Debug)]
pub struct Calibration {
    pub params: GappedParams,
    pub cells: Vec<CalibrationCell>,
    pub gaps: GapCosts,
    pub comp_bias: bool,
    pub seed: u64,
    pub samples_per_cell: usize,
    /// RMS of ln K across cells after the length fit (goodness of fit).
    pub ln_k_rms: f64,
}

/// Simulated optimal local scores of `samples` random background pairs of
/// lengths `(m, n)`. Deterministic for a given seed regardless of thread count.
pub fn simulate_scores(
    mat: &SubstitutionMatrix,
    gaps: GapCosts,
    comp_bias: bool,
    m: usize,
    n: usize,
    samples: usize,
    seed: u64,
    cell_id: u64,
) -> Vec<f64> {
    let bg = BackgroundSampler::new(mat);
    (0..samples)
        .into_par_iter()
        .map_init(SimdScratch::new, |scratch, i| {
            let mut rng = SplitMix64::new(stream_seed(seed, cell_id, i as u64));
            let q = bg.sequence(&mut rng, m);
            let t = bg.sequence(&mut rng, n);
            let sc = SeqScorer { query: &q, mat };
            let qp = if comp_bias {
                let bias = crate::compbias::local_bias_correction(&q, mat);
                QueryProfile::build(&ScoreTable::with_bias(&sc, Some(&bias)))
            } else {
                QueryProfile::build(&sc)
            };
            qp.sw_score_with(&t, gaps, scratch).0 as f64
        })
        .collect()
}

/// Fit λ, K and the length model from simulated cells.
///
/// * λ: joint ML over the cells with `min(m, n) ≥ lambda_min_len` (shared λ,
///   per-cell μ); short cells carry finite-size distortion of the shape.
/// * Integer scores: the continuous fit to integer data approximates
///   `P(S ≥ x) = 1 − F(x − ½)`, so each cell's location enters as `μ + ½`.
/// * Length model: for each candidate `(α/λ, β)` the per-cell
///   `ℓ_c = (α/λ)·λ·μ'_c + β` and `K_c = e^{λμ'_c}/((m−ℓ_c)(n−ℓ_c))`; the pair
///   minimising the variance of `ln K_c` is chosen and `K = exp(mean ln K_c)`.
pub fn fit_calibration(mut cells: Vec<(usize, usize, Vec<f64>)>, lambda_min_len: usize) -> Option<(GappedParams, Vec<CalibrationCell>, f64)> {
    cells.sort_by_key(|c| (c.0, c.1));
    // pooled λ by profile likelihood over the long cells
    let long: Vec<&(usize, usize, Vec<f64>)> = cells.iter().filter(|c| c.0.min(c.1) >= lambda_min_len).collect();
    let pool: Vec<&(usize, usize, Vec<f64>)> = if long.is_empty() { cells.iter().collect() } else { long };
    let loglik = |lam: f64| -> f64 {
        let mut ll = 0.0;
        for c in &pool {
            let mu = gumbel_mu_given_lambda(&c.2, lam);
            for &x in &c.2 {
                let z = lam * (x - mu);
                ll += lam.ln() - z - (-z).exp();
            }
        }
        ll
    };
    // golden-section search on λ ∈ [0.05, 1.5]
    let (mut a, mut b) = (0.05f64, 1.5f64);
    let gr = (5f64.sqrt() - 1.0) / 2.0;
    let mut c1 = b - gr * (b - a);
    let mut c2 = a + gr * (b - a);
    let (mut f1, mut f2) = (loglik(c1), loglik(c2));
    for _ in 0..80 {
        if f1 < f2 {
            a = c1;
            c1 = c2;
            f1 = f2;
            c2 = a + gr * (b - a);
            f2 = loglik(c2);
        } else {
            b = c2;
            c2 = c1;
            f2 = f1;
            c1 = b - gr * (b - a);
            f1 = loglik(c1);
        }
    }
    let lambda = 0.5 * (a + b);

    let mut out = Vec::new();
    for (m, n, xs) in &cells {
        let mean = xs.iter().sum::<f64>() / xs.len() as f64;
        let sd = (xs.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / (xs.len() as f64 - 1.0)).sqrt();
        let free = fit_gumbel_ml(xs)?;
        out.push(CalibrationCell {
            m: *m,
            n: *n,
            samples: xs.len(),
            mean,
            sd,
            lambda_ml: free.lambda,
            mu_ml: free.mu,
            mu: gumbel_mu_given_lambda(xs, lambda) + 0.5,
        });
    }
    // grid search for the length model
    let mut best: Option<(f64, f64, f64, f64)> = None; // (var, a, beta, mean lnK)
    let mut adl = 0.0f64;
    while adl <= 20.0 + 1e-9 {
        let mut beta = -80.0f64;
        while beta <= 20.0 + 1e-9 {
            let mut lnks = Vec::with_capacity(out.len());
            let mut ok = true;
            for c in &out {
                let ell = (adl * lambda * c.mu + beta).max(0.0);
                let me = c.m as f64 - ell;
                let ne = c.n as f64 - ell;
                if me < 1.0 || ne < 1.0 {
                    ok = false;
                    break;
                }
                lnks.push(lambda * c.mu - (me * ne).ln());
            }
            if ok && !lnks.is_empty() {
                let mean = lnks.iter().sum::<f64>() / lnks.len() as f64;
                let var = lnks.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / lnks.len() as f64;
                if best.map(|b| var < b.0).unwrap_or(true) {
                    best = Some((var, adl, beta, mean));
                }
            }
            beta += 1.0;
        }
        adl += 0.1;
    }
    let (var, adl, beta, mean_lnk) = best?;
    Some((
        GappedParams {
            lambda,
            k: mean_lnk.exp(),
            length: Some(LengthModel {
                alpha_d_lambda: (adl * 10.0).round() / 10.0,
                beta,
            }),
        },
        out,
        var.sqrt(),
    ))
}

/// Default calibration grid: equal query/target lengths.
pub const DEFAULT_CALIBRATION_LENGTHS: &[usize] = &[60, 100, 150, 220, 320, 480, 700, 1000];

/// Calibrate gapped sequence statistics by simulation.
pub fn calibrate_sequence(
    mat: &SubstitutionMatrix,
    gaps: GapCosts,
    comp_bias: bool,
    lengths: &[usize],
    samples: usize,
    seed: u64,
) -> Option<Calibration> {
    let mut cells = Vec::new();
    for (ci, &l) in lengths.iter().enumerate() {
        let xs = simulate_scores(mat, gaps, comp_bias, l, l, samples, seed, ci as u64);
        cells.push((l, l, xs));
    }
    let (params, fitted, rms) = fit_calibration(cells, 200)?;
    Some(Calibration {
        params,
        cells: fitted,
        gaps,
        comp_bias,
        seed,
        samples_per_cell: samples,
        ln_k_rms: rms,
    })
}

// ---------------------------------------------------------------------------
// Per-profile calibration
// ---------------------------------------------------------------------------

/// Per-PSSM statistics: random background targets whose lengths are drawn from
/// the target database's length distribution are aligned to the profile with
/// the same kernel and gap costs as the search, and a Gumbel is fitted by ML.
/// The resulting `E = N_targets · P(S ≥ x)` tail is expressed as
/// Karlin–Altschul `(λ, K)` over the global search space `m·N` (no separate
/// length model: the edge effects of this database's lengths are inside the
/// fit). `lengths` is a sample of target lengths (e.g. every k-th database
/// record).
pub fn calibrate_profile(
    table: &ScoreTable,
    mat: &SubstitutionMatrix,
    gaps: GapCosts,
    lengths: &[usize],
    samples: usize,
    seed: u64,
) -> Option<(GappedParams, GumbelFit)> {
    if lengths.is_empty() || table.is_empty() {
        return None;
    }
    let bg = BackgroundSampler::new(mat);
    let qp = QueryProfile::build(table);
    let xs: Vec<f64> = (0..samples)
        .into_par_iter()
        .map_init(SimdScratch::new, |scratch, i| {
            let mut rng = SplitMix64::new(stream_seed(seed, 0xBEEF, i as u64));
            let n = lengths[rng.below(lengths.len())].max(1);
            let t = bg.sequence(&mut rng, n);
            qp.sw_score_with(&t, gaps, scratch).0 as f64
        })
        .collect();
    let fit = fit_gumbel_ml(&xs)?;
    let mean_len = lengths.iter().sum::<usize>() as f64 / lengths.len() as f64;
    let m = table.len() as f64;
    // tail: E(S) = N_t·e^{−λ(S − μ − ½)} = K·m·(N_t·n̄)·e^{−λS}
    let k = (fit.lambda * (fit.mu + 0.5)).exp() / (m * mean_len);
    Some((
        GappedParams {
            lambda: fit.lambda,
            k,
            length: None,
        },
        fit,
    ))
}

/// Legacy single-length Monte-Carlo calibration (method of moments). Kept for
/// the `calibrate --legacy` report; superseded by [`calibrate_sequence`].
#[derive(Clone, Copy, Debug)]
pub struct GumbelParams {
    pub lambda: f64,
    pub k: f64,
    pub mu: f64,
    pub samples: usize,
    pub mean_score: f64,
    pub std_score: f64,
}

/// Legacy method-of-moments fit at one length (see [`GumbelParams`]).
pub fn calibrate_gapped(mat: &SubstitutionMatrix, gaps: GapCosts, seq_len: usize, num_pairs: usize, seed: u64) -> GumbelParams {
    let scores = simulate_scores(mat, gaps, false, seq_len, seq_len, num_pairs, seed, 0);
    let n = scores.len().max(1) as f64;
    let mean = scores.iter().sum::<f64>() / n;
    let var = if scores.len() > 1 {
        scores.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / (n - 1.0)
    } else {
        1.0
    };
    let std = var.sqrt().max(1e-9);
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gumbel_ml_recovers_known_parameters() {
        let mut rng = SplitMix64::new(1);
        let (lam, mu) = (0.27f64, 31.0f64);
        let xs: Vec<f64> = (0..20000)
            .map(|_| {
                let u = rng.next_f64().max(1e-300);
                mu - (-(u.ln())).ln() / lam
            })
            .collect();
        let f = fit_gumbel_ml(&xs).unwrap();
        assert!((f.lambda - lam).abs() / lam < 0.03, "lambda {}", f.lambda);
        assert!((f.mu - mu).abs() < 0.2, "mu {}", f.mu);
    }

    #[test]
    fn length_adjustment_is_a_fixed_point() {
        let (k, adl, beta) = (0.041f64, 7.1, -30.0);
        for &(m, n, ns) in &[(150.0, 3.5e6, 14000.0), (300.0, 1e9, 3e6), (1000.0, 1e7, 3e4)] {
            let ell = length_adjustment(k, adl, beta, m, n, ns);
            let rhs = adl * (k.ln() + ((m - ell) * (n - ns * ell)).ln()) + beta;
            assert!(ell >= 0.0 && ell < m);
            // ℓ = floor(fixed point): rhs(ℓ) ≥ ℓ and rhs(ℓ+1) < ℓ+1
            assert!(rhs >= ell - 1e-9, "m={} ell={} rhs={}", m, ell, rhs);
            let up = ell + 1.0;
            let rhs_up = adl * (k.ln() + ((m - up) * (n - ns * up)).ln()) + beta;
            assert!(rhs_up < up, "m={} ell={} not maximal", m, ell);
        }
        // tiny problems get no adjustment
        assert_eq!(length_adjustment(k, adl, beta, 5.0, 5.0, 1.0), 0.0);
    }

    #[test]
    fn evalues_are_monotone_and_bits_consistent() {
        let p = EValueParams::sequence(GapCosts::sequence_default(), false, 1_000_000, 4000).unwrap();
        let e_low = p.evalue(40, 300);
        let e_high = p.evalue(200, 300);
        assert!(e_high < e_low);
        assert!(p.bit_score(200) > p.bit_score(40));
        // E = space · 2^(−bits)
        let space = p.search_space(300);
        let e = p.evalue(55, 300);
        let via_bits = space * 2f64.powf(-p.bit_score(55));
        assert!((e / via_bits - 1.0).abs() < 1e-9);
        // min_score_for inverts evalue
        let s = p.min_score_for(1e-3, space);
        assert!(p.evalue_in_space(s, space) <= 1e-3);
        assert!(p.evalue_in_space(s - 1, space) > 1e-3);
    }

    #[test]
    fn v02_statistics_were_overconfident() {
        // The v0.2 defect, quantified: ungapped analytic λ of this integer
        // matrix (0.389) with BLAST's gapped K = 0.041 and no length
        // correction, versus the calibrated gapped model. At the score that the
        // calibrated model puts at E = 1e-3 (300-aa query, 1 M residues in 4,000
        // targets) v0.2 reported an E-value about an order of magnitude smaller.
        let m = SubstitutionMatrix::blosum62();
        let legacy = EValueParams::new(&m, 1_000_000);
        assert!((legacy.lambda - 0.389).abs() < 0.002, "ungapped lambda {}", legacy.lambda);
        let gapped = EValueParams::sequence(GapCosts::sequence_default(), false, 1_000_000, 4000).unwrap();
        assert!(gapped.lambda < legacy.lambda);
        let s = gapped.min_score_for(1e-3, gapped.search_space(300));
        let ratio = gapped.evalue(s, 300) / legacy.evalue(s, 300);
        // about an order of magnitude (9.5-10x for this setting)
        assert!(ratio > 8.0 && ratio < 12.0, "v0.2 E-values {}x too small", ratio);
    }

    #[test]
    fn calibration_is_reproducible_and_thread_independent() {
        let m = SubstitutionMatrix::blosum62();
        let a = simulate_scores(&m, GapCosts::sequence_default(), false, 120, 120, 300, 9, 3);
        let b = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap()
            .install(|| simulate_scores(&m, GapCosts::sequence_default(), false, 120, 120, 300, 9, 3));
        assert_eq!(a, b);
    }

    #[test]
    fn table_matches_fresh_calibration() {
        // an independent seed and fewer samples than the table; tolerances are
        // statistical (λ ±6 %, K within a factor 1.6, E at a typical threshold
        // within a factor 2)
        let m = SubstitutionMatrix::blosum62();
        for comp_bias in [false, true] {
            let gaps = GapCosts::sequence_default();
            let cal = calibrate_sequence(&m, gaps, comp_bias, &[100, 220, 480], 1500, 0xD1FF_5EED).unwrap();
            let t = table_params(gaps, comp_bias).unwrap().params;
            let rel = (cal.params.lambda - t.lambda).abs() / t.lambda;
            assert!(rel < 0.06, "cbc={} lambda fresh {} table {}", comp_bias, cal.params.lambda, t.lambda);
            for &(qlen, n_res, n_seqs) in &[(150usize, 3_500_000usize, 14_000usize), (400, 50_000_000, 150_000)] {
                let a = EValueParams::from_params(t, n_res, n_seqs, "table");
                let b = EValueParams::from_params(cal.params, n_res, n_seqs, "fresh");
                let s = a.min_score_for(1e-3, a.search_space(qlen));
                let r = a.evalue(s, qlen) / b.evalue(s, qlen);
                assert!(r > 0.5 && r < 2.0, "cbc={} q={} ratio {}", comp_bias, qlen, r);
            }
        }
    }

    #[test]
    fn profile_calibration_is_per_profile() {
        let m = SubstitutionMatrix::blosum62();
        let msa_a = ">a\nMKVLLACDEFGHIKLMNPQRSTVWYACDEFGHIK\n>b\nMKILLSCDEFGHLKLMNPQRSTVWFACDEYGHIK\n";
        let msa_b = ">a\nWWWWCCCCWWWWCCCCWWWWCCCCWWWWCCCCWW\n";
        let pa = crate::profile::Profile::from_msa_str(msa_a, &m, Default::default()).unwrap();
        let pb = crate::profile::Profile::from_msa_str(msa_b, &m, Default::default()).unwrap();
        let ta = ScoreTable::from_scorer(&crate::align::ProfileScorer { profile: &pa });
        let tb = ScoreTable::from_scorer(&crate::align::ProfileScorer { profile: &pb });
        let lens: Vec<usize> = (0..50).map(|i| 100 + 7 * i).collect();
        let gaps = GapCosts::profile_default();
        let (ga, _) = calibrate_profile(&ta, &m, gaps, &lens, 1500, 5).unwrap();
        let (gb, _) = calibrate_profile(&tb, &m, gaps, &lens, 1500, 5).unwrap();
        // gapped profile λ lies below the ungapped ln2/8
        assert!(ga.lambda < std::f64::consts::LN_2 / 8.0 * 1.05, "{}", ga.lambda);
        // different PSSMs get materially different statistics
        let ea = EValueParams::from_params(ga, 1_000_000, 4000, "a").evalue(300, pa.query_len);
        let eb = EValueParams::from_params(gb, 1_000_000, 4000, "b").evalue(300, pb.query_len);
        assert!((ea / eb).ln().abs() > 0.5, "ea {} eb {}", ea, eb);
    }
}
