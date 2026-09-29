//! Substitution matrix handling — a faithful port of MMseqs2's
//! `SubstitutionMatrix` / `BaseMatrix` reconstruction.
//!
//! The `.out` matrix files ship *log-odds scores* `S_ab` (in half-bits) plus an
//! optional precomputed background distribution and Karlin–Altschul `lambda` in
//! the header. From these, MMseqs2 reconstructs the joint probability matrix and
//! then derives everything else. Sabertooth reproduces the same chain so profile
//! pseudocounts and integer scores are numerically identical:
//!
//! 1. `P(a,b) = exp(lambda * S_ab) * pBack[a] * pBack[b]`      (joint prob)
//! 2. `R[a][b] = P(a,b) / pBack[b]              = P(a | b)`     (pseudocount matrix)
//! 3. `score[i][j] = round(bitFactor * log2(P/(pa*pb)) + bias)` (integer scores)
//!
//! `ANY_BACK = 1e-5` and the `xIsPositive == false` background rescale are
//! reproduced verbatim.

use crate::alphabet::{ALPHABET_SIZE, NUM2AA, PROFILE_AA_SIZE};

/// MMseqs2 `BaseMatrix::ANY_BACK`.
const ANY_BACK: f64 = 1e-5;

/// The BLOSUM62 matrix data, byte-for-byte from MMseqs2's `data/blosum62.out`.
pub const BLOSUM62_OUT: &str = include_str!("data/blosum62.out");

/// A fully reconstructed substitution matrix.
pub struct SubstitutionMatrix {
    /// Human-readable name (e.g. `blosum62`).
    pub name: String,
    /// ASCII byte -> internal residue index (256-entry table).
    pub aa2num: [u8; 256],
    /// Background frequencies, indexed by residue (length `ALPHABET_SIZE`).
    pub p_back: [f64; ALPHABET_SIZE],
    /// Karlin–Altschul lambda from the matrix header (natural-log scale).
    pub lambda: f64,
    /// Joint probability matrix `P(a,b)` (`ALPHABET_SIZE^2`).
    pub prob: [[f64; ALPHABET_SIZE]; ALPHABET_SIZE],
    /// Pseudocount matrix `R[a][b] = P(a|b)` over the 20 real amino acids.
    pub r: [[f32; PROFILE_AA_SIZE]; PROFILE_AA_SIZE],
    /// Integer score matrix used for alignment (`bit_factor`, `bias` applied).
    pub score: [[i8; ALPHABET_SIZE]; ALPHABET_SIZE],
    /// The bit factor used to build `score` (MMseqs2 default 2.0).
    pub bit_factor: f64,
    /// The score bias used to build `score` (MMseqs2 default -0.2).
    pub bias: f64,
}

impl SubstitutionMatrix {
    /// Load the embedded BLOSUM62 with MMseqs2's default search scaling
    /// (`bit_factor = 2.0`, `bias = -0.2`).
    pub fn blosum62() -> Self {
        Self::from_out_str("blosum62", BLOSUM62_OUT, 2.0, -0.2)
    }

    /// Parse a matrix in MMseqs2 `.out` format and reconstruct all derived
    /// quantities. `bit_factor` and `bias` control the integer `score` matrix.
    pub fn from_out_str(name: &str, data: &str, bit_factor: f64, bias: f64) -> Self {
        let aa2num = crate::alphabet::build_aa2num();

        // Column order as written in the file's header row of letters.
        let mut col_letters: Vec<u8> = Vec::new();
        // Raw file scores, keyed by (row_residue_idx, col_position).
        let mut raw: [[f64; ALPHABET_SIZE]; ALPHABET_SIZE] = [[0.0; ALPHABET_SIZE]; ALPHABET_SIZE];
        let mut p_back = [0.0f64; ALPHABET_SIZE];
        let mut lambda = 0.0f64;
        let mut has_lambda = false;
        let mut has_background = false;
        let mut matrix_started = false;

        for line in data.lines() {
            let trimmed = line.trim_end();
            if trimmed.is_empty() {
                continue;
            }
            if let Some(rest) = trimmed.strip_prefix("# Background (precomputed optional):") {
                for (i, tok) in rest.split_whitespace().enumerate() {
                    if i < ALPHABET_SIZE {
                        p_back[i] = tok.parse().unwrap_or(0.0);
                    }
                }
                has_background = true;
                continue;
            }
            if let Some(rest) = trimmed.strip_prefix("# Lambda     (precomputed optional):") {
                lambda = rest.trim().parse().unwrap_or(0.0);
                has_lambda = true;
                continue;
            }
            if trimmed.starts_with('#') {
                continue;
            }

            let tokens: Vec<&str> = trimmed.split_whitespace().collect();
            if !matrix_started {
                // The header row of column letters: all single alphabetic tokens.
                if tokens.iter().all(|t| t.len() == 1 && t.as_bytes()[0].is_ascii_alphabetic()) {
                    col_letters = tokens.iter().map(|t| aa2num[t.as_bytes()[0] as usize]).collect();
                    matrix_started = true;
                }
                continue;
            }

            // A score row: first token is the row residue letter.
            let row_byte = tokens[0].as_bytes()[0];
            if !row_byte.is_ascii_alphabetic() {
                continue;
            }
            let row = aa2num[row_byte as usize] as usize;
            for (j, tok) in tokens[1..].iter().enumerate() {
                if j < col_letters.len() {
                    let col = col_letters[j] as usize;
                    raw[row][col] = tok.parse().unwrap_or(0.0);
                }
            }
        }

        assert!(has_lambda && has_background,
            "matrix header must provide precomputed Lambda and Background (Sabertooth does not re-estimate them)");

        // Determine whether X is "positive" (MMseqs2 xIsPositive check).
        let x = ALPHABET_SIZE - 1; // index 20
        let mut x_is_positive = false;
        for j in 0..ALPHABET_SIZE {
            if raw[x][j] > 0.0 || raw[j][x] > 0.0 {
                x_is_positive = true;
                break;
            }
        }

        // Background rescale for the non-positive-X case (verbatim MMseqs2).
        if !x_is_positive {
            p_back[x] = ANY_BACK;
            for i in 0..PROFILE_AA_SIZE {
                p_back[i] *= 1.0 - p_back[x];
            }
        }

        // Reconstruct joint probability matrix: P(a,b) = exp(lambda*S) * pa * pb.
        let mut prob = [[0.0f64; ALPHABET_SIZE]; ALPHABET_SIZE];
        for i in 0..ALPHABET_SIZE {
            for j in 0..ALPHABET_SIZE {
                prob[i][j] = (lambda * raw[i][j]).exp() * p_back[i] * p_back[j];
            }
        }

        // R[a][b] = P(a,b) / pBack[b] = P(a|b), over the 20 real amino acids.
        let mut r = [[0.0f32; PROFILE_AA_SIZE]; PROFILE_AA_SIZE];
        for a in 0..PROFILE_AA_SIZE {
            for b in 0..PROFILE_AA_SIZE {
                r[a][b] = (prob[a][b] / p_back[b]) as f32;
            }
        }

        // Integer score matrix: round(bit_factor * log2(P/(pa*pb)) + bias).
        let mut score = [[0i8; ALPHABET_SIZE]; ALPHABET_SIZE];
        for i in 0..ALPHABET_SIZE {
            for j in 0..ALPHABET_SIZE {
                let s = prob[i][j] / (p_back[i] * p_back[j]);
                let val = bit_factor * s.log2() + bias;
                let rounded = if val < 0.0 { val - 0.5 } else { val + 0.5 };
                score[i][j] = rounded.clamp(-128.0, 127.0) as i8;
            }
        }

        SubstitutionMatrix {
            name: name.to_string(),
            aa2num,
            p_back,
            lambda,
            prob,
            r,
            score,
            bit_factor,
            bias,
        }
    }

    /// Score for a pair of internal residue indices.
    #[inline]
    pub fn score(&self, a: u8, b: u8) -> i32 {
        self.score[a as usize][b as usize] as i32
    }

    /// Background probability of a residue as `f32` (profile math is `f32`).
    #[inline]
    pub fn p_back_f32(&self, a: usize) -> f32 {
        self.p_back[a] as f32
    }

    /// Render the integer score matrix as a text table (used by `info`).
    pub fn score_table(&self) -> String {
        let mut out = String::from("     ");
        for &l in &NUM2AA {
            out.push_str(&format!("{:>4}", l as char));
        }
        out.push('\n');
        for i in 0..ALPHABET_SIZE {
            out.push_str(&format!("{:>3}  ", NUM2AA[i] as char));
            for j in 0..ALPHABET_SIZE {
                out.push_str(&format!("{:>4}", self.score[i][j]));
            }
            out.push('\n');
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::aa_to_num;

    #[test]
    fn blosum62_diagonal_is_positive_and_sane() {
        let m = SubstitutionMatrix::blosum62();
        // W-W is the strongest self-score in BLOSUM62.
        let w = aa_to_num(b'W', &m.aa2num);
        let a = aa_to_num(b'A', &m.aa2num);
        assert!(m.score(w, w) > m.score(a, a));
        // Self-scores are all positive.
        for i in 0..PROFILE_AA_SIZE as u8 {
            assert!(m.score(i, i) > 0, "residue {} self-score not positive", i);
        }
    }

    #[test]
    fn r_columns_are_conditional_distributions() {
        let m = SubstitutionMatrix::blosum62();
        // Sum_a P(a|b) should be ~1 for each b.
        for b in 0..PROFILE_AA_SIZE {
            let s: f32 = (0..PROFILE_AA_SIZE).map(|a| m.r[a][b]).sum();
            assert!((s - 1.0).abs() < 0.02, "column {} sums to {}", b, s);
        }
    }

    #[test]
    fn background_sums_to_one() {
        let m = SubstitutionMatrix::blosum62();
        let s: f64 = (0..PROFILE_AA_SIZE).map(|i| m.p_back[i]).sum();
        assert!((s - 1.0).abs() < 0.01, "background sums to {}", s);
    }
}
