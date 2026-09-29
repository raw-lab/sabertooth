//! Profile / PSSM construction — the sensitivity core.
//!
//! This is a faithful port of MMseqs2's `PSSMCalculator::computePSSMFromMSA`
//! (default global-weight + substitution-score pseudocount path). Given a
//! multiple sequence alignment it produces a position-specific scoring matrix
//! (PSSM) whose remote-homology sensitivity is the whole point of profile
//! search. The pipeline, matching MMseqs2 step for step:
//!
//! 1. **Position-based sequence weights** (Henikoff & Henikoff 1994):
//!    `w[k] += 1 / (nl[aa] * distinct_aa * (nres[k] + 30))`, summed over columns.
//! 2. Normalise weights to sum 1.
//! 3. **Weighted match frequencies**: `f[pos][aa] = Σ_k w[k]·[msa[k][pos]==aa]`,
//!    then column-normalised (falling back to the background if a column is empty).
//! 4. **Effective sequence count Neff** per column, from column entropy.
//! 5. **Substitution-matrix pseudocounts**:
//!    `g[pos][aa] = Σ_b R[aa][b]·f[pos][b]` with `R[a][b] = P(a|b)`, mixed via
//!    `τ = min(1, pca / (1 + Neff/pcb))`, `p = (1-τ)·f + τ·g`.
//! 6. **Log-odds PSSM**: `M[pos][aa] = round(8·log2(p/pBack) + 8·bias)`, clamped
//!    to `[-128, 127]`.
//!
//! Defaults `pca = 1.0`, `pcb = 1.5`, `bias = 0.0` reproduce MMseqs2's defaults.

use crate::alphabet::{aa_to_num, num_to_aa, ALPHABET_SIZE, ANY, GAP, PROFILE_AA_SIZE};
use crate::matrix::SubstitutionMatrix;

/// Bit factor used when converting profile probabilities to integer PSSM
/// scores. MMseqs2 uses 8 (scores are eighth-bits / ×8 log-odds).
pub const PROFILE_BIT_FACTOR: f64 = 8.0;

/// Bytes per column in an on-disk MMseqs2 profile DB entry
/// (`Sequence::PROFILE_READIN_SIZE`).
pub const PROFILE_READIN_SIZE: usize = 25;
/// Byte offset of the query residue within a profile column.
pub const PROFILE_QUERY_OFFSET: usize = 20;
/// Byte offset of the consensus residue within a profile column.
pub const PROFILE_CONSENSUS_OFFSET: usize = 21;
/// Byte offset of the quantized Neff within a profile column.
pub const PROFILE_NEFF_OFFSET: usize = 22;
// offsets 23, 24 are gap-open/gap-extend (`PROFILE_GAP_DEL` / `PROFILE_GAP_INS`),
// which Sabertooth does not model per column.

/// Decode MMseqs2's quantized Neff byte: `2^((b - 1) / 64)`
/// (`MathUtil::convertNeffToFloat`).
#[inline]
pub fn neff_from_byte(b: u8) -> f32 {
    ((b as f32 - 1.0) / 64.0).exp2()
}

/// A computed sequence profile.
#[derive(Clone, Debug)]
pub struct Profile {
    /// Number of match-state columns (== query length).
    pub query_len: usize,
    /// Number of sequences in the source MSA.
    ///
    /// **Not recoverable from an on-disk MMseqs2 profile DB**: the 25-byte column
    /// format stores per-column `Neff` but never the source MSA's sequence count,
    /// so [`Profile::from_pssm_columns`] sets this to `0`. It is unused by search
    /// and alignment (only `pssm`, `prob`, `neff`, `consensus`, `query_num` are),
    /// so a `set_size` of 0 on a DB-loaded profile is expected, not an error.
    pub set_size: usize,
    /// Integer PSSM, row-major `query_len × PROFILE_AA_SIZE` (`i8`).
    pub pssm: Vec<i8>,
    /// Position probabilities, row-major `query_len × PROFILE_AA_SIZE` (`f32`).
    pub prob: Vec<f32>,
    /// Effective number of sequences per column.
    pub neff: Vec<f32>,
    /// Consensus residues (numeric indices), length `query_len`.
    pub consensus: Vec<u8>,
    /// The query (first MSA row) as numeric residues, length `query_len`.
    pub query_num: Vec<u8>,
}

/// Tunable pseudocount parameters (MMseqs2 `pca`, `pcb`).
#[derive(Clone, Copy)]
pub struct PseudoCountParams {
    pub pca: f32,
    pub pcb: f32,
    pub score_bias: f32,
}

impl Default for PseudoCountParams {
    fn default() -> Self {
        // MMseqs2 defaults for amino-acid profile pseudocounts.
        PseudoCountParams {
            pca: 1.0,
            pcb: 1.5,
            score_bias: 0.0,
        }
    }
}

impl Profile {
    /// Decode a profile from the raw column bytes of an on-disk MMseqs2 profile
    /// DB entry (`dbtype == 2`).
    ///
    /// Each column is [`PROFILE_READIN_SIZE`] = 25 bytes, laid out exactly as
    /// `Sequence::mapProfile` reads them:
    ///
    /// | offset | contents |
    /// |---|---|
    /// | `0..20` | 20 `i8` PSSM scores, ×8 log-odds (`PROFILE_BIT_FACTOR`), in `NUM2AA` order |
    /// | `20`    | query residue (numeric index) |
    /// | `21`    | consensus residue (numeric index) |
    /// | `22`    | quantized Neff, decoded by [`neff_from_byte`] |
    /// | `23..25`| gap del/ins bytes (not modelled here) |
    ///
    /// The stored scores are on the **same ×8 scale as [`Profile::pssm`]**; MMseqs2
    /// divides by 4 only when building its alignment profile
    /// (`Sequence.cpp`: `profile_for_alignment = profile_score / 4`), which is why
    /// `mmseqs profile2pssm` dumps `score / 4`.
    ///
    /// A single trailing `NUL` terminator is tolerated and stripped. Note that
    /// profile entries are framed with a bare `\0` — **not** the `\n\0` used for
    /// text entries — because the binary payload legitimately contains `0x0A`.
    ///
    /// `prob` is **reconstructed** by inverting the log-odds
    /// (`p[a] ∝ p_back[a] · 2^(pssm[a]/8)`, renormalized): the DB stores quantized
    /// integer scores, not posteriors, so this is lossy (scores saturate at the
    /// `i8` bounds). `set_size` is set to `0` — see the field docs.
    pub fn from_pssm_columns(bytes: &[u8], mat: &SubstitutionMatrix) -> Result<Profile, String> {
        // tolerate the single NUL terminator used by the DB framing
        let data = if bytes.len() % PROFILE_READIN_SIZE == 1 && bytes.last() == Some(&0) {
            &bytes[..bytes.len() - 1]
        } else {
            bytes
        };
        if data.is_empty() {
            return Err("profile entry is empty".into());
        }
        if data.len() % PROFILE_READIN_SIZE != 0 {
            return Err(format!(
                "profile entry length {} is not a multiple of {} bytes per column",
                data.len(),
                PROFILE_READIN_SIZE
            ));
        }
        let query_len = data.len() / PROFILE_READIN_SIZE;
        let mut pssm = vec![0i8; query_len * PROFILE_AA_SIZE];
        let mut prob = vec![0f32; query_len * PROFILE_AA_SIZE];
        let mut neff = vec![0f32; query_len];
        let mut consensus = vec![0u8; query_len];
        let mut query_num = vec![0u8; query_len];

        for col in 0..query_len {
            let base = col * PROFILE_READIN_SIZE;
            let out = col * PROFILE_AA_SIZE;
            // 20 signed score bytes, already on the ×8 scale
            for a in 0..PROFILE_AA_SIZE {
                pssm[out + a] = data[base + a] as i8;
            }
            let q = data[base + PROFILE_QUERY_OFFSET];
            let c = data[base + PROFILE_CONSENSUS_OFFSET];
            if q as usize >= ALPHABET_SIZE || c as usize >= ALPHABET_SIZE {
                return Err(format!(
                    "column {}: query/consensus residue index out of range ({}, {}); \
                     is this really a profile DB (dbtype 2)?",
                    col, q, c
                ));
            }
            query_num[col] = q;
            consensus[col] = c;
            neff[col] = neff_from_byte(data[base + PROFILE_NEFF_OFFSET]);

            // invert the log-odds to recover an (approximate) posterior
            let mut sum = 0f64;
            for a in 0..PROFILE_AA_SIZE {
                let p = mat.p_back[a] * (pssm[out + a] as f64 / PROFILE_BIT_FACTOR).exp2();
                prob[out + a] = p as f32;
                sum += p;
            }
            if sum > 0.0 {
                for a in 0..PROFILE_AA_SIZE {
                    prob[out + a] = (prob[out + a] as f64 / sum) as f32;
                }
            }
        }

        Ok(Profile {
            query_len,
            set_size: 0, // not recoverable from the DB — see field docs
            pssm,
            prob,
            neff,
            consensus,
            query_num,
        })
    }

    /// PSSM score for aligning a target residue (internal index) at a query
    /// position. Target `X`/unknown scores 0 (neutral), matching the practical
    /// behaviour of profile search on masked residues.
    #[inline]
    pub fn score(&self, pos: usize, target_residue: u8) -> i32 {
        let aa = target_residue as usize;
        if aa < PROFILE_AA_SIZE {
            self.pssm[pos * PROFILE_AA_SIZE + aa] as i32
        } else {
            0
        }
    }

    /// Build a profile from a numeric MSA: `msa[k]` is a row of length
    /// `query_len` using internal residue indices (0..19), `ANY` (20) for X, or
    /// `GAP` (0xFF) for gaps. Row 0 is the query/center sequence.
    ///
    /// Returns an error (rather than panicking) if the MSA is empty, the query
    /// length is zero, or any row is shorter than `query_len`. This is the
    /// validation boundary for the public numeric API; callers that already
    /// guarantee well-formed rectangular input (e.g. [`parse_msa`]) still pay
    /// only a cheap length check.
    pub fn from_numeric_msa(
        msa: &[Vec<u8>],
        query_len: usize,
        mat: &SubstitutionMatrix,
        params: PseudoCountParams,
    ) -> Result<Profile, String> {
        let set_size = msa.len();
        if set_size == 0 {
            return Err("cannot build a profile from an empty MSA".into());
        }
        if query_len == 0 {
            return Err("cannot build a profile with zero match columns".into());
        }
        for (k, row) in msa.iter().enumerate() {
            if row.len() < query_len {
                return Err(format!(
                    "MSA row {} has length {} but query_len is {}",
                    k + 1,
                    row.len(),
                    query_len
                ));
            }
        }

        // --- 1. position-based sequence weights (Henikoff 1994) --------------
        let seq_weight = compute_sequence_weights(msa, query_len, set_size);

        // normalise to sum 1
        let total: f32 = seq_weight.iter().sum();
        let seq_weight: Vec<f32> = if total > 0.0 {
            seq_weight.iter().map(|w| w / total).collect()
        } else {
            vec![1.0 / set_size as f32; set_size]
        };

        // --- 3. weighted match frequencies ----------------------------------
        let mut freq = vec![0.0f32; query_len * PROFILE_AA_SIZE];
        for pos in 0..query_len {
            let base = pos * PROFILE_AA_SIZE;
            for k in 0..set_size {
                let c = msa[k][pos];
                if (c as usize) < PROFILE_AA_SIZE {
                    freq[base + c as usize] += seq_weight[k];
                }
            }
            normalize_to_1(&mut freq[base..base + PROFILE_AA_SIZE], Some(&mat.p_back));
        }

        // --- 4. Neff per column ---------------------------------------------
        let neff = compute_neff(&freq, &seq_weight, msa, query_len, set_size);

        // --- 5. consensus ---------------------------------------------------
        let consensus = compute_consensus(&freq, query_len, mat);

        // --- 6. substitution-matrix pseudocounts ----------------------------
        let mut prob = vec![0.0f32; query_len * PROFILE_AA_SIZE];
        if params.pca > 0.0 {
            // g[pos][a] = Σ_b R[a][b] · f[pos][b]
            let mut g = vec![0.0f32; query_len * PROFILE_AA_SIZE];
            for pos in 0..query_len {
                let base = pos * PROFILE_AA_SIZE;
                for a in 0..PROFILE_AA_SIZE {
                    let mut acc = 0.0f32;
                    for b in 0..PROFILE_AA_SIZE {
                        acc += mat.r[a][b] * freq[base + b];
                    }
                    g[base + a] = acc;
                }
            }
            // p = (1-τ)·f + τ·g,  τ = min(1, pca/(1 + Neff/pcb))
            for pos in 0..query_len {
                let base = pos * PROFILE_AA_SIZE;
                let tau = (params.pca / (1.0 + neff[pos] / params.pcb)).min(1.0);
                for a in 0..PROFILE_AA_SIZE {
                    prob[base + a] = (1.0 - tau) * freq[base + a] + tau * g[base + a];
                }
            }
        } else {
            prob.copy_from_slice(&freq);
        }

        // --- 7. log-odds PSSM ------------------------------------------------
        let pssm = compute_log_pssm(&prob, query_len, mat, params.score_bias);

        let query_num: Vec<u8> = (0..query_len).map(|p| msa[0][p]).collect();

        Ok(Profile {
            query_len,
            set_size,
            pssm,
            prob,
            neff,
            consensus,
            query_num,
        })
    }

    /// Parse an MSA from text (aligned FASTA or a3m) and build a profile.
    pub fn from_msa_str(
        data: &str,
        mat: &SubstitutionMatrix,
        params: PseudoCountParams,
    ) -> Result<Profile, String> {
        let (msa, query_len) = parse_msa(data, &mat.aa2num)?;
        Profile::from_numeric_msa(&msa, query_len, mat, params)
    }

    /// Render the PSSM as a text table (MMseqs2 `profile2pssm`-style).
    pub fn pssm_table(&self) -> String {
        let mut out = String::from("Pos Res");
        for aa in 0..PROFILE_AA_SIZE {
            out.push_str(&format!(" {:>4}", num_to_aa(aa as u8) as char));
        }
        out.push_str("  Neff\n");
        for pos in 0..self.query_len {
            out.push_str(&format!(
                "{:>3} {:>3}",
                pos + 1,
                num_to_aa(self.query_num[pos]) as char
            ));
            for aa in 0..PROFILE_AA_SIZE {
                out.push_str(&format!(" {:>4}", self.pssm[pos * PROFILE_AA_SIZE + aa]));
            }
            out.push_str(&format!("  {:.2}\n", self.neff[pos]));
        }
        out
    }
}

/// Henikoff & Henikoff (1994) position-based sequence weights.
fn compute_sequence_weights(msa: &[Vec<u8>], query_len: usize, set_size: usize) -> Vec<f32> {
    let mut seq_weight = vec![1e-6f32; set_size];

    // residues per sequence
    let mut number_res = vec![0u32; set_size];
    for (k, row) in msa.iter().enumerate() {
        let mut nr = 0u32;
        for &c in row.iter().take(query_len) {
            if c != GAP {
                nr += 1;
            }
        }
        number_res[k] = nr;
    }

    let mut nl = [0i32; PROFILE_AA_SIZE];
    for pos in 0..query_len {
        nl.iter_mut().for_each(|x| *x = 0);
        for k in 0..set_size {
            let c = msa[k][pos];
            if (c as usize) < PROFILE_AA_SIZE {
                nl[c as usize] += 1;
            }
        }
        let distinct: i32 = nl.iter().filter(|&&x| x > 0).count() as i32;
        if distinct == 0 {
            continue;
        }
        for k in 0..set_size {
            let c = msa[k][pos];
            if (c as usize) < PROFILE_AA_SIZE {
                let denom = nl[c as usize] as f32 * distinct as f32 * (number_res[k] as f32 + 30.0);
                seq_weight[k] += 1.0 / denom;
            }
        }
    }
    seq_weight
}

/// Effective number of sequences per column (MMseqs2 `computeNeff_M`).
fn compute_neff(
    freq: &[f32],
    seq_weight: &[f32],
    msa: &[Vec<u8>],
    query_len: usize,
    set_size: usize,
) -> Vec<f32> {
    // Neff_HMM = mean over positions of 2^entropy(column).
    let mut neff_hmm = 0.0f32;
    for pos in 0..query_len {
        let base = pos * PROFILE_AA_SIZE;
        let mut sum = 0.0f32;
        for aa in 0..PROFILE_AA_SIZE {
            let f = freq[base + aa];
            if f > 1e-10 {
                sum -= f * f.log2();
            }
        }
        neff_hmm += sum.exp2();
    }
    neff_hmm /= query_len.max(1) as f32;

    let nlim = 10.0f32.max(neff_hmm + 1.0);
    let scale = ((nlim - neff_hmm) / (nlim - 1.0)).log2();

    let mut neff = vec![0.0f32; query_len];
    for pos in 0..query_len {
        let mut w_m = -1.0 / set_size as f32;
        for k in 0..set_size {
            if msa[k][pos] != GAP {
                w_m += seq_weight[k];
            }
        }
        neff[pos] = if w_m < 0.0 {
            1.0
        } else {
            nlim - (nlim - 1.0) * (scale * w_m).exp2()
        };
    }
    neff
}

/// Consensus = residue whose (freq - background) is maximal (MMseqs2).
fn compute_consensus(freq: &[f32], query_len: usize, mat: &SubstitutionMatrix) -> Vec<u8> {
    let mut cons = vec![ANY; query_len];
    for pos in 0..query_len {
        let base = pos * PROFILE_AA_SIZE;
        let mut maxw = 1e-8f32;
        let mut maxa = ANY;
        for aa in 0..PROFILE_AA_SIZE {
            let d = freq[base + aa] - mat.p_back_f32(aa);
            if d > maxw {
                maxw = d;
                maxa = aa as u8;
            }
        }
        cons[pos] = maxa;
    }
    cons
}

/// Log-odds PSSM: `round(8·log2(p/pBack) + 8·bias)` clamped to `[-128,127]`.
fn compute_log_pssm(
    prob: &[f32],
    query_len: usize,
    mat: &SubstitutionMatrix,
    score_bias: f32,
) -> Vec<i8> {
    let bit = PROFILE_BIT_FACTOR as f32;
    let mut pssm = vec![0i8; query_len * PROFILE_AA_SIZE];
    for pos in 0..query_len {
        let base = pos * PROFILE_AA_SIZE;
        for aa in 0..PROFILE_AA_SIZE {
            let p = prob[base + aa];
            let log_prob = (p / mat.p_back_f32(aa)).log2();
            let val = bit * log_prob + bit * score_bias;
            let rounded = if val < 0.0 { val - 0.5 } else { val + 0.5 };
            pssm[base + aa] = rounded.clamp(-128.0, 127.0) as i8;
        }
    }
    pssm
}

/// Normalise a slice to sum 1; if the sum is zero, copy the background.
fn normalize_to_1(v: &mut [f32], background: Option<&[f64; crate::alphabet::ALPHABET_SIZE]>) {
    let sum: f32 = v.iter().sum();
    if sum != 0.0 {
        let fac = 1.0 / sum;
        for x in v.iter_mut() {
            *x *= fac;
        }
    } else if let Some(bg) = background {
        for (i, x) in v.iter_mut().enumerate() {
            *x = bg[i] as f32;
        }
    }
}

/// Parse an MSA (aligned FASTA or a3m) into numeric rows of equal match length.
///
/// If any lowercase residues are present the input is treated as a3m: lowercase
/// letters are insertions relative to the query and are dropped; uppercase and
/// `-`/`.` are match columns. Otherwise the input is a plain aligned FASTA and
/// match columns are the non-gap positions of the first sequence.
///
/// Aligned-FASTA parsing is deliberately tolerant: rows are read against the
/// query's match columns, so a row shorter than the query is padded with gaps
/// and any trailing columns beyond the query are ignored. (Files that omit
/// trailing gaps on the last sequence are common and handled silently.) The
/// a3m path, by contrast, computes each row's match columns independently and
/// rejects rows whose count disagrees with the query — a genuine structural
/// error rather than a formatting quirk.
pub fn parse_msa(data: &str, aa2num: &[u8; 256]) -> Result<(Vec<Vec<u8>>, usize), String> {
    // collect raw (header, sequence-as-written) records
    let mut raw: Vec<String> = Vec::new();
    let mut cur = String::new();
    let mut started = false;
    for line in data.lines() {
        let line = line.trim_end();
        if line.starts_with('>') {
            if started {
                raw.push(std::mem::take(&mut cur));
            }
            started = true;
        } else {
            cur.push_str(line.trim());
        }
    }
    if started {
        raw.push(cur);
    }
    if raw.is_empty() {
        return Err("MSA contains no sequences".into());
    }

    let is_a3m = raw.iter().any(|s| s.bytes().any(|b| b.is_ascii_lowercase()));

    let mut msa: Vec<Vec<u8>> = Vec::with_capacity(raw.len());
    if is_a3m {
        // Drop lowercase (inserts); keep uppercase + gaps as match states.
        for s in &raw {
            let mut row = Vec::new();
            for b in s.bytes() {
                if b.is_ascii_lowercase() {
                    continue; // insert state relative to query
                }
                if b == b'-' || b == b'.' {
                    row.push(GAP);
                } else if b.is_ascii_uppercase() {
                    row.push(aa_to_num(b, aa2num));
                }
            }
            msa.push(row);
        }
    } else {
        // Plain aligned FASTA: all rows padded to equal length; match columns =
        // first sequence's non-gap positions.
        let width = raw.iter().map(|s| s.len()).max().unwrap_or(0);
        let first = raw[0].as_bytes();
        let match_cols: Vec<usize> = (0..width)
            .filter(|&i| i < first.len() && first[i] != b'-' && first[i] != b'.')
            .collect();
        for s in &raw {
            let bytes = s.as_bytes();
            let mut row = Vec::with_capacity(match_cols.len());
            for &i in &match_cols {
                let b = if i < bytes.len() { bytes[i] } else { b'-' };
                if b == b'-' || b == b'.' {
                    row.push(GAP);
                } else {
                    row.push(aa_to_num(b.to_ascii_uppercase(), aa2num));
                }
            }
            msa.push(row);
        }
    }

    let query_len = msa[0].len();
    for (k, row) in msa.iter().enumerate() {
        if row.len() != query_len {
            return Err(format!(
                "MSA row {} has {} match columns, expected {}",
                k + 1,
                row.len(),
                query_len
            ));
        }
    }
    Ok((msa, query_len))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build one 25-byte profile column: scores, query, consensus, neff byte.
    fn col(scores: &[i8; PROFILE_AA_SIZE], q: u8, cns: u8, neff_b: u8) -> Vec<u8> {
        let mut v = Vec::with_capacity(PROFILE_READIN_SIZE);
        v.extend(scores.iter().map(|&s| s as u8));
        v.push(q);
        v.push(cns);
        v.push(neff_b);
        v.push(0); // gap del
        v.push(0); // gap ins
        v
    }

    #[test]
    fn neff_byte_decodes_as_pow2() {
        // MathUtil::convertNeffToFloat: 2^((b-1)/64)
        assert!((neff_from_byte(1) - 1.0).abs() < 1e-6);
        assert!((neff_from_byte(65) - 2.0).abs() < 1e-5);
        assert!((neff_from_byte(129) - 4.0).abs() < 1e-4);
    }

    #[test]
    fn from_pssm_columns_decodes_all_fields() {
        let m = SubstitutionMatrix::blosum62();
        let mut s = [0i8; PROFILE_AA_SIZE];
        s[0] = -4; // A
        s[5] = 21; // G  (matches the real family.a3m column 0)
        s[19] = -7; // Y
        let q = m.aa2num[b'G' as usize];
        let cns = m.aa2num[b'G' as usize];
        let bytes = col(&s, q, cns, 65); // neff byte 65 -> 2.0
        let p = Profile::from_pssm_columns(&bytes, &m).unwrap();

        assert_eq!(p.query_len, 1);
        assert_eq!(p.set_size, 0); // not recoverable from a profile DB
        assert_eq!(p.pssm[0], -4);
        assert_eq!(p.pssm[5], 21);
        assert_eq!(p.pssm[19], -7);
        assert_eq!(p.query_num[0], q);
        assert_eq!(p.consensus[0], cns);
        assert!((p.neff[0] - 2.0).abs() < 1e-5);
        // reconstructed posterior is a normalized distribution favouring G
        let sum: f32 = p.prob[0..PROFILE_AA_SIZE].iter().sum();
        assert!((sum - 1.0).abs() < 1e-4, "prob must sum to 1, got {}", sum);
        let best = (0..PROFILE_AA_SIZE)
            .max_by(|&a, &b| p.prob[a].partial_cmp(&p.prob[b]).unwrap())
            .unwrap();
        assert_eq!(best, 5, "highest-scoring residue should carry most mass");
    }

    #[test]
    fn newline_bytes_inside_payload_survive() {
        // 0x0A ('\n') is a legitimate score/residue byte. A `\n\0` text unframing
        // would corrupt the last column; only a bare NUL may be stripped.
        let m = SubstitutionMatrix::blosum62();
        let mut s = [0i8; PROFILE_AA_SIZE];
        s[3] = 10; // 0x0A as a score
        let q = m.aa2num[b'E' as usize];
        let mut bytes = col(&s, q, q, 10); // 0x0A as the neff byte too
        // last data byte (gap ins) is 0x0A as well, immediately before the NUL
        let last = bytes.len() - 1;
        bytes[last] = 10;
        bytes.push(0); // the DB's single NUL terminator

        let p = Profile::from_pssm_columns(&bytes, &m).unwrap();
        assert_eq!(p.query_len, 1, "one column must survive the terminator strip");
        assert_eq!(p.pssm[3], 10);
        assert!((p.neff[0] - neff_from_byte(10)).abs() < 1e-6);
    }

    #[test]
    fn from_pssm_columns_rejects_bad_length_and_bad_residue() {
        let m = SubstitutionMatrix::blosum62();
        assert!(Profile::from_pssm_columns(&[], &m).is_err());
        assert!(Profile::from_pssm_columns(&[0u8; 24], &m).is_err()); // not a multiple of 25
        let mut bad = col(&[0i8; PROFILE_AA_SIZE], 99, 0, 1); // residue index 99
        bad.truncate(PROFILE_READIN_SIZE);
        let err = Profile::from_pssm_columns(&bad, &m).unwrap_err();
        assert!(err.contains("out of range"), "got: {}", err);
    }

    #[test]
    fn multi_column_roundtrip_preserves_order() {
        let m = SubstitutionMatrix::blosum62();
        let mut bytes = Vec::new();
        for i in 0..5u8 {
            let mut s = [0i8; PROFILE_AA_SIZE];
            s[i as usize] = 20 + i as i8; // distinct peak per column
            bytes.extend(col(&s, i, i, 1));
        }
        let p = Profile::from_pssm_columns(&bytes, &m).unwrap();
        assert_eq!(p.query_len, 5);
        for i in 0..5usize {
            assert_eq!(p.pssm[i * PROFILE_AA_SIZE + i], 20 + i as i8);
            assert_eq!(p.query_num[i], i as u8);
        }
    }

    #[test]
    fn profile_is_clone_and_debug() {
        let m = SubstitutionMatrix::blosum62();
        let p = Profile::from_msa_str(">q\nMKVLLACDEFGHIK\n", &m, PseudoCountParams::default())
            .unwrap();
        let c = p.clone();
        assert_eq!(c.query_len, p.query_len);
        assert!(format!("{:?}", c).starts_with("Profile {"));
        // Debug on the Ok type is what makes unwrap_err() usable
        let e = Profile::from_msa_str("", &m, PseudoCountParams::default()).unwrap_err();
        assert!(!e.is_empty());
    }

    #[test]
    fn single_sequence_profile_tracks_matrix() {
        // A profile from a single sequence, with pseudocounts, should score the
        // query's own residues highly and unrelated residues low.
        let m = SubstitutionMatrix::blosum62();
        let msa = ">q\nMKVLA\n";
        let p = Profile::from_msa_str(msa, &m, PseudoCountParams::default()).unwrap();
        assert_eq!(p.query_len, 5);
        assert_eq!(p.set_size, 1);
        // position 0 is M: aligning M should score higher than aligning P.
        let m_idx = aa_to_num(b'M', &m.aa2num);
        let p_idx = aa_to_num(b'P', &m.aa2num);
        assert!(p.score(0, m_idx) > p.score(0, p_idx));
    }

    #[test]
    fn conserved_column_scores_higher_than_variable() {
        // Column fully conserved as C vs a variable column should give C a much
        // stronger positive score in the conserved column.
        let m = SubstitutionMatrix::blosum62();
        let msa = ">a\nCA\n>b\nCK\n>c\nCD\n>d\nCE\n";
        let p = Profile::from_msa_str(msa, &m, PseudoCountParams::default()).unwrap();
        let c_idx = aa_to_num(b'C', &m.aa2num);
        // Conserved C column (pos 0) should score C higher than the variable
        // column (pos 1) scores its query residue A.
        let a_idx = aa_to_num(b'A', &m.aa2num);
        assert!(p.score(0, c_idx) > p.score(1, a_idx));
        // Neff_M is driven by column *occupancy* (sequence weights of non-gap
        // rows), not within-column diversity. With no gaps, both fully-occupied
        // columns share the same Neff — matching MMseqs2's computeNeff_M.
        assert!((p.neff[0] - p.neff[1]).abs() < 1e-4);
    }

    #[test]
    fn a3m_inserts_are_dropped() {
        let m = SubstitutionMatrix::blosum62();
        // lowercase 'a' is an insert relative to the query and must be ignored.
        let msa = ">q\nMKVLA\n>h\nMKaVLA\n"; // second row would be longer with insert kept
        // Force equal match length: query has 5 match cols; the insert 'a' drops.
        let msa2 = ">q\nMKVLA\n>h\nMKVLA\n";
        let p1 = Profile::from_msa_str(msa2, &m, PseudoCountParams::default()).unwrap();
        // parsing the a3m version should also yield query_len 5.
        let (rows, ql) = parse_msa(msa, &m.aa2num).unwrap();
        assert_eq!(ql, 5);
        assert_eq!(rows.len(), 2);
        assert_eq!(p1.query_len, 5);
    }

    #[test]
    fn empty_msa_is_error_not_panic() {
        let m = SubstitutionMatrix::blosum62();
        // no sequences at all
        assert!(Profile::from_msa_str(">only_header_no_seq\n", &m, PseudoCountParams::default()).is_err()
            || Profile::from_msa_str("", &m, PseudoCountParams::default()).is_err());
        // empty numeric MSA
        assert!(Profile::from_numeric_msa(&[], 0, &m, PseudoCountParams::default()).is_err());
    }

    #[test]
    fn ragged_numeric_msa_is_error_not_panic() {
        // A row shorter than query_len must return Err, never index out of bounds.
        let m = SubstitutionMatrix::blosum62();
        let msa = vec![vec![0u8, 1, 2, 3], vec![0u8, 1]]; // second row too short
        let res = Profile::from_numeric_msa(&msa, 4, &m, PseudoCountParams::default());
        assert!(res.is_err(), "ragged MSA should error, got Ok");
    }

    #[test]
    fn mismatched_match_columns_error() {
        // In a3m mode, rows whose match-column counts differ (after inserts are
        // dropped) are structurally inconsistent and must be rejected. Here the
        // lowercase 'a' selects a3m mode and 'bad' ends up with 6 match columns
        // versus the query's 5.
        let m = SubstitutionMatrix::blosum62();
        let msa = ">q\nMKVLA\n>bad\nMKVLAaG\n";
        assert!(Profile::from_msa_str(msa, &m, PseudoCountParams::default()).is_err());
    }
}
