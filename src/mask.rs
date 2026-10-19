//! Low-complexity masking: **soft** for seeding (always), optionally **hard**
//! for scoring (`--mask-scoring`, on by default for sequence queries).
//!
//! Seeds (k-mer matches) that fall in low-complexity regions — homopolymer
//! runs, short tandem repeats, compositionally extreme stretches — hit huge
//! numbers of unrelated sequences, inflate prefilter work and feed false
//! positives into gapped alignment. Masked positions are excluded from seed
//! generation on the query side and from the seed index on the target side.
//! For seeding the original residues are untouched: gapped scoring still sees
//! the full sequences ("soft" masking). Soft masking alone does not stop a
//! candidate found through an ordinary seed from aligning its low-complexity
//! region against a target's (e.g. two `PEV` repeats); composition-bias
//! correction discounts but does not remove such scores. Search therefore
//! also hard-masks query positions in scoring with the stricter
//! [`MaskParams::scoring_default`] (their score rows become the `X` row);
//! `--mask-scoring 0` restores pure soft masking.
//!
//! The detector is a simplified SEG (Wootton & Federhen 1993): Shannon entropy
//! of residue composition in 12-residue windows; windows with entropy ≤ 2.2 bits
//! trigger, and the masked segment extends over neighbouring windows with
//! entropy ≤ 2.5 bits. Unlike NCBI SEG it does not optimise segment ends by
//! the compositional-probability criterion, so segment boundaries can differ by
//! a few residues.

/// SEG-like masking parameters.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "distributed", derive(serde::Serialize, serde::Deserialize))]
pub struct MaskParams {
    pub window: usize,
    pub locut: f64,
    pub hicut: f64,
}

impl Default for MaskParams {
    fn default() -> Self {
        MaskParams {
            window: 12,
            locut: 2.2,
            hicut: 2.5,
        }
    }
}

impl MaskParams {
    /// Stricter, tighter-bounded parameters for **hard** masking of query
    /// positions in scoring (`--mask-scoring`): 15-residue windows with entropy
    /// ≤ 2.3 bits, no extension beyond the triggering windows, so masked
    /// segments stop close to the low-complexity region instead of eating into
    /// flanking domain residues. Chosen on the tuning split: it removes almost
    /// all low-complexity-driven false positives at a 0.3% TP cost on clean
    /// queries (bench/README.md).
    pub fn scoring_default() -> Self {
        MaskParams {
            window: 15,
            locut: 2.3,
            hicut: 2.3,
        }
    }
}

/// Mark low-complexity positions (`true` = masked for seeding).
pub fn low_complexity_mask(seq: &[u8], p: MaskParams) -> Vec<bool> {
    let mut mask = vec![false; seq.len()];
    low_complexity_mask_into(seq, p, &mut mask);
    mask
}

/// As [`low_complexity_mask`], writing into `mask` (resized to `seq.len()`).
/// Returns the number of masked positions.
pub fn low_complexity_mask_into(seq: &[u8], p: MaskParams, mask: &mut Vec<bool>) -> usize {
    let n = seq.len();
    mask.clear();
    mask.resize(n, false);
    let w = p.window;
    if w < 2 || n < w {
        return 0;
    }
    // c·log2(c) table for counts 0..=w
    let mut clog = vec![0f64; w + 1];
    for (c, slot) in clog.iter_mut().enumerate().skip(1) {
        *slot = c as f64 * (c as f64).log2();
    }
    let log_w = (w as f64).log2();
    let nwin = n - w + 1;
    let mut ent = vec![0f64; nwin];
    let mut counts = [0usize; 32];
    let mut sum_clog = 0f64;
    let sym = |r: u8| -> usize { (r as usize).min(31) };
    for &r in &seq[..w] {
        let s = sym(r);
        sum_clog -= clog[counts[s]];
        counts[s] += 1;
        sum_clog += clog[counts[s]];
    }
    ent[0] = log_w - sum_clog / w as f64;
    for i in 1..nwin {
        let out = sym(seq[i - 1]);
        sum_clog -= clog[counts[out]];
        counts[out] -= 1;
        sum_clog += clog[counts[out]];
        let inn = sym(seq[i + w - 1]);
        sum_clog -= clog[counts[inn]];
        counts[inn] += 1;
        sum_clog += clog[counts[inn]];
        ent[i] = log_w - sum_clog / w as f64;
    }
    // trigger + extension over windows
    let mut i = 0;
    let mut masked_any = false;
    while i < nwin {
        if ent[i] <= p.locut + 1e-12 {
            let mut lo = i;
            while lo > 0 && ent[lo - 1] <= p.hicut + 1e-12 {
                lo -= 1;
            }
            let mut hi = i;
            while hi + 1 < nwin && ent[hi + 1] <= p.hicut + 1e-12 {
                hi += 1;
            }
            for m in mask.iter_mut().take(hi + w).skip(lo) {
                *m = true;
            }
            masked_any = true;
            i = hi + 1;
        } else {
            i += 1;
        }
    }
    if masked_any {
        mask.iter().filter(|&&b| b).count()
    } else {
        0
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
    fn homopolymers_and_short_repeats_are_masked() {
        let s = enc(&format!("MKVLLACDEFGHIK{}MKVLLACDEFGHIK", "Q".repeat(30)));
        let m = low_complexity_mask(&s, MaskParams::default());
        assert!(m[14..44].iter().all(|&b| b), "poly-Q not masked");
        assert!(!m[..8].iter().any(|&b| b), "diverse flank masked");
        let r = enc(&format!("MKVLLACDEFGHIKLMNPQRSTVWY{}ACDEFGHIKLMNPQRSTVWY", "GPP".repeat(10)));
        let m = low_complexity_mask(&r, MaskParams::default());
        let masked: usize = m.iter().filter(|&&b| b).count();
        assert!(masked >= 30, "collagen-like repeat masked only {} positions", masked);
    }

    #[test]
    fn ordinary_protein_is_mostly_unmasked() {
        // a real kinase fragment (demo/targets.fasta, true_homolog_diverged)
        let s = enc("MADENKLKLGSGSFGEVFLVKHKESGNIYAMKKIRLEVEEEGVPSTAIREISLLKELKHENIVRLYDVVHSDKKLTLVFEYLDKDLKKYMDSCPEFAKDPRLIKMFLYQILRGLAYCHRRKVLHRDLKPQNLLIDRRTNSLKLADFGLARAFGVPLRTYTHEVVTLWYRAPEILLG");
        let m = low_complexity_mask(&s, MaskParams::default());
        let frac = m.iter().filter(|&&b| b).count() as f64 / m.len() as f64;
        assert!(frac < 0.1, "masked fraction {}", frac);
    }

    #[test]
    fn short_sequences_are_untouched() {
        let s = enc("QQQQQ");
        assert!(low_complexity_mask(&s, MaskParams::default()).iter().all(|&b| !b));
    }
}
