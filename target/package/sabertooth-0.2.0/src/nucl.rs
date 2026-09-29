//! Nucleotide handling and translated search, mirroring MMseqs2's translated
//! modes (`extractorfs` + `translatenucs` + protein search).
//!
//! Protein search is far more sensitive than nucleotide search for remote
//! homology, so MMseqs2 (like BLASTx) translates nucleotide sequences in all six
//! reading frames, breaks them into ORFs at stop codons, and searches the ORF
//! **proteins** against a protein database. This module provides the standard
//! genetic code, reverse-complement, six-frame ORF extraction, and a helper that
//! packages the ORFs of a nucleotide database as a protein query `SeqDb` so the
//! existing protein search path can be reused unchanged.

use crate::fasta::{Record, SeqDb};

/// NCBI standard genetic code (transl_table = 1), amino-acid letters in codon
/// order with base index T=0, C=1, A=2, G=3 and `codon = 16·b1 + 4·b2 + b3`.
/// `*` marks a stop codon.
const AAS: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

/// Base index in the T,C,A,G ordering used by [`AAS`]; `None` for ambiguous.
#[inline]
fn base_idx(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'T' | b'U' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

/// True if the byte is a recognised nucleotide (incl. `U`, `N`, and IUPAC
/// ambiguity codes). Used to decide whether an input FASTA is nucleotide.
#[inline]
pub fn is_nucleotide_byte(b: u8) -> bool {
    matches!(
        b.to_ascii_uppercase(),
        b'A' | b'C' | b'G' | b'T' | b'U' | b'N' | b'R' | b'Y' | b'S' | b'W' | b'K'
            | b'M' | b'B' | b'D' | b'H' | b'V'
    )
}

/// Heuristic: does this raw sequence look like DNA/RNA? (≥ 90% ACGTUN.)
pub fn looks_like_nucleotide(seq: &[u8]) -> bool {
    if seq.is_empty() {
        return false;
    }
    let acgt = seq
        .iter()
        .filter(|&&b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'U' | b'N'))
        .count();
    acgt * 10 >= seq.len() * 9
}

/// Translate a single codon to an amino-acid letter, or `b'*'` for a stop and
/// `b'X'` if any base is ambiguous.
#[inline]
pub fn translate_codon(c0: u8, c1: u8, c2: u8) -> u8 {
    match (base_idx(c0), base_idx(c1), base_idx(c2)) {
        (Some(a), Some(b), Some(c)) => AAS[16 * a + 4 * b + c],
        _ => b'X',
    }
}

/// Reverse-complement a nucleotide sequence (IUPAC-aware for the common bases).
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' | b'U' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'N' => b'N',
            _ => b'N',
        })
        .collect()
}

/// An open reading frame extracted from a nucleotide sequence.
pub struct Orf {
    /// Frame: +1/+2/+3 forward, −1/−2/−3 reverse.
    pub frame: i8,
    /// 0-based nucleotide start of the ORF **on the forward strand**.
    pub nt_start: usize,
    /// Protein letters (uppercase; no stops).
    pub protein: Vec<u8>,
}

/// Translate one strand at a given frame offset and split into ORFs at stops.
/// `strand_len` is the length of the strand being translated (used to map
/// reverse-frame coordinates back onto the forward strand).
fn orfs_in_frame(
    strand: &[u8],
    offset: usize,
    frame: i8,
    forward_len: usize,
    min_aa: usize,
    out: &mut Vec<Orf>,
) {
    let mut protein: Vec<u8> = Vec::new();
    let mut seg_start_nt = offset; // nt index on `strand` where current ORF began
    let mut i = offset;
    while i + 3 <= strand.len() {
        let aa = translate_codon(strand[i], strand[i + 1], strand[i + 2]);
        if aa == b'*' {
            if protein.len() >= min_aa {
                push_orf(&protein, seg_start_nt, frame, forward_len, out);
            }
            protein.clear();
            seg_start_nt = i + 3;
        } else {
            protein.push(aa);
        }
        i += 3;
    }
    if protein.len() >= min_aa {
        push_orf(&protein, seg_start_nt, frame, forward_len, out);
    }
}

fn push_orf(protein: &[u8], seg_start_nt: usize, frame: i8, forward_len: usize, out: &mut Vec<Orf>) {
    // For reverse frames, seg_start_nt is on the reverse strand; map back to a
    // forward-strand coordinate for reporting.
    let nt_start = if frame > 0 {
        seg_start_nt
    } else {
        forward_len.saturating_sub(seg_start_nt + protein.len() * 3)
    };
    out.push(Orf {
        frame,
        nt_start,
        protein: protein.to_vec(),
    });
}

/// Extract all six-frame ORFs of length ≥ `min_aa` from a nucleotide sequence.
pub fn six_frame_orfs(nt: &[u8], min_aa: usize) -> Vec<Orf> {
    let mut out = Vec::new();
    let len = nt.len();
    for off in 0..3 {
        orfs_in_frame(nt, off, (off + 1) as i8, len, min_aa, &mut out);
    }
    let rc = revcomp(nt);
    for off in 0..3 {
        orfs_in_frame(&rc, off, -((off + 1) as i8), len, min_aa, &mut out);
    }
    out
}

/// Metadata linking an ORF query record back to its nucleotide source.
pub struct OrfMeta {
    pub source_id: String,
    pub frame: i8,
    pub nt_start: usize,
    pub aa_len: usize,
}

/// Build a protein query `SeqDb` from a nucleotide database by six-frame ORF
/// extraction, so the ordinary protein search path can run over it. Each ORF
/// record is named `"<source>_f<frame>_<nt_start>"`. Returns the db and parallel
/// metadata (same order as `db.records`).
pub fn orf_query_db(nt_db: &SeqDb, aa2num: &[u8; 256], min_aa: usize) -> (SeqDb, Vec<OrfMeta>) {
    let mut records = Vec::new();
    let mut meta = Vec::new();
    let mut total = 0usize;
    for rec in &nt_db.records {
        for orf in six_frame_orfs(&rec.seq, min_aa) {
            let num: Vec<u8> = orf.protein.iter().map(|&b| aa2num[b as usize]).collect();
            let id = format!("{}_f{}_{}", rec.id, orf.frame, orf.nt_start);
            total += num.len();
            meta.push(OrfMeta {
                source_id: rec.id.clone(),
                frame: orf.frame,
                nt_start: orf.nt_start,
                aa_len: num.len(),
            });
            records.push(Record {
                header: id.clone(),
                id,
                seq: orf.protein,
                num,
            });
        }
    }
    (
        SeqDb {
            records,
            total_residues: total,
        },
        meta,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::build_aa2num;

    #[test]
    fn standard_codons_translate() {
        assert_eq!(translate_codon(b'A', b'T', b'G'), b'M'); // Met / start
        assert_eq!(translate_codon(b'T', b'T', b'T'), b'F'); // Phe
        assert_eq!(translate_codon(b'A', b'T', b'A'), b'I'); // Ile (standard, not mito M)
        assert_eq!(translate_codon(b'T', b'G', b'G'), b'W'); // Trp
        assert_eq!(translate_codon(b'T', b'A', b'A'), b'*'); // ochre stop
        assert_eq!(translate_codon(b'T', b'G', b'A'), b'*'); // opal stop
        assert_eq!(translate_codon(b'G', b'G', b'G'), b'G'); // Gly
        assert_eq!(translate_codon(b'N', b'N', b'N'), b'X'); // ambiguous
    }

    #[test]
    fn revcomp_is_correct() {
        assert_eq!(revcomp(b"ATGC"), b"GCAT");
        assert_eq!(revcomp(b"AAAA"), b"TTTT");
    }

    #[test]
    fn forward_frame_translates_orf() {
        // ATG AAA TTT GGG TAA  -> M K F G (stop)
        let nt = b"ATGAAATTTGGGTAA";
        let orfs = six_frame_orfs(nt, 3);
        // the +1 frame ORF "MKFG" must be present
        let has = orfs.iter().any(|o| o.frame == 1 && o.protein == b"MKFG");
        assert!(has, "frames: {:?}", orfs.iter().map(|o| (o.frame, String::from_utf8_lossy(&o.protein).to_string())).collect::<Vec<_>>());
    }

    #[test]
    fn reverse_frame_is_found() {
        // revcomp(CATGGG...) etc. Build a sequence whose reverse strand has a
        // clean ORF: forward "TTA CCC AAA TTT CAT" -> reverse strand starts ATG.
        let fwd = b"ATGAAATTTGGGTAA";
        let rc = revcomp(fwd); // TTACCCAAATTTCAT
        let orfs = six_frame_orfs(&rc, 3);
        // some reverse-frame ORF should recover MKFG (since rc of rc = fwd)
        let has = orfs.iter().any(|o| o.frame < 0 && o.protein == b"MKFG");
        assert!(has, "no reverse ORF recovered MKFG");
    }

    #[test]
    fn orf_db_builds_protein_records() {
        let t = build_aa2num();
        let nt = SeqDb {
            records: vec![Record {
                header: "seq1".into(),
                id: "seq1".into(),
                seq: b"ATGAAATTTGGGTAA".to_vec(),
                num: vec![],
            }],
            total_residues: 0,
        };
        let (db, meta) = orf_query_db(&nt, &t, 3);
        assert!(!db.records.is_empty());
        assert_eq!(db.records.len(), meta.len());
        // the MKFG ORF should be encoded to protein numerics
        let mkfg: Vec<u8> = "MKFG".bytes().map(|b| t[b as usize]).collect();
        assert!(db.records.iter().any(|r| r.num == mkfg));
    }
}
