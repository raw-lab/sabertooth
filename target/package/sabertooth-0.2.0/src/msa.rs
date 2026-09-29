//! MSA generation from search results — the a3m producer behind MMseqs2's
//! `result2msa`.
//!
//! Given a query and its alignments to target sequences, we emit an **a3m** MSA:
//! the query is the first row and defines the *match columns*; every other row is
//! its target projected onto those columns. Match positions carry the target
//! residue (uppercase), query columns with no target residue (deletions) carry
//! `-`, and residues the target inserts relative to the query carry lowercase
//! letters between the match columns. This is exactly the a3m convention the rest
//! of Sabertooth already consumes in `parse_msa`, so a `result2msa` output can be
//! fed straight back into `msa2profile` / `profilesearch`.

use crate::align::Alignment;
use crate::alphabet::NUM2AA;

#[inline]
fn up(r: u8) -> char {
    if (r as usize) < NUM2AA.len() {
        NUM2AA[r as usize] as char
    } else {
        'X'
    }
}

#[inline]
fn lo(r: u8) -> char {
    up(r).to_ascii_lowercase()
}

/// Project a single hit's target sequence onto the query's match columns,
/// yielding one a3m row. The number of match-state characters (uppercase + `-`)
/// equals `query_len`; insertions add lowercase characters and do not count.
pub fn a3m_row(query_len: usize, target: &[u8], aln: &Alignment) -> String {
    let mut row = String::with_capacity(query_len + 8);
    for _ in 0..aln.query_start {
        row.push('-'); // query columns left of the alignment
    }
    let mut tpos = aln.target_start;
    let mut matched_qcols = 0usize;
    let mut num = 0usize;
    for ch in aln.cigar.chars() {
        if ch.is_ascii_digit() {
            num = num * 10 + (ch as usize - '0' as usize);
            continue;
        }
        let count = num.max(1);
        num = 0;
        match ch {
            'M' => {
                for _ in 0..count {
                    row.push(up(target[tpos]));
                    tpos += 1;
                    matched_qcols += 1;
                }
            }
            'D' => {
                // deletion in target: a query column with no target residue
                for _ in 0..count {
                    row.push('-');
                    matched_qcols += 1;
                }
            }
            'I' => {
                // insertion relative to the query: lowercase, no column consumed
                for _ in 0..count {
                    row.push(lo(target[tpos]));
                    tpos += 1;
                }
            }
            _ => {}
        }
    }
    let filled = aln.query_start + matched_qcols;
    for _ in filled..query_len {
        row.push('-'); // query columns right of the alignment
    }
    row
}

/// Assemble a full a3m MSA: the query as the first, all-uppercase row, then one
/// row per hit. `hits` are `(name, target_residues, alignment)` and should
/// already be filtered/sorted by the caller (typically by E-value).
pub fn build_a3m(query_name: &str, query: &[u8], hits: &[(String, Vec<u8>, Alignment)]) -> String {
    let query_len = query.len();
    let mut out = String::new();
    out.push('>');
    out.push_str(query_name);
    out.push('\n');
    for &r in query {
        out.push(up(r));
    }
    out.push('\n');
    for (name, target, aln) in hits {
        out.push('>');
        out.push_str(name);
        out.push('\n');
        out.push_str(&a3m_row(query_len, target, aln));
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{local_align, GapCosts, SeqScorer};
    use crate::alphabet::build_aa2num;
    use crate::matrix::SubstitutionMatrix;

    fn enc(s: &str) -> Vec<u8> {
        let t = build_aa2num();
        s.bytes().map(|b| t[b as usize]).collect()
    }

    #[test]
    fn identical_hit_is_all_uppercase_match() {
        let m = SubstitutionMatrix::blosum62();
        let q = enc("MKVLLACDEFGHIK");
        let t = enc("MKVLLACDEFGHIK");
        let sc = SeqScorer { query: &q, mat: &m };
        let aln = local_align(&sc, &t, GapCosts::sequence_default()).unwrap();
        let row = a3m_row(q.len(), &t, &aln);
        assert_eq!(row, "MKVLLACDEFGHIK"); // no gaps, no inserts, length == query_len
        // count of match-state chars equals query length
        assert_eq!(row.chars().filter(|c| *c == '-' || c.is_ascii_uppercase()).count(), q.len());
    }

    #[test]
    fn target_insertion_is_lowercase() {
        // target has an extra residue relative to the query -> lowercase insert
        let m = SubstitutionMatrix::blosum62();
        let q = enc("MKVLLACDEFGHIK");
        let t = enc("MKVLLACWDEFGHIK"); // inserted W after position 6
        let sc = SeqScorer { query: &q, mat: &m };
        let aln = local_align(&sc, &t, GapCosts::sequence_default()).unwrap();
        let row = a3m_row(q.len(), &t, &aln);
        // exactly query_len match-state chars; at least one lowercase insert
        let match_states = row.chars().filter(|c| *c == '-' || c.is_ascii_uppercase()).count();
        assert_eq!(match_states, q.len());
        assert!(row.chars().any(|c| c.is_ascii_lowercase()), "row = {}", row);
    }

    #[test]
    fn query_deletion_is_gap() {
        // target is missing a residue the query has -> '-' in the row
        let m = SubstitutionMatrix::blosum62();
        let q = enc("MKVLLACDEFGHIK");
        let t = enc("MKVLLADEFGHIK"); // dropped C at position 6
        let sc = SeqScorer { query: &q, mat: &m };
        let aln = local_align(&sc, &t, GapCosts::sequence_default()).unwrap();
        let row = a3m_row(q.len(), &t, &aln);
        assert_eq!(row.chars().filter(|c| *c == '-' || c.is_ascii_uppercase()).count(), q.len());
        assert!(row.contains('-'), "row = {}", row);
    }

    #[test]
    fn full_a3m_roundtrips_through_parser() {
        // build_a3m output must parse back with the same match-column count
        let m = SubstitutionMatrix::blosum62();
        let q = enc("MKVLLACDEFGHIKLMNPQR");
        let t1 = enc("MKILLSCDEFGHLKLMNPQR");
        let sc = SeqScorer { query: &q, mat: &m };
        let a1 = local_align(&sc, &t1, GapCosts::sequence_default()).unwrap();
        let a3m = build_a3m("query", &q, &[("hit1".into(), t1.clone(), a1)]);
        let (rows, cols) = crate::profile::parse_msa(&a3m, &m.aa2num).unwrap();
        assert_eq!(cols, q.len());
        assert_eq!(rows.len(), 2); // query + 1 hit
        assert!(rows.iter().all(|r| r.len() == q.len()));
    }
}
