//! Linear-time clustering in the spirit of MMseqs2 `linclust`.
//!
//! The linclust idea that makes clustering near-linear: instead of comparing all
//! pairs, reduce each sequence to a handful of **minimizer** k-mers (the lowest-
//! hashing distinct k-mers). Sequences that share a minimizer are grouped; within
//! a group the **longest** sequence is the centre, and every other member is
//! aligned only to that centre. Passing edges (identity + coverage above the
//! thresholds) feed a greedy, longest-first set-cover that emits clusters.
//!
//! This is the algorithmic core, not the full MMseqs2 cascade (no reduced
//! alphabet, no multi-round `clusthash`/`rescorediagonal` refinement). Output is
//! the same representative→member TSV MMseqs2 produces via `createtsv`.

use crate::align::{local_align, GapCosts, SeqScorer};
use crate::fasta::SeqDb;
use crate::matrix::SubstitutionMatrix;
use std::collections::HashMap;

const A: u64 = crate::alphabet::PROFILE_AA_SIZE as u64;

/// Clustering parameters.
#[derive(Clone, Copy)]
pub struct ClusterParams {
    /// K-mer length for minimizer selection.
    pub k: usize,
    /// Number of minimizer k-mers kept per sequence.
    pub num_minimizers: usize,
    /// Minimum sequence identity (0..1) of a member to its centre.
    pub min_seq_id: f64,
    /// Minimum coverage (aligned length / shorter sequence, 0..1).
    pub min_cov: f64,
}

impl Default for ClusterParams {
    fn default() -> Self {
        ClusterParams {
            k: 6,
            num_minimizers: 20,
            min_seq_id: 0.5,
            min_cov: 0.8,
        }
    }
}

/// One cluster: a representative sequence and the members assigned to it
/// (the representative is always the first member).
pub struct Cluster {
    pub representative: usize,
    pub members: Vec<usize>,
}

/// SplitMix64 finalizer — a good avalanche hash for k-mer codes.
#[inline]
fn hash64(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Select up to `num` distinct minimizer k-mers (smallest hash) from a sequence.
/// K-mers containing an unknown residue (≥ 20) are skipped.
pub fn select_minimizers(seq: &[u8], k: usize, num: usize) -> Vec<u64> {
    if seq.len() < k {
        return Vec::new();
    }
    let mut coded: Vec<(u64, u64)> = Vec::new(); // (hash, code)
    'windows: for start in 0..=(seq.len() - k) {
        let mut code = 0u64;
        for &r in &seq[start..start + k] {
            if r as u64 >= A {
                continue 'windows;
            }
            code = code * A + r as u64;
        }
        coded.push((hash64(code), code));
    }
    coded.sort_unstable();
    coded.dedup_by_key(|&mut (_, c)| c);
    coded.truncate(num);
    coded.into_iter().map(|(_, c)| c).collect()
}

/// Cluster a sequence database. Returns clusters; every sequence appears in
/// exactly one cluster (singletons included).
pub fn linclust(db: &SeqDb, mat: &SubstitutionMatrix, params: ClusterParams) -> Vec<Cluster> {
    let n = db.records.len();
    let gaps = GapCosts::sequence_default();

    // 1. minimizer index: k-mer code -> sequences containing it
    let mut kmer_map: HashMap<u64, Vec<usize>> = HashMap::new();
    for (sid, rec) in db.records.iter().enumerate() {
        for code in select_minimizers(&rec.num, params.k, params.num_minimizers) {
            kmer_map.entry(code).or_default().push(sid);
        }
    }

    // 2. per-group centre (longest member) → candidate member edges.
    //    neighbours[centre] = set of members that pass identity+coverage.
    let mut neighbours: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut seen_pairs: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
    for members in kmer_map.values() {
        if members.len() < 2 {
            continue;
        }
        // centre = longest (ties: smallest id)
        let centre = *members
            .iter()
            .max_by_key(|&&s| (db.records[s].num.len(), usize::MAX - s))
            .unwrap();
        let clen = db.records[centre].num.len();
        let sc = SeqScorer {
            query: &db.records[centre].num,
            mat,
        };
        for &mem in members {
            if mem == centre {
                continue;
            }
            // only evaluate each (centre, mem) once
            if !seen_pairs.insert((centre, mem)) {
                continue;
            }
            let mlen = db.records[mem].num.len();
            if let Some(aln) = local_align(&sc, &db.records[mem].num, gaps) {
                let ident = aln.pct_identity() / 100.0;
                let cov = aln.aln_len as f64 / clen.min(mlen).max(1) as f64;
                if ident >= params.min_seq_id && cov >= params.min_cov {
                    neighbours.entry(centre).or_default().push(mem);
                }
            }
        }
    }

    // 3. greedy set cover, longest-first: a sequence not yet assigned opens a
    //    cluster and absorbs its unassigned passing members.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&s| (usize::MAX - db.records[s].num.len(), s)); // longest first
    let mut assigned = vec![false; n];
    let mut clusters = Vec::new();
    for &s in &order {
        if assigned[s] {
            continue;
        }
        assigned[s] = true;
        let mut members = vec![s];
        if let Some(neigh) = neighbours.get(&s) {
            for &mem in neigh {
                if !assigned[mem] {
                    assigned[mem] = true;
                    members.push(mem);
                }
            }
        }
        clusters.push(Cluster {
            representative: s,
            members,
        });
    }
    clusters
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::build_aa2num;
    use std::io::Cursor;

    fn db_from(fasta: &str) -> SeqDb {
        let t = build_aa2num();
        SeqDb::from_reader(Cursor::new(fasta), &t).unwrap()
    }

    #[test]
    fn minimizers_are_stable_and_bounded() {
        let t = build_aa2num();
        let s: Vec<u8> = "MKVLLACDEFGHIKLMNPQRSTVWY".bytes().map(|b| t[b as usize]).collect();
        let m1 = select_minimizers(&s, 6, 5);
        let m2 = select_minimizers(&s, 6, 5);
        assert_eq!(m1, m2); // deterministic
        assert!(m1.len() <= 5);
        assert!(!m1.is_empty());
    }

    #[test]
    fn near_duplicates_cluster_together() {
        // three near-identical sequences + one unrelated
        let db = db_from(
            ">a\nMKVLLACDEFGHIKLMNPQRSTVWYACDEFGHIK\n\
             >b\nMKVLLACDEFGHLKLMNPQRSTVWYACDEFGHIK\n\
             >c\nMKVLLACDEFGHIKLMNPQRSTVWFACDEFGHIK\n\
             >z\nWYWYWYWYWYWYWYWYWYWYWYWYWYWYWYWYWY\n",
        );
        let m = SubstitutionMatrix::blosum62();
        let clusters = linclust(&db, &m, ClusterParams::default());
        // a/b/c should collapse into one cluster; z stands alone
        let sizes: Vec<usize> = clusters.iter().map(|c| c.members.len()).collect();
        assert!(sizes.contains(&3), "expected a size-3 cluster, got {:?}", sizes);
        // every sequence assigned exactly once
        let total: usize = clusters.iter().map(|c| c.members.len()).sum();
        assert_eq!(total, db.records.len());
        // z is a singleton
        assert!(clusters.iter().any(|c| c.members.len() == 1));
    }

    #[test]
    fn dissimilar_sequences_stay_separate() {
        let db = db_from(
            ">a\nMKVLLACDEFGHIKLMNPQR\n\
             >b\nWYFWYFWYFWYFWYFWYFWY\n\
             >c\nGGGGGGPPPPPPCCCCCCDD\n",
        );
        let m = SubstitutionMatrix::blosum62();
        let clusters = linclust(&db, &m, ClusterParams::default());
        assert_eq!(clusters.len(), 3, "all three should be singletons");
    }

    #[test]
    fn identity_threshold_is_enforced() {
        // two sequences sharing k-mers but below the identity cutoff shouldn't merge
        let db = db_from(
            ">a\nMKVLLACDEFGHIKLMNPQRSTVWYMKVLLACDEF\n\
             >b\nMKVLLACDEFGHIKLMNPQRWWWWWWWWWWWWWWWW\n",
        );
        let m = SubstitutionMatrix::blosum62();
        let mut p = ClusterParams::default();
        p.min_seq_id = 0.9; // strict
        p.min_cov = 0.9;
        let clusters = linclust(&db, &m, p);
        assert_eq!(clusters.len(), 2, "strict identity should keep them separate");
    }
}
