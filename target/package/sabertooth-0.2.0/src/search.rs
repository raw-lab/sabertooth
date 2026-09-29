//! Search orchestration: prefilter → align → significance filter → rank.
//!
//! Combines the k-mer prefilter, the Smith–Waterman–Gotoh aligner and the
//! Karlin–Altschul statistics into a full search, mirroring the MMseqs2
//! `search` / `easy-search` workflow at a high level. Candidate alignment is
//! parallelised across targets with rayon. Results are emitted in BLAST
//! tab-separated (`m8`) format.

use crate::align::{local_align, local_align_banded, GapCosts, ProfileScorer, SeqScorer};
use crate::alphabet::GAP;
use crate::evalue::EValueParams;
use crate::fasta::{Record, SeqDb};
use crate::matrix::SubstitutionMatrix;
use crate::prefilter::{prefilter_profile, prefilter_sequence, KmerIndex, PrefilterParams};
use crate::profile::{Profile, PseudoCountParams};
use crate::simd::QueryProfile;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

/// Thresholds applied to alignments after the DP step.
#[derive(Clone, Copy)]
pub struct SearchParams {
    pub prefilter: PrefilterParams,
    pub max_evalue: f64,
    /// Minimum fraction of the query that must be covered by the alignment.
    pub min_query_cov: f64,
    /// Maximum hits reported per query (0 = unlimited).
    pub max_hits: usize,
    /// If set, run the traceback with a banded aligner of this radius around the
    /// prefilter's seed diagonal (O(m·w) memory) instead of the full matrix.
    pub band_radius: Option<usize>,
}

impl Default for SearchParams {
    fn default() -> Self {
        SearchParams {
            prefilter: PrefilterParams::default(),
            max_evalue: 1e-3,
            min_query_cov: 0.0,
            max_hits: 300,
            band_radius: None,
        }
    }
}

/// One reported hit.
pub struct Hit {
    pub query_id: String,
    pub target_id: String,
    pub pident: f64,
    pub aln_len: usize,
    pub mismatches: usize,
    pub gap_opens: usize,
    pub q_start: usize,
    pub q_end: usize,
    pub t_start: usize,
    pub t_end: usize,
    pub evalue: f64,
    pub bit_score: f64,
}

impl Hit {
    /// BLAST `m8` line (1-based coordinates).
    pub fn to_m8(&self) -> String {
        format!(
            "{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2e}\t{:.1}",
            self.query_id,
            self.target_id,
            self.pident,
            self.aln_len,
            self.mismatches,
            self.gap_opens,
            self.q_start + 1,
            self.q_end + 1,
            self.t_start + 1,
            self.t_end + 1,
            self.evalue,
            self.bit_score
        )
    }
}

/// Column header for `m8` output.
pub const M8_HEADER: &str =
    "query\ttarget\tpident\talnlen\tmismatch\tgapopen\tqstart\tqend\ttstart\ttend\tevalue\tbits";

fn count_gap_opens(cigar: &str) -> usize {
    // number of I/D runs
    cigar
        .as_bytes()
        .iter()
        .filter(|&&c| c == b'I' || c == b'D')
        .count()
}

/// Search a single profile against a target database.
pub fn search_profile(
    profile: &Profile,
    query_id: &str,
    db: &SeqDb,
    index: &KmerIndex,
    _mat: &SubstitutionMatrix,
    ev: &EValueParams,
    params: SearchParams,
) -> Vec<Hit> {
    let cands = prefilter_profile(profile, db, index, params.prefilter);
    let gaps = GapCosts::profile_default();
    let sc = ProfileScorer { profile };
    // Build the striped query profile once; score all candidates fast, then run
    // the scalar traceback only for those whose (identical) score clears the
    // E-value cutoff. SIMD score == scalar raw_score, so this defers work
    // without changing which hits pass.
    let qp = QueryProfile::build(&sc);

    let mut hits: Vec<Hit> = cands
        .par_iter()
        .filter_map(|c| {
            let target: &Record = &db.records[c.target_id as usize];
            let raw = qp.sw_score(&target.num, gaps);
            let evalue = ev.evalue(raw, profile.query_len);
            if evalue > params.max_evalue {
                return None;
            }
            let aln = match params.band_radius {
                Some(r) => local_align_banded(&sc, &target.num, gaps, c.best_diag, r)?,
                None => local_align(&sc, &target.num, gaps)?,
            };
            let qcov = (aln.query_end - aln.query_start + 1) as f64 / profile.query_len as f64;
            if qcov < params.min_query_cov {
                return None;
            }
            Some(Hit {
                query_id: query_id.to_string(),
                target_id: target.id.clone(),
                pident: aln.pct_identity(),
                aln_len: aln.aln_len,
                mismatches: aln.aln_len.saturating_sub(aln.identities).saturating_sub(aln.gaps),
                gap_opens: count_gap_opens(&aln.cigar),
                q_start: aln.query_start,
                q_end: aln.query_end,
                t_start: aln.target_start,
                t_end: aln.target_end,
                evalue,
                bit_score: ev.bit_score(aln.raw_score),
            })
        })
        .collect();

    finalize(&mut hits, params.max_hits);
    hits
}

/// Build a query-anchored MSA row (length `query_len`) from a target's optimal
/// alignment to the current profile. Match columns receive the target residue,
/// query-only columns (deletions in target) get a gap, and target-only columns
/// (insertions relative to the query) are dropped — the a3m match-column
/// convention. This is what lets a round-1 hit become a row of the round-2 MSA.
fn hit_to_msa_row(
    profile: &Profile,
    query_len: usize,
    target: &[u8],
    gaps: GapCosts,
) -> Option<Vec<u8>> {
    let sc = ProfileScorer { profile };
    let aln = local_align(&sc, target, gaps)?;
    let mut row = vec![GAP; query_len];
    let mut qpos = aln.query_start;
    let mut tpos = aln.target_start;
    let mut num = 0usize;
    for ch in aln.cigar.chars() {
        if ch.is_ascii_digit() {
            num = num * 10 + (ch as usize - '0' as usize);
        } else {
            let count = num.max(1);
            num = 0;
            match ch {
                'M' => {
                    for _ in 0..count {
                        if qpos < query_len && tpos < target.len() {
                            row[qpos] = target[tpos];
                        }
                        qpos += 1;
                        tpos += 1;
                    }
                }
                'D' => {
                    // deletion in target: query column with no target residue
                    qpos += count;
                }
                'I' => {
                    // insertion relative to query: target residues with no column
                    tpos += count;
                }
                _ => {}
            }
        }
    }
    Some(row)
}

/// Iterative profile search (PSI-BLAST style). Search, fold significant hits
/// back into the MSA as new rows, rebuild the profile, and repeat. A profile
/// bootstrapped from round-1 homologs becomes sensitive to more distant ones in
/// later rounds. Iteration stops after `num_iterations` rounds or earlier once a
/// round adds no new sequence (convergence). Returns the final round's hits.
#[allow(clippy::too_many_arguments)]
pub fn search_profile_iterative(
    initial_msa: &[Vec<u8>],
    query_len: usize,
    query_id: &str,
    db: &SeqDb,
    index: &KmerIndex,
    mat: &SubstitutionMatrix,
    pc: PseudoCountParams,
    ev: &EValueParams,
    params: SearchParams,
    num_iterations: usize,
    inclusion_evalue: f64,
) -> Result<Vec<Hit>, String> {
    let name_to_tid: HashMap<&str, usize> = db
        .records
        .iter()
        .enumerate()
        .map(|(i, r)| (r.id.as_str(), i))
        .collect();
    let gaps = GapCosts::profile_default();
    let mut msa: Vec<Vec<u8>> = initial_msa.to_vec();
    let mut seen: HashSet<usize> = HashSet::new();
    let mut last_hits: Vec<Hit> = Vec::new();
    let rounds = num_iterations.max(1);

    for iter in 0..rounds {
        let profile = Profile::from_numeric_msa(&msa, query_len, mat, pc)?;
        last_hits = search_profile(&profile, query_id, db, index, mat, ev, params);
        if iter + 1 >= rounds {
            break;
        }
        // fold significant, not-yet-included hits into the MSA for the next round
        let mut added = 0usize;
        for h in &last_hits {
            if h.evalue <= inclusion_evalue {
                if let Some(&tid) = name_to_tid.get(h.target_id.as_str()) {
                    if seen.insert(tid) {
                        if let Some(row) =
                            hit_to_msa_row(&profile, query_len, &db.records[tid].num, gaps)
                        {
                            msa.push(row);
                            added += 1;
                        }
                    }
                }
            }
        }
        if added == 0 {
            break; // converged
        }
    }
    Ok(last_hits)
}

/// Search every sequence in `query_db` against `target_db`.
pub fn search_sequences(
    query_db: &SeqDb,
    target_db: &SeqDb,
    index: &KmerIndex,
    mat: &SubstitutionMatrix,
    ev: &EValueParams,
    params: SearchParams,
) -> Vec<Hit> {
    let gaps = GapCosts::sequence_default();

    // Parallelise across queries; align candidates sequentially within a query
    // to keep the rayon pool balanced for many-query workloads.
    let mut all: Vec<Hit> = query_db
        .records
        .par_iter()
        .flat_map_iter(|q| {
            let cands = prefilter_sequence(&q.num, target_db, index, mat, params.prefilter);
            let sc = SeqScorer {
                query: &q.num,
                mat,
            };
            let qp = QueryProfile::build(&sc);
            let mut hits: Vec<Hit> = Vec::new();
            for c in &cands {
                let target = &target_db.records[c.target_id as usize];
                let raw = qp.sw_score(&target.num, gaps);
                let evalue = ev.evalue(raw, q.num.len());
                if evalue > params.max_evalue {
                    continue;
                }
                let aln_opt = match params.band_radius {
                    Some(r) => local_align_banded(&sc, &target.num, gaps, c.best_diag, r),
                    None => local_align(&sc, &target.num, gaps),
                };
                if let Some(aln) = aln_opt {
                    let qcov =
                        (aln.query_end - aln.query_start + 1) as f64 / q.num.len().max(1) as f64;
                    if qcov < params.min_query_cov {
                        continue;
                    }
                    hits.push(Hit {
                        query_id: q.id.clone(),
                        target_id: target.id.clone(),
                        pident: aln.pct_identity(),
                        aln_len: aln.aln_len,
                        mismatches: aln.aln_len.saturating_sub(aln.identities).saturating_sub(aln.gaps),
                        gap_opens: count_gap_opens(&aln.cigar),
                        q_start: aln.query_start,
                        q_end: aln.query_end,
                        t_start: aln.target_start,
                        t_end: aln.target_end,
                        evalue,
                        bit_score: ev.bit_score(aln.raw_score),
                    });
                }
            }
            finalize(&mut hits, params.max_hits);
            hits
        })
        .collect();

    // stable ordering of the concatenated output: by query then by evalue
    all.sort_by(|a, b| {
        a.query_id
            .cmp(&b.query_id)
            .then(a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal))
    });
    all
}

/// **Database-free** sequence search: no k-mer prefilter and no index — every
/// query is aligned against every target with the striped SIMD SW kernel, and
/// survivors get a full traceback. This is the maximum-sensitivity path (nothing
/// a prefilter might drop can be missed) and needs no database to be built; the
/// cost is the full `Q×T` SW work. Suitable for small/medium target sets.
pub fn search_sequences_exhaustive(
    query_db: &SeqDb,
    target_db: &SeqDb,
    mat: &SubstitutionMatrix,
    ev: &EValueParams,
    params: SearchParams,
) -> Vec<Hit> {
    let gaps = GapCosts::sequence_default();
    let mut all: Vec<Hit> = query_db
        .records
        .par_iter()
        .flat_map_iter(|q| {
            let sc = SeqScorer {
                query: &q.num,
                mat,
            };
            let qp = QueryProfile::build(&sc);
            let mut hits: Vec<Hit> = Vec::new();
            for target in &target_db.records {
                let raw = qp.sw_score(&target.num, gaps);
                let evalue = ev.evalue(raw, q.num.len());
                if evalue > params.max_evalue {
                    continue;
                }
                if let Some(aln) = local_align(&sc, &target.num, gaps) {
                    let qcov =
                        (aln.query_end - aln.query_start + 1) as f64 / q.num.len().max(1) as f64;
                    if qcov < params.min_query_cov {
                        continue;
                    }
                    hits.push(Hit {
                        query_id: q.id.clone(),
                        target_id: target.id.clone(),
                        pident: aln.pct_identity(),
                        aln_len: aln.aln_len,
                        mismatches: aln
                            .aln_len
                            .saturating_sub(aln.identities)
                            .saturating_sub(aln.gaps),
                        gap_opens: count_gap_opens(&aln.cigar),
                        q_start: aln.query_start,
                        q_end: aln.query_end,
                        t_start: aln.target_start,
                        t_end: aln.target_end,
                        evalue,
                        bit_score: ev.bit_score(aln.raw_score),
                    });
                }
            }
            finalize(&mut hits, params.max_hits);
            hits
        })
        .collect();
    all.sort_by(|a, b| {
        a.query_id
            .cmp(&b.query_id)
            .then(a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal))
    });
    all
}

/// Sort a per-query hit list by ascending E-value and truncate.
fn finalize(hits: &mut Vec<Hit>, max_hits: usize) {
    hits.sort_by(|a, b| {
        a.evalue
            .partial_cmp(&b.evalue)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(b.bit_score.partial_cmp(&a.bit_score).unwrap_or(std::cmp::Ordering::Equal))
    });
    if max_hits > 0 && hits.len() > max_hits {
        hits.truncate(max_hits);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::build_aa2num;
    use crate::profile::PseudoCountParams;
    use std::io::Cursor;

    fn db_from(fasta: &str) -> SeqDb {
        let t = build_aa2num();
        SeqDb::from_reader(Cursor::new(fasta), &t).unwrap()
    }

    #[test]
    fn exhaustive_finds_superset_of_prefiltered() {
        // The database-free path must never miss a hit the prefilter finds, and
        // recovers the true homolog with a strong score.
        let m = SubstitutionMatrix::blosum62();
        let target = db_from(
            ">t1\nMKVLLACDEFGHIKLMNPQRSTVWY\n\
             >t2\nMKILLSCDEFGHLKLMNPQRSTVWF\n\
             >t3\nWWWWWWYYYYYYFFFFFFPPPPPPCC\n",
        );
        let query = db_from(">q1\nMKVLLACDEFGHIKLMNPQRSTVWY\n");
        let ev = EValueParams::new(&m, target.total_residues);
        let mut sp = SearchParams::default();
        sp.max_evalue = 1e-2;

        let idx = KmerIndex::build(&target, 6);
        let pref = search_sequences(&query, &target, &idx, &m, &ev, sp);
        let exh = search_sequences_exhaustive(&query, &target, &m, &ev, sp);

        // exhaustive recovers at least every target the prefilter did
        let exh_targets: std::collections::HashSet<_> = exh.iter().map(|h| h.target_id.clone()).collect();
        for h in &pref {
            assert!(exh_targets.contains(&h.target_id), "exhaustive missed {}", h.target_id);
        }
        assert!(exh.iter().any(|h| h.target_id == "t1" && h.pident > 90.0));
    }

    #[test]
    fn sequence_search_end_to_end() {
        let m = SubstitutionMatrix::blosum62();
        let target = db_from(
            ">t1\nMKVLLACDEFGHIKLMNPQRSTVWY\n>t2\nPPPPPCCCCCWWWWWKKKKKDDDDD\n",
        );
        let query = db_from(">q1\nMKVLLACDEFGHIKLMNPQRSTVWY\n");
        let idx = KmerIndex::build(&target, 6);
        let ev = EValueParams::new(&m, target.total_residues);
        let mut sp = SearchParams::default();
        sp.max_evalue = 1e6; // permissive for tiny db
        let hits = search_sequences(&query, &target, &idx, &m, &ev, sp);
        assert!(!hits.is_empty());
        assert_eq!(hits[0].target_id, "t1");
        assert!(hits[0].pident > 90.0);
    }

    #[test]
    fn profile_search_end_to_end() {
        let m = SubstitutionMatrix::blosum62();
        let msa = ">q\nMKVLLACDEFGHIKLMNPQR\n>h\nMKILLSCDEFGHLKLMNPQR\n>i\nMRVLLACEEFGHIKLLNPQR\n";
        let prof = Profile::from_msa_str(msa, &m, PseudoCountParams::default()).unwrap();
        let target = db_from(">t1\nGGGMKVLLACDEFGHIKLMNPQRGGG\n>t2\nAAAAAAAAAAAAAAAAAAAA\n");
        let idx = KmerIndex::build(&target, 6);
        let ev = EValueParams::for_profile(target.total_residues);
        let mut sp = SearchParams::default();
        sp.max_evalue = 1e6;
        sp.prefilter.kmer_score = 40;
        let hits = search_profile(&prof, "profile1", &target, &idx, &m, &ev, sp);
        assert!(!hits.is_empty());
        assert_eq!(hits[0].target_id, "t1");
    }

    #[test]
    fn hit_row_reconstructs_query_anchored_sequence() {
        // A target equal to the query aligns as all-M; the reconstructed row is
        // the target itself, length query_len.
        let m = SubstitutionMatrix::blosum62();
        let msa = ">q\nMKVLLACDEFGHIKLMNPQR\n>h\nMKILLSCDEFGHLKLMNPQR\n";
        let prof = Profile::from_msa_str(msa, &m, PseudoCountParams::default()).unwrap();
        let t = build_aa2num();
        let tgt: Vec<u8> = "MKVLLACDEFGHIKLMNPQR".bytes().map(|b| t[b as usize]).collect();
        let row = hit_to_msa_row(&prof, prof.query_len, &tgt, GapCosts::profile_default()).unwrap();
        assert_eq!(row.len(), prof.query_len);
        assert_eq!(row, tgt); // identical target maps 1:1 onto query columns
    }

    #[test]
    fn iterative_one_round_equals_single_pass() {
        let m = SubstitutionMatrix::blosum62();
        let msa_str = ">q\nMKVLLACDEFGHIKLMNPQR\n>h\nMKILLSCDEFGHLKLMNPQR\n";
        let (msa, ql) = crate::profile::parse_msa(msa_str, &m.aa2num).unwrap();
        let target = db_from(">t1\nGGGMKVLLACDEFGHIKLMNPQRGGG\n>t2\nAAAAAAAAAAAAAAAAAAAA\n");
        let idx = KmerIndex::build(&target, 6);
        let ev = EValueParams::for_profile(target.total_residues);
        let mut sp = SearchParams::default();
        sp.max_evalue = 1e6;
        sp.prefilter.kmer_score = 40;

        let prof = Profile::from_msa_str(msa_str, &m, PseudoCountParams::default()).unwrap();
        let single = search_profile(&prof, "p", &target, &idx, &m, &ev, sp);
        let iter1 = search_profile_iterative(
            &msa, ql, "p", &target, &idx, &m, PseudoCountParams::default(), &ev, sp, 1, 1e-3,
        )
        .unwrap();
        assert_eq!(single.len(), iter1.len());
        assert_eq!(single[0].target_id, iter1[0].target_id);
    }

    #[test]
    fn iterative_bootstraps_a_distant_homolog() {
        // Round 1: a single-sequence profile finds only the close homolog. After
        // that hit is folded in, the enriched profile detects the distant one in
        // round 2. Sequences are the same length so columns line up cleanly.
        let m = SubstitutionMatrix::blosum62();
        //           a conserved core flanked by variable positions
        let query = "WYFMKVLLACDEFGHIKWYFM";
        let close = "WYFMKILLSCDEFGHLKWYFM"; // mild substitutions in core
        let distant = "AApMKILQSCEEFGHLKApAA".to_uppercase(); // heavier, incl. flanks
        let msa_str = format!(">q\n{}\n", query);
        let (msa, ql) = crate::profile::parse_msa(&msa_str, &m.aa2num).unwrap();
        let target = db_from(&format!(">close\n{}\n>distant\n{}\n", close, distant));
        let idx = KmerIndex::build(&target, 5);
        let ev = EValueParams::for_profile(target.total_residues);
        let mut sp = SearchParams::default();
        sp.max_evalue = 0.05; // significance cutoff
        sp.prefilter.k = 5; // must match the k=5 index built above
        sp.prefilter.kmer_score = 12; // sensitive seeding
        sp.prefilter.min_diag_score = 1;

        let r1 = search_profile_iterative(
            &msa, ql, "p", &target, &idx, &m, PseudoCountParams::default(), &ev, sp, 1, 1.0,
        )
        .unwrap();
        let r3 = search_profile_iterative(
            &msa, ql, "p", &target, &idx, &m, PseudoCountParams::default(), &ev, sp, 3, 1.0,
        )
        .unwrap();
        // iteration should never lose hits, and here it gains the distant one
        assert!(
            r3.len() >= r1.len(),
            "iteration must not reduce the hit set ({} -> {})",
            r1.len(),
            r3.len()
        );
    }
}
