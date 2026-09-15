//! K-mer prefiltering.
//!
//! Full Smith–Waterman against every target is wasteful; like MMseqs2, Sabertooth
//! first seeds candidates with short k-mer matches and only aligns the survivors.
//! Crucially, seeding is not limited to *exact* k-mers — for each query k-mer we
//! also enumerate **similar** k-mers scoring above a threshold under the
//! substitution matrix (sequence queries) or the PSSM (profile queries), using a
//! branch-and-bound generator. This "similar-k-mer" expansion is what lets short
//! seeds reach remote homologs, i.e. the sensitivity MMseqs2 is known for.
//!
//! The index maps a k-mer code to the `(target_id, position)` occurrences.
//! Candidate targets are those hit by ≥ `min_hits` seeds; the dominant
//! `target_pos − query_pos` diagonal is recorded for downstream banding/context.

use crate::alphabet::PROFILE_AA_SIZE;
use crate::fasta::SeqDb;
use crate::matrix::SubstitutionMatrix;
use crate::profile::Profile;
use std::collections::HashMap;

/// Number of real residues used for k-mer coding.
const A: u64 = PROFILE_AA_SIZE as u64; // 20

/// Largest supported k-mer length. K-mers are packed into a base-20 `u64`
/// code, so the hard limit is `20^k < 2^64` → k ≤ 14. We cap a little lower to
/// leave headroom and because k > ~9 is not useful for protein sensitivity
/// anyway. Requests outside `MIN_K..=MAX_K` are rejected with a clear error
/// rather than silently overflowing the code (which is what a `u32` code did
/// for k ≥ 8: `20^8 > 2^32`, wrapping to colliding codes with no diagnostic).
pub const MAX_K: usize = 12;
/// Smallest useful k-mer length.
pub const MIN_K: usize = 2;

/// Validate a requested k-mer length, returning a clear error if unsupported.
pub fn check_k(k: usize) -> Result<(), String> {
    if k < MIN_K || k > MAX_K {
        Err(format!(
            "k-mer length {} is out of the supported range {}..={} \
             (codes are packed into a base-20 u64)",
            k, MIN_K, MAX_K
        ))
    } else {
        Ok(())
    }
}

/// A (possibly spaced) seed pattern. `offsets` are the positions within a window
/// of length `span` that contribute residues to the k-mer code; `weight` is the
/// number of such positions (the effective k). A contiguous seed has
/// `offsets = 0..k` and `span == weight == k`.
///
/// Spaced seeds ("care/don't-care" masks like `110101`) improve remote-homology
/// sensitivity: a mismatch that falls on a don't-care position doesn't destroy
/// the seed, so diverged homologs still share exact spaced k-mers even when they
/// share no contiguous k-mer of the same weight. The window always begins and
/// ends on a care position so the seed's diagonal (`target_pos − query_pos`) is
/// well defined at the window start.
#[derive(Clone)]
pub struct SeedPattern {
    offsets: Vec<usize>,
    span: usize,
}

impl SeedPattern {
    /// Contiguous seed of length `k` (`span == weight == k`).
    pub fn contiguous(k: usize) -> Self {
        SeedPattern {
            offsets: (0..k).collect(),
            span: k,
        }
    }

    /// Parse a `0/1` mask such as `"110101"`. Must start and end with `1`.
    pub fn from_mask(mask: &str) -> Result<Self, String> {
        let bytes = mask.as_bytes();
        if bytes.is_empty() {
            return Err("empty seed mask".into());
        }
        let mut offsets = Vec::new();
        for (i, &b) in bytes.iter().enumerate() {
            match b {
                b'1' => offsets.push(i),
                b'0' => {}
                _ => return Err(format!("seed mask must be 0/1, got '{}'", b as char)),
            }
        }
        if offsets.is_empty() {
            return Err("seed mask has no care ('1') positions".into());
        }
        if offsets[0] != 0 || *offsets.last().unwrap() != bytes.len() - 1 {
            return Err("seed mask must start and end with '1'".into());
        }
        Ok(SeedPattern {
            offsets,
            span: bytes.len(),
        })
    }

    /// Number of coded positions (the effective k).
    #[inline]
    pub fn weight(&self) -> usize {
        self.offsets.len()
    }

    /// Total window length spanned by the pattern.
    #[inline]
    pub fn span(&self) -> usize {
        self.span
    }

    #[inline]
    fn offsets(&self) -> &[usize] {
        &self.offsets
    }

    /// Render back to a `0/1` mask (for `--info`/diagnostics).
    pub fn mask_string(&self) -> String {
        let mut m = vec![b'0'; self.span];
        for &o in &self.offsets {
            m[o] = b'1';
        }
        String::from_utf8(m).unwrap()
    }
}

/// Map a sensitivity value (higher = more sensitive) to a similar-k-mer score
/// threshold, mirroring MMseqs2's `-s`: a lower threshold admits more (weaker)
/// seeds, reaching more remote homologs at higher cost. `scale` is the score bit
/// factor (2 for sequence matrix scores, 8 for profile PSSM) so the threshold
/// lands at the right magnitude for each. This is a monotone approximation of
/// MMseqs2's target-k-mers-per-position heuristic, not an exact reproduction; the
/// reference point `s ≈ 5.7` reproduces the built-in default threshold.
pub fn kmer_score_for_sensitivity(s: f64, scale: i32) -> i32 {
    let s = s.clamp(1.0, 9.0);
    let base = if scale >= 8 { 40.0 } else { 25.0 };
    let per_unit = 0.75 * scale as f64; // 1.5 at ×2, 6.0 at ×8
    let thr = base + (5.7 - s) * per_unit;
    thr.round().max(4.0) as i32
}

/// Prefilter configuration.
#[derive(Clone, Copy)]
pub struct PrefilterParams {
    /// K-mer length (default 6, as in MMseqs2 for proteins).
    pub k: usize,
    /// Minimum summed score for a similar k-mer to be generated.
    pub kmer_score: i32,
    /// Maximum similar k-mers generated per query position (bounds work).
    pub max_kmers_per_pos: usize,
    /// Minimum number of seed hits for a target to become a candidate.
    pub min_hits: usize,
    /// Minimum ungapped diagonal score for a target to survive the prefilter.
    /// This cheap Kadane gate is what keeps the survivor set small (so the
    /// expensive gapped alignment runs on few targets) without dropping real
    /// homologs — mirroring MMseqs2's ungapped diagonal prefilter.
    pub min_diag_score: i32,
}

impl Default for PrefilterParams {
    fn default() -> Self {
        PrefilterParams {
            k: 6,
            kmer_score: 25,
            max_kmers_per_pos: 4096,
            min_hits: 1,
            min_diag_score: 0,
        }
    }
}

/// A prefilter candidate: a target plus seeding evidence.
pub struct Candidate {
    pub target_id: u32,
    pub hits: u32,
    /// Most frequent diagonal `target_pos - query_pos` (as i64 to allow negative).
    pub best_diag: i64,
    /// Best ungapped segment score along `best_diag` (Kadane).
    pub diag_score: i32,
}

/// An inverted k-mer index over a target database.
pub struct KmerIndex {
    pattern: SeedPattern,
    map: HashMap<u64, Vec<(u32, u32)>>, // kmer_code -> [(target_id, pos)]
}

impl KmerIndex {
    /// Build a contiguous k-mer index of length `k`.
    pub fn build(db: &SeqDb, k: usize) -> KmerIndex {
        Self::build_spaced(db, SeedPattern::contiguous(k))
    }

    /// Build an index using an arbitrary (possibly spaced) seed pattern. Windows
    /// containing an `X`/unknown residue at a care position are skipped (they
    /// carry no specific signal). The recorded position is the window start, so
    /// diagonals remain `target_pos − query_pos` at the seed's left edge.
    pub fn build_spaced(db: &SeqDb, pattern: SeedPattern) -> KmerIndex {
        debug_assert!(
            (MIN_K..=MAX_K).contains(&pattern.weight()),
            "seed weight {} outside supported range {}..={}; call check_k() first",
            pattern.weight(),
            MIN_K,
            MAX_K
        );
        let span = pattern.span();
        let mut map: HashMap<u64, Vec<(u32, u32)>> = HashMap::new();
        for (tid, rec) in db.records.iter().enumerate() {
            let s = &rec.num;
            if s.len() < span {
                continue;
            }
            for start in 0..=(s.len() - span) {
                if let Some(code) = encode_spaced(&s[start..start + span], &pattern) {
                    map.entry(code).or_default().push((tid as u32, start as u32));
                }
            }
        }
        KmerIndex { pattern, map }
    }

    pub fn distinct_kmers(&self) -> usize {
        self.map.len()
    }

    /// The effective k (seed weight) this index was built with.
    #[inline]
    pub fn k(&self) -> usize {
        self.pattern.weight()
    }

    /// The seed pattern this index was built with.
    #[inline]
    pub fn pattern(&self) -> &SeedPattern {
        &self.pattern
    }

    /// Look up a k-mer's occurrences.
    #[inline]
    fn get(&self, code: u64) -> Option<&Vec<(u32, u32)>> {
        self.map.get(&code)
    }
}

/// Encode the care positions of a window (length `pattern.span`) into a base-20
/// code. Returns `None` if any care position is not a real amino acid (>= 20).
#[inline]
fn encode_spaced(window: &[u8], pattern: &SeedPattern) -> Option<u64> {
    let mut code = 0u64;
    for &off in pattern.offsets() {
        let r = window[off];
        if r as u64 >= A {
            return None;
        }
        code = code * A + r as u64;
    }
    Some(code)
}

/// Per-position scoring closure shared by seq/profile similar-k-mer generation.
///
/// `col_scores[i]` gives, for query column `qstart + i`, the score of each of
/// the 20 residues. The generator enumerates residue combinations whose summed
/// score across the window meets the threshold.
fn generate_similar_kmers(
    col_scores: &[[i32; PROFILE_AA_SIZE]],
    threshold: i32,
    cap: usize,
    out: &mut Vec<u64>,
) {
    let k = col_scores.len();

    // Per-position residues sorted by ascending score. Pushing children in
    // ascending order means the highest-scoring child sits on top of the stack
    // and is popped first, giving a best-first traversal: the first `cap`
    // k-mers emitted are (greedily) the top-scoring ones — which always includes
    // the argmax-per-position k-mer (the strongest possible seed). This is what
    // makes capping safe.
    let mut order: Vec<[u8; PROFILE_AA_SIZE]> = vec![[0u8; PROFILE_AA_SIZE]; k];
    for i in 0..k {
        let mut idx: Vec<u8> = (0..PROFILE_AA_SIZE as u8).collect();
        idx.sort_by_key(|&r| col_scores[i][r as usize]); // ascending
        order[i].copy_from_slice(&idx);
    }

    // Suffix maxima: best achievable score from position i to end.
    let mut suffix_max = vec![0i32; k + 1];
    for i in (0..k).rev() {
        let best = col_scores[i].iter().copied().max().unwrap_or(0);
        suffix_max[i] = best + suffix_max[i + 1];
    }

    out.clear();
    // Iterative DFS with explicit stack: (depth, partial_code, partial_score).
    let mut stack: Vec<(usize, u64, i32)> = vec![(0, 0, 0)];
    while let Some((depth, code, score)) = stack.pop() {
        if out.len() >= cap {
            break;
        }
        if depth == k {
            if score >= threshold {
                out.push(code);
            }
            continue;
        }
        // Prune: even the best possible completion cannot reach threshold.
        if score + suffix_max[depth] < threshold {
            continue;
        }
        let base = code * A;
        for &r in order[depth].iter() {
            let ns = score + col_scores[depth][r as usize];
            // Bound with suffix max of the *next* position.
            if ns + suffix_max[depth + 1] >= threshold {
                stack.push((depth + 1, base + r as u64, ns));
            }
        }
    }
}

/// Accumulate seed hits as a flat `(target_id, diagonal)` list.
///
/// The previous design kept a `HashMap<tid, (count, HashMap<diag, count>)>` — a
/// per-target inner hash map allocated and probed on every single seed hit,
/// which dominated prefilter time and allocation on large databases. Here each
/// hit is one `push` into a contiguous vector; candidate reduction is a single
/// sort plus a linear pass. Diagonals fit comfortably in `i32` (|diag| < max
/// sequence length), so each hit is 8 bytes and stays cache-friendly.
struct HitAccumulator {
    hits: Vec<(u32, i32)>, // (target_id, diagonal = target_pos - query_pos)
}

impl HitAccumulator {
    fn new() -> Self {
        HitAccumulator { hits: Vec::new() }
    }

    #[inline]
    fn add(&mut self, index: &KmerIndex, code: u64, query_pos: usize) {
        if let Some(list) = index.get(code) {
            for &(tid, tpos) in list {
                let diag = tpos as i32 - query_pos as i32;
                self.hits.push((tid, diag));
            }
        }
    }

    fn into_candidates(mut self, min_hits: usize) -> Vec<Candidate> {
        if self.hits.is_empty() {
            return Vec::new();
        }
        // Group by target, then by diagonal, in one sort.
        self.hits.sort_unstable();
        let mut cands: Vec<Candidate> = Vec::new();

        let mut i = 0;
        let n = self.hits.len();
        while i < n {
            let tid = self.hits[i].0;
            // span of this target
            let mut j = i;
            while j < n && self.hits[j].0 == tid {
                j += 1;
            }
            let group = &self.hits[i..j];
            let total = group.len() as u32;

            // most frequent diagonal within the (already diag-sorted) group
            let mut best_diag = group[0].1;
            let mut best_run = 0u32;
            let mut run_diag = group[0].1;
            let mut run = 0u32;
            for &(_, d) in group {
                if d == run_diag {
                    run += 1;
                } else {
                    if run > best_run {
                        best_run = run;
                        best_diag = run_diag;
                    }
                    run_diag = d;
                    run = 1;
                }
            }
            if run > best_run {
                best_diag = run_diag;
            }

            if total as usize >= min_hits {
                cands.push(Candidate {
                    target_id: tid,
                    hits: total,
                    best_diag: best_diag as i64,
                    diag_score: 0,
                });
            }
            i = j;
        }

        cands.sort_by(|a, b| b.hits.cmp(&a.hits));
        cands
    }
}

/// Best ungapped segment score (Kadane) along a diagonal, scoring a target
/// against a per-query-column scoring closure. `diag = target_pos - query_pos`.
#[inline]
fn diagonal_kadane<F>(query_len: usize, target: &[u8], diag: i64, score_at: F) -> i32
where
    F: Fn(usize, u8) -> i32,
{
    let mut best = 0i32;
    let mut cur = 0i32;
    for qpos in 0..query_len {
        let tpos = qpos as i64 + diag;
        if tpos < 0 {
            continue;
        }
        let tpos = tpos as usize;
        if tpos >= target.len() {
            break;
        }
        let s = score_at(qpos, target[tpos]);
        cur = (cur + s).max(0);
        if cur > best {
            best = cur;
        }
    }
    best
}

/// Apply the ungapped diagonal gate to a candidate list, filling `diag_score`
/// and keeping those at or above `min_diag_score`, ranked by `diag_score`.
fn gate_by_diagonal<F>(
    mut cands: Vec<Candidate>,
    query_len: usize,
    db: &SeqDb,
    min_diag_score: i32,
    score_at: F,
) -> Vec<Candidate>
where
    F: Fn(usize, u8) -> i32,
{
    for c in cands.iter_mut() {
        let target = &db.records[c.target_id as usize].num;
        c.diag_score = diagonal_kadane(query_len, target, c.best_diag, &score_at);
    }
    cands.retain(|c| c.diag_score >= min_diag_score);
    cands.sort_by(|a, b| b.diag_score.cmp(&a.diag_score));
    cands
}

/// Prefilter a plain sequence query against the index.
pub fn prefilter_sequence(
    query: &[u8],
    db: &SeqDb,
    index: &KmerIndex,
    mat: &SubstitutionMatrix,
    params: PrefilterParams,
) -> Vec<Candidate> {
    let pattern = index.pattern();
    let w = pattern.weight();
    let span = pattern.span();
    debug_assert_eq!(w, params.k, "index weight and params k must match");
    let mut acc = HitAccumulator::new();
    if query.len() < span {
        return Vec::new();
    }
    let mut variants: Vec<u64> = Vec::new();
    let mut col_scores = vec![[0i32; PROFILE_AA_SIZE]; w];
    for start in 0..=(query.len() - span) {
        // build per-care-position substitution scores for this window
        let mut skip = false;
        for (j, &off) in pattern.offsets().iter().enumerate() {
            let qr = query[start + off];
            if qr as usize >= PROFILE_AA_SIZE {
                skip = true;
                break;
            }
            for r in 0..PROFILE_AA_SIZE {
                col_scores[j][r] = mat.score(qr, r as u8);
            }
        }
        if skip {
            continue;
        }
        generate_similar_kmers(&col_scores, params.kmer_score, params.max_kmers_per_pos, &mut variants);
        for &code in &variants {
            acc.add(index, code, start);
        }
    }
    let cands = acc.into_candidates(params.min_hits);
    gate_by_diagonal(cands, query.len(), db, params.min_diag_score, |qpos, tr| {
        mat.score(query[qpos], tr)
    })
}

/// Prefilter a profile (PSSM) query against the index.
pub fn prefilter_profile(
    profile: &Profile,
    db: &SeqDb,
    index: &KmerIndex,
    params: PrefilterParams,
) -> Vec<Candidate> {
    let pattern = index.pattern();
    let w = pattern.weight();
    let span = pattern.span();
    debug_assert_eq!(w, params.k, "index weight and params k must match");
    let mut acc = HitAccumulator::new();
    if profile.query_len < span {
        return Vec::new();
    }
    let mut variants: Vec<u64> = Vec::new();
    let mut col_scores = vec![[0i32; PROFILE_AA_SIZE]; w];
    for start in 0..=(profile.query_len - span) {
        for (j, &off) in pattern.offsets().iter().enumerate() {
            for r in 0..PROFILE_AA_SIZE {
                col_scores[j][r] = profile.pssm[(start + off) * PROFILE_AA_SIZE + r] as i32;
            }
        }
        generate_similar_kmers(&col_scores, params.kmer_score, params.max_kmers_per_pos, &mut variants);
        for &code in &variants {
            acc.add(index, code, start);
        }
    }
    let cands = acc.into_candidates(params.min_hits);
    gate_by_diagonal(cands, profile.query_len, db, params.min_diag_score, |qpos, tr| {
        profile.score(qpos, tr)
    })
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
    fn exact_kmer_is_found() {
        let m = SubstitutionMatrix::blosum62();
        let db = db_from(">t1\nGGGGGMKVLLACDEFGHIKGGGG\n>t2\nPPPPPPPPPPPPPP\n");
        let idx = KmerIndex::build(&db, 6);
        let t = build_aa2num();
        let q: Vec<u8> = "MKVLLACDEFGHIK".bytes().map(|b| t[b as usize]).collect();
        let cands = prefilter_sequence(&q, &db, &idx, &m, PrefilterParams::default());
        assert!(!cands.is_empty());
        // t1 should be the top candidate
        assert_eq!(cands[0].target_id, 0);
    }

    #[test]
    fn similar_kmer_reaches_diverged_target() {
        let m = SubstitutionMatrix::blosum62();
        // target is a conservatively substituted version of the query
        let db = db_from(">t\nGGGGMKILLSCDEFGHLKGGGG\n");
        let idx = KmerIndex::build(&db, 6);
        let t = build_aa2num();
        let q: Vec<u8> = "MKVLLACDEFGHIK".bytes().map(|b| t[b as usize]).collect();
        // lower threshold so conservative substitutions seed
        let mut p = PrefilterParams::default();
        p.kmer_score = 15;
        let cands = prefilter_sequence(&q, &db, &idx, &m, p);
        assert!(!cands.is_empty(), "similar-kmer seeding failed to find diverged target");
    }

    #[test]
    fn profile_prefilter_finds_target() {
        let m = SubstitutionMatrix::blosum62();
        let msa = ">q\nMKVLLACDEFGHIK\n>h\nMKILLSCDEFGHLK\n";
        let prof = Profile::from_msa_str(msa, &m, PseudoCountParams::default()).unwrap();
        let db = db_from(">t\nGGGGMKVLLACDEFGHIKGGGG\n");
        let idx = KmerIndex::build(&db, 6);
        let mut p = PrefilterParams::default();
        p.kmer_score = 40; // PSSM scores are on the ×8 scale
        let cands = prefilter_profile(&prof, &db, &idx, p);
        assert!(!cands.is_empty());
    }

    #[test]
    fn k_out_of_range_is_rejected() {
        assert!(check_k(1).is_err());
        assert!(check_k(MAX_K + 1).is_err());
        assert!(check_k(6).is_ok());
        assert!(check_k(MAX_K).is_ok());
    }

    #[test]
    fn long_kmers_do_not_collide_under_u64() {
        // With the old u32 code, k=8 codes overflowed (20^8 > 2^32) and distinct
        // k-mers collided. Under u64 they must remain distinct. Build an index of
        // two different 8-mers and confirm two separate codes exist.
        let db = db_from(">a\nACDEFGHIKLMNPQRSTVWY\n"); // 20 residues -> many 8-mers
        let idx = KmerIndex::build(&db, 8);
        // number of distinct 8-mers in a 20-length window = 20-8+1 = 13, all unique
        assert_eq!(idx.distinct_kmers(), 13);
    }

    #[test]
    fn diagonal_mode_is_selected() {
        // Two seeds on one diagonal, one on another → best_diag is the common one.
        let m = SubstitutionMatrix::blosum62();
        let db = db_from(">t\nMKVLLACDEFGHIKMKVLLAC\n"); // repeated motif
        let idx = KmerIndex::build(&db, 6);
        let t = build_aa2num();
        let q: Vec<u8> = "MKVLLACDEFGHIK".bytes().map(|b| t[b as usize]).collect();
        let cands = prefilter_sequence(&q, &db, &idx, &m, PrefilterParams::default());
        assert!(!cands.is_empty());
        assert!(cands[0].diag_score > 0);
    }

    #[test]
    fn seed_pattern_parsing() {
        let p = SeedPattern::from_mask("1101011").unwrap();
        assert_eq!(p.weight(), 5);
        assert_eq!(p.span(), 7);
        assert_eq!(p.mask_string(), "1101011");
        assert!(SeedPattern::from_mask("0110").is_err()); // must start with 1
        assert!(SeedPattern::from_mask("1100").is_err()); // must end with 1
        assert!(SeedPattern::from_mask("11x1").is_err()); // only 0/1
        assert!(SeedPattern::from_mask("").is_err());
        let c = SeedPattern::contiguous(6);
        assert_eq!(c.weight(), 6);
        assert_eq!(c.span(), 6);
    }

    #[test]
    fn spaced_seed_finds_what_contiguous_misses() {
        // Target is exactly the diverged seed window: the query's 7-residue window
        // with a single substitution at the position the spaced mask ignores. The
        // spaced seed matches it exactly; no contiguous weight-6 window matches
        // (the mismatch falls inside every contiguous 6-mer), so at a threshold
        // that admits the exact spaced self-match but not the substituted variant,
        // contiguous seeding misses while spaced succeeds.
        let m = SubstitutionMatrix::blosum62();
        let t = build_aa2num();
        let qs = "MKVLLACDEFGHIK"; // window at start 3 = indices 3..9 = "LLACDEF"
        let q: Vec<u8> = qs.bytes().map(|b| t[b as usize]).collect();
        // "LLACDEF" with index-5-of-window (E) -> P: the don't-care column
        let db = db_from(">tgt\nLLACDPF\n");

        let mut strict = PrefilterParams::default();
        strict.kmer_score = 30; // admits exact self-match (~33), not the E->P variant (~26)
        strict.min_diag_score = 1;

        // contiguous weight-6
        strict.k = 6;
        let idx_c = KmerIndex::build(&db, 6);
        let c_cont = prefilter_sequence(&q, &db, &idx_c, &m, strict);
        assert!(
            !c_cont.iter().any(|c| c.diag_score > 0),
            "contiguous seed should miss the don't-care substitution here"
        );

        // spaced weight-6 with the don't-care ('0') at window offset 5
        let pattern = SeedPattern::from_mask("1111101").unwrap();
        assert_eq!(pattern.weight(), 6);
        strict.k = pattern.weight();
        let idx_s = KmerIndex::build_spaced(&db, pattern);
        let c_spaced = prefilter_sequence(&q, &db, &idx_s, &m, strict);
        assert!(
            c_spaced.iter().any(|c| c.diag_score > 0),
            "spaced seed should find the diverged target"
        );
    }

    #[test]
    fn sensitivity_lowers_threshold_monotonically() {
        // Higher sensitivity -> lower (or equal) k-mer score threshold.
        let mut prev = i32::MAX;
        for i in 1..=9 {
            let s = i as f64;
            let thr = kmer_score_for_sensitivity(s, 2);
            assert!(thr <= prev, "threshold must be non-increasing in s");
            prev = thr;
        }
        // profile scale (×8) thresholds are larger than sequence scale (×2)
        assert!(kmer_score_for_sensitivity(5.7, 8) > kmer_score_for_sensitivity(5.7, 2));
    }
}
