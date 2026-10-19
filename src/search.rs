//! Search orchestration: staged, bounded and deterministic.
//!
//! Per query:
//!
//! 1. **Seed evidence** — [`crate::prefilter`] (or every target in exhaustive
//!    mode) yields ranked candidates.
//! 2. **Gapped score** — the striped SIMD kernels (u8 first pass, i16 when its
//!    headroom is exceeded) score each candidate (whole target, or a
//!    seed-guided window of a long target) and compare it with the query's
//!    minimum score for `--evalue`, derived once from the calibrated statistics
//!    (no per-candidate E-value computation).
//! 3. **Deferred traceback** — survivors are sorted by score and aligned with
//!    [`crate::align::align_canonical`] in that order only until `max_hits` hits
//!    are *accepted*; when a traceback-dependent filter (coverage) rejects one,
//!    the next survivor is aligned. The traceback is exact (three-state), so the
//!    reported raw score, bit score, E-value and coordinates all describe the
//!    same alignment; a disagreement with the SIMD score is counted
//!    (`score_mismatches`, expected 0).
//!
//! Low-complexity query positions can be hard-masked in scoring
//! (`mask_scoring`, default on for sequence queries): their score rows become
//! the matrix `X` row before steps 2–3.
//!
//! Results per query are bounded by `max_hits`; ordering is by query ordinal,
//! then E-value, bit score and target ordinal, so output is deterministic and
//! duplicate query identifiers stay separate records.

use crate::align::{
    align_canonical, local_align, AlignScratch, Alignment, GapCosts, ProfileScorer, DEFAULT_BAND_RADIUS,
    ScoreTable, SeqScorer,
};
use crate::alphabet::{ALPHABET_SIZE, GAP};
use crate::evalue::{calibrate_profile, EValueParams};
use crate::fasta::SeqDb;
use crate::index::SeedIndex;
use crate::mask::low_complexity_mask_into;
use crate::matrix::SubstitutionMatrix;
use crate::prefilter::{prefilter_query, Candidate, PrefilterParams, PrefilterQuery, PrefilterScratch, PrefilterStats};
use crate::profile::{Profile, PseudoCountParams};
use crate::simd::{Kernel, QueryProfile, SimdScratch};
use rayon::prelude::*;
use std::collections::HashSet;
use std::time::Instant;

/// Seed-guided target windows for long targets (a **heuristic**: an alignment
/// entirely outside the window cannot be found; one that touches the window
/// edge triggers re-alignment on a doubled window).
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "distributed", derive(serde::Serialize, serde::Deserialize))]
pub struct WindowParams {
    pub enabled: bool,
    /// Only targets at least this long are windowed.
    pub min_target_len: usize,
    /// Margin added on each side beyond `query_len`.
    pub margin: usize,
}

impl Default for WindowParams {
    fn default() -> Self {
        WindowParams {
            enabled: false,
            min_target_len: 3000,
            margin: 200,
        }
    }
}

/// Per-worker scratch limits (the CLI derives them from `--memory-limit`).
/// They bound memory only; results do not depend on them. (The traceback
/// method switch is a fixed problem-size constant,
/// [`crate::align::FULL_TRACEBACK_MAX_CELLS`], for the same reason.)
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "distributed", derive(serde::Serialize, serde::Deserialize))]
pub struct ScratchLimits {
    /// Posting ranges held by the per-worker variant cache.
    pub variant_cache_ranges: usize,
}

impl Default for ScratchLimits {
    fn default() -> Self {
        ScratchLimits {
            variant_cache_ranges: 1 << 19,
        }
    }
}

/// Complete, validated search configuration (shared by CLI, library and the
/// distributed backend).
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "distributed", derive(serde::Serialize, serde::Deserialize))]
pub struct SearchParams {
    pub prefilter: PrefilterParams,
    /// Report hits with E-value ≤ this.
    pub max_evalue: f64,
    /// Minimum fraction of the query covered by the alignment (0 = off).
    pub min_query_cov: f64,
    /// Minimum fraction of the target covered by the alignment (0 = off).
    pub min_target_cov: f64,
    /// Maximum hits reported per query (0 = unlimited). An output limit only.
    pub max_hits: usize,
    /// Initial traceback band radius (exact: widened until the optimum is
    /// reproduced). `None` = automatic (full traceback for small problems).
    pub band_radius: Option<usize>,
    /// Gap costs (open includes the first position).
    pub gaps: GapCosts,
    /// MMseqs2-style local composition-bias correction (sequence queries).
    pub comp_bias_corr: bool,
    /// Clip positive composition corrections to 0 (experimental).
    pub comp_bias_penalize_only: bool,
    /// Apply low-complexity masking to query seed windows (the index side is
    /// fixed at index build time).
    pub mask_query: bool,
    /// Also hard-mask low-complexity query positions in **scoring** (their
    /// score row becomes the matrix `X` row), so similarity confined to
    /// low-complexity regions cannot produce a hit (sequence queries).
    pub mask_scoring: Option<crate::mask::MaskParams>,
    /// Seed-guided windows for long targets.
    pub window: WindowParams,
    /// Skip the prefilter: score every target (easysearch).
    pub exhaustive: bool,
    /// Random sequences per profile for per-PSSM calibration.
    pub profile_calibration_samples: usize,
    /// Seed for per-profile calibration.
    pub calibration_seed: u64,
    /// Per-worker memory limits.
    pub scratch: ScratchLimits,
}

impl SearchParams {
    pub fn sequence_default() -> Self {
        SearchParams {
            prefilter: PrefilterParams::sequence_default(),
            max_evalue: 1e-3,
            min_query_cov: 0.0,
            min_target_cov: 0.0,
            max_hits: 300,
            band_radius: None,
            gaps: GapCosts::sequence_default(),
            comp_bias_corr: true,
            comp_bias_penalize_only: false,
            mask_query: true,
            mask_scoring: Some(crate::mask::MaskParams::scoring_default()),
            window: WindowParams::default(),
            exhaustive: false,
            profile_calibration_samples: 2000,
            calibration_seed: 0x5EED,
            scratch: ScratchLimits::default(),
        }
    }

    pub fn profile_default() -> Self {
        SearchParams {
            prefilter: PrefilterParams::profile_default(),
            gaps: GapCosts::profile_default(),
            comp_bias_corr: false,
            mask_scoring: None,
            ..Self::sequence_default()
        }
    }

    pub fn validate(&self) -> Result<(), String> {
        self.prefilter.validate()?;
        if self.max_evalue.is_nan() || self.max_evalue <= 0.0 || !self.max_evalue.is_finite() {
            return Err("--evalue must be a positive number".into());
        }
        for (v, name) in [(self.min_query_cov, "--min-cov/--min-qcov"), (self.min_target_cov, "--min-tcov")] {
            if !(0.0..=1.0).contains(&v) {
                return Err(format!("{} must be within 0..1", name));
            }
        }
        if self.gaps.open < 0 || self.gaps.extend < 0 {
            return Err("gap costs must be non-negative".into());
        }
        Ok(())
    }
}

impl Default for SearchParams {
    fn default() -> Self {
        Self::sequence_default()
    }
}

/// One reported hit.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "distributed", derive(serde::Serialize, serde::Deserialize))]
pub struct Hit {
    /// Ordinal of the query in its input (stable; duplicates stay distinct).
    pub query_index: usize,
    /// Ordinal of the target in the database.
    pub target_index: usize,
    pub query_id: String,
    pub target_id: String,
    pub raw_score: i32,
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
    pub qlen: usize,
    pub tlen: usize,
    pub cigar: String,
}

/// The order of one query's hits: E-value ascending, bit score descending,
/// target ordinal ascending (a total order: the target ordinal is unique per
/// query). Shared by the local engine and the distributed merge.
pub fn hit_order(a: &Hit, b: &Hit) -> std::cmp::Ordering {
    a.evalue
        .partial_cmp(&b.evalue)
        .unwrap_or(std::cmp::Ordering::Equal)
        .then(b.bit_score.partial_cmp(&a.bit_score).unwrap_or(std::cmp::Ordering::Equal))
        .then(a.target_index.cmp(&b.target_index))
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

    /// Extra columns: raw score, query/target length, query/target coverage, CIGAR.
    pub fn extra_cols(&self) -> String {
        format!(
            "{}\t{}\t{}\t{:.3}\t{:.3}\t{}",
            self.raw_score,
            self.qlen,
            self.tlen,
            (self.q_end + 1 - self.q_start) as f64 / self.qlen.max(1) as f64,
            (self.t_end + 1 - self.t_start) as f64 / self.tlen.max(1) as f64,
            self.cigar
        )
    }
}

/// Column header for `m8` output.
pub const M8_HEADER: &str = "query\ttarget\tpident\talnlen\tmismatch\tgapopen\tqstart\tqend\ttstart\ttend\tevalue\tbits";
/// Header of [`Hit::extra_cols`].
pub const EXTRA_HEADER: &str = "raw\tqlen\ttlen\tqcov\ttcov\tcigar";

/// Read access to target sequences (in-memory database or mapped index).
pub trait Targets: Sync {
    fn count(&self) -> usize;
    fn residues_of(&self, i: usize) -> &[u8];
    fn ident(&self, i: usize) -> &str;
    fn residue_total(&self) -> usize;
}

impl Targets for SeqDb {
    fn count(&self) -> usize {
        self.len()
    }
    fn residues_of(&self, i: usize) -> &[u8] {
        self.num(i)
    }
    fn ident(&self, i: usize) -> &str {
        self.id(i)
    }
    fn residue_total(&self) -> usize {
        self.total_residues
    }
}

impl Targets for SeedIndex {
    fn count(&self) -> usize {
        self.len()
    }
    fn residues_of(&self, i: usize) -> &[u8] {
        self.seq(i)
    }
    fn ident(&self, i: usize) -> &str {
        self.id(i)
    }
    fn residue_total(&self) -> usize {
        self.total_residues()
    }
}

/// Counters and stage timings (summable across queries and workers).
#[derive(Clone, Copy, Debug, Default)]
#[cfg_attr(feature = "distributed", derive(serde::Serialize, serde::Deserialize))]
pub struct SearchStats {
    pub prefilter: PrefilterStats,
    pub queries: u64,
    pub query_residues: u64,
    pub candidates: u64,
    pub gapped_scored: u64,
    pub gapped_cells: u64,
    pub scalar_fallbacks: u64,
    pub windowed: u64,
    pub window_expansions: u64,
    pub passed_score: u64,
    pub tracebacks: u64,
    pub rejected_coverage: u64,
    pub hits: u64,
    pub score_mismatches: u64,
    pub profile_calibrations: u64,
    pub masked_scoring_positions: u64,
    pub ns_setup: u64,
    pub ns_prefilter: u64,
    pub ns_gapped: u64,
    pub ns_traceback: u64,
}

impl SearchStats {
    pub fn add(&mut self, o: &SearchStats) {
        self.prefilter.add(&o.prefilter);
        macro_rules! sum { ($($f:ident),*) => { $( self.$f += o.$f; )* } }
        sum!(
            queries, query_residues, candidates, gapped_scored, gapped_cells, scalar_fallbacks, windowed,
            window_expansions, passed_score, tracebacks, rejected_coverage, hits, score_mismatches,
            profile_calibrations, masked_scoring_positions, ns_setup, ns_prefilter, ns_gapped, ns_traceback
        );
    }
}

/// Per-worker reusable memory.
#[derive(Default)]
pub struct WorkerScratch {
    pub prefilter: PrefilterScratch,
    pub simd: SimdScratch,
    pub align: AlignScratch,
    mask: Vec<bool>,
}

impl WorkerScratch {
    pub fn new() -> Self {
        Self::default()
    }
    /// Scratch bounded by `limits`.
    pub fn with_limits(limits: &ScratchLimits) -> Self {
        let mut s = WorkerScratch::default();
        s.prefilter.cache_limit = limits.variant_cache_ranges;
        s
    }
    pub fn heap_bytes(&self) -> usize {
        self.prefilter.heap_bytes() + self.align.heap_bytes() + self.mask.capacity()
    }
}

/// Statistics source.
#[derive(Clone, Debug)]
pub enum Statistics {
    /// Fixed parameters (sequence search: the calibrated table).
    Fixed(EValueParams),
    /// Calibrate each profile by simulation; `db_residues`/`db_seqs` are the
    /// global database size.
    PerProfile { db_residues: usize, db_seqs: usize, length_sample: Vec<usize> },
}

/// A query to search.
pub enum QueryKind<'a> {
    Sequence(&'a [u8]),
    Profile(&'a Profile),
}

/// Everything shared by all queries of a search.
pub struct Engine<'a, T: Targets> {
    pub mat: &'a SubstitutionMatrix,
    pub targets: &'a T,
    /// Seed index for prefiltering (required unless `params.exhaustive`).
    pub index: Option<&'a SeedIndex>,
    pub params: SearchParams,
    pub stats: Statistics,
}

struct Pass {
    score: i32,
    tid: u32,
    w0: u32,
    w1: u32,
}

impl<'a, T: Targets> Engine<'a, T> {
    fn window_for(&self, c: &Candidate, m: usize, n: usize) -> (usize, usize) {
        let wp = self.params.window;
        if !wp.enabled || n < wp.min_target_len {
            return (0, n);
        }
        let ext = m + wp.margin;
        if n <= 2 * ext + (c.t_hi as usize - c.t_lo as usize) {
            return (0, n);
        }
        let w0 = (c.t_lo as usize).saturating_sub(ext);
        let w1 = (c.t_hi as usize + 1 + ext).min(n);
        (w0, w1)
    }

    /// Search one query. `parallel` uses rayon across candidates (for
    /// one/few-query workloads; callers must not also parallelise across
    /// queries).
    pub fn search_query(
        &self,
        qi: usize,
        qid: &str,
        query: QueryKind,
        scratch: &mut WorkerScratch,
        stats: &mut SearchStats,
        parallel: bool,
    ) -> Vec<(Hit, Alignment)> {
        let p = &self.params;
        let t_setup = Instant::now();
        // ---- score tables --------------------------------------------------
        // low-complexity mask of a sequence query for seeding
        let mut have_mask = false;
        if let (Some(mp), QueryKind::Sequence(q)) = (p.prefilter.mask, &query) {
            if p.mask_query {
                low_complexity_mask_into(q, mp, &mut scratch.mask);
                have_mask = true;
            }
        }
        let (seed_rows, align_table, residues): (Vec<[i32; ALPHABET_SIZE]>, ScoreTable, Option<&[u8]>) = match query {
            QueryKind::Sequence(q) => {
                let sc = SeqScorer { query: q, mat: self.mat };
                let plain = ScoreTable::from_scorer(&sc);
                let mut table = if p.comp_bias_corr {
                    let bias = crate::compbias::local_bias_correction_mode(q, self.mat, p.comp_bias_penalize_only);
                    ScoreTable::with_bias(&sc, Some(&bias))
                } else {
                    plain.clone()
                };
                if let Some(msp) = p.mask_scoring {
                    let x = crate::alphabet::ALPHABET_SIZE as u8 - 1;
                    let mut xrow = [0i32; ALPHABET_SIZE];
                    for (r, v) in xrow.iter_mut().enumerate() {
                        *v = self.mat.score(x, r as u8);
                    }
                    let masked: Vec<usize> = if have_mask && Some(msp) == p.prefilter.mask {
                        (0..q.len()).filter(|&i| scratch.mask[i]).collect()
                    } else {
                        let mk = crate::mask::low_complexity_mask(q, msp);
                        (0..q.len()).filter(|&i| mk[i]).collect()
                    };
                    stats.masked_scoring_positions += masked.len() as u64;
                    table.set_rows(&masked, &xrow);
                }
                (plain.rows().to_vec(), table, Some(q))
            }
            QueryKind::Profile(pr) => {
                let t = ScoreTable::from_scorer(&ProfileScorer { profile: pr });
                (t.rows().to_vec(), t, None)
            }
        };
        let m = align_table.len();
        stats.queries += 1;
        stats.query_residues += m as u64;
        let ev = match &self.stats {
            Statistics::Fixed(e) => e.clone(),
            Statistics::PerProfile { db_residues, db_seqs, length_sample } => {
                stats.profile_calibrations += 1;
                match calibrate_profile(
                    &align_table,
                    self.mat,
                    p.gaps,
                    length_sample,
                    p.profile_calibration_samples.max(100),
                    p.calibration_seed,
                ) {
                    Some((gp, _)) => EValueParams::from_params(gp, *db_residues, *db_seqs, "per-profile calibrated"),
                    None => EValueParams::for_profile(*db_residues),
                }
            }
        };
        let space = ev.search_space(m);
        let min_score = ev.min_score_for(p.max_evalue, space);
        let qp = QueryProfile::from_rows(align_table.rows().to_vec());
        stats.ns_setup += t_setup.elapsed().as_nanos() as u64;

        // ---- candidates ------------------------------------------------------
        let t_pf = Instant::now();
        let cands: Vec<Candidate> = match self.index {
            Some(idx) if !p.exhaustive => {
                let mask_ref: Option<&[bool]> = if have_mask && p.mask_query && residues.is_some() {
                    Some(&scratch.mask)
                } else {
                    None
                };
                let pq = PrefilterQuery {
                    seed_rows: &seed_rows,
                    residues,
                    ext_rows: align_table.rows(),
                    seed_mask: mask_ref,
                };
                prefilter_query(&pq, idx, &p.prefilter, p.gaps, &mut scratch.prefilter, &mut stats.prefilter)
            }
            // exhaustive (or no index): every target is a candidate
            _ => (0..self.targets.count())
                .map(|t| Candidate {
                    target_id: t as u32,
                    hits: 0,
                    best_diag: 0,
                    diag_score: 0,
                    t_lo: 0,
                    t_hi: self.targets.residues_of(t).len().saturating_sub(1) as u32,
                })
                .collect(),
        };
        stats.candidates += cands.len() as u64;
        stats.ns_prefilter += t_pf.elapsed().as_nanos() as u64;

        // ---- gapped scores -------------------------------------------------
        let t_g = Instant::now();
        let score_one = |c: &Candidate, simd: &mut SimdScratch| -> (Option<Pass>, u64, bool, bool) {
            let t = self.targets.residues_of(c.target_id as usize);
            let (w0, w1) = self.window_for(c, m, t.len());
            let (s, k) = qp.sw_score_with(&t[w0..w1], p.gaps, simd);
            let cells = (m * (w1 - w0)) as u64;
            let pass = if s >= min_score && s > 0 {
                Some(Pass {
                    score: s,
                    tid: c.target_id,
                    w0: w0 as u32,
                    w1: w1 as u32,
                })
            } else {
                None
            };
            (pass, cells, k == Kernel::Scalar, (w0, w1) != (0, t.len()))
        };
        let mut passes: Vec<Pass> = Vec::new();
        if parallel && cands.len() > 64 {
            let res: Vec<(Option<Pass>, u64, bool, bool)> = cands
                .par_iter()
                .map_init(SimdScratch::new, |s, c| score_one(c, s))
                .collect();
            for (pass, cells, scalar, win) in res {
                stats.gapped_scored += 1;
                stats.gapped_cells += cells;
                stats.scalar_fallbacks += scalar as u64;
                stats.windowed += win as u64;
                if let Some(x) = pass {
                    passes.push(x);
                }
            }
        } else {
            for c in &cands {
                let (pass, cells, scalar, win) = score_one(c, &mut scratch.simd);
                stats.gapped_scored += 1;
                stats.gapped_cells += cells;
                stats.scalar_fallbacks += scalar as u64;
                stats.windowed += win as u64;
                if let Some(x) = pass {
                    passes.push(x);
                }
            }
        }
        stats.passed_score += passes.len() as u64;
        passes.sort_unstable_by(|a, b| b.score.cmp(&a.score).then(a.tid.cmp(&b.tid)));
        stats.ns_gapped += t_g.elapsed().as_nanos() as u64;

        // ---- deferred traceback --------------------------------------------
        let t_tb = Instant::now();
        let mut out: Vec<(Hit, Alignment)> = Vec::new();
        for ps in &passes {
            if p.max_hits > 0 && out.len() >= p.max_hits {
                break;
            }
            let tid = ps.tid as usize;
            let t = self.targets.residues_of(tid);
            let n = t.len();
            let (mut w0, mut w1) = (ps.w0 as usize, ps.w1 as usize);
            let mut expected = ps.score;
            let aln = loop {
                stats.tracebacks += 1;
                let tw = &t[w0..w1];
                let a = align_canonical(
                    &qp,
                    &align_table,
                    tw,
                    p.gaps,
                    p.band_radius.unwrap_or(DEFAULT_BAND_RADIUS),
                    &mut scratch.align,
                );
                let mut a = match a {
                    Some(a) => a,
                    None => break None,
                };
                a.shift_target(w0);
                let touches = (w0 > 0 && a.target_start <= w0) || (w1 < n && a.target_end + 1 >= w1);
                if touches {
                    // the window clipped the alignment: widen and redo
                    stats.window_expansions += 1;
                    let ext = (w1 - w0).max(m);
                    w0 = w0.saturating_sub(ext);
                    w1 = (w1 + ext).min(n);
                    expected = qp.sw_score_with(&t[w0..w1], p.gaps, &mut scratch.simd).0;
                    continue;
                }
                break Some(a);
            };
            let aln = match aln {
                Some(a) => a,
                None => continue,
            };
            if aln.raw_score != expected {
                stats.score_mismatches += 1;
            }
            if p.min_query_cov > 0.0 && aln.query_cov(m) < p.min_query_cov {
                stats.rejected_coverage += 1;
                continue;
            }
            if p.min_target_cov > 0.0 && aln.target_cov(n) < p.min_target_cov {
                stats.rejected_coverage += 1;
                continue;
            }
            let evalue = ev.evalue_in_space(aln.raw_score, space);
            if evalue > p.max_evalue {
                continue;
            }
            let hit = Hit {
                query_index: qi,
                target_index: tid,
                query_id: qid.to_string(),
                target_id: self.targets.ident(tid).to_string(),
                raw_score: aln.raw_score,
                pident: aln.pct_identity(),
                aln_len: aln.aln_len,
                mismatches: aln.mismatches,
                gap_opens: aln.gap_opens,
                q_start: aln.query_start,
                q_end: aln.query_end,
                t_start: aln.target_start,
                t_end: aln.target_end,
                evalue,
                bit_score: ev.bit_score(aln.raw_score),
                qlen: m,
                tlen: n,
                cigar: aln.cigar.clone(),
            };
            out.push((hit, aln));
        }
        // final order: E-value, bits, target ordinal
        out.sort_by(|a, b| hit_order(&a.0, &b.0));
        stats.hits += out.len() as u64;
        stats.ns_traceback += t_tb.elapsed().as_nanos() as u64;
        out
    }

    /// Search a batch of sequence queries on the calling thread only (no
    /// rayon): used where an outer scheduler (the distributed backend) already
    /// runs one task per CPU slot, so nothing is nested.
    pub fn search_batch_sequential(&self, queries: &SeqDb, first_ordinal: usize, stats: &mut SearchStats) -> Vec<Vec<(Hit, Alignment)>> {
        let mut scratch = WorkerScratch::with_limits(&self.params.scratch);
        (0..queries.len())
            .map(|i| {
                self.search_query(
                    first_ordinal + i,
                    queries.id(i),
                    QueryKind::Sequence(queries.num(i)),
                    &mut scratch,
                    stats,
                    false,
                )
            })
            .collect()
    }

    /// Search a batch of sequence queries (`first_ordinal` = ordinal of
    /// `queries[0]`). Parallel across queries, or across candidates when the
    /// batch has fewer queries than threads (never both).
    pub fn search_batch(&self, queries: &SeqDb, first_ordinal: usize, stats: &mut SearchStats) -> Vec<Vec<(Hit, Alignment)>> {
        let threads = rayon::current_num_threads().max(1);
        if queries.len() < threads {
            let mut scratch = WorkerScratch::with_limits(&self.params.scratch);
            (0..queries.len())
                .map(|i| {
                    self.search_query(
                        first_ordinal + i,
                        queries.id(i),
                        QueryKind::Sequence(queries.num(i)),
                        &mut scratch,
                        stats,
                        true,
                    )
                })
                .collect()
        } else {
            let per: Vec<(Vec<(Hit, Alignment)>, SearchStats)> = (0..queries.len())
                .into_par_iter()
                .map_init(|| WorkerScratch::with_limits(&self.params.scratch), |scratch, i| {
                    let mut st = SearchStats::default();
                    let r = self.search_query(
                        first_ordinal + i,
                        queries.id(i),
                        QueryKind::Sequence(queries.num(i)),
                        scratch,
                        &mut st,
                        false,
                    );
                    (r, st)
                })
                .collect();
            per.into_iter()
                .map(|(r, st)| {
                    stats.add(&st);
                    r
                })
                .collect()
        }
    }
}

fn check_index_matches(target_db: &SeqDb, index: &SeedIndex) {
    assert!(
        index.len() == target_db.len() && index.total_residues() == target_db.total_residues,
        "seed index ({} seqs, {} residues) was not built from this target database ({} seqs, {} residues)",
        index.len(),
        index.total_residues(),
        target_db.len(),
        target_db.total_residues
    );
}

fn synced(mut params: SearchParams, index: &SeedIndex) -> SearchParams {
    params.prefilter.k = index.k();
    params
}

/// Search every sequence in `query_db` against `target_db` (library API).
/// `index` must have been built from `target_db`; `ev` supplies the statistics.
pub fn search_sequences(
    query_db: &SeqDb,
    target_db: &SeqDb,
    index: &SeedIndex,
    mat: &SubstitutionMatrix,
    ev: &EValueParams,
    params: SearchParams,
) -> Vec<Hit> {
    check_index_matches(target_db, index);
    let eng = Engine {
        mat,
        targets: index,
        index: Some(index),
        params: synced(params, index),
        stats: Statistics::Fixed(ev.clone()),
    };
    let mut st = SearchStats::default();
    eng.search_batch(query_db, 0, &mut st)
        .into_iter()
        .flatten()
        .map(|(h, _)| h)
        .collect()
}

/// **Database-free** exhaustive search: every query against every target with
/// the striped SIMD kernel (no prefilter, no index).
pub fn search_sequences_exhaustive(
    query_db: &SeqDb,
    target_db: &SeqDb,
    mat: &SubstitutionMatrix,
    ev: &EValueParams,
    mut params: SearchParams,
) -> Vec<Hit> {
    params.exhaustive = true;
    let eng = Engine {
        mat,
        targets: target_db,
        index: None,
        params,
        stats: Statistics::Fixed(ev.clone()),
    };
    let mut st = SearchStats::default();
    eng.search_batch(query_db, 0, &mut st)
        .into_iter()
        .flatten()
        .map(|(h, _)| h)
        .collect()
}

/// Profile search directly on a seed index (the index holds the targets).
/// `db_size` = (residues, sequences) of the *global* database (equal to the
/// index unless the search is sharded). With
/// `params.profile_calibration_samples > 0` the PSSM is calibrated by
/// simulation; otherwise `fixed` statistics are used.
pub fn profile_hits(
    profile: &Profile,
    query_id: &str,
    index: &SeedIndex,
    mat: &SubstitutionMatrix,
    params: SearchParams,
    db_size: (usize, usize),
    fixed: Option<&EValueParams>,
    stats_out: Option<&mut SearchStats>,
) -> Vec<Hit> {
    let stats = match (params.profile_calibration_samples > 0, fixed) {
        (false, Some(ev)) => Statistics::Fixed(ev.clone()),
        _ => Statistics::PerProfile {
            db_residues: db_size.0,
            db_seqs: db_size.1,
            length_sample: index.length_sample(2000),
        },
    };
    let eng = Engine {
        mat,
        targets: index,
        index: Some(index),
        params: synced(params, index),
        stats,
    };
    let mut local = SearchStats::default();
    let st = match stats_out {
        Some(s) => s,
        None => &mut local,
    };
    let mut scratch = WorkerScratch::with_limits(&eng.params.scratch);
    eng.search_query(0, query_id, QueryKind::Profile(profile), &mut scratch, st, true)
        .into_iter()
        .map(|(h, _)| h)
        .collect()
}

/// Search a single profile against a target database (library API). The
/// profile is calibrated by simulation when `params.profile_calibration_samples
/// > 0`; `ev` then only supplies the database size.
pub fn search_profile(
    profile: &Profile,
    query_id: &str,
    db: &SeqDb,
    index: &SeedIndex,
    mat: &SubstitutionMatrix,
    ev: &EValueParams,
    params: SearchParams,
) -> Vec<Hit> {
    check_index_matches(db, index);
    profile_hits(
        profile,
        query_id,
        index,
        mat,
        params,
        (ev.db_residues as usize, (ev.db_seqs as usize).max(db.len())),
        Some(ev),
        None,
    )
}

/// Build a query-anchored MSA row (length `query_len`) from a target's optimal
/// alignment to the current profile (a3m match-column convention).
fn hit_to_msa_row(profile: &Profile, query_len: usize, target: &[u8], gaps: GapCosts) -> Option<Vec<u8>> {
    let sc = ProfileScorer { profile };
    let aln = local_align(&sc, target, gaps)?;
    let mut row = vec![GAP; query_len];
    let mut qpos = aln.query_start;
    let mut tpos = aln.target_start;
    for (op, count) in crate::align::cigar_runs(&aln.cigar).ok()? {
        match op {
            b'M' => {
                for _ in 0..count {
                    if qpos < query_len && tpos < target.len() {
                        row[qpos] = target[tpos];
                    }
                    qpos += 1;
                    tpos += 1;
                }
            }
            b'D' => qpos += count,
            _ => tpos += count,
        }
    }
    Some(row)
}

/// Iterative profile search (PSI-BLAST style) on a seed index: search, fold
/// significant hits into the MSA, rebuild the profile, repeat until
/// `num_iterations` rounds or no new sequence is added. Returns the final
/// round's hits.
#[allow(clippy::too_many_arguments)]
pub fn iterate_profile(
    initial_msa: &[Vec<u8>],
    query_len: usize,
    query_id: &str,
    index: &SeedIndex,
    mat: &SubstitutionMatrix,
    pc: PseudoCountParams,
    params: SearchParams,
    db_size: (usize, usize),
    fixed: Option<&EValueParams>,
    num_iterations: usize,
    inclusion_evalue: f64,
    mut stats_out: Option<&mut SearchStats>,
) -> Result<Vec<Hit>, String> {
    let gaps = params.gaps;
    let mut msa: Vec<Vec<u8>> = initial_msa.to_vec();
    let mut seen: HashSet<usize> = HashSet::new();
    let mut last_hits: Vec<Hit> = Vec::new();
    let rounds = num_iterations.max(1);
    for iter in 0..rounds {
        let profile = Profile::from_numeric_msa(&msa, query_len, mat, pc)?;
        last_hits = profile_hits(&profile, query_id, index, mat, params, db_size, fixed, stats_out.as_deref_mut());
        if iter + 1 >= rounds {
            break;
        }
        let mut added = 0usize;
        for h in &last_hits {
            if h.evalue <= inclusion_evalue && seen.insert(h.target_index) {
                if let Some(row) = hit_to_msa_row(&profile, query_len, index.seq(h.target_index), gaps) {
                    msa.push(row);
                    added += 1;
                }
            }
        }
        if added == 0 {
            break;
        }
    }
    Ok(last_hits)
}

/// [`iterate_profile`] over an in-memory database (library API).
#[allow(clippy::too_many_arguments)]
pub fn search_profile_iterative(
    initial_msa: &[Vec<u8>],
    query_len: usize,
    query_id: &str,
    db: &SeqDb,
    index: &SeedIndex,
    mat: &SubstitutionMatrix,
    pc: PseudoCountParams,
    ev: &EValueParams,
    params: SearchParams,
    num_iterations: usize,
    inclusion_evalue: f64,
) -> Result<Vec<Hit>, String> {
    check_index_matches(db, index);
    iterate_profile(
        initial_msa,
        query_len,
        query_id,
        index,
        mat,
        pc,
        params,
        (ev.db_residues as usize, (ev.db_seqs as usize).max(db.len())),
        Some(ev),
        num_iterations,
        inclusion_evalue,
        None,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::build_aa2num;
    use crate::compbias::local_bias_correction;
    use crate::prefilter::KmerIndex;
    use crate::profile::PseudoCountParams;
    use std::io::Cursor;

    fn db_from(fasta: &str) -> SeqDb {
        let t = build_aa2num();
        SeqDb::from_reader(Cursor::new(fasta), &t).unwrap()
    }

    fn ev_for(db: &SeqDb) -> EValueParams {
        EValueParams::sequence(GapCosts::sequence_default(), true, db.total_residues, db.len()).unwrap()
    }

    #[test]
    fn exhaustive_finds_superset_of_prefiltered() {
        let m = SubstitutionMatrix::blosum62();
        let target = db_from(
            ">t1\nMKVLLACDEFGHIKLMNPQRSTVWY\n>t2\nMKILLSCDEFGHLKLMNPQRSTVWF\n>t3\nWWWWWWYYYYYYFFFFFFPPPPPPCC\n",
        );
        let query = db_from(">q1\nMKVLLACDEFGHIKLMNPQRSTVWY\n");
        let ev = ev_for(&target);
        let mut sp = SearchParams::default();
        sp.max_evalue = 1e-2;
        let idx = KmerIndex::build(&target, 6);
        let pref = search_sequences(&query, &target, &idx, &m, &ev, sp);
        let exh = search_sequences_exhaustive(&query, &target, &m, &ev, sp);
        let exh_targets: HashSet<_> = exh.iter().map(|h| h.target_id.clone()).collect();
        for h in &pref {
            assert!(exh_targets.contains(&h.target_id));
        }
        assert!(exh.iter().any(|h| h.target_id == "t1" && h.pident > 90.0));
    }

    #[test]
    fn sequence_search_end_to_end() {
        let m = SubstitutionMatrix::blosum62();
        let target = db_from(">t1\nMKVLLACDEFGHIKLMNPQRSTVWY\n>t2\nPPPPPCCCCCWWWWWKKKKKDDDDD\n");
        let query = db_from(">q1\nMKVLLACDEFGHIKLMNPQRSTVWY\n");
        let idx = KmerIndex::build(&target, 6);
        let mut sp = SearchParams::default();
        sp.max_evalue = 1e6;
        let hits = search_sequences(&query, &target, &idx, &m, &ev_for(&target), sp);
        assert!(!hits.is_empty());
        assert_eq!(hits[0].target_id, "t1");
        assert!(hits[0].pident > 90.0);
    }

    #[test]
    fn reported_numbers_describe_the_same_alignment() {
        // raw score, CIGAR, coordinates, bits and E-value are mutually consistent
        let m = SubstitutionMatrix::blosum62();
        let mut rng = crate::evalue::SplitMix64::new(77);
        let aa = b"ACDEFGHIKLMNPQRSTVWY";
        let base: String = (0..200).map(|_| aa[rng.below(20)] as char).collect();
        let mut fasta = String::new();
        for i in 0..40 {
            let mut s: Vec<u8> = base.bytes().collect();
            for _ in 0..(i * 2) {
                let p = rng.below(s.len());
                s[p] = aa[rng.below(20)];
            }
            if i % 3 == 0 {
                s.splice(50..50, b"GGGGGG".iter().cloned());
            }
            if i % 4 == 0 {
                s.drain(120..127);
            }
            fasta.push_str(&format!(">t{}\n{}\n", i, String::from_utf8(s).unwrap()));
        }
        let target = db_from(&fasta);
        let query = db_from(&format!(">q\n{}\n", base));
        let idx = KmerIndex::build(&target, 5);
        let ev = ev_for(&target);
        for band in [None, Some(2usize)] {
            let mut sp = SearchParams::default();
            sp.prefilter.k = 5;
            sp.max_evalue = 10.0;
            sp.band_radius = band;
            let eng = Engine { mat: &m, targets: &idx, index: Some(&idx), params: sp, stats: Statistics::Fixed(ev.clone()) };
            let mut st = SearchStats::default();
            let res = eng.search_batch(&query, 0, &mut st);
            assert!(res[0].len() > 10);
            assert_eq!(st.score_mismatches, 0);
            let table = {
                let sc = SeqScorer { query: query.num(0), mat: &m };
                ScoreTable::with_bias(&sc, Some(&local_bias_correction(query.num(0), &m)))
            };
            for (h, a) in &res[0] {
                let t = idx.seq(h.target_index);
                let rescored = crate::align::rescore_alignment(&table, t, sp.gaps, a).unwrap();
                assert_eq!(rescored, h.raw_score);
                assert_eq!((h.q_start, h.q_end, h.t_start, h.t_end), (a.query_start, a.query_end, a.target_start, a.target_end));
                assert!((h.bit_score - ev.bit_score(h.raw_score)).abs() < 1e-9);
                assert!((h.evalue / ev.evalue(h.raw_score, h.qlen) - 1.0).abs() < 1e-9);
                assert_eq!(h.mismatches, a.mismatches);
                assert_eq!(h.gap_opens, a.gap_opens);
            }
        }
    }

    #[test]
    fn coverage_filters_continue_past_rejected_hits() {
        // the best hit is a short local match rejected by --min-qcov; the
        // deferred traceback must keep going and report the full-length hit
        let m = SubstitutionMatrix::blosum62();
        let q = "MKVLLACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY";
        let target = db_from(&format!(
            ">short\nWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW{}WWWW\n>full\n{}\n",
            &q[..20],
            q.replace('K', "R")
        ));
        let query = db_from(&format!(">q\n{}\n", q));
        let idx = KmerIndex::build(&target, 5);
        let mut sp = SearchParams::default();
        sp.prefilter.k = 5;
        sp.max_evalue = 1e3;
        sp.max_hits = 1;
        sp.min_query_cov = 0.9;
        let hits = search_sequences(&query, &target, &idx, &m, &ev_for(&target), sp);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].target_id, "full");
        sp.min_query_cov = 0.0;
        sp.min_target_cov = 0.9;
        let hits = search_sequences(&query, &target, &idx, &m, &ev_for(&target), sp);
        assert!(hits.iter().all(|h| h.target_id == "full"));
    }

    #[test]
    fn windows_preserve_coordinates_and_expand_when_clipped() {
        let m = SubstitutionMatrix::blosum62();
        let mut rng = crate::evalue::SplitMix64::new(3);
        let aa = b"ACDEFGHIKLMNPQRSTVWY";
        let q: String = (0..120).map(|_| aa[rng.below(20)] as char).collect();
        let pre: String = (0..5000).map(|_| aa[rng.below(20)] as char).collect();
        let post: String = (0..3000).map(|_| aa[rng.below(20)] as char).collect();
        let target = db_from(&format!(">long\n{}{}{}\n", pre, q, post));
        let query = db_from(&format!(">q\n{}\n", q));
        let idx = KmerIndex::build(&target, 6);
        let ev = ev_for(&target);
        let mut sp = SearchParams::default();
        let full = search_sequences(&query, &target, &idx, &m, &ev, sp);
        sp.window = WindowParams { enabled: true, min_target_len: 1000, margin: 10 };
        let eng = Engine { mat: &m, targets: &idx, index: Some(&idx), params: sp, stats: Statistics::Fixed(ev.clone()) };
        let mut st = SearchStats::default();
        let win = eng.search_batch(&query, 0, &mut st);
        assert!(st.windowed > 0);
        assert_eq!(win[0].len(), 1);
        let h = &win[0][0].0;
        assert_eq!(h.t_start, 5000);
        assert_eq!(h, &full[0]);
    }

    #[test]
    fn profile_search_end_to_end() {
        let m = SubstitutionMatrix::blosum62();
        let msa = ">q\nMKVLLACDEFGHIKLMNPQR\n>h\nMKILLSCDEFGHLKLMNPQR\n>i\nMRVLLACEEFGHIKLLNPQR\n";
        let prof = Profile::from_msa_str(msa, &m, PseudoCountParams::default()).unwrap();
        let target = db_from(">t1\nGGGMKVLLACDEFGHIKLMNPQRGGG\n>t2\nAAAAAAAAAAAAAAAAAAAA\n");
        let idx = KmerIndex::build(&target, 5);
        let ev = EValueParams::for_profile(target.total_residues);
        let mut sp = SearchParams::profile_default();
        sp.max_evalue = 1e6;
        let hits = search_profile(&prof, "profile1", &target, &idx, &m, &ev, sp);
        assert!(!hits.is_empty());
        assert_eq!(hits[0].target_id, "t1");
    }

    #[test]
    fn hit_row_reconstructs_query_anchored_sequence() {
        let m = SubstitutionMatrix::blosum62();
        let msa = ">q\nMKVLLACDEFGHIKLMNPQR\n>h\nMKILLSCDEFGHLKLMNPQR\n";
        let prof = Profile::from_msa_str(msa, &m, PseudoCountParams::default()).unwrap();
        let t = build_aa2num();
        let tgt: Vec<u8> = "MKVLLACDEFGHIKLMNPQR".bytes().map(|b| t[b as usize]).collect();
        let row = hit_to_msa_row(&prof, prof.query_len, &tgt, GapCosts::profile_default()).unwrap();
        assert_eq!(row, tgt);
    }

    #[test]
    fn iterative_one_round_equals_single_pass() {
        let m = SubstitutionMatrix::blosum62();
        let msa_str = ">q\nMKVLLACDEFGHIKLMNPQR\n>h\nMKILLSCDEFGHLKLMNPQR\n";
        let (msa, ql) = crate::profile::parse_msa(msa_str, &m.aa2num).unwrap();
        let target = db_from(">t1\nGGGMKVLLACDEFGHIKLMNPQRGGG\n>t2\nAAAAAAAAAAAAAAAAAAAA\n");
        let idx = KmerIndex::build(&target, 5);
        let ev = EValueParams::for_profile(target.total_residues);
        let mut sp = SearchParams::profile_default();
        sp.max_evalue = 1e6;
        let prof = Profile::from_msa_str(msa_str, &m, PseudoCountParams::default()).unwrap();
        let single = search_profile(&prof, "p", &target, &idx, &m, &ev, sp);
        let iter1 = search_profile_iterative(&msa, ql, "p", &target, &idx, &m, PseudoCountParams::default(), &ev, sp, 1, 1e-3).unwrap();
        assert_eq!(single, iter1);
    }

    #[test]
    fn iterative_never_loses_hits() {
        let m = SubstitutionMatrix::blosum62();
        let query = "WYFMKVLLACDEFGHIKWYFM";
        let close = "WYFMKILLSCDEFGHLKWYFM";
        let distant = "AAPMKILQSCEEFGHLKAPAA";
        let msa_str = format!(">q\n{}\n", query);
        let (msa, ql) = crate::profile::parse_msa(&msa_str, &m.aa2num).unwrap();
        let target = db_from(&format!(">close\n{}\n>distant\n{}\n", close, distant));
        let idx = KmerIndex::build(&target, 5);
        let ev = EValueParams::for_profile(target.total_residues);
        let mut sp = SearchParams::profile_default();
        sp.max_evalue = 0.05;
        sp.prefilter.kmer_score = 12;
        sp.prefilter.min_diag_score = 1;
        let r1 = search_profile_iterative(&msa, ql, "p", &target, &idx, &m, PseudoCountParams::default(), &ev, sp, 1, 1.0).unwrap();
        let r3 = search_profile_iterative(&msa, ql, "p", &target, &idx, &m, PseudoCountParams::default(), &ev, sp, 3, 1.0).unwrap();
        assert!(r3.len() >= r1.len());
    }

    #[test]
    fn single_query_candidate_parallelism_matches_sequential() {
        let m = SubstitutionMatrix::blosum62();
        let mut rng = crate::evalue::SplitMix64::new(8);
        let aa = b"ACDEFGHIKLMNPQRSTVWY";
        let base: String = (0..150).map(|_| aa[rng.below(20)] as char).collect();
        let mut fasta = String::new();
        for i in 0..300 {
            let mut s: Vec<u8> = base.bytes().collect();
            for _ in 0..(i % 60) {
                let p = rng.below(s.len());
                s[p] = aa[rng.below(20)];
            }
            fasta.push_str(&format!(">t{}\n{}\n", i, String::from_utf8(s).unwrap()));
        }
        let target = db_from(&fasta);
        let idx = KmerIndex::build(&target, 6);
        let query = db_from(&format!(">q\n{}\n", base));
        let mut sp = SearchParams::default();
        sp.max_evalue = 1.0;
        let eng = Engine { mat: &m, targets: &idx, index: Some(&idx), params: sp, stats: Statistics::Fixed(ev_for(&target)) };
        let mut s1 = SearchStats::default();
        let mut s2 = SearchStats::default();
        let a = eng.search_query(0, "q", QueryKind::Sequence(query.num(0)), &mut WorkerScratch::new(), &mut s1, true);
        let b = eng.search_query(0, "q", QueryKind::Sequence(query.num(0)), &mut WorkerScratch::new(), &mut s2, false);
        assert_eq!(a, b);
        assert!(a.len() > 50);
    }
}
