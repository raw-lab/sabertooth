//! Distributed search on **HydraMPP** (RAW-lab's Ray-like scheduler), compiled
//! only with `--features distributed`.
//!
//! ## Design (and what v0.2 got wrong)
//!
//! v0.2 shipped every shard's raw target sequences *and the whole query set* in
//! each task, rebuilt a hash-map index per task, re-created `SearchParams`
//! from a hand-picked subset of fields (so prefilter, masking, coverage and
//! gap settings silently fell back to defaults) and used ungapped statistics.
//! Here:
//!
//! * **Persistent, residue-balanced shards** ([`crate::shard`], shared with
//!   the out-of-core local search). The target's seed index is cut into
//!   contiguous ordinal ranges holding ≈ equal residues; each shard gets its
//!   own persistent seed index (same seed pattern and masking as the parent)
//!   under a shard directory with a manifest. Later runs reuse the shards when
//!   the manifest still matches the parent (no rebuild, no rescan). Workers
//!   open shard files by path — on a multi-node cluster the shard directory
//!   must be on shared storage (or replicated to the same path).
//! * **Complete configuration.** The validated [`SearchParams`] travels as one
//!   serde value, so local and distributed runs use the same struct, defaults
//!   and option parsing; nothing is re-derived on the worker.
//! * **Bounded query batches.** Queries are streamed in the same bounded
//!   batches as a local run; a batch is sent to every shard, never the whole
//!   query file.
//! * **Global statistics.** Each worker computes E-values from the calibrated
//!   gapped parameters and the **whole** database size (residues and
//!   sequences), so a hit's E-value does not depend on its shard.
//! * **Deterministic merge.** Hits carry the query input ordinal and the
//!   global target ordinal; per query they are merged with the local engine's
//!   own order ([`crate::search::hit_order`]) and cut at `--max-hits`. Every
//!   stage before that is per target (seeds are generated from the query
//!   alone; evidence, gapped scores and tracebacks depend only on the pair),
//!   so with the default unlimited `--max-candidates` / `--max-seed-occ` the
//!   output is **identical** to a single-process search — tests check every
//!   output field. (Those two caps are per-shard when set.)
//! * **No nested parallelism.** Each task reserves one CPU slot and searches
//!   its batch on its own thread; HydraMPP runs as many tasks as there are
//!   slots. Use at least as many shards as total CPU slots.
//! * **GPUs.** `--gpus N` reserves N GPU slots per task (scheduling and
//!   `CUDA_VISIBLE_DEVICES` pinning, reported per shard). There is no GPU
//!   kernel: all alignment runs on the CPU.

use hydra_mpp_core::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, OnceLock};

use crate::evalue::EValueParams;
use crate::fasta::{LoadOptions, SeqDb};
use crate::index::SeedIndex;
use crate::matrix::SubstitutionMatrix;
use crate::search::{Engine, Hit, SearchParams, SearchStats, Statistics};
use crate::shard::{merge_shard_hits, open_shard};
pub use crate::shard::{balance, prepare_shards, Shard, ShardSet};

/// HydraMPP runtime configuration (`Config::from_args()` reads `--hydra-*`).
pub use hydra_mpp_core::Config;

/// Registered task name (versioned: workers must run the same build).
pub const TASK: &str = "sabertooth_search_shard_v3";

/// A unit of work: one query batch against one shard.
#[derive(Serialize, Deserialize)]
pub struct ShardTask {
    pub shard: Shard,
    /// The complete, validated configuration of the run.
    pub params: SearchParams,
    /// Whole-database size for the statistics.
    pub global_residues: u64,
    pub global_seqs: usize,
    /// Input ordinal of `queries[0]`.
    pub first_ordinal: usize,
    /// `(identifier, numeric residues)` of this batch only.
    pub queries: Vec<(String, Vec<u8>)>,
}

/// A task's result: hits per batch query (global ordinals), counters, pinned GPUs.
#[derive(Serialize, Deserialize)]
pub struct ShardOutput {
    pub shard: usize,
    pub hits: Vec<Vec<Hit>>,
    pub stats: SearchStats,
    pub gpu: String,
    pub error: Option<String>,
}

/// Worker-side cache of opened shard indexes, keyed by path and revalidated
/// by (size, mtime) so a rebuilt shard is reopened, never served stale.
fn open_cached(shard: &Shard, k: usize) -> Result<Arc<SeedIndex>, String> {
    type Entry = ((u64, u128), Arc<SeedIndex>);
    static CACHE: OnceLock<Mutex<HashMap<PathBuf, Entry>>> = OnceLock::new();
    let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    let path: &Path = &shard.path;
    let md = fs::metadata(path).map_err(|e| format!("cannot open shard index {}: {}", path.display(), e))?;
    let stamp = (
        md.len(),
        md.modified().ok().and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok()).map(|d| d.as_nanos()).unwrap_or(0),
    );
    if let Some((st, ix)) = cache.lock().unwrap().get(path) {
        if *st == stamp && ix.len() == shard.end - shard.start && ix.total_residues() as u64 == shard.residues && ix.k() == k {
            return Ok(ix.clone());
        }
    }
    let ix = Arc::new(open_shard(shard, k)?);
    cache.lock().unwrap().insert(path.to_path_buf(), (stamp, ix.clone()));
    Ok(ix)
}

fn run_inner(t: &ShardTask) -> Result<(Vec<Vec<Hit>>, SearchStats), String> {
    let idx = open_cached(&t.shard, t.params.prefilter.k)?;
    let m = SubstitutionMatrix::blosum62();
    let ev = EValueParams::sequence(t.params.gaps, t.params.comp_bias_corr, t.global_residues as usize, t.global_seqs)?;
    let mut q = SeqDb::new(LoadOptions::SEARCH);
    for (id, num) in &t.queries {
        q.push_encoded(id, None, num, None);
    }
    let eng = Engine {
        mat: &m,
        targets: idx.as_ref(),
        index: Some(idx.as_ref()),
        params: t.params,
        stats: Statistics::Fixed(ev),
    };
    let mut st = SearchStats::default();
    let res = eng.search_batch_sequential(&q, t.first_ordinal, &mut st);
    let hits = res
        .into_iter()
        .map(|v| {
            v.into_iter()
                .map(|(mut h, _)| {
                    h.target_index += t.shard.start;
                    h
                })
                .collect()
        })
        .collect();
    Ok((hits, st))
}

/// The registered task (a pure function of its input, safe on any node).
pub fn run_shard_task(t: ShardTask) -> ShardOutput {
    let gpu = hydra_mpp_core::cuda_visible_devices();
    match run_inner(&t) {
        Ok((hits, stats)) => ShardOutput {
            shard: t.shard.id,
            hits,
            stats,
            gpu,
            error: None,
        },
        Err(e) => ShardOutput {
            shard: t.shard.id,
            hits: Vec::new(),
            stats: SearchStats::default(),
            gpu,
            error: Some(e),
        },
    }
}

/// Merged result of one query batch.
pub struct BatchResult {
    /// Per batch query, merged and capped, in the local engine's order.
    pub hits: Vec<Vec<Hit>>,
    /// Counters summed over shards.
    pub stats: SearchStats,
    /// GPU devices pinned to the batch's tasks (empty without `--gpus`).
    pub gpu_devices: Vec<String>,
}

/// Head-side driver: owns the HydraMPP runtime and the shard set.
pub struct DistributedSearcher {
    hydra: Hydra,
    shards: Vec<Shard>,
    params: SearchParams,
    global_residues: u64,
    global_seqs: usize,
    gpus_per_shard: usize,
}

impl DistributedSearcher {
    /// Start HydraMPP (`cfg`: e.g. `Config::from_args()` for `--hydra-*`
    /// flags, or `Config::local()`) and register the shard task.
    pub fn start(set: &ShardSet, params: SearchParams, gpus_per_shard: usize, cfg: Config) -> Result<Self, String> {
        let hydra = Hydra::new();
        hydra.register(TASK, run_shard_task);
        hydra.init(cfg).map_err(|e| format!("HydraMPP init failed: {e}"))?;
        Ok(DistributedSearcher {
            hydra,
            shards: set.shards.clone(),
            params,
            global_residues: set.total_residues,
            global_seqs: set.total_seqs,
            gpus_per_shard,
        })
    }

    /// Search one query batch on every shard and merge deterministically.
    pub fn search_batch(&self, batch: &SeqDb, first_ordinal: usize) -> Result<BatchResult, String> {
        let queries: Vec<(String, Vec<u8>)> = (0..batch.len()).map(|i| (batch.id(i).to_string(), batch.num(i).to_vec())).collect();
        let mut ids = Vec::with_capacity(self.shards.len());
        for sh in &self.shards {
            let task = ShardTask {
                shard: sh.clone(),
                params: self.params,
                global_residues: self.global_residues,
                global_seqs: self.global_seqs,
                first_ordinal,
                queries: queries.clone(),
            };
            let id = self
                .hydra
                .task(TASK)
                .cpus(1)
                .gpus(self.gpus_per_shard)
                .submit(&task)
                .map_err(|e| format!("HydraMPP submit failed: {e}"))?;
            ids.push(id);
        }
        let mut merged: Vec<Vec<Hit>> = vec![Vec::new(); batch.len()];
        let mut stats = SearchStats::default();
        let mut gpu_devices = Vec::new();
        for id in ids {
            let out: ShardOutput = self.hydra.get_typed(id).map_err(|e| {
                if self.gpus_per_shard > 0 {
                    format!(
                        "HydraMPP could not run a shard task reserving {} GPU(s): {e} \
                         (declare or simulate devices with --hydra-gpus N; GPUs are only reserved \
                         and pinned, all alignment runs on the CPU)",
                        self.gpus_per_shard
                    )
                } else {
                    format!("HydraMPP get failed: {e}")
                }
            })?;
            if let Some(e) = out.error {
                return Err(format!("shard {}: {}", out.shard, e));
            }
            if out.hits.len() != batch.len() {
                return Err(format!("shard {} returned {} result lists for {} queries", out.shard, out.hits.len(), batch.len()));
            }
            for (slot, hs) in merged.iter_mut().zip(out.hits) {
                slot.extend(hs);
            }
            stats.add(&out.stats);
            if !out.gpu.is_empty() {
                gpu_devices.push(format!("shard {}: {}", out.shard, out.gpu));
            }
        }
        stats.hits = merge_shard_hits(&mut merged, self.params.max_hits);
        Ok(BatchResult {
            hits: merged,
            stats,
            gpu_devices,
        })
    }

    /// Stop the runtime.
    pub fn shutdown(self) {
        self.hydra.shutdown();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::BuildOptions;
    use crate::prefilter::SeedPattern;
    use crate::shard::tests::{local, synthetic};

    #[test]
    fn params_round_trip_completely() {
        let mut sp = SearchParams::sequence_default();
        sp.max_evalue = 0.123;
        sp.min_query_cov = 0.4;
        sp.min_target_cov = 0.2;
        sp.max_hits = 17;
        sp.band_radius = Some(9);
        sp.prefilter.kmer_score = 19;
        sp.prefilter.max_candidates = 77;
        sp.window.enabled = true;
        sp.mask_scoring = None;
        let bytes = bincode::serialize(&sp).unwrap();
        let back: SearchParams = bincode::deserialize(&bytes).unwrap();
        assert_eq!(back, sp);
        assert_eq!(format!("{:?}", back), format!("{:?}", sp));
    }

    #[test]
    fn distributed_output_is_identical_to_single_process() {
        let (queries, targets) = synthetic(7);
        let dir = std::env::temp_dir().join(format!("sabertooth_dist_test_{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        let mut opts = BuildOptions::new(SeedPattern::contiguous(6));
        opts.mask = Some(crate::mask::MaskParams::default());
        opts.tmp_dir = std::env::temp_dir();
        let (parent, _) = SeedIndex::build_transient(&targets, &opts).unwrap();
        let mut sp = SearchParams::sequence_default();
        sp.max_evalue = 10.0;
        sp.max_hits = 25;
        sp.prefilter.mask = Some(crate::mask::MaskParams::default());
        let want = local(&queries, &parent, sp);
        assert!(want.iter().map(|v| v.len()).sum::<usize>() > 100, "test data too easy/hard");
        // the cap must bind for some query so the global re-cap is exercised
        assert!(want.iter().any(|v| v.len() == 25));
        for shards in [1usize, 2, 3, 7] {
            let sdir = dir.join(format!("s{}", shards));
            let set = prepare_shards(&parent, "test-parent", shards, &sdir, &opts).unwrap();
            assert!(!set.reused);
            assert_eq!(set.shards.len(), shards);
            let ds = DistributedSearcher::start(&set, sp, 0, Config::local().quiet(true).port(0)).unwrap();
            // two batches with an ordinal offset, as the CLI streams them
            let split = queries.len() / 2;
            let a = ds.search_batch(&queries.subset(0..split), 0).unwrap();
            let b = ds.search_batch(&queries.subset(split..queries.len()), split).unwrap();
            ds.shutdown();
            let got: Vec<Vec<Hit>> = a.hits.into_iter().chain(b.hits).collect();
            assert_eq!(got.len(), want.len());
            for (qi, (g, w)) in got.iter().zip(&want).enumerate() {
                assert_eq!(g, w, "query {} differs with {} shards", qi, shards);
            }
            // reuse: same manifest → no rebuild
            let again = prepare_shards(&parent, "test-parent", shards, &sdir, &opts).unwrap();
            assert!(again.reused);
            // a different parent identity invalidates the shards
            let other = prepare_shards(&parent, "another-parent", shards, &sdir, &opts).unwrap();
            assert!(!other.reused);
        }
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn impossible_gpu_reservation_fails_with_a_hint() {
        let (queries, targets) = synthetic(5);
        let dir = std::env::temp_dir().join(format!("sabertooth_dist_gpu_{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        let mut opts = BuildOptions::new(SeedPattern::contiguous(6));
        opts.tmp_dir = std::env::temp_dir();
        let (parent, _) = SeedIndex::build_transient(&targets, &opts).unwrap();
        let mut sp = SearchParams::sequence_default();
        sp.prefilter.mask = None;
        let set = prepare_shards(&parent, "gpu", 2, &dir, &opts).unwrap();
        // no GPUs offered (num_gpus 0): a task reserving one can never be placed
        let ds = DistributedSearcher::start(&set, sp, 1, Config::local().quiet(true).port(0).num_gpus(0)).unwrap();
        let r = ds.search_batch(&queries, 0);
        ds.shutdown();
        let e = r.err().expect("an unplaceable GPU reservation must fail");
        assert!(e.contains("--hydra-gpus"), "{}", e);
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn corrupt_shard_is_reported_not_silently_empty() {
        let (queries, targets) = synthetic(3);
        let dir = std::env::temp_dir().join(format!("sabertooth_dist_bad_{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        let mut opts = BuildOptions::new(SeedPattern::contiguous(6));
        opts.tmp_dir = std::env::temp_dir();
        let (parent, _) = SeedIndex::build_transient(&targets, &opts).unwrap();
        let mut sp = SearchParams::sequence_default();
        sp.prefilter.mask = None;
        let set = prepare_shards(&parent, "p", 2, &dir, &opts).unwrap();
        // truncate shard 1 behind the manifest's back
        let p1 = &set.shards[1].path;
        let len = fs::metadata(p1).unwrap().len();
        let f = fs::OpenOptions::new().write(true).open(p1).unwrap();
        f.set_len(len / 2).unwrap();
        let ds = DistributedSearcher::start(&set, sp, 0, Config::local().quiet(true).port(0)).unwrap();
        let r = ds.search_batch(&queries, 0);
        ds.shutdown();
        let e = r.err().expect("a corrupt shard must fail the search");
        assert!(e.contains("shard 1"), "{}", e);
        // and prepare_shards notices and rebuilds it
        let again = prepare_shards(&parent, "p", 2, &dir, &opts).unwrap();
        assert!(!again.reused);
        let _ = fs::remove_dir_all(&dir);
    }
}
