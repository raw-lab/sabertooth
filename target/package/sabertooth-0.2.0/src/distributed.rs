//! Optional distributed search backend built on **HydraMPP** (RAW-lab's Ray-like
//! scheduler). Compiled only under `--features distributed`.
//!
//! The target database is split into shards; each shard is an independent,
//! fully self-contained task (`ShardInput` carries the shard's target sequences,
//! the broadcast query set, and the Karlin–Altschul parameters — nothing is
//! captured), so the same code runs on one workstation's cores or across a
//! multi-node cluster with `--hydra-host` / `--hydra-client`. The substitution
//! matrix is deterministic and rebuilt on each worker rather than shipped.
//!
//! E-values use the *global* database residue count passed in the input, so a
//! sharded run yields the same significance as a single-process search.

use hydra_mpp_core::prelude::*;
use serde::{Deserialize, Serialize};

use crate::evalue::EValueParams;
use crate::fasta::{Record, SeqDb};
use crate::matrix::SubstitutionMatrix;
use crate::prefilter::{KmerIndex, PrefilterParams};
use crate::search::{search_sequences, Hit, SearchParams};

/// A self-contained unit of search work: one target shard + the full query set.
#[derive(Serialize, Deserialize)]
struct ShardInput {
    /// (id, numeric-residue sequence) for each query (broadcast to every shard).
    queries: Vec<(String, Vec<u8>)>,
    /// (id, numeric-residue sequence) for the targets in this shard.
    targets: Vec<(String, Vec<u8>)>,
    k: usize,
    max_evalue: f64,
    min_query_cov: f64,
    max_hits: usize,
    kmer_score: i32,
    lambda: f64,
    ka_k: f64,
    total_db_residues: usize,
}

/// A serializable hit (mirrors `search::Hit`, which is not itself serde-derived
/// to keep the default build dependency-free).
#[derive(Serialize, Deserialize)]
struct HitWire {
    query_id: String,
    target_id: String,
    pident: f64,
    aln_len: usize,
    mismatches: usize,
    gap_opens: usize,
    q_start: usize,
    q_end: usize,
    t_start: usize,
    t_end: usize,
    evalue: f64,
    bit_score: f64,
}

/// A shard's output plus the GPU device(s) HydraMPP pinned to the task while it
/// ran (empty when no GPU was reserved). The device string is where a real GPU
/// kernel would be handed `CUDA_VISIBLE_DEVICES`.
#[derive(Serialize, Deserialize)]
struct ShardResult {
    gpu: String,
    hits: Vec<HitWire>,
}

fn db_from_pairs(pairs: &[(String, Vec<u8>)]) -> SeqDb {
    let mut total = 0usize;
    let records = pairs
        .iter()
        .map(|(id, num)| {
            total += num.len();
            Record {
                header: id.clone(),
                id: id.clone(),
                seq: Vec::new(), // search uses `num`; letters not needed here
                num: num.clone(),
            }
        })
        .collect();
    SeqDb {
        records,
        total_residues: total,
    }
}

/// The registered task: search one shard and return its hits together with the
/// GPU device HydraMPP pinned to it. Pure function of its input — safe to run on
/// any node. When a GPU is reserved (`--gpus`), `hydra_mpp_core::current_gpus()`
/// reports the pinned device; a CUDA kernel would run here on that device.
fn search_shard(input: ShardInput) -> ShardResult {
    let gpu = hydra_mpp_core::cuda_visible_devices();
    let m = SubstitutionMatrix::blosum62();
    let query = db_from_pairs(&input.queries);
    let target = db_from_pairs(&input.targets);
    let idx = KmerIndex::build(&target, input.k);
    let ev = EValueParams::from_raw(input.lambda, input.ka_k, input.total_db_residues);
    let mut sp = SearchParams::default();
    sp.prefilter = PrefilterParams::default();
    sp.prefilter.k = input.k;
    sp.prefilter.kmer_score = input.kmer_score;
    sp.max_evalue = input.max_evalue;
    sp.min_query_cov = input.min_query_cov;
    sp.max_hits = input.max_hits;

    let hits = search_sequences(&query, &target, &idx, &m, &ev, sp)
        .into_iter()
        .map(|h| HitWire {
            query_id: h.query_id,
            target_id: h.target_id,
            pident: h.pident,
            aln_len: h.aln_len,
            mismatches: h.mismatches,
            gap_opens: h.gap_opens,
            q_start: h.q_start,
            q_end: h.q_end,
            t_start: h.t_start,
            t_end: h.t_end,
            evalue: h.evalue,
            bit_score: h.bit_score,
        })
        .collect();
    ShardResult { gpu, hits }
}

/// Result of a distributed run: merged hits plus, for the GPU path, the device
/// string HydraMPP pinned to each shard task (for reporting).
pub struct DistOutcome {
    pub hits: Vec<Hit>,
    pub gpu_devices: Vec<String>,
    pub shards: usize,
}

/// Run a sequence search across HydraMPP workers. `num_shards` target shards are
/// dispatched; when `gpus_per_shard > 0` each shard task reserves that many GPU
/// devices (GPU-aware scheduling + device pinning) instead of a CPU slot, and the
/// pinned device is reported. Results are merged, sorted by (query, E-value), and
/// capped per query at `params.max_hits`.
pub fn search_sequences_distributed(
    query_db: &SeqDb,
    target_db: &SeqDb,
    ev: &EValueParams,
    params: SearchParams,
    num_shards: usize,
    gpus_per_shard: usize,
) -> Result<DistOutcome, String> {
    let hydra = Hydra::new();
    hydra.register("sabertooth_search_shard", search_shard);
    // Config::from_args() consumes any --hydra-* flags; with none this is local
    // multi-core mode. On a cluster the head runs --hydra-host and workers join;
    // --hydra-gpus N advertises (or simulates) GPU devices on a node.
    hydra
        .init(Config::from_args())
        .map_err(|e| format!("HydraMPP init failed: {e}"))?;

    let n = target_db.records.len();
    let shards = num_shards.max(1).min(n.max(1));
    let shard_size = (n + shards - 1) / shards.max(1);

    let queries: Vec<(String, Vec<u8>)> = query_db
        .records
        .iter()
        .map(|r| (r.id.clone(), r.num.clone()))
        .collect();

    let mut inputs: Vec<ShardInput> = Vec::new();
    let mut start = 0usize;
    while start < n {
        let end = (start + shard_size).min(n);
        inputs.push(ShardInput {
            queries: queries.clone(),
            targets: target_db.records[start..end]
                .iter()
                .map(|r| (r.id.clone(), r.num.clone()))
                .collect(),
            k: params.prefilter.k,
            max_evalue: params.max_evalue,
            min_query_cov: params.min_query_cov,
            max_hits: params.max_hits,
            kmer_score: params.prefilter.kmer_score,
            lambda: ev.lambda,
            ka_k: ev.k,
            total_db_residues: target_db.total_residues,
        });
        start = end;
    }
    let launched = inputs.len();

    // GPU-scheduled path: reserve devices per shard (cpus(0) so scheduling is
    // gated on GPU slots, not CPU slots). CPU path: the batch `map`.
    let shard_results: Vec<ShardResult> = if gpus_per_shard > 0 {
        let mut ids = Vec::with_capacity(inputs.len());
        for input in &inputs {
            let id = hydra
                .task("sabertooth_search_shard")
                .cpus(0)
                .gpus(gpus_per_shard)
                .submit(input)
                .map_err(|e| format!("HydraMPP GPU submit failed: {e}"))?;
            ids.push(id);
        }
        let mut out = Vec::with_capacity(ids.len());
        for id in ids {
            out.push(
                hydra
                    .get_typed::<ShardResult>(id)
                    .map_err(|e| format!("HydraMPP get failed: {e}"))?,
            );
        }
        out
    } else {
        hydra
            .map("sabertooth_search_shard", &inputs)
            .map_err(|e| format!("HydraMPP map failed: {e}"))?
    };
    hydra.shutdown();

    let mut gpu_devices: Vec<String> = Vec::new();
    let mut hits: Vec<Hit> = Vec::new();
    for sr in shard_results {
        if !sr.gpu.is_empty() {
            gpu_devices.push(sr.gpu);
        }
        for w in sr.hits {
            hits.push(Hit {
                query_id: w.query_id,
                target_id: w.target_id,
                pident: w.pident,
                aln_len: w.aln_len,
                mismatches: w.mismatches,
                gap_opens: w.gap_opens,
                q_start: w.q_start,
                q_end: w.q_end,
                t_start: w.t_start,
                t_end: w.t_end,
                evalue: w.evalue,
                bit_score: w.bit_score,
            });
        }
    }

    // global ordering + per-query cap (shards each capped locally; a hit could be
    // bumped from one shard but survive globally, so re-cap after the merge)
    hits.sort_by(|a, b| {
        a.query_id
            .cmp(&b.query_id)
            .then(a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal))
    });
    if params.max_hits > 0 {
        let mut capped: Vec<Hit> = Vec::with_capacity(hits.len());
        let mut cur = String::new();
        let mut count = 0usize;
        for h in hits {
            if h.query_id != cur {
                cur = h.query_id.clone();
                count = 0;
            }
            if count < params.max_hits {
                capped.push(h);
                count += 1;
            }
        }
        hits = capped;
    }
    Ok(DistOutcome {
        hits,
        gpu_devices,
        shards: launched,
    })
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
    fn distributed_matches_single_process() {
        let m = SubstitutionMatrix::blosum62();
        let target = db_from(
            ">t1\nMKVLLACDEFGHIKLMNPQRSTVWY\n\
             >t2\nMKILLSCDEFGHLKLMNPQRSTVWF\n\
             >t3\nWWWWWWYYYYYYFFFFFFPPPPPPCC\n\
             >t4\nMADENKLKLGSGSFGEVFLVKHKESG\n",
        );
        let query = db_from(">q1\nMKVLLACDEFGHIKLMNPQRSTVWY\n>q2\nMADENKLKLGSGSFGEVFLVKHKESG\n");
        let ev = EValueParams::new(&m, target.total_residues);
        let mut sp = SearchParams::default();
        sp.max_evalue = 1e-2;

        let idx = KmerIndex::build(&target, sp.prefilter.k);
        let single = search_sequences(&query, &target, &idx, &m, &ev, sp);
        let dist = search_sequences_distributed(&query, &target, &ev, sp, 3, 0).unwrap();

        // same (query,target) hit set, regardless of sharding
        let key = |h: &Hit| (h.query_id.clone(), h.target_id.clone());
        let mut a: Vec<_> = single.iter().map(key).collect();
        let mut b: Vec<_> = dist.hits.iter().map(key).collect();
        a.sort();
        b.sort();
        assert_eq!(a, b, "distributed hit set differs from single-process");
    }
}
