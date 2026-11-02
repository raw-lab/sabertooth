//! Residue-balanced target shards with persistent seed indexes, and the
//! **out-of-core** search that uses them.
//!
//! A database's seed index is cut into contiguous ordinal ranges holding
//! ≈ equal residues; each shard gets its own persistent index (same seed
//! pattern and masking as the parent) in a directory with a manifest. Later
//! runs reuse the shards while the manifest still matches the parent (no
//! rebuild, no rescan); a shard whose file no longer matches the manifest is
//! rebuilt.
//!
//! Shards serve two callers:
//!
//! * **Out-of-core local search** (`search --shards N`, or automatically when
//!   `--memory-limit` is smaller than the index): every query batch is searched
//!   against one shard at a time ([`search_batch_sharded`]), so only one
//!   shard's pages need to be resident instead of the whole index. On a
//!   machine with less RAM than the index this replaces page-cache thrashing
//!   by one sequential pass over the shards per query batch.
//! * The distributed backend (`--features distributed`), which sends query
//!   batches to shard tasks on HydraMPP workers.
//!
//! Both merge with [`merge_shard_hits`]: hits carry the global target ordinal,
//! and per query they are sorted with the local engine's own order
//! ([`crate::search::hit_order`]) and cut at `max_hits`. Every stage before
//! that is per target (seeds come from the query alone; evidence, gapped
//! scores and tracebacks depend only on the pair) and E-values use the
//! whole-database size, so with the default unlimited `--max-candidates` /
//! `--max-seed-occ` the output is **identical** to an unsharded search (tests
//! compare every output field). Those two caps apply per shard when set.

use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::evalue::EValueParams;
use crate::fasta::SeqDb;
use crate::index::{build_index, BuildOptions, Hasher64, RecordSource, SeedIndex, SourceFingerprint};
use crate::matrix::SubstitutionMatrix;
use crate::search::{hit_order, Engine, Hit, SearchParams, SearchStats, Statistics};

const MANIFEST: &str = "MANIFEST";
const MANIFEST_MAGIC: &str = "sabertooth-shards 1";

/// One shard: target ordinals `start..end` with their own seed index.
#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "distributed", derive(serde::Serialize, serde::Deserialize))]
pub struct Shard {
    pub id: usize,
    pub start: usize,
    pub end: usize,
    pub residues: u64,
    pub path: PathBuf,
}

/// The shards of one database.
#[derive(Clone, Debug)]
pub struct ShardSet {
    pub dir: PathBuf,
    pub shards: Vec<Shard>,
    pub total_seqs: usize,
    pub total_residues: u64,
    /// `true` if existing shard files were reused.
    pub reused: bool,
    pub seconds: f64,
}

/// Residue-balanced contiguous split of `lens` into at most `n` non-empty parts.
pub fn balance(lens: &[usize], n: usize) -> Vec<(usize, usize)> {
    let total: u64 = lens.iter().map(|&l| l as u64).sum();
    let n = n.clamp(1, lens.len().max(1));
    let mut out = Vec::with_capacity(n);
    let mut start = 0usize;
    let mut acc = 0u64;
    for s in 0..n {
        if start >= lens.len() {
            break;
        }
        let remaining_parts = n - s;
        if remaining_parts == 1 {
            out.push((start, lens.len()));
            break;
        }
        let target = total * (s as u64 + 1) / n as u64;
        let mut end = start;
        // take at least one sequence; leave at least one per remaining part
        let max_end = lens.len() - (remaining_parts - 1);
        while end < max_end && (end == start || acc + (lens[end] as u64) / 2 < target) {
            acc += lens[end] as u64;
            end += 1;
        }
        out.push((start, end));
        start = end;
    }
    out
}

struct RangeSource<'a> {
    idx: &'a SeedIndex,
    start: usize,
    end: usize,
}

impl<'a> RecordSource for RangeSource<'a> {
    fn for_each(&mut self, f: &mut dyn FnMut(&[u8], &[u8]) -> io::Result<()>) -> io::Result<()> {
        for i in self.start..self.end {
            f(self.idx.id(i).as_bytes(), self.idx.seq(i))?;
        }
        Ok(())
    }
}

fn manifest_text(parent_identity: &str, parent: &SeedIndex, ranges: &[(usize, usize)], residues: &[u64]) -> String {
    let m = parent.meta();
    let mut s = format!(
        "{}\nparent {}\nsequences {}\nresidues {}\npattern {}\nmask {}\nshards {}\n",
        MANIFEST_MAGIC,
        parent_identity,
        parent.len(),
        parent.total_residues(),
        m.pattern.mask_string(),
        if m.masked {
            format!("{},{},{}", m.mask.window, m.mask.locut, m.mask.hicut)
        } else {
            "none".into()
        },
        ranges.len()
    );
    for (i, (&(a, b), r)) in ranges.iter().zip(residues).enumerate() {
        s.push_str(&format!("shard {} {} {} {} shard_{}.stidx\n", i, a, b, r, i));
    }
    s
}

/// A string that changes whenever the content of `parent` can change: the
/// index path, size and modification time for a persistent index, or the
/// fingerprint of `source` (the FASTA / database it was built from) for a
/// temporary one.
pub fn parent_identity(parent: &SeedIndex, source: &str) -> io::Result<String> {
    if parent.is_transient() {
        let fp = SourceFingerprint::of(source)?;
        Ok(format!("source {:?}", fp))
    } else {
        let p = parent.path();
        let md = fs::metadata(p)?;
        let mt = md
            .modified()
            .ok()
            .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        Ok(format!("index {} {} {}", p.display(), md.len(), mt))
    }
}

/// Default shard directory: `<index>.shards-<n>` next to a persistent index,
/// or a directory in `tmp` named after the identity of a temporary one.
pub fn default_shard_dir(parent: &SeedIndex, identity: &str, n: usize, tmp: &Path) -> PathBuf {
    if parent.is_transient() {
        let mut h = Hasher64::default();
        h.update(identity.as_bytes());
        tmp.join(format!("sabertooth-shards-{:016x}-{}", h.finish(), n))
    } else {
        PathBuf::from(format!("{}.shards-{}", parent.path().display(), n))
    }
}

/// Create (or reuse) the residue-balanced shards of `parent` in `dir`.
/// `parent_identity` must change whenever the parent's content can change
/// (see [`parent_identity`]).
pub fn prepare_shards(parent: &SeedIndex, parent_identity: &str, n: usize, dir: &Path, opts: &BuildOptions) -> io::Result<ShardSet> {
    let t0 = Instant::now();
    if parent.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "cannot shard an empty database"));
    }
    let lens: Vec<usize> = (0..parent.len()).map(|i| parent.seq_len(i)).collect();
    let ranges = balance(&lens, n);
    let residues: Vec<u64> = ranges.iter().map(|&(a, b)| lens[a..b].iter().map(|&l| l as u64).sum()).collect();
    let want = manifest_text(parent_identity, parent, &ranges, &residues);
    let shards: Vec<Shard> = ranges
        .iter()
        .zip(&residues)
        .enumerate()
        .map(|(i, (&(a, b), &r))| Shard {
            id: i,
            start: a,
            end: b,
            residues: r,
            path: dir.join(format!("shard_{}.stidx", i)),
        })
        .collect();
    let mpath = dir.join(MANIFEST);
    let reusable = fs::read_to_string(&mpath).map(|have| have == want).unwrap_or(false)
        && shards.iter().all(|sh| match SeedIndex::open(&sh.path) {
            Ok(ix) => ix.len() == sh.end - sh.start && ix.total_residues() as u64 == sh.residues,
            Err(_) => false,
        });
    if !reusable {
        fs::create_dir_all(dir)?;
        let _ = fs::remove_file(&mpath); // invalidate before rewriting shards
        let mut o = opts.clone();
        o.pattern = parent.pattern().clone();
        o.mask = if parent.meta().masked { Some(parent.meta().mask) } else { None };
        for sh in &shards {
            let mut src = RangeSource {
                idx: parent,
                start: sh.start,
                end: sh.end,
            };
            build_index(&mut src, SourceFingerprint::transient(), &sh.path, &o)?;
        }
        let tmp = dir.join(format!("{}.tmp-{}", MANIFEST, std::process::id()));
        fs::write(&tmp, &want)?;
        fs::rename(&tmp, &mpath)?;
    }
    Ok(ShardSet {
        dir: dir.to_path_buf(),
        shards,
        total_seqs: parent.len(),
        total_residues: parent.total_residues() as u64,
        reused: reusable,
        seconds: t0.elapsed().as_secs_f64(),
    })
}

/// Open a shard index and check it against the manifest and the search.
pub fn open_shard(sh: &Shard, k: usize) -> Result<SeedIndex, String> {
    let idx = SeedIndex::open(&sh.path).map_err(|e| format!("cannot open shard index {}: {}", sh.path.display(), e))?;
    if idx.len() != sh.end - sh.start || idx.total_residues() as u64 != sh.residues {
        return Err(format!(
            "shard index {} holds {} sequences / {} residues, the manifest expects {} / {} (rebuild the shards)",
            sh.path.display(),
            idx.len(),
            idx.total_residues(),
            sh.end - sh.start,
            sh.residues
        ));
    }
    if idx.k() != k {
        return Err(format!("shard index seed weight {} != search k {}", idx.k(), k));
    }
    Ok(idx)
}

/// Sort each query's hits (pooled from several shards, global target
/// ordinals) in the engine's order and cut at `max_hits` (0 = no limit).
/// Returns the number of hits kept.
pub fn merge_shard_hits(merged: &mut [Vec<Hit>], max_hits: usize) -> u64 {
    let mut kept = 0u64;
    for v in merged.iter_mut() {
        v.sort_by(hit_order);
        if max_hits > 0 {
            v.truncate(max_hits);
        }
        kept += v.len() as u64;
    }
    kept
}

/// Search one query batch against every shard in turn and merge the hits.
///
/// Each shard index is mapped only while it is searched (and its presence
/// bitmap freed afterwards), so the resident set is one shard plus the
/// per-worker scratch. `ev` must describe the **whole** database. Counters in
/// `stats` are summed over shards (work actually done, e.g. seeds are
/// generated once per shard), except `queries`, `query_residues` and `hits`,
/// which count the batch once.
pub fn search_batch_sharded(
    set: &ShardSet,
    mat: &SubstitutionMatrix,
    params: SearchParams,
    ev: &EValueParams,
    batch: &SeqDb,
    first_ordinal: usize,
    stats: &mut SearchStats,
) -> Result<Vec<Vec<Hit>>, String> {
    let mut merged: Vec<Vec<Hit>> = vec![Vec::new(); batch.len()];
    let mut sum = SearchStats::default();
    let mut once: Option<(u64, u64)> = None;
    for sh in &set.shards {
        let idx = open_shard(sh, params.prefilter.k).map_err(|e| format!("shard {}: {}", sh.id, e))?;
        let eng = Engine {
            mat,
            targets: &idx,
            index: Some(&idx),
            params,
            stats: Statistics::Fixed(ev.clone()),
        };
        let mut st = SearchStats::default();
        let res = eng.search_batch(batch, first_ordinal, &mut st);
        if res.len() != batch.len() {
            return Err(format!("shard {} returned {} result lists for {} queries", sh.id, res.len(), batch.len()));
        }
        for (slot, hs) in merged.iter_mut().zip(res) {
            slot.extend(hs.into_iter().map(|(mut h, _)| {
                h.target_index += sh.start;
                h
            }));
        }
        once.get_or_insert((st.queries, st.query_residues));
        sum.add(&st);
    }
    let kept = merge_shard_hits(&mut merged, params.max_hits);
    if let Some((q, r)) = once {
        sum.queries = q;
        sum.query_residues = r;
    }
    sum.hits = kept;
    stats.add(&sum);
    Ok(merged)
}

/// Smallest shard worth cutting: indexes up to this size are never sharded,
/// and shards are never planned smaller than this.
pub const MIN_SHARD_BYTES: usize = 16 << 20;

/// Number of shards so that one shard plus `other_bytes` of (estimated,
/// worst-case) working memory fits in `budget` bytes, leaving a quarter of the
/// budget for page-cache turnover. `index_bytes` is the size of the whole
/// index file. Returns 1 when the index fits or is smaller than
/// [`MIN_SHARD_BYTES`]; when the budget leaves less than that per shard,
/// shards of `MIN_SHARD_BYTES` are used (the budget is then tight).
pub fn shards_for_budget(budget: usize, index_bytes: usize, other_bytes: usize, max_shards: usize) -> usize {
    let usable = (budget / 4) * 3;
    let room = usable.saturating_sub(other_bytes).max(MIN_SHARD_BYTES);
    if index_bytes <= room {
        return 1;
    }
    index_bytes.div_ceil(room).clamp(1, max_shards.max(1))
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::alphabet::build_aa2num;
    use crate::prefilter::SeedPattern;
    use std::io::Cursor;

    fn db_from(fasta: &str) -> SeqDb {
        SeqDb::from_reader(Cursor::new(fasta), &build_aa2num()).unwrap()
    }

    /// Related families with indels, repeats, low-complexity tails and
    /// duplicate query identifiers.
    pub(crate) fn synthetic(seed: u64) -> (SeqDb, SeqDb) {
        let mut rng = crate::evalue::SplitMix64::new(seed);
        let aa = b"ACDEFGHIKLMNPQRSTVWY";
        let mut fam = Vec::new();
        for _ in 0..12 {
            let l = 60 + rng.below(220);
            fam.push((0..l).map(|_| aa[rng.below(20)]).collect::<Vec<u8>>());
        }
        let mutate = |rng: &mut crate::evalue::SplitMix64, s: &[u8], rate: usize| -> Vec<u8> {
            let mut o = Vec::new();
            for &c in s {
                let r = rng.below(100);
                if r < 3 {
                    continue;
                }
                if r < 6 {
                    o.push(aa[rng.below(20)]);
                }
                o.push(if rng.below(100) < rate { aa[rng.below(20)] } else { c });
            }
            o
        };
        let mut t = String::new();
        for i in 0..300 {
            let f = &fam[i % fam.len()];
            let mut s = mutate(&mut rng, f, 20 + (i % 50));
            if i % 9 == 0 {
                s.extend(std::iter::repeat_n(b'Q', 30));
            }
            let pre: Vec<u8> = (0..rng.below(80)).map(|_| aa[rng.below(20)]).collect();
            t.push_str(&format!(">t{} family {}\n{}{}\n", i, i % fam.len(), String::from_utf8(pre).unwrap(), String::from_utf8(s).unwrap()));
        }
        let mut q = String::new();
        for (i, f) in fam.iter().enumerate() {
            let s = mutate(&mut rng, f, 25);
            q.push_str(&format!(">q{}\n{}\n", i, String::from_utf8(s).unwrap()));
        }
        q.push_str(&format!(">q0\n{}\n", String::from_utf8(fam[3].clone()).unwrap())); // duplicate id
        (db_from(&q), db_from(&t))
    }

    pub(crate) fn local(queries: &SeqDb, idx: &SeedIndex, sp: SearchParams) -> Vec<Vec<Hit>> {
        let m = SubstitutionMatrix::blosum62();
        let ev = EValueParams::sequence(sp.gaps, sp.comp_bias_corr, idx.total_residues(), idx.len()).unwrap();
        let eng = Engine {
            mat: &m,
            targets: idx,
            index: Some(idx),
            params: sp,
            stats: Statistics::Fixed(ev),
        };
        let mut st = SearchStats::default();
        eng.search_batch(queries, 0, &mut st)
            .into_iter()
            .map(|v| v.into_iter().map(|(h, _)| h).collect())
            .collect()
    }

    #[test]
    fn balance_is_contiguous_nonempty_and_residue_balanced() {
        let lens: Vec<usize> = (0..1000).map(|i| 50 + (i * 37) % 400).collect();
        for n in [1usize, 2, 3, 7, 16, 999, 1000, 5000] {
            let r = balance(&lens, n);
            assert_eq!(r.len(), n.min(lens.len()));
            assert_eq!(r[0].0, 0);
            assert_eq!(r.last().unwrap().1, lens.len());
            for w in r.windows(2) {
                assert_eq!(w[0].1, w[1].0);
            }
            assert!(r.iter().all(|&(a, b)| b > a));
            if n <= 16 {
                let total: usize = lens.iter().sum();
                let ideal = total / n;
                for &(a, b) in &r {
                    let s: usize = lens[a..b].iter().sum();
                    assert!(s.abs_diff(ideal) <= 450, "n={} part {} vs ideal {}", n, s, ideal);
                }
            }
        }
    }

    #[test]
    fn sharded_search_is_identical_to_unsharded() {
        let (queries, targets) = synthetic(7);
        let dir = std::env::temp_dir().join(format!("sabertooth_shard_test_{}", std::process::id()));
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
        // the cap must bind for some query so the merge's re-cap is exercised
        assert!(want.iter().any(|v| v.len() == 25));
        let m = SubstitutionMatrix::blosum62();
        let ev = EValueParams::sequence(sp.gaps, sp.comp_bias_corr, parent.total_residues(), parent.len()).unwrap();
        let mut unsharded_stats = SearchStats::default();
        {
            let eng = Engine { mat: &m, targets: &parent, index: Some(&parent), params: sp, stats: Statistics::Fixed(ev.clone()) };
            eng.search_batch(&queries, 0, &mut unsharded_stats);
        }
        for shards in [1usize, 2, 3, 7] {
            let sdir = dir.join(format!("s{}", shards));
            let set = prepare_shards(&parent, "test-parent", shards, &sdir, &opts).unwrap();
            assert!(!set.reused);
            assert_eq!(set.shards.len(), shards);
            // two batches with an ordinal offset, as the CLI streams them
            let split = queries.len() / 2;
            let mut st = SearchStats::default();
            let a = search_batch_sharded(&set, &m, sp, &ev, &queries.subset(0..split), 0, &mut st).unwrap();
            let b = search_batch_sharded(&set, &m, sp, &ev, &queries.subset(split..queries.len()), split, &mut st).unwrap();
            let got: Vec<Vec<Hit>> = a.into_iter().chain(b).collect();
            assert_eq!(got.len(), want.len());
            for (qi, (g, w)) in got.iter().zip(&want).enumerate() {
                assert_eq!(g, w, "query {} differs with {} shards", qi, shards);
            }
            assert_eq!(st.queries, queries.len() as u64);
            assert_eq!(st.query_residues, unsharded_stats.query_residues);
            assert_eq!(st.hits, want.iter().map(|v| v.len() as u64).sum::<u64>());
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
    fn corrupt_shard_fails_loudly_and_is_rebuilt() {
        let (queries, targets) = synthetic(3);
        let dir = std::env::temp_dir().join(format!("sabertooth_shard_bad_{}", std::process::id()));
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
        let m = SubstitutionMatrix::blosum62();
        let ev = EValueParams::sequence(sp.gaps, sp.comp_bias_corr, parent.total_residues(), parent.len()).unwrap();
        let mut st = SearchStats::default();
        let e = search_batch_sharded(&set, &m, sp, &ev, &queries, 0, &mut st).expect_err("a corrupt shard must fail the search");
        assert!(e.contains("shard 1"), "{}", e);
        // and prepare_shards notices and rebuilds it
        let again = prepare_shards(&parent, "p", 2, &dir, &opts).unwrap();
        assert!(!again.reused);
        assert!(search_batch_sharded(&again, &m, sp, &ev, &queries, 0, &mut st).is_ok());
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn shard_count_follows_the_budget() {
        let mib = 1usize << 20;
        // fits: one shard
        assert_eq!(shards_for_budget(1024 * mib, 500 * mib, 100 * mib, 64), 1);
        // 873 MiB index, 60 MiB working memory, 512 MiB budget: 384 - 60 = 324 MiB per shard
        assert_eq!(shards_for_budget(512 * mib, 873 * mib, 60 * mib, 64), 3);
        // 256 MiB budget: 192 - 60 = 132 MiB per shard → 7 shards
        assert_eq!(shards_for_budget(256 * mib, 873 * mib, 60 * mib, 64), 7);
        // capped by the number of sequences
        assert_eq!(shards_for_budget(256 * mib, 873 * mib, 60 * mib, 4), 4);
        // no room left by the working-memory estimate: minimum-size shards
        assert_eq!(shards_for_budget(64 * mib, 873 * mib, 60 * mib, 64), 55);
        // a small index is never sharded, whatever the budget
        assert_eq!(shards_for_budget(64 * mib, 10 * mib, 60 * mib, 64), 1);
    }
}
