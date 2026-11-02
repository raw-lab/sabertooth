//! Sabertooth command-line interface.
//!
//! Dependency-light hand-rolled argument parsing. Every sequence argument goes
//! through the shared format-aware loader ([`sabertooth::seqio`]); searches run
//! on a seed index (persistent `<db>.stidx` when present and fresh, otherwise a
//! temporary one), stream queries in bounded batches and stream results through
//! a buffered writer.

use sabertooth::align::{local_align, GapCosts, SeqScorer};
use sabertooth::alphabet::{num_to_aa, PROFILE_AA_SIZE};
use sabertooth::banner;
use sabertooth::evalue::EValueParams;
use sabertooth::fasta::{read_to_string, FastaReader, LoadOptions, SeqDb};
use sabertooth::index::{self, BuildOptions, BuildStats, SeedIndex, SourceFingerprint, Unusable};
use sabertooth::mask::MaskParams;
use sabertooth::matrix::SubstitutionMatrix;
use sabertooth::metrics::{parse_bytes, proc_stats, Json};
use sabertooth::prefilter::{check_k, SeedPattern};
use sabertooth::profile::{Profile, PseudoCountParams};
use sabertooth::search::{Engine, Hit, SearchParams, SearchStats, Statistics, WindowParams, EXTRA_HEADER, M8_HEADER};
use sabertooth::seqio::{self, DynReader, SeqKind};
use sabertooth::VERSION;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::process::ExitCode;
use std::time::Instant;

// ---------------------------------------------------------------------------
// argument parsing
// ---------------------------------------------------------------------------

/// Parsed command line: positionals + `--flag value` options + bare flags.
struct Args {
    positionals: Vec<String>,
    options: HashMap<String, String>,
}

impl Args {
    fn parse(argv: &[String]) -> Args {
        let mut positionals = Vec::new();
        let mut options = HashMap::new();
        let mut i = 0;
        while i < argv.len() {
            let a = &argv[i];
            if a == "-s" {
                if i + 1 < argv.len() && !argv[i + 1].starts_with('-') {
                    options.insert("sensitivity".to_string(), argv[i + 1].clone());
                    i += 2;
                } else {
                    options.insert("sensitivity".to_string(), String::new());
                    i += 1;
                }
            } else if let Some(key) = a.strip_prefix("--") {
                if let Some((k, v)) = key.split_once('=') {
                    options.insert(k.to_string(), v.to_string());
                    i += 1;
                } else if i + 1 < argv.len() && !argv[i + 1].starts_with("--") {
                    options.insert(key.to_string(), argv[i + 1].clone());
                    i += 2;
                } else {
                    options.insert(key.to_string(), String::new());
                    i += 1;
                }
            } else {
                positionals.push(a.clone());
                i += 1;
            }
        }
        Args { positionals, options }
    }

    fn has(&self, key: &str) -> bool {
        self.options.contains_key(key)
    }

    fn opt_str(&self, key: &str) -> Option<&str> {
        self.options.get(key).map(|s| s.as_str())
    }

    fn parse_opt<T: std::str::FromStr>(&self, key: &str, default: T) -> Result<T, String> {
        match self.options.get(key) {
            None => Ok(default),
            Some(v) if v.is_empty() => Err(format!("--{} needs a value", key)),
            Some(v) => v
                .parse::<T>()
                .map_err(|_| format!("invalid value for --{}: '{}'", key, v)),
        }
    }

    fn opt_f64(&self, key: &str, default: f64) -> Result<f64, String> {
        self.parse_opt(key, default)
    }
    fn opt_usize(&self, key: &str, default: usize) -> Result<usize, String> {
        self.parse_opt(key, default)
    }
    fn opt_i32(&self, key: &str, default: i32) -> Result<i32, String> {
        self.parse_opt(key, default)
    }
    fn opt_bool(&self, key: &str, default: bool) -> Result<bool, String> {
        match self.options.get(key).map(|s| s.as_str()) {
            None => Ok(default),
            Some("") | Some("1") | Some("true") | Some("yes") | Some("on") => Ok(true),
            Some("0") | Some("false") | Some("no") | Some("off") => Ok(false),
            Some(v) => Err(format!("invalid value for --{}: '{}' (use 0 or 1)", key, v)),
        }
    }
    fn opt_bytes(&self, key: &str) -> Result<Option<usize>, String> {
        match self.options.get(key) {
            None => Ok(None),
            Some(v) => parse_bytes(v).map(Some).map_err(|e| format!("--{}: {}", key, e)),
        }
    }

    /// Reject options this command does not know: a misspelled option (say
    /// `--evalu 1e-5`) must not silently run with the default (v0.3 only
    /// warned).
    fn check_known(&self, cmd: &str, known: &[&str]) -> Result<(), String> {
        let mut bad: Vec<&str> = self
            .options
            .keys()
            .map(|k| k.as_str())
            .filter(|k| !(k.starts_with("hydra-") || known.contains(k) || COMMON_OPTS.contains(k)))
            .collect();
        if bad.is_empty() {
            return Ok(());
        }
        bad.sort_unstable();
        Err(format!(
            "unknown option{} for `{}`: {} (see `sabertooth {} --help`)",
            if bad.len() > 1 { "s" } else { "" },
            cmd,
            bad.iter().map(|k| format!("--{}", k)).collect::<Vec<_>>().join(", "),
            cmd
        ))
    }
}

const COMMON_OPTS: &[&str] = &["out", "threads", "stats", "stats-json", "memory-limit", "tmp-dir", "quiet"];

fn main() -> ExitCode {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    if argv.is_empty() {
        print_help();
        return ExitCode::SUCCESS;
    }
    let cmd = argv[0].clone();
    // `<command> --help` / `-h`: the full option list
    if argv[1..].iter().any(|a| a == "--help" || a == "-h") {
        print_help();
        return ExitCode::SUCCESS;
    }
    let args = Args::parse(&argv[1..]);

    let result = match cmd.as_str() {
        "version" | "--version" | "-V" => {
            print_version();
            Ok(())
        }
        "info" => {
            print_info();
            Ok(())
        }
        "doctor" => run_doctor(),
        "help" | "--help" | "-h" => {
            print_help();
            Ok(())
        }
        "createdb" => run_createdb(&args),
        "msa2profile" => run_msa2profile(&args),
        "search" => with_pool(&args, || run_search(&args)),
        "profilesearch" => with_pool(&args, || run_profilesearch(&args)),
        "align" => run_align(&args),
        "result2msa" => with_pool(&args, || run_result2msa(&args)),
        "translatesearch" => with_pool(&args, || run_translatesearch(&args)),
        "profileprofile" => run_profileprofile(&args),
        "cluster" => with_pool(&args, || run_cluster(&args)),
        "easysearch" => with_pool(&args, || run_easysearch(&args)),
        "makedb" => run_makedb(&args),
        "createindex" => with_pool(&args, || run_createindex(&args)),
        "indexinfo" => run_indexinfo(&args),
        "convert2fasta" => run_convert2fasta(&args),
        "taxonomy" => with_pool(&args, || run_taxonomy(&args)),
        "calibrate" => with_pool(&args, || run_calibrate(&args)),
        other => {
            eprintln!("sabertooth: unknown command '{}'\n", other);
            print_help();
            return ExitCode::from(2);
        }
    };

    sabertooth::search::flush_candidate_trace();
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("sabertooth: error: {}", e);
            ExitCode::FAILURE
        }
    }
}

/// Run `f` inside a rayon pool sized by `--threads` and the memory plan.
fn with_pool<F: FnOnce() -> Result<(), String> + Send>(args: &Args, f: F) -> Result<(), String> {
    let threads = plan_threads(args)?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| format!("cannot create thread pool: {}", e))?;
    pool.install(f)
}

fn cpus() -> usize {
    std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
}

// ---------------------------------------------------------------------------
// memory plan
// ---------------------------------------------------------------------------

/// How a memory budget is spent (all values in bytes).
#[derive(Clone, Debug)]
struct Plan {
    budget: Option<usize>,
    threads: usize,
    hit_block: usize,
    tb_cells: usize,
    batch_seqs: usize,
    batch_residues: usize,
    build_budget: usize,
    per_worker: usize,
    notes: Vec<String>,
}

const FIXED_OVERHEAD: usize = 24 << 20;

/// Per-worker scratch estimate: target-block hit buffer and its radix-sort
/// copy (16 B/hit), segment-evidence chunk and SIMD buffers (~1 MiB),
/// traceback matrix (1 B/cell) and fixed buffers. Seeds of the query being
/// processed (24 B per seed with postings) come on top and scale with query
/// length × `--max-kmers`.
fn per_worker_bytes(hit_block: usize, tb_cells: usize) -> usize {
    hit_block * 16 + tb_cells + (3 << 20)
}

fn plan_threads(args: &Args) -> Result<usize, String> {
    let requested = args.opt_usize("threads", cpus())?.max(1);
    Ok(requested)
}

fn make_plan(args: &Args, sp: &mut SearchParams) -> Result<Plan, String> {
    let requested = args.opt_usize("threads", cpus())?.max(1);
    let budget = args.opt_bytes("memory-limit")?;
    let mut hit_block = sp.prefilter.hit_block;
    let tb_cells = sabertooth::align::FULL_TRACEBACK_MAX_CELLS;
    let mut batch_seqs = args.opt_usize("query-batch", 1024)?.max(1);
    let mut batch_residues = 1usize << 26;
    let mut threads = requested;
    let mut build_budget = 1usize << 30;
    let mut notes = Vec::new();
    if let Some(b) = budget {
        if b < FIXED_OVERHEAD + (8 << 20) {
            return Err(format!(
                "--memory-limit {} is too small (need at least {} MiB)",
                b,
                (FIXED_OVERHEAD + (8 << 20)) >> 20
            ));
        }
        let results_per_query = sp.max_hits.clamp(10, 1000) * 256;
        // query batch + its results: at most 1/4 of the budget
        let batch_bytes = (b / 4).min(256 << 20);
        batch_seqs = batch_seqs.min((batch_bytes / (results_per_query + 1024)).max(16));
        batch_residues = (batch_bytes / 2).max(1 << 16);
        let avail = b.saturating_sub(FIXED_OVERHEAD + batch_bytes);
        // shrink per-worker scratch until at least one worker fits (the
        // traceback matrix limit is fixed so alignments never depend on it)
        while per_worker_bytes(hit_block, tb_cells) > avail && hit_block > 1 << 14 {
            hit_block = (hit_block / 2).max(1 << 14);
        }
        let pw = per_worker_bytes(hit_block, tb_cells);
        let fit = (avail / pw).max(1);
        if fit < threads {
            notes.push(format!("memory limit allows {} worker(s) of {} requested", fit, threads));
            threads = fit;
        }
        build_budget = (b / 2).max(16 << 20);
    }
    sp.prefilter.hit_block = hit_block;
    Ok(Plan {
        budget,
        threads,
        hit_block,
        tb_cells,
        batch_seqs,
        batch_residues,
        build_budget,
        per_worker: per_worker_bytes(hit_block, tb_cells),
        notes,
    })
}

fn plan_json(p: &Plan) -> Json {
    Json::new()
        .raw("budget_bytes", &p.budget.map(|b| b.to_string()).unwrap_or_else(|| "null".into()))
        .num("threads", p.threads)
        .num("hit_block", p.hit_block)
        .num("full_traceback_cells", p.tb_cells)
        .num("per_worker_scratch_estimate_bytes", p.per_worker)
        .num("query_batch_seqs", p.batch_seqs)
        .num("query_batch_residues", p.batch_residues)
        .num("index_build_budget_bytes", p.build_budget)
}

// ---------------------------------------------------------------------------
// banners
// ---------------------------------------------------------------------------

fn print_version() {
    print!("{}", banner::version_banner(VERSION));
}

fn print_info() {
    print!("{}", banner::INFO_BANNER);
    let m = SubstitutionMatrix::blosum62();
    println!("  Matrix        : {}", m.name);
    println!("  Alphabet      : ACDEFGHIKLMNPQRSTVWY + X  (20 aa, MMseqs2 order)");
    println!("  Bit factor    : {}  (search scores)", m.bit_factor);
    println!("  Score bias    : {}", m.bias);
    println!("  Lambda (hdr)  : {:.5}  (real-valued half-bit matrix)", m.lambda);
    let ug = EValueParams::new(&m, 1);
    println!("  Lambda (int)  : {:.4}  (ungapped, integer matrix; diagnostic only)", ug.lambda);
    for e in sabertooth::evalue::SEQUENCE_TABLE {
        let lm = e.params.length.unwrap();
        println!(
            "  Gapped stats  : gaps {}/{} comp-bias {}: lambda {:.4} K {:.4} alpha/lambda {:.1} beta {:.0}",
            e.gaps.open, e.gaps.extend, e.comp_bias as u8, e.params.lambda, e.params.k, lm.alpha_d_lambda, lm.beta
        );
    }
    print!("  Background    :");
    for i in 0..PROFILE_AA_SIZE {
        if i % 10 == 0 {
            print!("\n     ");
        }
        print!(" {}:{:.4}", num_to_aa(i as u8) as char, m.p_back[i]);
    }
    println!("\n");
    println!("  Profile PSSM  : bit factor 8, pca=1.0 pcb=1.5 (MMseqs2 defaults); statistics calibrated per profile");
    println!("  Pipeline      : seed index -> k-mer prefilter -> SIMD gapped score -> exact traceback -> calibrated E-values");
    println!(
        "  SIMD kernels  : {} (SABERTOOTH_SIMD / SABERTOOTH_SIMD_U8 override)",
        sabertooth::simd::kernel_name(sabertooth::simd::kernel_choice())
    );
    println!();
    println!("  Integer substitution scores (BLOSUM62, bit_factor={}, bias={}):", m.bit_factor, m.bias);
    print!("{}", indent(&m.score_table(), "  "));
}

fn indent(s: &str, pad: &str) -> String {
    s.lines().map(|l| format!("{}{}\n", pad, l)).collect()
}

fn print_help() {
    print!("{}", banner::HELP_BANNER);
    println!(
        r#"
USAGE:
  sabertooth <command> [args] [options]

DATABASES (two separate steps):
  makedb <in.fasta[.gz]> <db>            Sequence database: write an MMseqs2-format DB
                                         (<db>, .index, .dbtype, _h, ...; readable by mmseqs)
  createindex <db|fasta> [out.stidx]     Seed index: build a persistent, memory-mapped
                                         k-mer index (default <db>.stidx). Searches use it
                                         automatically when it is present, up to date and
                                         built with the search's seed (default 111001101;
                                         --k / --spaced / --mask as for search;
                                         --index-layout auto|sparse|direct, default auto:
                                         direct k-mer table once the sparse one would
                                         exceed 64 MiB)
  indexinfo <file.stidx> [--verify]      Show index metadata; --verify checks all checksums
  convert2fasta <db> [out.fasta]         Read a database (or FASTA) back to FASTA
  createdb <in> <out.fasta>              Validate + normalise sequences to FASTA

SEARCH:
  search <query> <target>                Sequence-vs-sequence search
  easysearch <query> <target>            Exhaustive SW against every target (no prefilter)
  profilesearch <msa> <target>           Profile (PSSM) search
  profileprofile <query.msa> <t.msa...>  Profile-vs-profile search (co-emission)
  translatesearch <nucl_q> <prot_target> Six-frame translated nucleotide search
  result2msa <query> <target>            Search, then emit an a3m MSA per query
  align <query> <target>                 Pairwise local alignments (query 1 vs targets)
  taxonomy <query> <target>              LCA taxonomic assignment from hits
  cluster <in>                           Linclust-style clustering (rep→member TSV)

  <query>/<target>/<in> may be FASTA, gzip FASTA, '-' (stdin), an MMseqs2 database
  prefix, or a .stidx seed index (identifiers only, no header descriptions).

OTHER:
  msa2profile <msa>                      Build a PSSM from an MSA (a3m/aligned FASTA)
  calibrate                              Fit gapped statistics (lambda, K, length model)
  info | doctor | version | help

COMMON OPTIONS:
  --out <path>            Write results to a file instead of stdout
  --threads <n>           Worker threads (default: all cores)
  --memory-limit <size>   Memory budget, e.g. 2G: sizes threads, batches and buffers;
                          `search` also counts the mapped index: if it does not fit next
                          to the working memory, the index is searched as residue-balanced
                          shards one at a time (out-of-core; identical output)
  --tmp-dir <dir>         Where temporary indexes go (default: system temp)
  --stats                 Print stage timings and counters to stderr
  --stats-json <file>     Write timings, counters, peak RSS and page faults as JSON

SEARCH OPTIONS:
  --evalue <float>        Max E-value to report (default 1e-3). E-values are calibrated
                          expected chance hits in this database, not probabilities
  --max-hits <int>        Max hits per query, 0 = unlimited (default 300; output limit)
  --min-cov, --min-qcov <f>  Min fraction of the query in the alignment (default 0)
  --min-tcov <f>          Min fraction of the target in the alignment (default 0)
  -s <float>              Sensitivity preset 1..9 (default 5.7). Sequence search: sets
                          --kmer-score/--max-kmers together (1: 32/64, 2: 30/128,
                          3: 28/128, 4: 26/256, 5: 24/256, 5.7: 22/512, 6.5: 20/512,
                          7.5: 18/1024, 8.5: 18/2048, 9: 16/4096; interpolated between).
                          Explicit options below override it
  --spaced <mask>         Seed pattern, e.g. 110101 (must start/end with 1). Default for
                          sequence search and createindex: 111001101 (weight 6, span 9)
  --k <int>               Contiguous k-mer seed instead (profile search default: 5;
                          `--k 6` = the seed of v0.3 and its indexes)
  --kmer-score <int>      Similar-k-mer score threshold (seq 22, profile 40)
  --max-kmers <int>       Similar k-mers per window: the exact top-N (default 512)
  --min-diag-score <int>  Ungapped prefilter gate, raw score. Default (sequence search):
                          42 for a 10^8-residue database, 2.5 lower per e-fold below and
                          1.0 higher per e-fold above (33 at 2.6 M, 38 at 17 M, 41 at
                          68 M, 43 at 265 M residues); profile search 30
  --ungapped-xdrop <int>  Ungapped evidence. 0 (default): best segment within
                          --diag-window of a diagonal's seed hits, scored in SIMD
                          batches; x > 0: v0.3's X-drop extension (drop-off x)
  --diag-window <int>     Residues beyond the seed hits searched for that segment
                          (default 48; 0 = the whole diagonal)
  --min-hits <int>        Min seed hits per target (default 1)
  --diag-hits <int>       Min seed hits on one diagonal (default 1; 2 = double hit)
  --prefilter-diags <int> Diagonals per target scored/chained (default 4)
  --max-candidates <int>  Candidate budget per query, best first (default 0 = all)
  --max-seed-occ <int>    Skip k-mers with longer posting lists (default 0 = never)
  --mask 0|1              Soft low-complexity masking of seeds (default 1)
  --mask-scoring 0|1      Also hard-mask low-complexity query positions in scoring
                          (default 1 for sequence queries; SEG 15/2.3/2.3)
  --comp-bias-corr 0|1    Composition-bias correction of scores (default 1)
  --band <radius>         Initial traceback band radius (default 16; exact: widened
                          until the optimal score is reproduced)
  --window 0|1            Seed-guided windows for targets >= --window-min-len (default 0;
                          heuristic: see README)
  --hit-block <int>       Seed hits buffered per target block (default 1048576; reduced
                          by --memory-limit; results do not depend on it)
  --query-batch <int>     Queries per batch (default 1024)
  --shards <int>          Search the index as N residue-balanced shards, one at a time
                          (default: from --memory-limit, else 1). Shards are persistent
                          (<index>.shards-<N>, or --shard-dir) and reused; output is
                          identical (--max-candidates / --max-seed-occ apply per shard)
  --shard-dir <dir>       Where shards live (default <index>.shards-<N>, or a hashed
                          directory under --tmp-dir for FASTA targets)
  --require-index         Fail instead of building a temporary index
  --extra-cols            Append raw, qlen, tlen, qcov, tcov, cigar columns

PROFILE OPTIONS:
  --pca <float> / --pcb <float> / --bias <float>   Pseudocounts (1.0 / 1.5 / 0.0)
  --num-iterations <n>    PSI-BLAST-style rounds (default 1)
  --inclusion-evalue <f>  Max E-value to fold a hit into the next round (default 1e-3)
  --profile-calibration-samples <n>  Random targets per PSSM fit (default 2000)

TRANSLATED / CLUSTER / TAXONOMY OPTIONS:
  --min-orf <int>         Min ORF length in aa for translatesearch (default 20)
  --min-seq-id <float>    Min identity for cluster membership (default 0.5)
  --min-cov <float>       Min coverage for cluster membership (default 0.8)
  --nodes / --names / --seqmap   NCBI taxonomy inputs (taxonomy)
  --lca-frac <float>      Keep hits within this score fraction of best; 0 = all

DISTRIBUTED OPTIONS (build with `--features distributed`):
  --distributed           Run search across HydraMPP workers/nodes; output is identical
                          to a single-process search (same parameters, global statistics)
  --shards <int>          Residue-balanced target shards, each with its own persistent
                          seed index (default 4 with --distributed; use >= total CPU slots)
  --shard-dir <dir>       As for local shards; reused while the manifest matches. Must
                          be shared storage on a cluster
  --gpus <int>            GPUs reserved per shard (scheduling/pinning only: the
                          alignment kernels run on CPU; no GPU compute is used)
  --hydra-host / --hydra-client <ip> / --hydra-gpus <n> / --hydra-port <n> /
  --hydra-quiet           (consumed by HydraMPP)

ENVIRONMENT:
  SABERTOOTH_SIMD=avx2|sse2|scalar   Force a kernel family (default: best available)
  SABERTOOTH_SIMD_U8=0               Disable the u8 first scoring pass
                                     (results are identical for every choice)

EXAMPLES:
  sabertooth makedb proteins.fasta.gz protDB
  sabertooth createindex protDB --threads 8 --memory-limit 4G
  sabertooth search queries.fasta protDB --out hits.m8 --stats-json stats.json
  sabertooth profilesearch family.a3m protDB --evalue 1e-5
  sabertooth easysearch queries.fasta small_db.fasta
"#
    );
}

// ---------------------------------------------------------------------------
// doctor
// ---------------------------------------------------------------------------

fn run_doctor() -> Result<(), String> {
    print!("{}", banner::DOCTOR_BANNER);
    println!();
    let mut ok = true;
    let m = SubstitutionMatrix::blosum62();
    let ww = m.score(m.aa2num[b'W' as usize], m.aa2num[b'W' as usize]);
    let check_matrix = ww > 0 && m.lambda > 0.0;
    report("substitution matrix (BLOSUM62) reconstructed", check_matrix);
    ok &= check_matrix;
    let bsum: f64 = (0..PROFILE_AA_SIZE).map(|i| m.p_back[i]).sum();
    let check_bg = (bsum - 1.0).abs() < 0.01;
    report(&format!("background distribution sums to 1 ({:.5})", bsum), check_bg);
    ok &= check_bg;
    let mut check_r = true;
    for b in 0..PROFILE_AA_SIZE {
        let s: f32 = (0..PROFILE_AA_SIZE).map(|a| m.r[a][b]).sum();
        if (s - 1.0).abs() > 0.03 {
            check_r = false;
        }
    }
    report("pseudocount matrix R columns normalised", check_r);
    ok &= check_r;
    let ev = EValueParams::new(&m, 1_000_000);
    let check_lambda = ev.lambda > 0.25 && ev.lambda < 0.45;
    report(&format!("ungapped Karlin-Altschul lambda in range ({:.4})", ev.lambda), check_lambda);
    ok &= check_lambda;
    let gapped = EValueParams::sequence(GapCosts::sequence_default(), true, 1_000_000, 1000);
    let check_table = gapped.as_ref().map(|g| g.lambda > 0.2 && g.lambda < ev.lambda).unwrap_or(false);
    report("calibrated gapped statistics table present", check_table);
    ok &= check_table;
    let msa = ">q\nMKVLLACDEFGHIK\n>h\nMKILLSCDEFGHLK\n>i\nMRVLLACEEFGHIK\n";
    let prof = Profile::from_msa_str(msa, &m, PseudoCountParams::default());
    let check_profile = match &prof {
        Ok(p) => p.query_len == 14 && p.set_size == 3,
        Err(_) => false,
    };
    report("profile/PSSM construction from MSA", check_profile);
    ok &= check_profile;
    let check_align = if let Ok(p) = &prof {
        let sc = sabertooth::align::ProfileScorer { profile: p };
        let t: Vec<u8> = "MKVLLACDEFGHIK".bytes().map(|b| m.aa2num[b as usize]).collect();
        local_align(&sc, &t, GapCosts::profile_default())
            .map(|a| a.raw_score > 0 && sabertooth::align::rescore_alignment(&sc, &t, GapCosts::profile_default(), &a).is_ok())
            .unwrap_or(false)
    } else {
        false
    };
    report("Smith-Waterman-Gotoh alignment (traceback rescored)", check_align);
    ok &= check_align;
    let t = sabertooth::alphabet::build_aa2num();
    let db = SeqDb::from_reader(std::io::Cursor::new(">t\nGGGMKVLLACDEFGHIKGGG\n"), &t).map_err(|e| e.to_string())?;
    let check_index = match SeedIndex::build_transient(&db, &BuildOptions::new(SeedPattern::sequence_default())) {
        Ok((idx, _)) => idx.distinct_kmers() > 0,
        Err(_) => false,
    };
    report("seed index build + memory map", check_index);
    ok &= check_index;
    println!();
    println!("  environment:");
    println!("    sabertooth version : {}", VERSION);
    println!("    available cores    : {}", cpus());
    println!("    pointer width      : {} bit", std::mem::size_of::<usize>() * 8);
    println!();
    if ok {
        println!("  RESULT: all checks passed — fangs are sharp. \\V/");
        Ok(())
    } else {
        Err("one or more self-checks failed".into())
    }
}

fn report(label: &str, ok: bool) {
    let mark = if ok { "[ ok ]" } else { "[FAIL]" };
    println!("  {} {}", mark, label);
}

// ---------------------------------------------------------------------------
// loading helpers
// ---------------------------------------------------------------------------

/// Load a sequence argument through the shared format-aware loader.
fn load_seqs(path: &str, opts: LoadOptions, kind: SeqKind, role: &str) -> Result<SeqDb, String> {
    let aa2num = sabertooth::alphabet::build_aa2num();
    let loaded = seqio::load_sequences(path, opts, &aa2num, kind)
        .map_err(|e| format!("cannot read {} '{}': {}", role, path, e))?;
    for w in &loaded.warnings {
        eprintln!("sabertooth: warning: {}", w);
    }
    Ok(loaded.db)
}

/// Queries, streamed in bounded batches.
enum QueryStream {
    Fasta(Box<FastaReader<DynReader>>, String),
    Db(SeqDb, usize),
}

impl QueryStream {
    fn open(path: &str, kind: SeqKind) -> Result<QueryStream, String> {
        let fmt = seqio::detect_format(path).map_err(|e| format!("cannot read query '{}': {}", path, e))?;
        match fmt {
            seqio::InputFormat::Fasta { .. } => {
                let fr = seqio::open_fasta(path).map_err(|e| format!("cannot read query '{}': {}", path, e))?;
                Ok(QueryStream::Fasta(Box::new(fr), path.to_string()))
            }
            seqio::InputFormat::MmseqsDb | seqio::InputFormat::SeedIndex => {
                Ok(QueryStream::Db(load_seqs(path, LoadOptions::SEARCH, kind, "query")?, 0))
            }
        }
    }

    /// Next batch of at most `max_seqs` sequences / `max_residues` residues.
    fn next_batch(&mut self, max_seqs: usize, max_residues: usize, aa2num: &[u8; 256]) -> Result<Option<SeqDb>, String> {
        match self {
            QueryStream::Fasta(fr, path) => {
                let mut db = SeqDb::new(LoadOptions::SEARCH);
                let (mut h, mut s) = (Vec::new(), Vec::new());
                while db.len() < max_seqs && db.total_residues < max_residues {
                    let more = fr.next_record(&mut h, &mut s).map_err(|e| format!("cannot read query '{}': {}", path, e))?;
                    if !more {
                        break;
                    }
                    db.push_record(&h, &s, aa2num);
                }
                if db.is_empty() {
                    if fr.records_read() == 0 {
                        return Err(format!("cannot read query '{}': no FASTA records found (the input is empty)", path));
                    }
                    return Ok(None);
                }
                Ok(Some(db))
            }
            QueryStream::Db(all, pos) => {
                if *pos >= all.len() {
                    return Ok(None);
                }
                let mut end = *pos;
                let mut res = 0usize;
                while end < all.len() && end - *pos < max_seqs && res < max_residues {
                    res += all.seq_len(end);
                    end += 1;
                }
                let b = all.subset(*pos..end);
                *pos = end;
                Ok(Some(b))
            }
        }
    }
}

/// Buffered output (file or stdout).
fn open_output(args: &Args) -> Result<BufWriter<Box<dyn Write>>, String> {
    let w: Box<dyn Write> = match args.opt_str("out") {
        Some(p) if !p.is_empty() => Box::new(File::create(p).map_err(|e| format!("cannot create '{}': {}", p, e))?),
        _ => Box::new(io::stdout()),
    };
    Ok(BufWriter::with_capacity(1 << 20, w))
}

fn write_hits<W: Write>(out: &mut W, hits: &[Hit], extra: bool) -> Result<(), String> {
    for h in hits {
        let line = if extra {
            format!("{}\t{}\n", h.to_m8(), h.extra_cols())
        } else {
            format!("{}\n", h.to_m8())
        };
        out.write_all(line.as_bytes()).map_err(|e| format!("write failed: {}", e))?;
    }
    Ok(())
}

fn write_header<W: Write>(out: &mut W, extra: bool) -> Result<(), String> {
    let h = if extra {
        format!("{}\t{}\n", M8_HEADER, EXTRA_HEADER)
    } else {
        format!("{}\n", M8_HEADER)
    };
    out.write_all(h.as_bytes()).map_err(|e| format!("write failed: {}", e))
}

// ---------------------------------------------------------------------------
// search parameters
// ---------------------------------------------------------------------------

fn pc_params(args: &Args) -> Result<PseudoCountParams, String> {
    Ok(PseudoCountParams {
        pca: args.opt_f64("pca", 1.0)? as f32,
        pcb: args.opt_f64("pcb", 1.5)? as f32,
        score_bias: args.opt_f64("bias", 0.0)? as f32,
    })
}

const SEARCH_OPTS: &[&str] = &[
    "evalue", "max-hits", "min-cov", "min-qcov", "min-tcov", "k", "spaced", "kmer-score", "sensitivity",
    "max-kmers", "min-diag-score", "min-hits", "diag-hits", "prefilter-diags", "max-candidates",
    "max-seed-occ", "exact-seeds", "hit-block", "ungapped-xdrop", "diag-window", "mask", "comp-bias-corr", "band", "window",
    "window-min-len", "window-margin", "query-batch", "require-index", "extra-cols", "mask-query", "mask-scoring", "mask-scoring-params",
    "comp-bias-mode",
    "distributed", "shards", "gpus", "shard-dir",
    // consumed by HydraMPP's own parser (Config::from_args)
    "hydra-host", "hydra-client", "hydra-slurm", "hydra-cpus", "hydra-gpus", "hydra-port", "hydra-pool", "hydra-stack",
    "hydra-quiet",
];

/// Requested seed pattern: `--spaced` mask, else contiguous `--k`, else the
/// command's default pattern.
fn requested_pattern(args: &Args, default: &SeedPattern) -> Result<SeedPattern, String> {
    if args.opt_str("spaced").is_some_and(|m| !m.is_empty()) && args.has("k") {
        return Err("--spaced and --k are mutually exclusive (the mask sets the seed weight)".into());
    }
    match args.opt_str("spaced").filter(|m| !m.is_empty()) {
        Some(mask) => {
            let p = SeedPattern::from_mask(mask)?;
            check_k(p.weight())?;
            Ok(p)
        }
        None if args.has("k") => {
            let k = args.opt_usize("k", default.weight())?;
            check_k(k)?;
            Ok(SeedPattern::contiguous(k))
        }
        None => Ok(default.clone()),
    }
}

/// Scale the ungapped gate with the size of the searched database (unless
/// `--min-diag-score` fixes it): see [`sabertooth::prefilter::size_adjusted_gate`].
fn apply_db_size_gate(args: &Args, sp: &mut SearchParams, db_residues: usize) {
    if !args.has("min-diag-score") {
        sp.prefilter.min_diag_score = sabertooth::prefilter::size_adjusted_gate(sp.prefilter.min_diag_score, db_residues);
    }
}

/// `W,LO,HI` → SEG-like mask parameters.
fn parse_mask_params(v: &str) -> Result<MaskParams, String> {
    let p: Vec<&str> = v.split(',').collect();
    let bad = || format!("invalid mask parameters '{}' (expected WINDOW,LOCUT,HICUT, e.g. 12,2.2,2.5)", v);
    if p.len() != 3 {
        return Err(bad());
    }
    let window: usize = p[0].trim().parse().map_err(|_| bad())?;
    let locut: f64 = p[1].trim().parse().map_err(|_| bad())?;
    let hicut: f64 = p[2].trim().parse().map_err(|_| bad())?;
    if window < 2 || locut > hicut || locut < 0.0 {
        return Err(bad());
    }
    Ok(MaskParams { window, locut, hicut })
}

fn search_params(args: &Args, profile_mode: bool) -> Result<SearchParams, String> {
    let mut sp = if profile_mode {
        SearchParams::profile_default()
    } else {
        SearchParams::sequence_default()
    };
    let scale = if profile_mode { 8 } else { 2 };
    let pf = &mut sp.prefilter;
    // -s preset first; explicit --kmer-score / --max-kmers / --min-diag-score win
    if let Some(v) = args.opt_str("sensitivity") {
        let s: f64 = v.parse().map_err(|_| format!("invalid value for -s: '{}'", v))?;
        if !(1.0..=9.0).contains(&s) {
            return Err(format!("-s must be between 1 and 9, got {}", s));
        }
        if profile_mode {
            pf.kmer_score = sabertooth::prefilter::kmer_score_for_sensitivity(s, scale);
        } else {
            let (k, c, g) = sabertooth::prefilter::sequence_preset(s);
            pf.kmer_score = k;
            pf.max_kmers_per_pos = c;
            pf.min_diag_score = g;
        }
    }
    pf.kmer_score = args.opt_i32("kmer-score", pf.kmer_score)?;
    let default_seed = if profile_mode { SeedPattern::profile_default() } else { SeedPattern::sequence_default() };
    pf.k = requested_pattern(args, &default_seed)?.weight();
    pf.max_kmers_per_pos = args.opt_usize("max-kmers", pf.max_kmers_per_pos)?;
    pf.min_hits = args.opt_usize("min-hits", pf.min_hits)?;
    pf.min_diag_score = args.opt_i32("min-diag-score", pf.min_diag_score)?;
    pf.diag_hits = args.opt_usize("diag-hits", pf.diag_hits)?;
    pf.max_diagonals = args.opt_usize("prefilter-diags", pf.max_diagonals)?;
    pf.max_candidates = args.opt_usize("max-candidates", pf.max_candidates)?;
    pf.max_seed_occ = args.opt_usize("max-seed-occ", pf.max_seed_occ)?;
    pf.exact_seeds = args.opt_bool("exact-seeds", pf.exact_seeds)?;
    pf.hit_block = args.opt_usize("hit-block", pf.hit_block)?;
    pf.xdrop = args.opt_i32("ungapped-xdrop", pf.xdrop)?;
    if pf.xdrop < 0 {
        return Err("--ungapped-xdrop must be >= 0 (0 = whole-diagonal scoring)".into());
    }
    pf.diag_window = args.opt_usize("diag-window", pf.diag_window)?;
    let mask = args.opt_bool("mask", true)?;
    pf.mask = if mask { Some(MaskParams::default()) } else { None };
    sp.max_evalue = args.opt_f64("evalue", sp.max_evalue)?;
    let qcov = args.opt_f64("min-cov", 0.0)?;
    sp.min_query_cov = args.opt_f64("min-qcov", qcov)?;
    sp.min_target_cov = args.opt_f64("min-tcov", 0.0)?;
    sp.max_hits = args.opt_usize("max-hits", sp.max_hits)?;
    sp.band_radius = match args.opt_str("band") {
        Some(_) => Some(args.opt_usize("band", 32)?),
        None => None,
    };
    sp.mask_query = args.opt_bool("mask-query", mask)?;
    let mask_scoring = args.opt_bool("mask-scoring", sp.mask_scoring.is_some())? && !profile_mode;
    sp.mask_scoring = if mask_scoring {
        Some(match args.opt_str("mask-scoring-params") {
            Some(v) if !v.is_empty() => parse_mask_params(v)?,
            _ => MaskParams::scoring_default(),
        })
    } else {
        None
    };
    if !profile_mode {
        sp.comp_bias_corr = args.opt_bool("comp-bias-corr", true)?;
        sp.comp_bias_penalize_only = match args.opt_str("comp-bias-mode") {
            None | Some("full") => false,
            Some("neg") => true,
            Some(v) => return Err(format!("invalid --comp-bias-mode '{}' (full|neg)", v)),
        };
    } else if args.has("comp-bias-corr") {
        eprintln!("sabertooth: note: --comp-bias-corr applies to sequence queries; profiles are calibrated per PSSM");
    }
    sp.window = WindowParams {
        enabled: args.opt_bool("window", false)?,
        min_target_len: args.opt_usize("window-min-len", WindowParams::default().min_target_len)?,
        margin: args.opt_usize("window-margin", WindowParams::default().margin)?,
    };
    sp.profile_calibration_samples = args.opt_usize("profile-calibration-samples", sp.profile_calibration_samples)?;
    sp.calibration_seed = args.opt_usize("calibration-seed", sp.calibration_seed as usize)? as u64;
    sp.validate()?;
    Ok(sp)
}

// ---------------------------------------------------------------------------
// target resolution
// ---------------------------------------------------------------------------

struct TargetHandle {
    index: SeedIndex,
    how: String,
    build: Option<BuildStats>,
    seconds: f64,
}

fn tmp_dir(args: &Args) -> std::path::PathBuf {
    args.opt_str("tmp-dir")
        .filter(|s| !s.is_empty())
        .map(std::path::PathBuf::from)
        .unwrap_or_else(std::env::temp_dir)
}

/// Find or build the seed index for a target argument, syncing
/// `sp.prefilter.k`/`mask` with the index used.
///
/// The search always uses the command's own seed pattern and masking (its
/// defaults or the explicit `--k` / `--spaced` / `--mask`): an index built
/// with different settings is not adopted (a temporary index is built and a
/// warning says why), so results never depend on which indexes exist.
fn resolve_target(path: &str, args: &Args, sp: &mut SearchParams, plan: &Plan, default_seed: &SeedPattern, quiet: bool) -> Result<TargetHandle, String> {
    let t0 = Instant::now();
    let want_pattern = requested_pattern(args, default_seed)?;
    let want_mask = sp.prefilter.mask;
    let mut opts = BuildOptions::new(want_pattern.clone());
    opts.mask = want_mask;
    opts.memory_budget = plan.build_budget;
    opts.tmp_dir = tmp_dir(args);
    let require = args.has("require-index");
    let adopt = |idx: &SeedIndex, sp: &mut SearchParams| -> Result<(), Unusable> {
        idx.usable_for(&want_pattern, want_mask, None)?;
        sp.prefilter.k = want_pattern.weight();
        sp.prefilter.mask = want_mask;
        Ok(())
    };
    let hint = |e: &Unusable| -> &'static str {
        match e {
            Unusable::Incompatible(_) => " (pass the index's --k / --spaced / --mask to use it)",
            Unusable::Stale(_) => " (rebuild it with `sabertooth createindex`)",
        }
    };
    // 1) the argument is an index file
    if index::is_index_file(path) {
        let idx = SeedIndex::open(std::path::Path::new(path)).map_err(|e| e.to_string())?;
        match adopt(&idx, sp) {
            Ok(()) => {
                return Ok(TargetHandle {
                    how: format!("persistent index {}", path),
                    index: idx,
                    build: None,
                    seconds: t0.elapsed().as_secs_f64(),
                })
            }
            Err(e) => {
                if require {
                    return Err(format!("{}: {}", path, e));
                }
                if !quiet {
                    eprintln!("sabertooth: warning: {}: {}; re-indexing its sequences in a temporary index{}", path, e, hint(&e));
                }
                let (tidx, st) = SeedIndex::build_transient_from(&mut index::IndexSource(&idx), &opts).map_err(|e| e.to_string())?;
                sp.prefilter.k = want_pattern.weight();
                sp.prefilter.mask = want_mask;
                return Ok(TargetHandle {
                    index: tidx,
                    how: "temporary index (re-indexed)".into(),
                    build: Some(st),
                    seconds: t0.elapsed().as_secs_f64(),
                });
            }
        }
    }
    // 2) a sibling persistent index
    let sib = index::default_index_path(path);
    if sib.is_file() {
        match SeedIndex::open(&sib) {
            Ok(idx) => {
                let fresh = idx.usable_for(idx.pattern(), if idx.meta().masked { Some(idx.meta().mask) } else { None }, Some(path));
                match fresh.and_then(|_| adopt(&idx, sp)) {
                    Ok(()) => {
                        return Ok(TargetHandle {
                            how: format!("persistent index {}", sib.display()),
                            index: idx,
                            build: None,
                            seconds: t0.elapsed().as_secs_f64(),
                        })
                    }
                    Err(e) => {
                        if require {
                            return Err(format!("{}: {}", sib.display(), e));
                        }
                        if !quiet {
                            eprintln!("sabertooth: warning: {}: {}; using a temporary index{}", sib.display(), e, hint(&e));
                        }
                    }
                }
            }
            Err(e) => {
                if require {
                    return Err(e.to_string());
                }
                if !quiet {
                    eprintln!("sabertooth: warning: ignoring unreadable index {}: {}", sib.display(), e);
                }
            }
        }
    } else if require {
        return Err(format!(
            "no index for '{}' (expected {}); build one with `sabertooth createindex {}`",
            path,
            sib.display(),
            path
        ));
    }
    // 3) temporary index, streamed from the source
    if !quiet && !sib.is_file() {
        eprintln!(
            "sabertooth: no persistent seed index for '{}'; building a temporary one (keep one with `sabertooth createindex {}`)",
            path, path
        );
    }
    let mut src = index::source_for(path).map_err(|e| format!("cannot read target '{}': {}", path, e))?;
    let (idx, st) = SeedIndex::build_transient_from(src.as_mut(), &opts).map_err(|e| format!("cannot read target '{}': {}", path, e))?;
    sp.prefilter.k = want_pattern.weight();
    sp.prefilter.mask = want_mask;
    Ok(TargetHandle {
        index: idx,
        how: "temporary index".into(),
        build: Some(st),
        seconds: t0.elapsed().as_secs_f64(),
    })
}

fn build_stats_json(b: &BuildStats) -> Json {
    Json::new()
        .num("sequences", b.num_seqs)
        .num("residues", b.total_residues)
        .num("postings", b.num_postings)
        .num("distinct_kmers", b.num_kmers)
        .num("masked_residues", b.masked_residues)
        .num("partitions", b.partitions)
        .num("partition_buffer_bytes", b.partition_buffer_bytes)
        .num("bucket_table_bytes", b.bucket_table_bytes)
        .str("layout", b.layout)
        .num("posting_width_bytes", b.posting_width)
        .num("lookup_bytes", b.lookup_bytes)
        .num("file_bytes", b.file_bytes)
        .float("seconds_sequences", b.seconds_sequences)
        .float("seconds_count", b.seconds_count)
        .float("seconds_postings", b.seconds_postings)
        .float("seconds_finish", b.seconds_finish)
}

fn index_json(t: &TargetHandle) -> Json {
    let m = t.index.meta();
    let mut j = Json::new()
        .str("source", &t.how)
        .str("seed_pattern", &m.pattern.mask_string())
        .bool("masked", m.masked)
        .num("sequences", m.num_seqs)
        .num("residues", m.total_residues)
        .num("postings", m.num_postings)
        .num("distinct_kmers", m.num_kmers)
        .num("mapped_bytes", t.index.mapped_bytes())
        .float("seconds_open_or_build", t.seconds);
    if let Some(b) = &t.build {
        j = j.obj("build", build_stats_json(b));
    }
    j
}

fn stats_json(st: &SearchStats) -> Json {
    let p = &st.prefilter;
    let s = |ns: u64| ns as f64 / 1e9;
    Json::new()
        .num("queries", st.queries)
        .num("query_residues", st.query_residues)
        .num("windows", p.windows)
        .num("windows_skipped_mask_or_x", p.windows_skipped)
        .num("kmers_generated", p.kmers_generated)
        .num("windows_at_kmer_cap", p.windows_at_cap)
        .num("exact_seeds_added", p.exact_seeds_added)
        .num("index_lookups", p.lookups)
        .num("index_prefix_probes", p.prefix_probes)
        .num("presence_bitmap_rejects", p.presence_rejects)
        .num("seeds_with_postings", p.seeds)
        .num("postings_visited", p.postings_visited)
        .num("seeds_over_occ_cap", p.seeds_over_occ_cap)
        .num("postings_skipped_occ_cap", p.postings_skipped_occ)
        .num("target_blocks", p.blocks)
        .num("hit_buffer_peak_entries", p.hit_buffer_peak)
        .num("hit_buffer_peak_bytes", p.hit_buffer_peak * 16)
        .num("hit_buffer_overruns", p.hit_buffer_overruns)
        .num("targets_hit", p.targets_hit)
        .num("targets_evaluated", p.targets_evaluated)
        .num("diagonals_scored", p.diagonals_scored)
        .num("xdrop_extensions", p.xdrop_extensions)
        .num("chains_used", p.chains_used)
        .num("candidates", st.candidates)
        .num("candidates_capped", p.candidates_capped)
        .num("gapped_scored", st.gapped_scored)
        .num("gapped_cells", st.gapped_cells)
        .num("scalar_fallbacks", st.scalar_fallbacks)
        .num("windowed", st.windowed)
        .num("window_expansions", st.window_expansions)
        .num("passed_score", st.passed_score)
        .num("tracebacks", st.tracebacks)
        .num("rejected_coverage", st.rejected_coverage)
        .num("hits", st.hits)
        .num("score_mismatches", st.score_mismatches)
        .num("profile_calibrations", st.profile_calibrations)
        .num("masked_scoring_positions", st.masked_scoring_positions)
        .float("cpu_s_setup", s(st.ns_setup))
        .float("cpu_s_prefilter", s(st.ns_prefilter))
        .float("cpu_s_seed_generation", s(p.ns_seeds))
        .float("cpu_s_hit_collection_sort", s(p.ns_hits))
        .float("cpu_s_diagonal_evidence", s(p.ns_evidence))
        .float("cpu_s_gapped_scoring", s(st.ns_gapped))
        .float("cpu_s_traceback", s(st.ns_traceback))
}

struct Report<'a> {
    command: &'a str,
    plan: &'a Plan,
    target: Option<&'a TargetHandle>,
    stats: &'a SearchStats,
    ev_model: String,
    wall_load: f64,
    wall_search: f64,
    wall_output: f64,
    wall_total: f64,
    /// Out-of-core / distributed shards used (None: the whole index at once).
    shards: Option<&'a sabertooth::shard::ShardSet>,
    /// Effective prefilter settings (None: no prefilter, e.g. easysearch).
    prefilter: Option<sabertooth::prefilter::PrefilterParams>,
}

fn prefilter_json(p: &sabertooth::prefilter::PrefilterParams) -> Json {
    Json::new()
        .num("k", p.k)
        .num("kmer_score", p.kmer_score)
        .num("max_kmers_per_window", p.max_kmers_per_pos)
        .num("min_diag_score", p.min_diag_score)
        .num("ungapped_xdrop", p.xdrop)
        .num("diag_window", p.diag_window)
        .num("max_diagonals", p.max_diagonals)
}

fn shards_json(set: &sabertooth::shard::ShardSet) -> Json {
    Json::new()
        .num("count", set.shards.len())
        .str("dir", &set.dir.display().to_string())
        .bool("reused", set.reused)
        .float("seconds_prepare", set.seconds)
        .num("largest_shard_residues", set.shards.iter().map(|s| s.residues).max().unwrap_or(0))
}

fn emit_report(args: &Args, r: &Report) -> Result<(), String> {
    let ps = proc_stats();
    let st = r.stats;
    let nq = st.queries.max(1) as f64;
    if args.has("stats") {
        let p = &st.prefilter;
        eprintln!("--- {} stats ---", r.command);
        eprintln!("  threads {} | memory plan: {}", r.plan.threads, r.plan.budget.map(|b| format!("{} MiB", b >> 20)).unwrap_or("unlimited".into()));
        for n in &r.plan.notes {
            eprintln!("  note: {}", n);
        }
        if let Some(t) = r.target {
            eprintln!(
                "  target: {} ({} seqs, {} residues, seed {}, masked {}) in {:.2}s",
                t.how,
                t.index.len(),
                t.index.total_residues(),
                t.index.pattern().mask_string(),
                t.index.meta().masked,
                t.seconds
            );
        }
        if let Some(set) = r.shards {
            eprintln!(
                "  shards: {} in {} ({} in {:.2}s)",
                set.shards.len(),
                set.dir.display(),
                if set.reused { "reused" } else { "built" },
                set.seconds
            );
        }
        if let Some(pf) = &r.prefilter {
            eprintln!(
                "  prefilter: k-mer score {}  max k-mers/window {}  ungapped gate {}  evidence {}",
                pf.kmer_score,
                pf.max_kmers_per_pos,
                pf.min_diag_score,
                if pf.xdrop > 0 { format!("X-drop {}", pf.xdrop) } else if pf.diag_window > 0 { format!("segment within {} of seeds", pf.diag_window) } else { "whole diagonal".to_string() }
            );
        }
        eprintln!("  statistics: {}", r.ev_model);
        eprintln!(
            "  wall: load {:.2}s  search {:.2}s  output {:.2}s  total {:.2}s  ({:.1} queries/s)",
            r.wall_load,
            r.wall_search,
            r.wall_output,
            r.wall_total,
            st.queries as f64 / r.wall_search.max(1e-9)
        );
        eprintln!(
            "  per query: postings {:.0}  candidates {:.1}  gapped {:.1}  tracebacks {:.1}  hits {:.1}",
            p.postings_visited as f64 / nq,
            st.candidates as f64 / nq,
            st.gapped_scored as f64 / nq,
            st.tracebacks as f64 / nq,
            st.hits as f64 / nq
        );
        eprintln!(
            "  cpu s: seeds {:.2}  hits {:.2}  evidence {:.2}  gapped {:.2}  traceback {:.2}",
            p.ns_seeds as f64 / 1e9,
            p.ns_hits as f64 / 1e9,
            p.ns_evidence as f64 / 1e9,
            st.ns_gapped as f64 / 1e9,
            st.ns_traceback as f64 / 1e9
        );
        eprintln!(
            "  caps: kmer-cap windows {}  occ-capped seeds {}  capped candidates {}  hit-buffer peak {} ({} overruns)",
            p.windows_at_cap, p.seeds_over_occ_cap, p.candidates_capped, p.hit_buffer_peak, p.hit_buffer_overruns
        );
        eprintln!(
            "  memory: peak RSS {:.1} MiB  page faults minor {} major {}",
            ps.peak_rss as f64 / 1048576.0,
            ps.minflt,
            ps.majflt
        );
        if st.score_mismatches > 0 {
            eprintln!("  WARNING: {} traceback/SIMD score mismatches", st.score_mismatches);
        }
    }
    if let Some(path) = args.opt_str("stats-json") {
        let mut j = Json::new()
            .str("command", r.command)
            .str("version", VERSION)
            .obj("plan", plan_json(r.plan))
            .str("statistics_model", &r.ev_model)
            .obj(
                "wall_s",
                Json::new()
                    .float("load", r.wall_load)
                    .float("search", r.wall_search)
                    .float("output", r.wall_output)
                    .float("total", r.wall_total),
            )
            .float("queries_per_s", st.queries as f64 / r.wall_search.max(1e-9))
            .obj(
                "per_query",
                Json::new()
                    .float("postings_visited", st.prefilter.postings_visited as f64 / nq)
                    .float("candidates", st.candidates as f64 / nq)
                    .float("gapped_alignments", st.gapped_scored as f64 / nq)
                    .float("tracebacks", st.tracebacks as f64 / nq)
                    .float("hits", st.hits as f64 / nq),
            )
            .obj("counters", stats_json(st))
            .obj(
                "memory",
                Json::new()
                    .num("peak_rss_bytes", ps.peak_rss)
                    .num("seed_hit_buffer_peak_bytes", st.prefilter.hit_buffer_peak * 16)
                    .num("minor_page_faults", ps.minflt)
                    .num("major_page_faults", ps.majflt)
                    .float("user_cpu_s", ps.utime)
                    .float("system_cpu_s", ps.stime),
            );
        if let Some(t) = r.target {
            j = j.obj("index", index_json(t));
        }
        if let Some(set) = r.shards {
            j = j.obj("shards", shards_json(set));
        }
        if let Some(pf) = &r.prefilter {
            j = j.obj("prefilter", prefilter_json(pf));
        }
        std::fs::write(path, j.finish() + "\n").map_err(|e| format!("cannot write --stats-json '{}': {}", path, e))?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// search
// ---------------------------------------------------------------------------

fn run_search(args: &Args) -> Result<(), String> {
    if args.positionals.len() < 2 {
        return Err("usage: sabertooth search <query> <target> [options]".into());
    }
    args.check_known("search", SEARCH_OPTS)?;
    let m = SubstitutionMatrix::blosum62();
    let mut sp = search_params(args, false)?;

    #[cfg(feature = "distributed")]
    if args.has("distributed") {
        return run_search_distributed(args, &m, sp);
    }
    #[cfg(not(feature = "distributed"))]
    if args.has("distributed") {
        return Err("this build has no distributed backend; rebuild with `--features distributed`".into());
    }

    let plan = make_plan(args, &mut sp)?;
    let t_total = Instant::now();
    let target = resolve_target(&args.positionals[1], args, &mut sp, &plan, &SeedPattern::sequence_default(), false)?;
    apply_db_size_gate(args, &mut sp, target.index.total_residues());
    let ev = EValueParams::sequence(sp.gaps, sp.comp_bias_corr, target.index.total_residues(), target.index.len())?;
    let nshards = local_shard_count(args, &plan, &target)?;
    if nshards > 1 {
        return run_search_sharded(args, &m, sp, &plan, &target, &ev, nshards, t_total);
    }
    let eng = Engine {
        mat: &m,
        targets: &target.index,
        index: Some(&target.index),
        params: sp,
        stats: Statistics::Fixed(ev.clone()),
    };
    let extra = args.has("extra-cols");
    let mut out = open_output(args)?;
    write_header(&mut out, extra)?;
    let mut qs = QueryStream::open(&args.positionals[0], SeqKind::Protein)?;
    let mut stats = SearchStats::default();
    let (mut wall_search, mut wall_output, mut wall_load) = (0f64, 0f64, target.seconds);
    let mut ordinal = 0usize;
    let mut nhits = 0usize;
    loop {
        let tl = Instant::now();
        let batch = match qs.next_batch(plan.batch_seqs, plan.batch_residues, &m.aa2num)? {
            Some(b) => b,
            None => break,
        };
        wall_load += tl.elapsed().as_secs_f64();
        let ts = Instant::now();
        let res = eng.search_batch(&batch, ordinal, &mut stats);
        wall_search += ts.elapsed().as_secs_f64();
        let to = Instant::now();
        for q in &res {
            let hits: Vec<Hit> = q.iter().map(|(h, _)| h.clone()).collect();
            nhits += hits.len();
            write_hits(&mut out, &hits, extra)?;
        }
        wall_output += to.elapsed().as_secs_f64();
        ordinal += batch.len();
    }
    out.flush().map_err(|e| format!("write failed: {}", e))?;
    eprintln!(
        "search: {} queries x {} targets, {} hits in {:.2?}",
        ordinal,
        target.index.len(),
        nhits,
        t_total.elapsed()
    );
    emit_report(
        args,
        &Report {
            command: "search",
            plan: &plan,
            target: Some(&target),
            stats: &stats,
            ev_model: format!("{} (lambda {:.4}, K {:.4})", ev.model, ev.lambda, ev.k),
            wall_load,
            wall_search,
            wall_output,
            wall_total: t_total.elapsed().as_secs_f64(),
            shards: None,
            prefilter: Some(sp.prefilter),
        },
    )
}

/// Shards for an out-of-core local search: `--shards N` (1 = never shard), or,
/// with `--memory-limit`, enough residue-balanced shards for one shard plus
/// the planned working memory to fit in the budget (1 when the whole index
/// fits).
fn local_shard_count(args: &Args, plan: &Plan, target: &TargetHandle) -> Result<usize, String> {
    let max = target.index.len().clamp(1, 4096);
    if args.has("shards") {
        return Ok(args.opt_usize("shards", 1)?.clamp(1, max));
    }
    let Some(budget) = plan.budget else { return Ok(1) };
    // working memory next to the mapped index: fixed overhead, every worker's
    // scratch and the k-mer presence bitmap (built per opened index)
    let w = target.index.pattern().weight().min(6) as u32;
    let presence = 20usize.pow(w) / 8;
    let working = FIXED_OVERHEAD + plan.threads * plan.per_worker + presence;
    Ok(sabertooth::shard::shards_for_budget(budget, target.index.mapped_bytes(), working, max.min(256)))
}

/// Build (or reuse) the shards of `target` for `n` shards.
fn shard_set_for(args: &Args, plan: &Plan, target: &TargetHandle, n: usize) -> Result<sabertooth::shard::ShardSet, String> {
    use sabertooth::shard::{default_shard_dir, parent_identity, prepare_shards};
    let tpath = &args.positionals[1];
    let identity = parent_identity(&target.index, tpath).map_err(|e| format!("cannot identify '{}': {}", tpath, e))?;
    let dir = match args.opt_str("shard-dir").filter(|s| !s.is_empty()) {
        Some(d) => std::path::PathBuf::from(d),
        None => default_shard_dir(&target.index, &identity, n, &tmp_dir(args)),
    };
    let mut opts = BuildOptions::new(target.index.pattern().clone());
    opts.memory_budget = plan.build_budget;
    opts.tmp_dir = tmp_dir(args);
    prepare_shards(&target.index, &identity, n, &dir, &opts).map_err(|e| format!("cannot prepare shards in {}: {}", dir.display(), e))
}

/// Out-of-core search: every query batch is searched against one shard at a
/// time (only that shard's index pages need to be resident) and the hits are
/// merged per query — output identical to the unsharded search.
#[allow(clippy::too_many_arguments)]
fn run_search_sharded(
    args: &Args,
    m: &SubstitutionMatrix,
    sp: SearchParams,
    plan: &Plan,
    target: &TargetHandle,
    ev: &EValueParams,
    n: usize,
    t_total: Instant,
) -> Result<(), String> {
    let set = shard_set_for(args, plan, target, n)?;
    eprintln!(
        "search: index of {} MiB searched as {} residue-balanced shard(s), one at a time ({} in {} in {:.2}s){}",
        target.index.mapped_bytes() >> 20,
        set.shards.len(),
        if set.reused { "reused" } else { "built" },
        set.dir.display(),
        set.seconds,
        if args.has("shards") { "" } else { "; --memory-limit is smaller than the index" }
    );
    let extra = args.has("extra-cols");
    let mut out = open_output(args)?;
    write_header(&mut out, extra)?;
    let mut qs = QueryStream::open(&args.positionals[0], SeqKind::Protein)?;
    let mut stats = SearchStats::default();
    let (mut wall_search, mut wall_output, mut wall_load) = (0f64, 0f64, target.seconds + set.seconds);
    let mut ordinal = 0usize;
    let mut nhits = 0usize;
    loop {
        let tl = Instant::now();
        let batch = match qs.next_batch(plan.batch_seqs, plan.batch_residues, &m.aa2num)? {
            Some(b) => b,
            None => break,
        };
        wall_load += tl.elapsed().as_secs_f64();
        let ts = Instant::now();
        let res = sabertooth::shard::search_batch_sharded(&set, m, sp, ev, &batch, ordinal, &mut stats)?;
        wall_search += ts.elapsed().as_secs_f64();
        let to = Instant::now();
        for q in &res {
            nhits += q.len();
            write_hits(&mut out, q, extra)?;
        }
        wall_output += to.elapsed().as_secs_f64();
        ordinal += batch.len();
    }
    out.flush().map_err(|e| format!("write failed: {}", e))?;
    eprintln!(
        "search: {} queries x {} targets ({} shards), {} hits in {:.2?}",
        ordinal,
        target.index.len(),
        set.shards.len(),
        nhits,
        t_total.elapsed()
    );
    emit_report(
        args,
        &Report {
            command: "search",
            plan,
            target: Some(target),
            stats: &stats,
            ev_model: format!("{} (lambda {:.4}, K {:.4}; whole-database size)", ev.model, ev.lambda, ev.k),
            wall_load,
            wall_search,
            wall_output,
            wall_total: t_total.elapsed().as_secs_f64(),
            shards: Some(&set),
            prefilter: Some(sp.prefilter),
        },
    )
}

#[cfg(feature = "distributed")]
fn run_search_distributed(args: &Args, m: &SubstitutionMatrix, mut sp: SearchParams) -> Result<(), String> {
    use sabertooth::distributed::DistributedSearcher;
    let shards = args.opt_usize("shards", 4)?.max(1);
    let gpus_per_shard = args.opt_usize("gpus", 0)?;
    let plan = make_plan(args, &mut sp)?;
    let t_total = Instant::now();
    // the same target resolution (persistent or temporary index) and the same
    // validated parameters as a local search
    let target = resolve_target(&args.positionals[1], args, &mut sp, &plan, &SeedPattern::sequence_default(), false)?;
    apply_db_size_gate(args, &mut sp, target.index.total_residues());
    let set = shard_set_for(args, &plan, &target, shards)?;
    eprintln!(
        "distributed: {} residue-balanced shard(s) in {} ({} in {:.2}s)",
        set.shards.len(),
        set.dir.display(),
        if set.reused { "reused" } else { "built" },
        set.seconds
    );
    let searcher = DistributedSearcher::start(&set, sp, gpus_per_shard, sabertooth::distributed::Config::from_args())?;
    let extra = args.has("extra-cols");
    let mut out = open_output(args)?;
    write_header(&mut out, extra)?;
    let mut qs = QueryStream::open(&args.positionals[0], SeqKind::Protein)?;
    let mut stats = SearchStats::default();
    let (mut wall_search, mut wall_output, mut wall_load) = (0f64, 0f64, target.seconds + set.seconds);
    let mut ordinal = 0usize;
    let mut nhits = 0usize;
    let mut gpu_devices: Vec<String> = Vec::new();
    let res: Result<(), String> = (|| {
        loop {
            let tl = Instant::now();
            let batch = match qs.next_batch(plan.batch_seqs, plan.batch_residues, &m.aa2num)? {
                Some(b) => b,
                None => break,
            };
            wall_load += tl.elapsed().as_secs_f64();
            let ts = Instant::now();
            let r = searcher.search_batch(&batch, ordinal)?;
            wall_search += ts.elapsed().as_secs_f64();
            stats.add(&r.stats);
            for d in r.gpu_devices {
                if !gpu_devices.contains(&d) {
                    gpu_devices.push(d);
                }
            }
            let to = Instant::now();
            for q in &r.hits {
                nhits += q.len();
                write_hits(&mut out, q, extra)?;
            }
            wall_output += to.elapsed().as_secs_f64();
            ordinal += batch.len();
        }
        Ok(())
    })();
    searcher.shutdown();
    res?;
    out.flush().map_err(|e| format!("write failed: {}", e))?;
    eprintln!(
        "search (HydraMPP, {} shards): {} queries x {} targets, {} hits in {:.2?}",
        set.shards.len(),
        ordinal,
        target.index.len(),
        nhits,
        t_total.elapsed()
    );
    if gpus_per_shard > 0 {
        if gpu_devices.is_empty() {
            eprintln!("  no GPU devices were available/reserved (run with `--hydra-gpus N` to simulate)");
        } else {
            eprintln!("  GPU devices reserved/pinned: {}", gpu_devices.join(" | "));
        }
        eprintln!("  note: --gpus is GPU-aware scheduling only; all alignment ran on the CPU (no GPU kernel exists)");
    }
    let ev = EValueParams::sequence(sp.gaps, sp.comp_bias_corr, target.index.total_residues(), target.index.len())?;
    emit_report(
        args,
        &Report {
            command: "search --distributed",
            plan: &plan,
            target: Some(&target),
            stats: &stats,
            ev_model: format!("{} (lambda {:.4}, K {:.4}; global database size)", ev.model, ev.lambda, ev.k),
            wall_load,
            wall_search,
            wall_output,
            wall_total: t_total.elapsed().as_secs_f64(),
            shards: Some(&set),
            prefilter: Some(sp.prefilter),
        },
    )
}

fn run_profilesearch(args: &Args) -> Result<(), String> {
    if args.positionals.len() < 2 {
        return Err("usage: sabertooth profilesearch <msa> <target> [options]".into());
    }
    let mut known = SEARCH_OPTS.to_vec();
    known.extend_from_slice(&["pca", "pcb", "bias", "num-iterations", "inclusion-evalue", "profile-calibration-samples", "calibration-seed"]);
    args.check_known("profilesearch", &known)?;
    let m = SubstitutionMatrix::blosum62();
    let msa_data = read_to_string(&args.positionals[0]).map_err(|e| e.to_string())?;
    let pc = pc_params(args)?;
    let prof = Profile::from_msa_str(&msa_data, &m, pc)?;
    let query_id = std::path::Path::new(&args.positionals[0])
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("profile")
        .to_string();
    let mut sp = search_params(args, true)?;
    let plan = make_plan(args, &mut sp)?;
    let t_total = Instant::now();
    let target = resolve_target(&args.positionals[1], args, &mut sp, &plan, &SeedPattern::profile_default(), false)?;
    let db_size = (target.index.total_residues(), target.index.len());
    let num_iter = args.opt_usize("num-iterations", 1)?;
    let mut stats = SearchStats::default();
    let ts = Instant::now();
    let hits = if num_iter > 1 {
        let (msa, query_len) = sabertooth::profile::parse_msa(&msa_data, &m.aa2num)?;
        let incl = args.opt_f64("inclusion-evalue", 1e-3)?;
        eprintln!("profilesearch: iterative ({} rounds, inclusion E <= {:.0e})", num_iter, incl);
        sabertooth::search::iterate_profile(&msa, query_len, &query_id, &target.index, &m, pc, sp, db_size, None, num_iter, incl, Some(&mut stats))?
    } else {
        sabertooth::search::profile_hits(&prof, &query_id, &target.index, &m, sp, db_size, None, Some(&mut stats))
    };
    let wall_search = ts.elapsed().as_secs_f64();
    let extra = args.has("extra-cols");
    let mut out = open_output(args)?;
    write_header(&mut out, extra)?;
    write_hits(&mut out, &hits, extra)?;
    out.flush().map_err(|e| e.to_string())?;
    eprintln!(
        "profilesearch: profile ({} cols, {} seqs) x {} targets, {} hits in {:.2?}",
        prof.query_len,
        prof.set_size,
        target.index.len(),
        hits.len(),
        t_total.elapsed()
    );
    emit_report(
        args,
        &Report {
            command: "profilesearch",
            plan: &plan,
            target: Some(&target),
            stats: &stats,
            ev_model: "per-profile calibrated gapped Gumbel".into(),
            wall_load: target.seconds,
            wall_search,
            wall_output: 0.0,
            wall_total: t_total.elapsed().as_secs_f64(),
            shards: None,
            prefilter: Some(sp.prefilter),
        },
    )
}

fn run_easysearch(args: &Args) -> Result<(), String> {
    if args.positionals.len() < 2 {
        return Err("usage: sabertooth easysearch <query> <target> [options]".into());
    }
    args.check_known("easysearch", SEARCH_OPTS)?;
    let m = SubstitutionMatrix::blosum62();
    let mut sp = search_params(args, false)?;
    sp.exhaustive = true;
    let plan = make_plan(args, &mut sp)?;
    let t_total = Instant::now();
    let target = load_seqs(&args.positionals[1], LoadOptions::SEARCH, SeqKind::Protein, "target")?;
    let wall_load = t_total.elapsed().as_secs_f64();
    let ev = EValueParams::sequence(sp.gaps, sp.comp_bias_corr, target.total_residues, target.len())?;
    let eng = Engine {
        mat: &m,
        targets: &target,
        index: None,
        params: sp,
        stats: Statistics::Fixed(ev.clone()),
    };
    let extra = args.has("extra-cols");
    let mut out = open_output(args)?;
    write_header(&mut out, extra)?;
    let mut qs = QueryStream::open(&args.positionals[0], SeqKind::Protein)?;
    let mut stats = SearchStats::default();
    let mut ordinal = 0usize;
    let mut nhits = 0;
    let ts = Instant::now();
    while let Some(batch) = qs.next_batch(plan.batch_seqs, plan.batch_residues, &m.aa2num)? {
        let res = eng.search_batch(&batch, ordinal, &mut stats);
        for q in &res {
            let hits: Vec<Hit> = q.iter().map(|(h, _)| h.clone()).collect();
            nhits += hits.len();
            write_hits(&mut out, &hits, extra)?;
        }
        ordinal += batch.len();
    }
    let wall_search = ts.elapsed().as_secs_f64();
    out.flush().map_err(|e| e.to_string())?;
    eprintln!(
        "easysearch (database-free, no prefilter index): {} queries x {} targets = {} full SW alignments, {} hits in {:.2?}",
        ordinal,
        target.len(),
        ordinal * target.len(),
        nhits,
        t_total.elapsed()
    );
    emit_report(
        args,
        &Report {
            command: "easysearch",
            plan: &plan,
            target: None,
            stats: &stats,
            ev_model: format!("{} (lambda {:.4}, K {:.4})", ev.model, ev.lambda, ev.k),
            wall_load,
            wall_search,
            wall_output: 0.0,
            wall_total: t_total.elapsed().as_secs_f64(),
            shards: None,
            prefilter: None,
        },
    )
}

// ---------------------------------------------------------------------------
// createdb / msa2profile / align
// ---------------------------------------------------------------------------

fn run_createdb(args: &Args) -> Result<(), String> {
    args.check_known("createdb", &[])?;
    if args.positionals.len() < 2 {
        return Err("usage: sabertooth createdb <in> <out.fasta>".into());
    }
    let db = load_seqs(&args.positionals[0], LoadOptions::FULL, SeqKind::Any, "input")?;
    db.write_fasta(&args.positionals[1]).map_err(|e| e.to_string())?;
    println!(
        "createdb: {} sequences, {} residues -> {}",
        db.len(),
        db.total_residues,
        args.positionals[1]
    );
    Ok(())
}

fn run_msa2profile(args: &Args) -> Result<(), String> {
    args.check_known("msa2profile", &["pca", "pcb", "bias"])?;
    if args.positionals.is_empty() {
        return Err("usage: sabertooth msa2profile <msa> [--out file]".into());
    }
    let m = SubstitutionMatrix::blosum62();
    let data = read_to_string(&args.positionals[0]).map_err(|e| e.to_string())?;
    let prof = Profile::from_msa_str(&data, &m, pc_params(args)?)?;
    let table = prof.pssm_table();
    write_out(args.opt_str("out"), &table)?;
    eprintln!("msa2profile: {} sequences, {} match columns", prof.set_size, prof.query_len);
    Ok(())
}

fn run_align(args: &Args) -> Result<(), String> {
    args.check_known("align", &[])?;
    if args.positionals.len() < 2 {
        return Err("usage: sabertooth align <query> <target>".into());
    }
    let m = SubstitutionMatrix::blosum62();
    let query = load_seqs(&args.positionals[0], LoadOptions::SEARCH, SeqKind::Protein, "query")?;
    let target = load_seqs(&args.positionals[1], LoadOptions::SEARCH, SeqKind::Protein, "target")?;
    let qnum = query.num(0);
    let gaps = GapCosts::sequence_default();
    let ev = EValueParams::sequence(gaps, false, target.total_residues, target.len())?;
    let sc = SeqScorer { query: qnum, mat: &m };
    let mut out = String::new();
    for ti in 0..target.len() {
        if let Some(aln) = local_align(&sc, target.num(ti), gaps) {
            out.push_str(&format!(
                "{} vs {}: score={} bits={:.1} evalue={:.2e} id={:.1}% len={} cigar={}\n  Q[{}-{}] T[{}-{}]\n",
                query.id(0),
                target.id(ti),
                aln.raw_score,
                ev.bit_score(aln.raw_score),
                ev.evalue(aln.raw_score, qnum.len()),
                aln.pct_identity(),
                aln.aln_len,
                aln.cigar,
                aln.query_start + 1,
                aln.query_end + 1,
                aln.target_start + 1,
                aln.target_end + 1,
            ));
        }
    }
    if out.is_empty() {
        out.push_str("(no local alignments with positive score)\n");
    }
    write_out(args.opt_str("out"), &out)?;
    Ok(())
}

// ---------------------------------------------------------------------------
// result2msa / translatesearch / profileprofile / taxonomy
// ---------------------------------------------------------------------------

/// Run a full in-memory sequence search through the engine (helper for the
/// commands that post-process hits).
#[allow(clippy::type_complexity)]
fn engine_search(args: &Args, query: &SeqDb, target_path: &str, sp: &mut SearchParams) -> Result<(Vec<Vec<(Hit, sabertooth::align::Alignment)>>, TargetHandle), String> {
    let m = SubstitutionMatrix::blosum62();
    let plan = make_plan(args, sp)?;
    let target = resolve_target(target_path, args, sp, &plan, &SeedPattern::sequence_default(), false)?;
    apply_db_size_gate(args, sp, target.index.total_residues());
    let ev = EValueParams::sequence(sp.gaps, sp.comp_bias_corr, target.index.total_residues(), target.index.len())?;
    let eng = Engine {
        mat: &m,
        targets: &target.index,
        index: Some(&target.index),
        params: *sp,
        stats: Statistics::Fixed(ev),
    };
    let mut stats = SearchStats::default();
    let mut all = Vec::with_capacity(query.len());
    let step = plan.batch_seqs.max(1);
    let mut i = 0;
    while i < query.len() {
        let j = (i + step).min(query.len());
        let b = query.subset(i..j);
        all.extend(eng.search_batch(&b, i, &mut stats));
        i = j;
    }
    drop(eng);
    Ok((all, target))
}

fn run_result2msa(args: &Args) -> Result<(), String> {
    if args.positionals.len() < 2 {
        return Err("usage: sabertooth result2msa <query> <target> [options]".into());
    }
    args.check_known("result2msa", SEARCH_OPTS)?;
    let query = load_seqs(&args.positionals[0], LoadOptions::SEARCH, SeqKind::Protein, "query")?;
    let mut sp = search_params(args, false)?;
    // the MSA is built from the plain-matrix alignment of each hit
    let start = Instant::now();
    let (res, target) = engine_search(args, &query, &args.positionals[1], &mut sp)?;
    let mut out = String::new();
    let mut total_rows = 0usize;
    for (qi, hits) in res.iter().enumerate() {
        let rows: Vec<(String, Vec<u8>, sabertooth::align::Alignment)> = hits
            .iter()
            .map(|(h, a)| (h.target_id.clone(), target.index.seq(h.target_index).to_vec(), a.clone()))
            .collect();
        total_rows += rows.len();
        out.push_str(&sabertooth::msa::build_a3m(query.id(qi), query.num(qi), &rows));
    }
    write_out(args.opt_str("out"), &out)?;
    eprintln!(
        "result2msa: {} queries, {} aligned rows total in {:.2?}",
        query.len(),
        total_rows,
        start.elapsed()
    );
    Ok(())
}

fn run_translatesearch(args: &Args) -> Result<(), String> {
    if args.positionals.len() < 2 {
        return Err("usage: sabertooth translatesearch <nucl_query> <protein_target> [options]".into());
    }
    let mut known = SEARCH_OPTS.to_vec();
    known.push("min-orf");
    args.check_known("translatesearch", &known)?;
    let m = SubstitutionMatrix::blosum62();
    let nt_query = load_seqs(&args.positionals[0], LoadOptions::LETTERS, SeqKind::Nucleotide, "nucleotide query")?;
    if !seqio::looks_nucleotide_db(&nt_query) {
        eprintln!("translatesearch: warning — query does not look like nucleotide sequence");
    }
    let min_aa = args.opt_usize("min-orf", 20)?;
    let (orf_db, _meta) = sabertooth::nucl::orf_query_db(&nt_query, &m.aa2num, min_aa);
    if orf_db.is_empty() {
        return Err(format!("no ORFs of length >= {} found in the query", min_aa));
    }
    let mut sp = search_params(args, false)?;
    let start = Instant::now();
    let (res, target) = engine_search(args, &orf_db, &args.positionals[1], &mut sp)?;
    let extra = args.has("extra-cols");
    let mut out = open_output(args)?;
    write_header(&mut out, extra)?;
    let mut n = 0;
    for q in &res {
        let hits: Vec<Hit> = q.iter().map(|(h, _)| h.clone()).collect();
        n += hits.len();
        write_hits(&mut out, &hits, extra)?;
    }
    out.flush().map_err(|e| e.to_string())?;
    eprintln!(
        "translatesearch: {} nucleotide queries -> {} ORFs (>= {} aa) x {} targets, {} hits in {:.2?}",
        nt_query.len(),
        orf_db.len(),
        min_aa,
        target.index.len(),
        n,
        start.elapsed()
    );
    eprintln!("  (query IDs encode source_frame_ntstart, e.g. contig1_f-2_57)");
    Ok(())
}

fn run_profileprofile(args: &Args) -> Result<(), String> {
    args.check_known("profileprofile", &["pca", "pcb", "bias"])?;
    if args.positionals.len() < 2 {
        return Err("usage: sabertooth profileprofile <query.msa> <target.msa> [more_target.msa ...]".into());
    }
    let m = SubstitutionMatrix::blosum62();
    let query_data = read_to_string(&args.positionals[0]).map_err(|e| e.to_string())?;
    let pc = pc_params(args)?;
    let qprof = Profile::from_msa_str(&query_data, &m, pc)?;
    let pback: Vec<f64> = (0..PROFILE_AA_SIZE).map(|i| m.p_back[i]).collect();
    let gaps = GapCosts::profile_default();
    let start = Instant::now();
    let mut results: Vec<(i32, String, sabertooth::align::Alignment)> = Vec::new();
    for path in &args.positionals[1..] {
        let data = read_to_string(path).map_err(|e| e.to_string())?;
        let tprof = Profile::from_msa_str(&data, &m, pc)?;
        let name = std::path::Path::new(path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or(path)
            .to_string();
        if let Some(aln) = sabertooth::align::align_profile_profile(&qprof, &tprof, gaps, &pback) {
            results.push((aln.raw_score, name, aln));
        }
    }
    results.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
    // profile–profile significance needs its own calibration; only raw scores
    // (and their ×8-scale bit equivalents) are reported. E-values are omitted.
    let mut out = String::new();
    out.push_str("target\tscore\tbits\tcons_id%\tq_start\tq_end\tt_start\tt_end\tcigar\n");
    for (score, name, aln) in &results {
        out.push_str(&format!(
            "{}\t{}\t{:.1}\t{:.1}\t{}\t{}\t{}\t{}\t{}\n",
            name,
            score,
            *score as f64 / 8.0,
            aln.pct_identity(),
            aln.query_start + 1,
            aln.query_end + 1,
            aln.target_start + 1,
            aln.target_end + 1,
            aln.cigar
        ));
    }
    write_out(args.opt_str("out"), &out)?;
    eprintln!(
        "profileprofile: query ({} cols) vs {} target profiles, {} aligned in {:.2?}",
        qprof.query_len,
        args.positionals.len() - 1,
        results.len(),
        start.elapsed()
    );
    Ok(())
}

fn run_taxonomy(args: &Args) -> Result<(), String> {
    if args.positionals.len() < 2 {
        return Err("usage: sabertooth taxonomy <query> <target> --nodes nodes.dmp --seqmap seqid2taxid.tsv [--names names.dmp] [--lca-frac f]".into());
    }
    let mut known = SEARCH_OPTS.to_vec();
    known.extend_from_slice(&["nodes", "names", "seqmap", "lca-frac"]);
    args.check_known("taxonomy", &known)?;
    let nodes_path = args.opt_str("nodes").ok_or("taxonomy requires --nodes <nodes.dmp>")?;
    let seqmap_path = args.opt_str("seqmap").ok_or("taxonomy requires --seqmap <seqid2taxid.tsv>")?;
    let query = load_seqs(&args.positionals[0], LoadOptions::SEARCH, SeqKind::Protein, "query")?;
    let mut sp = search_params(args, false)?;
    let (res, _target) = engine_search(args, &query, &args.positionals[1], &mut sp)?;
    let mut tree = sabertooth::taxonomy::TaxTree::parse_nodes(&read_to_string(nodes_path).map_err(|e| e.to_string())?);
    if let Some(names) = args.opt_str("names") {
        tree.load_names(&read_to_string(names).map_err(|e| e.to_string())?);
    }
    let seqmap = sabertooth::taxonomy::parse_seqmap(&read_to_string(seqmap_path).map_err(|e| e.to_string())?);
    let frac = args.opt_f64("lca-frac", 0.0)?;
    let mut by_query: Vec<Vec<(u32, f64)>> = vec![Vec::new(); query.len()];
    let mut mapped = 0usize;
    for hits in &res {
        for (h, _) in hits {
            if let Some(&taxid) = seqmap.get(&h.target_id) {
                mapped += 1;
                by_query[h.query_index].push((taxid, h.bit_score));
            }
        }
    }
    let mut out = String::new();
    out.push_str("query\ttaxid\trank\tname\tassigned_hits\n");
    let mut classified = 0usize;
    for (qi, scored) in by_query.iter().enumerate() {
        let qid = query.id(qi);
        let assignment = if scored.is_empty() {
            None
        } else {
            let taxa = sabertooth::taxonomy::weighted_taxa(scored, frac);
            tree.lca(&taxa).map(|taxid| (taxid, taxa.len()))
        };
        match assignment {
            Some((taxid, n)) => {
                classified += 1;
                out.push_str(&format!("{}\t{}\t{}\t{}\t{}\n", qid, taxid, tree.rank_of(taxid), tree.name_of(taxid), n));
            }
            None => out.push_str(&format!("{}\t0\tno rank\tunclassified\t0\n", qid)),
        }
    }
    write_out(args.opt_str("out"), &out)?;
    eprintln!(
        "taxonomy: {} queries, {} taxon-mapped hits, {} classified (weighted LCA frac={})",
        query.len(),
        mapped,
        classified,
        frac
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// makedb / createindex / indexinfo / convert2fasta
// ---------------------------------------------------------------------------

fn run_makedb(args: &Args) -> Result<(), String> {
    args.check_known("makedb", &[])?;
    if args.positionals.len() < 2 {
        return Err("usage: sabertooth makedb <in.fasta[.gz]|db_prefix> <db_prefix>".into());
    }
    let input = &args.positionals[0];
    let out = std::path::Path::new(&args.positionals[1]);
    let fmt = seqio::detect_format(input).map_err(|e| format!("cannot read input '{}': {}", input, e))?;
    let (n, nucleotide) = match fmt {
        seqio::InputFormat::Fasta { .. } => {
            let mut fr = seqio::open_fasta(input).map_err(|e| format!("cannot read input '{}': {}", input, e))?;
            let (mut h, mut s) = (Vec::new(), Vec::new());
            if !fr.next_record(&mut h, &mut s).map_err(|e| format!("cannot read input '{}': {}", input, e))? {
                return Err(format!("cannot read input '{}': no FASTA records found (the input is empty)", input));
            }
            let nucleotide = sabertooth::nucl::looks_like_nucleotide(&s);
            let mut w = sabertooth::mmdb::DbWriter::create(out, nucleotide).map_err(|e| e.to_string())?;
            loop {
                w.push(&h, &s).map_err(|e| e.to_string())?;
                if !fr.next_record(&mut h, &mut s).map_err(|e| format!("cannot read input '{}': {}", input, e))? {
                    break;
                }
            }
            if fr.gap_chars_removed > 0 {
                eprintln!("sabertooth: warning: {}: removed {} alignment gap characters", input, fr.gap_chars_removed);
            }
            let n = w.len();
            w.finish(input).map_err(|e| e.to_string())?;
            (n, nucleotide)
        }
        seqio::InputFormat::MmseqsDb | seqio::InputFormat::SeedIndex => {
            let db = load_seqs(input, LoadOptions::FULL, SeqKind::Any, "input")?;
            let nucleotide = seqio::looks_nucleotide_db(&db);
            sabertooth::mmdb::write_db(out, &db, nucleotide, input).map_err(|e| e.to_string())?;
            (db.len(), nucleotide)
        }
    };
    eprintln!(
        "makedb: wrote {} sequences to {}{{,.index,.dbtype,_h,_h.index,.lookup,.source}} ({})",
        n,
        args.positionals[1],
        if nucleotide { "nucleotide" } else { "amino acid" }
    );
    eprintln!("  (sequence database only; build the seed index with `sabertooth createindex {}`)", args.positionals[1]);
    Ok(())
}

fn run_createindex(args: &Args) -> Result<(), String> {
    if args.positionals.is_empty() {
        return Err("usage: sabertooth createindex <db_prefix|fasta> [out.stidx] [--spaced mask (default 111001101) | --k n] [--mask 0|1] [--index-layout auto|sparse|direct] [--memory-limit size]".into());
    }
    args.check_known("createindex", &["k", "spaced", "mask", "profile", "index-layout"])?;
    let input = &args.positionals[0];
    let out = args
        .positionals
        .get(1)
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| index::default_index_path(input));
    let default_seed = if args.has("profile") { SeedPattern::profile_default() } else { SeedPattern::sequence_default() };
    let pattern = requested_pattern(args, &default_seed)?;
    let mut opts = BuildOptions::new(pattern);
    opts.mask = if args.opt_bool("mask", true)? { Some(MaskParams::default()) } else { None };
    opts.layout = match args.opt_str("index-layout").unwrap_or("auto") {
        "auto" => index::LayoutChoice::Auto,
        "sparse" => index::LayoutChoice::Sparse,
        "direct" => index::LayoutChoice::Direct,
        other => return Err(format!("--index-layout must be auto, sparse or direct (got '{}')", other)),
    };
    opts.memory_budget = match args.opt_bytes("memory-limit")? {
        Some(b) => (b / 2).max(16 << 20),
        None => 1 << 30,
    };
    opts.tmp_dir = out
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| std::path::PathBuf::from("."));
    if let Some(t) = args.opt_str("tmp-dir") {
        opts.tmp_dir = std::path::PathBuf::from(t);
    }
    let fp = SourceFingerprint::of(input).map_err(|e| format!("cannot read '{}': {}", input, e))?;
    let start = Instant::now();
    let mut src = index::source_for(input).map_err(|e| format!("cannot read '{}': {}", input, e))?;
    let st = index::build_index(src.as_mut(), fp, &out, &opts).map_err(|e| format!("index build failed: {}", e))?;
    let ps = proc_stats();
    eprintln!(
        "createindex: {} sequences, {} residues -> {} ({} postings of {} bytes, {} distinct {}-mers, {} lookup, {} masked residues, {} partition(s), {:.1} MiB) in {:.2?}; peak RSS {:.1} MiB",
        st.num_seqs,
        st.total_residues,
        out.display(),
        st.num_postings,
        st.posting_width,
        st.num_kmers,
        opts.pattern.weight(),
        st.layout,
        st.masked_residues,
        st.partitions,
        st.file_bytes as f64 / 1048576.0,
        start.elapsed(),
        ps.peak_rss as f64 / 1048576.0
    );
    if let Some(path) = args.opt_str("stats-json") {
        let j = Json::new()
            .str("command", "createindex")
            .str("version", VERSION)
            .str("output", &out.display().to_string())
            .str("seed_pattern", &opts.pattern.mask_string())
            .bool("masked", opts.mask.is_some())
            .num("memory_budget_bytes", opts.memory_budget)
            .num("threads", rayon::current_num_threads())
            .float("wall_s", start.elapsed().as_secs_f64())
            .obj("build", build_stats_json(&st))
            .obj(
                "memory",
                Json::new()
                    .num("peak_rss_bytes", ps.peak_rss)
                    .num("minor_page_faults", ps.minflt)
                    .num("major_page_faults", ps.majflt),
            );
        std::fs::write(path, j.finish() + "\n").map_err(|e| e.to_string())?;
    }
    Ok(())
}

fn run_indexinfo(args: &Args) -> Result<(), String> {
    args.check_known("indexinfo", &["verify"])?;
    if args.positionals.is_empty() {
        return Err("usage: sabertooth indexinfo <file.stidx> [--verify]".into());
    }
    let p = std::path::Path::new(&args.positionals[0]);
    let idx = SeedIndex::open(p).map_err(|e| e.to_string())?;
    let m = idx.meta();
    println!("index          : {}", p.display());
    println!("format version : {}", m.version);
    println!("built by       : sabertooth {} (unix time {})", m.tool_version, m.created_unix);
    println!("seed pattern   : {} (k = {})", m.pattern.mask_string(), m.pattern.weight());
    println!("seed masking   : {}", if m.masked { format!("SEG-like w={} {}/{}", m.mask.window, m.mask.locut, m.mask.hicut) } else { "off".into() });
    println!("sequences      : {}", m.num_seqs);
    println!("residues       : {} ({} masked for seeding)", m.total_residues, m.masked_residues);
    println!("longest        : {}", m.max_seq_len);
    println!("postings       : {} ({} bytes each)", m.num_postings, m.posting_width);
    println!("distinct k-mers: {}", m.num_kmers);
    println!(
        "k-mer lookup   : {}",
        match m.layout {
            index::IndexLayout::Sparse => "sparse (bucket tables + one entry per distinct k-mer)".to_string(),
            index::IndexLayout::Direct if m.offset_width == 2 => format!(
                "direct (16-bit offsets per possible k-mer in groups of {}, with group bases)",
                index::OFFSET_GROUP
            ),
            index::IndexLayout::Direct => format!("direct (one {}-byte offset per possible k-mer)", m.offset_width),
        }
    );
    println!("file size      : {:.1} MiB", idx.mapped_bytes() as f64 / 1048576.0);
    println!(
        "source         : {} ({})",
        if m.source.path.is_empty() { "<transient>" } else { &m.source.path },
        match m.source.kind {
            1 => "FASTA",
            2 => "MMseqs2 database",
            _ => "none",
        }
    );
    if !m.source.path.is_empty() {
        match idx.usable_for(idx.pattern(), if m.masked { Some(m.mask) } else { None }, Some(&m.source.path)) {
            Ok(()) => println!("freshness      : up to date with source"),
            Err(e) => println!("freshness      : {}", e),
        }
    }
    if args.has("verify") {
        let t = Instant::now();
        idx.verify().map_err(|e| e.to_string())?;
        println!("verify         : all section checksums OK ({:.2?})", t.elapsed());
    }
    Ok(())
}

fn run_convert2fasta(args: &Args) -> Result<(), String> {
    args.check_known("convert2fasta", &[])?;
    if args.positionals.is_empty() {
        return Err("usage: sabertooth convert2fasta <db_prefix> [out.fasta]".into());
    }
    let db = load_seqs(&args.positionals[0], LoadOptions::FULL, SeqKind::Any, "database")?;
    let dest = args.positionals.get(1).map(|s| s.as_str()).or_else(|| args.opt_str("out"));
    let w: Box<dyn Write> = match dest {
        Some(p) => Box::new(File::create(p).map_err(|e| format!("cannot create '{}': {}", p, e))?),
        None => Box::new(io::stdout()),
    };
    let mut out = BufWriter::with_capacity(1 << 20, w);
    for rec in db.iter() {
        out.write_all(b">").map_err(|e| e.to_string())?;
        out.write_all(rec.header.as_bytes()).map_err(|e| e.to_string())?;
        out.write_all(b"\n").map_err(|e| e.to_string())?;
        out.write_all(&rec.seq()).map_err(|e| e.to_string())?;
        out.write_all(b"\n").map_err(|e| e.to_string())?;
    }
    out.flush().map_err(|e| e.to_string())?;
    eprintln!("convert2fasta: {} sequences read from {}", db.len(), args.positionals[0]);
    Ok(())
}

// ---------------------------------------------------------------------------
// cluster
// ---------------------------------------------------------------------------

fn run_cluster(args: &Args) -> Result<(), String> {
    args.check_known("cluster", &["k", "min-count", "min-cov", "min-seq-id"])?;
    if args.positionals.is_empty() {
        return Err("usage: sabertooth cluster <in> [--min-seq-id f] [--min-cov f] [--k n]".into());
    }
    let m = SubstitutionMatrix::blosum62();
    let db = load_seqs(&args.positionals[0], LoadOptions::SEARCH, SeqKind::Protein, "input")?;
    let mut p = sabertooth::cluster::ClusterParams::default();
    p.k = args.opt_usize("k", p.k)?;
    p.num_minimizers = args.opt_usize("min-count", p.num_minimizers)?;
    p.min_seq_id = args.opt_f64("min-seq-id", p.min_seq_id)?;
    p.min_cov = args.opt_f64("min-cov", p.min_cov)?;
    let start = Instant::now();
    let clusters = sabertooth::cluster::linclust(&db, &m, p);
    let elapsed = start.elapsed();
    let mut out = String::new();
    for c in &clusters {
        let rep = db.id(c.representative);
        for &mem in &c.members {
            out.push_str(rep);
            out.push('\t');
            out.push_str(db.id(mem));
            out.push('\n');
        }
    }
    write_out(args.opt_str("out"), &out)?;
    let largest = clusters.iter().map(|c| c.members.len()).max().unwrap_or(0);
    let singletons = clusters.iter().filter(|c| c.members.len() == 1).count();
    eprintln!(
        "cluster: {} sequences -> {} clusters (largest {}, {} singletons) at id>={:.0}% cov>={:.0}% in {:.2?}",
        db.len(),
        clusters.len(),
        largest,
        singletons,
        p.min_seq_id * 100.0,
        p.min_cov * 100.0,
        elapsed
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// calibrate
// ---------------------------------------------------------------------------

fn run_calibrate(args: &Args) -> Result<(), String> {
    args.check_known("calibrate", &["comp-bias-corr", "gap-extend", "gap-open", "json", "legacy", "lengths", "samples", "seed", "table", "length", "pairs"])?;
    let m = SubstitutionMatrix::blosum62();
    if args.has("legacy") {
        return run_calibrate_legacy(args, &m);
    }
    let gaps = GapCosts {
        open: args.opt_i32("gap-open", 11)?,
        extend: args.opt_i32("gap-extend", 1)?,
    };
    let samples = args.opt_usize("samples", 4000)?;
    let seed = args.opt_usize("seed", 0x5EED)? as u64;
    let lengths: Vec<usize> = match args.opt_str("lengths") {
        Some(l) if !l.is_empty() => l
            .split(',')
            .map(|x| x.trim().parse::<usize>().map_err(|_| format!("bad --lengths entry '{}'", x)))
            .collect::<Result<_, _>>()?,
        _ => sabertooth::evalue::DEFAULT_CALIBRATION_LENGTHS.to_vec(),
    };
    let modes: Vec<bool> = match args.opt_str("comp-bias-corr") {
        Some("0") => vec![false],
        Some("1") => vec![true],
        _ => vec![false, true],
    };
    let mut configs = Vec::new();
    let mut table_src = String::new();
    for &cbc in &modes {
        let start = Instant::now();
        let cal = sabertooth::evalue::calibrate_sequence(&m, gaps, cbc, &lengths, samples, seed)
            .ok_or("calibration failed (degenerate score distribution)")?;
        let p = cal.params;
        let lm = p.length.unwrap();
        println!(
            "gaps {}/{}  comp-bias-corr {}  ({} samples x {} lengths, seed {}, {:.1?})",
            gaps.open,
            gaps.extend,
            cbc as u8,
            samples,
            lengths.len(),
            seed,
            start.elapsed()
        );
        println!(
            "  lambda = {:.4}   K = {:.4}   alpha/lambda = {:.1}   beta = {:.0}   (ln K rms across lengths {:.3})",
            p.lambda, p.k, lm.alpha_d_lambda, lm.beta, cal.ln_k_rms
        );
        println!("  {:>6} {:>6} {:>8} {:>8} {:>10} {:>9}", "m", "n", "mean", "sd", "lambda_ml", "mu");
        let mut cells = Vec::new();
        for c in &cal.cells {
            println!("  {:>6} {:>6} {:>8.2} {:>8.2} {:>10.4} {:>9.2}", c.m, c.n, c.mean, c.sd, c.lambda_ml, c.mu);
            cells.push(
                Json::new()
                    .num("m", c.m)
                    .num("n", c.n)
                    .num("samples", c.samples)
                    .float("mean", c.mean)
                    .float("sd", c.sd)
                    .float("lambda_ml", c.lambda_ml)
                    .float("mu_ml", c.mu_ml)
                    .float("mu_pooled", c.mu)
                    .finish(),
            );
        }
        match sabertooth::evalue::table_params(gaps, cbc) {
            Some(t) => println!("  table: lambda {:.4} K {:.4} ({})", t.params.lambda, t.params.k, t.provenance),
            None => println!("  table: no entry for this configuration"),
        }
        println!();
        configs.push(
            Json::new()
                .bool("comp_bias_corr", cbc)
                .float("lambda", p.lambda)
                .float("k", p.k)
                .float("alpha_d_lambda", lm.alpha_d_lambda)
                .float("beta", lm.beta)
                .float("ln_k_rms", cal.ln_k_rms)
                .raw("cells", &format!("[{}]", cells.join(", ")))
                .finish(),
        );
        table_src.push_str(&format!(
            "    TableEntry {{\n        gaps: GapCosts {{ open: {}, extend: {} }},\n        comp_bias: {},\n        params: GappedParams {{\n            lambda: {:.4},\n            k: {:.4},\n            length: Some(LengthModel {{\n                alpha_d_lambda: {:.1},\n                beta: {:.1},\n            }}),\n        }},\n        provenance: \"sabertooth {} calibrate --table --samples {} --seed {} (lengths {:?})\",\n    }},\n",
            gaps.open, gaps.extend, cbc, p.lambda, p.k, lm.alpha_d_lambda, lm.beta, VERSION, samples, seed, lengths
        ));
    }
    if let Some(path) = args.opt_str("json") {
        let j = Json::new()
            .str("tool", "sabertooth calibrate")
            .str("version", VERSION)
            .str("matrix", "BLOSUM62 bit_factor 2 bias -0.2")
            .num("gap_open", gaps.open)
            .num("gap_extend", gaps.extend)
            .num("samples_per_cell", samples)
            .num("seed", seed)
            .raw("configs", &format!("[{}]", configs.join(", ")));
        std::fs::write(path, j.finish() + "\n").map_err(|e| e.to_string())?;
    }
    if args.has("table") {
        println!("// paste into SEQUENCE_TABLE in src/evalue.rs");
        print!("{}", table_src);
    }
    Ok(())
}

fn run_calibrate_legacy(args: &Args, m: &SubstitutionMatrix) -> Result<(), String> {
    let seq_len = args.opt_usize("length", 250)?;
    let num_pairs = args.opt_usize("pairs", 3000)?;
    let seed = args.opt_usize("seed", 0x5EED)? as u64;
    let gaps = GapCosts::sequence_default();
    let start = Instant::now();
    let g = sabertooth::evalue::calibrate_gapped(m, gaps, seq_len, num_pairs, seed);
    let elapsed = start.elapsed();
    let ungapped = EValueParams::new(m, 1).lambda;
    println!("Sabertooth legacy single-length calibration (method of moments)");
    println!(
        "  simulated {} random pairs of length {} (gap open {}, extend {}) in {:.2?}",
        g.samples, seq_len, gaps.open, gaps.extend, elapsed
    );
    println!("  optimal local score: mean {:.2}, std {:.2}", g.mean_score, g.std_score);
    println!("  analytic ungapped lambda = {:.4}", ungapped);
    println!("  moments gapped    lambda = {:.4}   K = {:.4}   (mu = {:.2})", g.lambda, g.k, g.mu);
    Ok(())
}

// ---------------------------------------------------------------------------
// output helper
// ---------------------------------------------------------------------------

fn write_out(path: Option<&str>, content: &str) -> Result<(), String> {
    match path {
        Some(p) => {
            let mut f = File::create(p).map_err(|e| format!("cannot create '{}': {}", p, e))?;
            f.write_all(content.as_bytes()).map_err(|e| e.to_string())?;
            Ok(())
        }
        None => {
            let stdout = io::stdout();
            let mut lock = stdout.lock();
            lock.write_all(content.as_bytes()).map_err(|e| e.to_string())?;
            Ok(())
        }
    }
}
