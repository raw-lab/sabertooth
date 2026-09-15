//! Sabertooth command-line interface.
//!
//! Dependency-light hand-rolled argument parsing (only `rayon` is pulled in for
//! parallelism). Subcommands: `createdb`, `msa2profile`, `search`,
//! `profilesearch`, `align`, plus the four flavoured banners `doctor`,
//! `version`, `info`, and `help`.

use sabertooth::align::{local_align, GapCosts, SeqScorer};
use sabertooth::alphabet::{num_to_aa, PROFILE_AA_SIZE};
use sabertooth::banner;
use sabertooth::evalue::EValueParams;
use sabertooth::fasta::{read_to_string, SeqDb};
use sabertooth::matrix::SubstitutionMatrix;
use sabertooth::prefilter::{check_k, KmerIndex, PrefilterParams, SeedPattern};
use sabertooth::profile::{Profile, PseudoCountParams};
use sabertooth::search::{
    search_profile, search_sequences, search_sequences_exhaustive, SearchParams, M8_HEADER,
};
use sabertooth::VERSION;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Write};
use std::process::ExitCode;
use std::time::Instant;

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
                // short sensitivity flag, mirrors MMseqs2's -s
                if i + 1 < argv.len() && !argv[i + 1].starts_with('-') {
                    options.insert("sensitivity".to_string(), argv[i + 1].clone());
                    i += 2;
                } else {
                    i += 1;
                }
            } else if let Some(key) = a.strip_prefix("--") {
                // is the next token a value or another flag?
                if i + 1 < argv.len() && !argv[i + 1].starts_with("--") {
                    options.insert(key.to_string(), argv[i + 1].clone());
                    i += 2;
                } else {
                    // bare flag with no value; recorded as an empty-valued option
                    options.insert(key.to_string(), String::new());
                    i += 1;
                }
            } else {
                positionals.push(a.clone());
                i += 1;
            }
        }
        Args {
            positionals,
            options,
        }
    }

    fn opt_f64(&self, key: &str, default: f64) -> f64 {
        self.options
            .get(key)
            .and_then(|v| v.parse().ok())
            .unwrap_or(default)
    }
    fn opt_usize(&self, key: &str, default: usize) -> usize {
        self.options
            .get(key)
            .and_then(|v| v.parse().ok())
            .unwrap_or(default)
    }
    fn opt_i32(&self, key: &str, default: i32) -> i32 {
        self.options
            .get(key)
            .and_then(|v| v.parse().ok())
            .unwrap_or(default)
    }
    fn opt_str(&self, key: &str) -> Option<&str> {
        self.options.get(key).map(|s| s.as_str())
    }
}

fn main() -> ExitCode {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    if argv.is_empty() {
        print_help();
        return ExitCode::SUCCESS;
    }

    let cmd = argv[0].clone();
    let rest = &argv[1..];
    let args = Args::parse(rest);

    // configure rayon thread pool if requested
    if let Some(t) = args.opt_str("threads").and_then(|v| v.parse::<usize>().ok()) {
        let _ = rayon::ThreadPoolBuilder::new().num_threads(t).build_global();
    }

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
        "search" => run_search(&args),
        "profilesearch" => run_profilesearch(&args),
        "align" => run_align(&args),
        "result2msa" => run_result2msa(&args),
        "translatesearch" => run_translatesearch(&args),
        "profileprofile" => run_profileprofile(&args),
        "cluster" => run_cluster(&args),
        "easysearch" => run_easysearch(&args),
        "makedb" => run_makedb(&args),
        "convert2fasta" => run_convert2fasta(&args),
        "taxonomy" => run_taxonomy(&args),
        "calibrate" => run_calibrate(&args),
        other => {
            eprintln!("sabertooth: unknown command '{}'\n", other);
            print_help();
            return ExitCode::from(2);
        }
    };

    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("sabertooth: error: {}", e);
            ExitCode::FAILURE
        }
    }
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
    println!("  Lambda (hdr)  : {:.5}", m.lambda);
    print!("  Background    :");
    for i in 0..PROFILE_AA_SIZE {
        if i % 10 == 0 {
            print!("\n     ");
        }
        print!(" {}:{:.4}", num_to_aa(i as u8) as char, m.p_back[i]);
    }
    println!("\n");
    println!("  Profile PSSM  : bit factor 8, pca=1.0 pcb=1.5 (MMseqs2 defaults)");
    println!("  Pipeline      : kmer prefilter -> Smith-Waterman-Gotoh -> Karlin-Altschul");
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

COMMANDS:
  createdb <in.fasta> <out.fasta>        Validate + normalise a sequence FASTA
  makedb <in.fasta> <db_prefix>          Write an MMseqs2-format on-disk DB
  convert2fasta <db_prefix> [out.fasta]  Read an MMseqs2-format DB back to FASTA
  msa2profile <msa>                      Build a PSSM from an MSA (a3m/aligned FASTA)
  result2msa <query.fasta> <target>      Search, then emit an a3m MSA per query
  search <query.fasta> <target.fasta>    Sequence-vs-sequence search
  easysearch <query.fasta> <target>      Database-free exhaustive SW (no prefilter)
  profilesearch <msa> <target.fasta>     Profile (PSSM) search vs a target DB
  profileprofile <query.msa> <t.msa...>  Profile-vs-profile search (co-emission)
  translatesearch <nucl_q> <prot_target> Six-frame translated nucleotide search
  cluster <in.fasta>                     Linclust-style clustering (rep→member TSV)
  taxonomy <query> <target>              LCA taxonomic assignment from hits
  align <query.fasta> <target.fasta>     Pairwise local alignments (query 1 vs targets)
  calibrate                              Monte-Carlo gapped Gumbel statistics (λ, K)
  info                                   Matrix / statistics report
  doctor                                 Self-check and environment report
  version                                Version banner
  help                                   This help

COMMON OPTIONS:
  --out <path>          Write results to a file instead of stdout
  --threads <n>         Number of worker threads (default: all cores)

SEARCH OPTIONS:
  --k <int>             K-mer length (default: seq 6, profile 5)
  --kmer-score <int>    Similar-k-mer score threshold (seq default 25, profile 40)
  -s <float>            Sensitivity 1..9 (higher = more sensitive; overrides --kmer-score)
  --spaced <mask>       Spaced seed mask, e.g. 110101 (must start/end with 1)
  --band <radius>       Banded traceback of this radius around the seed diagonal
                        (O(m·w) memory; use for very long sequences)
  --evalue <float>      Max E-value to report (default 1e-3)
  --min-cov <float>     Min query coverage fraction 0..1 (default 0)
  --max-hits <int>      Max hits per query, 0 = unlimited (default 300)

TRANSLATED / CLUSTER / TAXONOMY OPTIONS:
  --min-orf <int>       Min ORF length in aa for translatesearch (default 20)
  --min-seq-id <float>  Min identity for cluster membership (default 0.5)
  --min-cov <float>     Min coverage for cluster membership (default 0.8)
  --nodes <nodes.dmp>   NCBI taxonomy tree            (taxonomy, required)
  --names <names.dmp>   NCBI scientific names         (taxonomy, optional)
  --seqmap <tsv>        seqid<TAB>taxid map            (taxonomy, required)
  --lca-frac <float>    Keep hits within this score fraction of best; 0 = all (default 0)

DISTRIBUTED OPTIONS (build with `--features distributed`):
  --distributed         Run search across HydraMPP workers/nodes
  --shards <int>        Number of target shards (default 4)
  --gpus <int>          GPUs reserved per shard (GPU-aware scheduling + pinning)
  --hydra-host          Run as the HydraMPP head node (consumed by HydraMPP)
  --hydra-client <ip>   Join a HydraMPP head as a worker (consumed by HydraMPP)
  --hydra-gpus <int>    Advertise/simulate N GPUs on this node (consumed by HydraMPP)

PROFILE OPTIONS:
  --pca <float>            Pseudocount admixture a (default 1.0)
  --pcb <float>            Pseudocount admixture b (default 1.5)
  --bias <float>           PSSM score bias (default 0.0)
  --num-iterations <n>     Iterative (PSI-BLAST-style) rounds (default 1)
  --inclusion-evalue <f>   Max E-value to fold a hit into the next round (default 1e-3)

CALIBRATE OPTIONS:
  --length <int>        Random sequence length (default 250)
  --pairs <int>         Number of random pairs to align (default 3000)
  --seed <int>          PRNG seed (default 24301)

EXAMPLES:
  sabertooth msa2profile family.a3m --out family.pssm
  sabertooth profilesearch family.a3m proteins.fasta --out hits.m8 --evalue 1e-5
  sabertooth result2msa queries.fasta proteins.fasta --out msa.a3m
  sabertooth translatesearch contigs.fasta proteins.fasta --min-orf 30
  sabertooth profileprofile familyA.a3m familyB.a3m familyC.a3m
  sabertooth cluster proteins.fasta --min-seq-id 0.5 --out clusters.tsv
  sabertooth easysearch queries.fasta small_db.fasta
  sabertooth makedb proteins.fasta protDB   # interoperable with `mmseqs`
  sabertooth taxonomy queries.fasta proteins.fasta --nodes nodes.dmp --names names.dmp --seqmap map.tsv
  sabertooth search q.fasta db.fasta --distributed --shards 8 --gpus 1 --hydra-gpus 4
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

    // 1. matrix reconstruction
    let m = SubstitutionMatrix::blosum62();
    let ww = m.score(
        m.aa2num[b'W' as usize],
        m.aa2num[b'W' as usize],
    );
    let check_matrix = ww > 0 && m.lambda > 0.0;
    report("substitution matrix (BLOSUM62) reconstructed", check_matrix);
    ok &= check_matrix;

    // 2. background sums to ~1
    let bsum: f64 = (0..PROFILE_AA_SIZE).map(|i| m.p_back[i]).sum();
    let check_bg = (bsum - 1.0).abs() < 0.01;
    report(&format!("background distribution sums to 1 ({:.5})", bsum), check_bg);
    ok &= check_bg;

    // 3. R columns are conditional distributions
    let mut check_r = true;
    for b in 0..PROFILE_AA_SIZE {
        let s: f32 = (0..PROFILE_AA_SIZE).map(|a| m.r[a][b]).sum();
        if (s - 1.0).abs() > 0.03 {
            check_r = false;
        }
    }
    report("pseudocount matrix R columns normalised", check_r);
    ok &= check_r;

    // 4. E-value lambda in range
    let ev = EValueParams::new(&m, 1_000_000);
    let check_lambda = ev.lambda > 0.25 && ev.lambda < 0.45;
    report(&format!("Karlin-Altschul lambda in range ({:.4})", ev.lambda), check_lambda);
    ok &= check_lambda;

    // 5. profile round-trip on a toy MSA
    let msa = ">q\nMKVLLACDEFGHIK\n>h\nMKILLSCDEFGHLK\n>i\nMRVLLACEEFGHIK\n";
    let prof = Profile::from_msa_str(msa, &m, PseudoCountParams::default());
    let check_profile = match &prof {
        Ok(p) => p.query_len == 14 && p.set_size == 3,
        Err(_) => false,
    };
    report("profile/PSSM construction from MSA", check_profile);
    ok &= check_profile;

    // 6. alignment sanity
    let check_align = if let Ok(p) = &prof {
        let sc = sabertooth::align::ProfileScorer { profile: p };
        let t: Vec<u8> = "MKVLLACDEFGHIK".bytes().map(|b| m.aa2num[b as usize]).collect();
        local_align(&sc, &t, GapCosts::profile_default())
            .map(|a| a.raw_score > 0)
            .unwrap_or(false)
    } else {
        false
    };
    report("Smith-Waterman-Gotoh alignment", check_align);
    ok &= check_align;

    // 7. prefilter index
    let t = sabertooth::alphabet::build_aa2num();
    let db = SeqDb::from_reader(
        std::io::Cursor::new(">t\nGGGMKVLLACDEFGHIKGGG\n"),
        &t,
    )
    .map_err(|e| e.to_string())?;
    let idx = KmerIndex::build(&db, 6);
    let check_index = idx.distinct_kmers() > 0;
    report(&format!("k-mer index build ({} distinct 6-mers)", idx.distinct_kmers()), check_index);
    ok &= check_index;

    // environment
    println!();
    println!("  environment:");
    println!("    sabertooth version : {}", VERSION);
    println!("    worker threads     : {}", rayon::current_num_threads());
    println!(
        "    pointer width      : {} bit",
        std::mem::size_of::<usize>() * 8
    );
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
// createdb
// ---------------------------------------------------------------------------

fn run_createdb(args: &Args) -> Result<(), String> {
    if args.positionals.len() < 2 {
        return Err("usage: sabertooth createdb <in.fasta> <out.fasta>".into());
    }
    let t = sabertooth::alphabet::build_aa2num();
    let db = SeqDb::read_fasta(&args.positionals[0], &t).map_err(|e| e.to_string())?;
    db.write_fasta(&args.positionals[1]).map_err(|e| e.to_string())?;
    println!(
        "createdb: {} sequences, {} residues -> {}",
        db.len(),
        db.total_residues,
        args.positionals[1]
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// msa2profile
// ---------------------------------------------------------------------------

fn pc_params(args: &Args) -> PseudoCountParams {
    PseudoCountParams {
        pca: args.opt_f64("pca", 1.0) as f32,
        pcb: args.opt_f64("pcb", 1.5) as f32,
        score_bias: args.opt_f64("bias", 0.0) as f32,
    }
}

fn run_msa2profile(args: &Args) -> Result<(), String> {
    if args.positionals.is_empty() {
        return Err("usage: sabertooth msa2profile <msa> [--out file]".into());
    }
    let m = SubstitutionMatrix::blosum62();
    let data = read_to_string(&args.positionals[0]).map_err(|e| e.to_string())?;
    let prof = Profile::from_msa_str(&data, &m, pc_params(args))?;
    let table = prof.pssm_table();
    write_out(args.opt_str("out"), &table)?;
    eprintln!(
        "msa2profile: {} sequences, {} match columns",
        prof.set_size, prof.query_len
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// search
// ---------------------------------------------------------------------------

/// Build the prefilter index, honoring an optional `--spaced <mask>`. Adjusts
/// `sp.prefilter.k` to the seed weight so downstream assertions hold.
fn build_prefilter_index(
    args: &Args,
    target: &sabertooth::fasta::SeqDb,
    sp: &mut SearchParams,
) -> Result<KmerIndex, String> {
    match args.opt_str("spaced").filter(|m| !m.is_empty()) {
        Some(mask) => {
            let pattern = SeedPattern::from_mask(mask)?;
            check_k(pattern.weight())?;
            eprintln!(
                "seed: spaced mask {} (weight {}, span {})",
                pattern.mask_string(),
                pattern.weight(),
                pattern.span()
            );
            sp.prefilter.k = pattern.weight();
            Ok(KmerIndex::build_spaced(target, pattern))
        }
        None => {
            check_k(sp.prefilter.k)?;
            Ok(KmerIndex::build(target, sp.prefilter.k))
        }
    }
}

fn search_params(args: &Args, profile_mode: bool) -> SearchParams {
    let mut sp = SearchParams::default();
    let scale = if profile_mode { 8 } else { 2 };
    // -s / --sensitivity overrides the explicit --kmer-score when given.
    let kmer_score = match args
        .opt_str("sensitivity")
        .and_then(|v| v.parse::<f64>().ok())
    {
        Some(s) => sabertooth::prefilter::kmer_score_for_sensitivity(s, scale),
        None => args.opt_i32("kmer-score", if profile_mode { 40 } else { 25 }),
    };
    sp.prefilter = PrefilterParams {
        k: args.opt_usize("k", if profile_mode { 5 } else { 6 }),
        kmer_score,
        max_kmers_per_pos: args.opt_usize("max-kmers", 256),
        min_hits: args.opt_usize("min-hits", 1),
        min_diag_score: args.opt_i32("min-diag-score", if profile_mode { 30 } else { 15 }),
    };
    sp.max_evalue = args.opt_f64("evalue", 1e-3);
    sp.min_query_cov = args.opt_f64("min-cov", 0.0);
    sp.max_hits = args.opt_usize("max-hits", 300);
    sp.band_radius = args.opt_str("band").and_then(|v| v.parse::<usize>().ok());
    sp
}

fn run_search(args: &Args) -> Result<(), String> {
    if args.positionals.len() < 2 {
        return Err("usage: sabertooth search <query.fasta> <target.fasta> [options]".into());
    }
    let m = SubstitutionMatrix::blosum62();
    let t = &m.aa2num;
    let query = SeqDb::read_fasta(&args.positionals[0], t).map_err(|e| e.to_string())?;
    let target = SeqDb::read_fasta(&args.positionals[1], t).map_err(|e| e.to_string())?;

    let mut sp = search_params(args, false);
    if sp.max_evalue <= 0.0 {
        return Err("--evalue must be positive".into());
    }

    // Optional HydraMPP-distributed path (only when built with --features distributed).
    #[cfg(feature = "distributed")]
    if args.options.contains_key("distributed") {
        let shards = args.opt_usize("shards", 4);
        let gpus_per_shard = args.opt_usize("gpus", 0);
        let ev = EValueParams::new(&m, target.total_residues);
        let start = Instant::now();
        let outcome = sabertooth::distributed::search_sequences_distributed(
            &query,
            &target,
            &ev,
            sp,
            shards,
            gpus_per_shard,
        )?;
        let elapsed = start.elapsed();
        let mut out = String::new();
        out.push_str(M8_HEADER);
        out.push('\n');
        for h in &outcome.hits {
            out.push_str(&h.to_m8());
            out.push('\n');
        }
        write_out(args.opt_str("out"), &out)?;
        if gpus_per_shard > 0 {
            eprintln!(
                "search (HydraMPP, {} shards, {} GPU/shard): {} queries x {} targets, {} hits in {:.2?}",
                outcome.shards, gpus_per_shard, query.len(), target.len(), outcome.hits.len(), elapsed
            );
            if outcome.gpu_devices.is_empty() {
                eprintln!("  no GPU devices were available/reserved (run with `-- --hydra-gpus N` to simulate)");
            } else {
                eprintln!("  GPU devices pinned to shards: {}", outcome.gpu_devices.join(" | "));
                eprintln!("  (compute ran on the SIMD CPU kernel; this is the pin-point where a CUDA kernel would execute)");
            }
        } else {
            eprintln!(
                "search (HydraMPP, {} shards): {} queries x {} targets, {} hits in {:.2?}",
                outcome.shards, query.len(), target.len(), outcome.hits.len(), elapsed
            );
        }
        return Ok(());
    }

    let start = Instant::now();
    let idx = build_prefilter_index(args, &target, &mut sp)?;
    let ev = EValueParams::new(&m, target.total_residues);
    let hits = search_sequences(&query, &target, &idx, &m, &ev, sp);
    let elapsed = start.elapsed();

    let mut out = String::new();
    out.push_str(M8_HEADER);
    out.push('\n');
    for h in &hits {
        out.push_str(&h.to_m8());
        out.push('\n');
    }
    write_out(args.opt_str("out"), &out)?;
    eprintln!(
        "search: {} queries x {} targets, {} hits in {:.2?}",
        query.len(),
        target.len(),
        hits.len(),
        elapsed
    );
    Ok(())
}

fn run_profilesearch(args: &Args) -> Result<(), String> {
    if args.positionals.len() < 2 {
        return Err("usage: sabertooth profilesearch <msa> <target.fasta> [options]".into());
    }
    let m = SubstitutionMatrix::blosum62();
    let t = &m.aa2num;
    let msa_data = read_to_string(&args.positionals[0]).map_err(|e| e.to_string())?;
    let prof = Profile::from_msa_str(&msa_data, &m, pc_params(args))?;
    let target = SeqDb::read_fasta(&args.positionals[1], t).map_err(|e| e.to_string())?;

    // profile id = basename of the MSA file
    let query_id = std::path::Path::new(&args.positionals[0])
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("profile")
        .to_string();

    let mut sp = search_params(args, true);
    if sp.max_evalue <= 0.0 {
        return Err("--evalue must be positive".into());
    }
    let start = Instant::now();
    let idx = build_prefilter_index(args, &target, &mut sp)?;
    let ev = EValueParams::for_profile(target.total_residues);
    let num_iter = args.opt_usize("num-iterations", 1);
    let hits = if num_iter > 1 {
        let (msa, query_len) = sabertooth::profile::parse_msa(&msa_data, &m.aa2num)?;
        let incl = args.opt_f64("inclusion-evalue", 1e-3);
        eprintln!(
            "profilesearch: iterative ({} rounds, inclusion E <= {:.0e})",
            num_iter, incl
        );
        sabertooth::search::search_profile_iterative(
            &msa,
            query_len,
            &query_id,
            &target,
            &idx,
            &m,
            pc_params(args),
            &ev,
            sp,
            num_iter,
            incl,
        )?
    } else {
        search_profile(&prof, &query_id, &target, &idx, &m, &ev, sp)
    };
    let elapsed = start.elapsed();

    let mut out = String::new();
    out.push_str(M8_HEADER);
    out.push('\n');
    for h in &hits {
        out.push_str(&h.to_m8());
        out.push('\n');
    }
    write_out(args.opt_str("out"), &out)?;
    eprintln!(
        "profilesearch: profile ({} cols, {} seqs) x {} targets, {} hits in {:.2?}",
        prof.query_len,
        prof.set_size,
        target.len(),
        hits.len(),
        elapsed
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// align (pairwise)
// ---------------------------------------------------------------------------

fn run_align(args: &Args) -> Result<(), String> {
    if args.positionals.len() < 2 {
        return Err("usage: sabertooth align <query.fasta> <target.fasta>".into());
    }
    let m = SubstitutionMatrix::blosum62();
    let t = &m.aa2num;
    let query = SeqDb::read_fasta(&args.positionals[0], t).map_err(|e| e.to_string())?;
    let target = SeqDb::read_fasta(&args.positionals[1], t).map_err(|e| e.to_string())?;
    if query.is_empty() {
        return Err("query FASTA is empty".into());
    }
    let q = &query.records[0];
    let ev = EValueParams::new(&m, target.total_residues);
    let sc = SeqScorer {
        query: &q.num,
        mat: &m,
    };
    let mut out = String::new();
    for trec in &target.records {
        if let Some(aln) = local_align(&sc, &trec.num, GapCosts::sequence_default()) {
            out.push_str(&format!(
                "{} vs {}: score={} bits={:.1} evalue={:.2e} id={:.1}% len={} cigar={}\n  Q[{}-{}] T[{}-{}]\n",
                q.id,
                trec.id,
                aln.raw_score,
                ev.bit_score(aln.raw_score),
                ev.evalue(aln.raw_score, q.num.len()),
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
// result2msa (a3m MSA generation)
// ---------------------------------------------------------------------------

fn run_result2msa(args: &Args) -> Result<(), String> {
    if args.positionals.len() < 2 {
        return Err("usage: sabertooth result2msa <query.fasta> <target.fasta> [options]".into());
    }
    let m = SubstitutionMatrix::blosum62();
    let t = &m.aa2num;
    let query = SeqDb::read_fasta(&args.positionals[0], t).map_err(|e| e.to_string())?;
    let target = SeqDb::read_fasta(&args.positionals[1], t).map_err(|e| e.to_string())?;
    if query.is_empty() {
        return Err("query FASTA is empty".into());
    }
    let mut sp = search_params(args, false);
    let idx = build_prefilter_index(args, &target, &mut sp)?;
    let ev = EValueParams::new(&m, target.total_residues);
    let gaps = GapCosts::sequence_default();
    let start = Instant::now();

    let mut out = String::new();
    let mut total_rows = 0usize;
    for q in &query.records {
        let cands =
            sabertooth::prefilter::prefilter_sequence(&q.num, &target, &idx, &m, sp.prefilter);
        let sc = SeqScorer { query: &q.num, mat: &m };
        let mut hits: Vec<(f64, String, Vec<u8>, sabertooth::align::Alignment)> = Vec::new();
        for c in &cands {
            let tgt = &target.records[c.target_id as usize];
            let aln_opt = match sp.band_radius {
                Some(r) => sabertooth::align::local_align_banded(&sc, &tgt.num, gaps, c.best_diag, r),
                None => local_align(&sc, &tgt.num, gaps),
            };
            if let Some(aln) = aln_opt {
                let e = ev.evalue(aln.raw_score, q.num.len());
                if e <= sp.max_evalue {
                    hits.push((e, tgt.id.clone(), tgt.num.clone(), aln));
                }
            }
        }
        hits.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        if sp.max_hits > 0 && hits.len() > sp.max_hits {
            hits.truncate(sp.max_hits);
        }
        let rows: Vec<(String, Vec<u8>, sabertooth::align::Alignment)> =
            hits.into_iter().map(|(_, n, s, a)| (n, s, a)).collect();
        total_rows += rows.len();
        out.push_str(&sabertooth::msa::build_a3m(&q.id, &q.num, &rows));
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

// ---------------------------------------------------------------------------
// translatesearch (six-frame translated nucleotide search)
// ---------------------------------------------------------------------------

fn run_translatesearch(args: &Args) -> Result<(), String> {
    if args.positionals.len() < 2 {
        return Err(
            "usage: sabertooth translatesearch <nucl_query.fasta> <protein_target.fasta> [options]"
                .into(),
        );
    }
    let m = SubstitutionMatrix::blosum62();
    let t = &m.aa2num;
    let nt_query = SeqDb::read_fasta(&args.positionals[0], t).map_err(|e| e.to_string())?;
    let target = SeqDb::read_fasta(&args.positionals[1], t).map_err(|e| e.to_string())?;
    if nt_query.is_empty() {
        return Err("nucleotide query FASTA is empty".into());
    }
    // sanity: warn if the "nucleotide" input doesn't look like DNA/RNA
    if let Some(first) = nt_query.records.first() {
        if !sabertooth::nucl::looks_like_nucleotide(&first.seq) {
            eprintln!("translatesearch: warning — query does not look like nucleotide sequence");
        }
    }
    let min_aa = args.opt_usize("min-orf", 20);
    let (orf_db, _meta) = sabertooth::nucl::orf_query_db(&nt_query, t, min_aa);
    if orf_db.records.is_empty() {
        return Err(format!("no ORFs of length >= {} found in the query", min_aa));
    }

    let mut sp = search_params(args, false);
    let idx = build_prefilter_index(args, &target, &mut sp)?;
    let ev = EValueParams::new(&m, target.total_residues);
    let start = Instant::now();
    let hits = search_sequences(&orf_db, &target, &idx, &m, &ev, sp);
    let elapsed = start.elapsed();

    let mut out = String::new();
    out.push_str(M8_HEADER);
    out.push('\n');
    for h in &hits {
        out.push_str(&h.to_m8());
        out.push('\n');
    }
    write_out(args.opt_str("out"), &out)?;
    eprintln!(
        "translatesearch: {} nucleotide queries -> {} ORFs (>= {} aa) x {} targets, {} hits in {:.2?}",
        nt_query.len(),
        orf_db.records.len(),
        min_aa,
        target.len(),
        hits.len(),
        elapsed
    );
    eprintln!("  (query IDs encode source_frame_ntstart, e.g. contig1_f-2_57)");
    Ok(())
}

// ---------------------------------------------------------------------------
// profileprofile (profile-vs-profile search)
// ---------------------------------------------------------------------------

fn run_profileprofile(args: &Args) -> Result<(), String> {
    if args.positionals.len() < 2 {
        return Err(
            "usage: sabertooth profileprofile <query.msa> <target.msa> [more_target.msa ...]".into(),
        );
    }
    let m = SubstitutionMatrix::blosum62();
    let query_data = read_to_string(&args.positionals[0]).map_err(|e| e.to_string())?;
    let qprof = Profile::from_msa_str(&query_data, &m, pc_params(args))?;
    let pback: Vec<f64> = (0..sabertooth::alphabet::PROFILE_AA_SIZE)
        .map(|i| m.p_back[i])
        .collect();
    let gaps = GapCosts::profile_default();
    let start = Instant::now();

    let mut results: Vec<(i32, String, sabertooth::align::Alignment)> = Vec::new();
    for path in &args.positionals[1..] {
        let data = read_to_string(path).map_err(|e| e.to_string())?;
        let tprof = Profile::from_msa_str(&data, &m, pc_params(args))?;
        let name = std::path::Path::new(path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or(path)
            .to_string();
        if let Some(aln) = sabertooth::align::align_profile_profile(&qprof, &tprof, gaps, &pback) {
            results.push((aln.raw_score, name, aln));
        }
    }
    results.sort_by(|a, b| b.0.cmp(&a.0));

    // profile–profile significance needs its own calibration; we report the raw
    // co-emission score and its bit-equivalent (÷8, the PSSM scale), plus the
    // consensus-identity and aligned ranges. E-values are intentionally omitted.
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

// ---------------------------------------------------------------------------
// taxonomy (LCA assignment)
// ---------------------------------------------------------------------------

fn run_taxonomy(args: &Args) -> Result<(), String> {
    if args.positionals.len() < 2 {
        return Err("usage: sabertooth taxonomy <query.fasta> <target.fasta> --nodes nodes.dmp --seqmap seqid2taxid.tsv [--names names.dmp] [--lca-frac f]".into());
    }
    let nodes_path = args
        .opt_str("nodes")
        .ok_or("taxonomy requires --nodes <nodes.dmp>")?;
    let seqmap_path = args
        .opt_str("seqmap")
        .ok_or("taxonomy requires --seqmap <seqid2taxid.tsv>")?;

    let m = SubstitutionMatrix::blosum62();
    let t = &m.aa2num;
    let query = SeqDb::read_fasta(&args.positionals[0], t).map_err(|e| e.to_string())?;
    let target = SeqDb::read_fasta(&args.positionals[1], t).map_err(|e| e.to_string())?;
    if query.is_empty() {
        return Err("query FASTA is empty".into());
    }

    let mut sp = search_params(args, false);
    let idx = build_prefilter_index(args, &target, &mut sp)?;
    let ev = EValueParams::new(&m, target.total_residues);
    let hits = search_sequences(&query, &target, &idx, &m, &ev, sp);

    let mut tree =
        sabertooth::taxonomy::TaxTree::parse_nodes(&read_to_string(nodes_path).map_err(|e| e.to_string())?);
    if let Some(names) = args.opt_str("names") {
        tree.load_names(&read_to_string(names).map_err(|e| e.to_string())?);
    }
    let seqmap = sabertooth::taxonomy::parse_seqmap(&read_to_string(seqmap_path).map_err(|e| e.to_string())?);
    let frac = args.opt_f64("lca-frac", 0.0);

    // group hits per query as (taxid, bit_score) via the seq->taxid map
    let mut by_query: std::collections::HashMap<String, Vec<(u32, f64)>> =
        std::collections::HashMap::new();
    let mut mapped = 0usize;
    for h in &hits {
        if let Some(&taxid) = seqmap.get(&h.target_id) {
            mapped += 1;
            by_query
                .entry(h.query_id.clone())
                .or_default()
                .push((taxid, h.bit_score));
        }
    }

    let mut out = String::new();
    out.push_str("query\ttaxid\trank\tname\tassigned_hits\n");
    let mut classified = 0usize;
    for q in &query.records {
        let assignment = by_query.get(&q.id).and_then(|scored| {
            let taxa = sabertooth::taxonomy::weighted_taxa(scored, frac);
            tree.lca(&taxa).map(|taxid| (taxid, taxa.len()))
        });
        match assignment {
            Some((taxid, n)) => {
                classified += 1;
                out.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{}\n",
                    q.id,
                    taxid,
                    tree.rank_of(taxid),
                    tree.name_of(taxid),
                    n
                ));
            }
            None => out.push_str(&format!("{}\t0\tno rank\tunclassified\t0\n", q.id)),
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
// makedb / convert2fasta (MMseqs2 on-disk DB format)
// ---------------------------------------------------------------------------

fn run_makedb(args: &Args) -> Result<(), String> {
    if args.positionals.len() < 2 {
        return Err("usage: sabertooth makedb <in.fasta> <db_prefix>".into());
    }
    let m = SubstitutionMatrix::blosum62();
    let db = SeqDb::read_fasta(&args.positionals[0], &m.aa2num).map_err(|e| e.to_string())?;
    if db.is_empty() {
        return Err("input FASTA is empty".into());
    }
    let nucleotide = db
        .records
        .first()
        .map(|r| sabertooth::nucl::looks_like_nucleotide(&r.seq))
        .unwrap_or(false);
    sabertooth::mmdb::write_db(
        std::path::Path::new(&args.positionals[1]),
        &db,
        nucleotide,
        &args.positionals[0],
    )
    .map_err(|e| e.to_string())?;
    eprintln!(
        "makedb: wrote {} sequences to {}{{,.index,.dbtype,_h,_h.index,.lookup,.source}} ({})",
        db.len(),
        args.positionals[1],
        if nucleotide { "nucleotide" } else { "amino acid" }
    );
    Ok(())
}

fn run_convert2fasta(args: &Args) -> Result<(), String> {
    if args.positionals.is_empty() {
        return Err("usage: sabertooth convert2fasta <db_prefix> [out.fasta]".into());
    }
    let m = SubstitutionMatrix::blosum62();
    let db = sabertooth::mmdb::read_db(std::path::Path::new(&args.positionals[0]), &m.aa2num)
        .map_err(|e| e.to_string())?;
    let mut out = String::new();
    for rec in &db.records {
        out.push('>');
        out.push_str(&rec.header);
        out.push('\n');
        out.push_str(&String::from_utf8_lossy(&rec.seq));
        out.push('\n');
    }
    // second positional or --out selects the destination
    let dest = args
        .positionals
        .get(1)
        .map(|s| s.as_str())
        .or_else(|| args.opt_str("out"));
    write_out(dest, &out)?;
    eprintln!("convert2fasta: {} sequences read from {}", db.len(), args.positionals[0]);
    Ok(())
}

// ---------------------------------------------------------------------------
// easysearch (database-free exhaustive SW search)
// ---------------------------------------------------------------------------

fn run_easysearch(args: &Args) -> Result<(), String> {
    if args.positionals.len() < 2 {
        return Err("usage: sabertooth easysearch <query.fasta> <target.fasta> [options]".into());
    }
    let m = SubstitutionMatrix::blosum62();
    let t = &m.aa2num;
    let query = SeqDb::read_fasta(&args.positionals[0], t).map_err(|e| e.to_string())?;
    let target = SeqDb::read_fasta(&args.positionals[1], t).map_err(|e| e.to_string())?;
    if query.is_empty() || target.is_empty() {
        return Err("query and target FASTA must both be non-empty".into());
    }
    // no prefilter index is built — that is the point of this mode
    let sp = search_params(args, false);
    let ev = EValueParams::new(&m, target.total_residues);
    let start = Instant::now();
    let hits = search_sequences_exhaustive(&query, &target, &m, &ev, sp);
    let elapsed = start.elapsed();

    let mut out = String::new();
    out.push_str(M8_HEADER);
    out.push('\n');
    for h in &hits {
        out.push_str(&h.to_m8());
        out.push('\n');
    }
    write_out(args.opt_str("out"), &out)?;
    eprintln!(
        "easysearch (database-free, no prefilter index): {} queries x {} targets = {} full SW alignments, {} hits in {:.2?}",
        query.len(),
        target.len(),
        query.len() * target.len(),
        hits.len(),
        elapsed
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// cluster (linclust-style clustering)
// ---------------------------------------------------------------------------

fn run_cluster(args: &Args) -> Result<(), String> {
    if args.positionals.is_empty() {
        return Err("usage: sabertooth cluster <in.fasta> [--min-seq-id f] [--min-cov f] [--k n]".into());
    }
    let m = SubstitutionMatrix::blosum62();
    let db = SeqDb::read_fasta(&args.positionals[0], &m.aa2num).map_err(|e| e.to_string())?;
    if db.is_empty() {
        return Err("input FASTA is empty".into());
    }
    let mut p = sabertooth::cluster::ClusterParams::default();
    p.k = args.opt_usize("k", p.k);
    p.num_minimizers = args.opt_usize("min-count", p.num_minimizers);
    p.min_seq_id = args.opt_f64("min-seq-id", p.min_seq_id);
    p.min_cov = args.opt_f64("min-cov", p.min_cov);
    let start = Instant::now();
    let clusters = sabertooth::cluster::linclust(&db, &m, p);
    let elapsed = start.elapsed();

    // MMseqs2 createtsv format: representative<TAB>member, one line per member.
    let mut out = String::new();
    for c in &clusters {
        let rep = &db.records[c.representative].id;
        for &mem in &c.members {
            out.push_str(rep);
            out.push('\t');
            out.push_str(&db.records[mem].id);
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
// calibrate (Monte-Carlo gapped statistics)
// ---------------------------------------------------------------------------

fn run_calibrate(args: &Args) -> Result<(), String> {
    let m = SubstitutionMatrix::blosum62();
    let seq_len = args.opt_usize("length", 250);
    let num_pairs = args.opt_usize("pairs", 3000);
    let seed = args.opt_usize("seed", 0x5EED) as u64;
    let gaps = GapCosts::sequence_default();

    let start = Instant::now();
    let g = sabertooth::evalue::calibrate_gapped(&m, gaps, seq_len, num_pairs, seed);
    let elapsed = start.elapsed();

    let ungapped = EValueParams::new(&m, 1).lambda;
    // Effect on a sample E-value: score 60 (×2 scale) against a 250-residue query
    // and a 100 Mbp search space, ungapped-λ vs calibrated-λ.
    let raw = 60i32;
    let (m_len, search_space) = (seq_len as f64, 1.0e8f64);
    let e_ungapped = search_space * (-ungapped * raw as f64).exp();
    let e_gapped = g.k * m_len * search_space * (-g.lambda * raw as f64).exp();

    println!("Sabertooth gapped-statistics calibration");
    println!("----------------------------------------");
    println!(
        "  simulated {} random pairs of length {} (gap open {}, extend {}) in {:.2?}",
        g.samples, seq_len, gaps.open, gaps.extend, elapsed
    );
    println!(
        "  optimal local score: mean {:.2}, std {:.2}",
        g.mean_score, g.std_score
    );
    println!();
    println!("  analytic ungapped   lambda = {:.4}", ungapped);
    println!(
        "  calibrated gapped    lambda = {:.4}   K = {:.4}   (mu = {:.2})",
        g.lambda, g.k, g.mu
    );
    println!(
        "  -> gaps lower lambda by {:.1}% (wider score distribution)",
        100.0 * (ungapped - g.lambda) / ungapped
    );
    println!();
    println!(
        "  sample E-value for raw score {} (len {}, {:.0e} search space):",
        raw, seq_len, search_space
    );
    println!("      ungapped-lambda estimate : {:.2e}", e_ungapped);
    println!("      calibrated-gapped estimate: {:.2e}", e_gapped);
    Ok(())
}

// ---------------------------------------------------------------------------
// output helper
// ---------------------------------------------------------------------------

fn write_out(path: Option<&str>, content: &str) -> Result<(), String> {
    match path {
        Some(p) => {
            let mut f = File::create(p).map_err(|e| e.to_string())?;
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
