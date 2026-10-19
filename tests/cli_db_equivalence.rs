//! Integration test: equivalent FASTA and database inputs must produce
//! equivalent hits, identifiers, coordinates and scores, for every command that
//! reads sequences. Runs the real CLI binary.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_sabertooth")
}

fn tmpdir(name: &str) -> PathBuf {
    let mut d = std::env::temp_dir();
    d.push(format!("sabertooth_it_{}_{}", std::process::id(), name));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

const AA: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

/// Families of mutated sequences (substitutions + indels) plus unrelated ones,
/// with descriptions, so search output is non-trivial.
fn synthetic_fasta(seed: u64, families: usize, members: usize) -> String {
    let mut rng = Rng(seed);
    let mut out = String::new();
    for f in 0..families {
        let len = 60 + rng.below(200);
        let anc: Vec<u8> = (0..len).map(|_| AA[rng.below(20)]).collect();
        for m in 0..members {
            let mut s = Vec::new();
            let mut i = 0;
            while i < anc.len() {
                let r = rng.below(100);
                if r < 3 {
                    i += 1 + rng.below(4);
                    continue;
                }
                if r < 6 {
                    for _ in 0..1 + rng.below(4) {
                        s.push(AA[rng.below(20)]);
                    }
                }
                if r < 30 {
                    s.push(AA[rng.below(20)]);
                } else {
                    s.push(anc[i]);
                }
                i += 1;
            }
            out.push_str(&format!(">fam{}_m{} family {} member {}\n", f, m, f, m));
            for chunk in s.chunks(60) {
                out.push_str(std::str::from_utf8(chunk).unwrap());
                out.push('\n');
            }
        }
    }
    for u in 0..5 {
        let len = 50 + rng.below(100);
        let s: String = (0..len).map(|_| AA[rng.below(20)] as char).collect();
        out.push_str(&format!(">unrelated{}\n{}\n", u, s));
    }
    out
}

fn run(args: &[&str]) -> (String, String) {
    let o = Command::new(bin()).args(args).output().expect("spawn sabertooth");
    let stdout = String::from_utf8_lossy(&o.stdout).to_string();
    let stderr = String::from_utf8_lossy(&o.stderr).to_string();
    assert!(
        o.status.success(),
        "sabertooth {:?} failed: {}\n{}",
        args,
        stderr,
        stdout
    );
    (stdout, stderr)
}

fn gzip(src: &Path, dst: &Path) {
    let data = std::fs::read(src).unwrap();
    let f = std::fs::File::create(dst).unwrap();
    let mut gz = flate2::write::GzEncoder::new(f, flate2::Compression::default());
    gz.write_all(&data).unwrap();
    gz.finish().unwrap();
}

#[test]
fn fasta_and_database_inputs_give_identical_results() {
    let d = tmpdir("equiv");
    let fa = d.join("t.fasta");
    std::fs::write(&fa, synthetic_fasta(0xC0FFEE, 6, 5)).unwrap();
    let gz = d.join("t.fasta.gz");
    gzip(&fa, &gz);
    let db = d.join("tDB");
    run(&["makedb", fa.to_str().unwrap(), db.to_str().unwrap()]);
    // a DB rebuilt from the DB (DB input path of makedb)
    let db2 = d.join("tDB2");
    run(&["makedb", db.to_str().unwrap(), db2.to_str().unwrap()]);

    let inputs = [
        fa.to_str().unwrap().to_string(),
        gz.to_str().unwrap().to_string(),
        db.to_str().unwrap().to_string(),
        db2.to_str().unwrap().to_string(),
    ];

    // sequence search: every query/target input combination
    let reference = run(&["search", &inputs[0], &inputs[0], "--evalue", "10"]).0;
    let lines = reference.lines().count();
    assert!(lines > 30, "test data should give many hits, got {}", lines - 1);
    for q in &inputs {
        for t in &inputs {
            let out = run(&["search", q, t, "--evalue", "10"]).0;
            assert_eq!(out, reference, "search {} vs {} differs", q, t);
        }
    }

    // other commands that read sequences (profile built from one query's MSA)
    let text = std::fs::read_to_string(&fa).unwrap();
    let first = format!(">{}", text.split('>').nth(1).unwrap());
    let q1 = d.join("q1.fasta");
    std::fs::write(&q1, first).unwrap();
    let msa = d.join("fam.a3m");
    run(&["result2msa", q1.to_str().unwrap(), &inputs[0], "--out", msa.to_str().unwrap()]);
    let cmds: Vec<Vec<String>> = vec![
        vec!["easysearch".into(), "{q}".into(), "{t}".into(), "--evalue".into(), "10".into()],
        vec!["result2msa".into(), "{q}".into(), "{t}".into()],
        vec!["align".into(), "{q}".into(), "{t}".into()],
        vec![
            "profilesearch".into(),
            msa.to_str().unwrap().into(),
            "{t}".into(),
            "--evalue".into(),
            "10".into(),
        ],
        vec!["cluster".into(), "{t}".into()],
        vec!["convert2fasta".into(), "{t}".into()],
    ];
    for c in &cmds {
        let render = |q: &str, t: &str| -> Vec<String> {
            c.iter()
                .map(|a| a.replace("{q}", q).replace("{t}", t))
                .collect()
        };
        let ref_args = render(&inputs[0], &inputs[0]);
        let ref_out = run(&ref_args.iter().map(|s| s.as_str()).collect::<Vec<_>>()).0;
        assert!(!ref_out.is_empty(), "{:?} produced nothing", ref_args);
        for i in &inputs[1..] {
            let a = render(i, i);
            let out = run(&a.iter().map(|s| s.as_str()).collect::<Vec<_>>()).0;
            assert_eq!(out, ref_out, "{:?} differs from FASTA input", a);
        }
    }
}

fn run_stdin(args: &[&str], input: &[u8]) -> String {
    use std::process::Stdio;
    let mut child = Command::new(bin())
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn sabertooth");
    child.stdin.take().unwrap().write_all(input).unwrap();
    let o = child.wait_with_output().unwrap();
    assert!(o.status.success(), "sabertooth {:?} (stdin) failed: {}", args, String::from_utf8_lossy(&o.stderr));
    String::from_utf8_lossy(&o.stdout).to_string()
}

/// Back-translate a protein record with one fixed codon per residue.
fn back_translate(protein: &str) -> String {
    protein
        .chars()
        .map(|c| match c {
            'A' => "GCT", 'C' => "TGT", 'D' => "GAT", 'E' => "GAA", 'F' => "TTT", 'G' => "GGT", 'H' => "CAT",
            'I' => "ATT", 'K' => "AAA", 'L' => "CTG", 'M' => "ATG", 'N' => "AAT", 'P' => "CCG", 'Q' => "CAG",
            'R' => "CGT", 'S' => "TCT", 'T' => "ACT", 'V' => "GTT", 'W' => "TGG", 'Y' => "TAT",
            _ => "NNN",
        })
        .collect()
}

/// The commands not covered above read targets through the same loader:
/// translated search, taxonomy, and queries streamed from stdin (plain and
/// gzip) give identical output for every target input kind.
#[test]
fn translated_taxonomy_and_stdin_inputs_give_identical_results() {
    let d = tmpdir("equiv2");
    let fa = d.join("t.fasta");
    let text = synthetic_fasta(0xFEED, 6, 5);
    std::fs::write(&fa, &text).unwrap();
    let gz = d.join("t.fasta.gz");
    gzip(&fa, &gz);
    let db = d.join("tDB");
    run(&["makedb", fa.to_str().unwrap(), db.to_str().unwrap()]);
    let dbi = d.join("tDBi");
    run(&["makedb", fa.to_str().unwrap(), dbi.to_str().unwrap()]);
    run(&["createindex", dbi.to_str().unwrap()]);
    let targets = [
        fa.to_str().unwrap().to_string(),
        gz.to_str().unwrap().to_string(),
        db.to_str().unwrap().to_string(),
        dbi.to_str().unwrap().to_string(),
        format!("{}.stidx", dbi.to_str().unwrap()),
    ];
    let recs: Vec<&str> = text.split('>').filter(|s| !s.is_empty()).collect();
    let seq = |r: &str| r.lines().skip(1).collect::<String>();

    // queries from stdin: plain and gzip FASTA
    let q = d.join("q.fasta");
    std::fs::write(&q, every_nth_record(&text, 4)).unwrap();
    let qgz = d.join("q.fasta.gz");
    gzip(&q, &qgz);
    let reference = run(&["search", q.to_str().unwrap(), &targets[0], "--evalue", "10"]).0;
    assert!(reference.lines().count() > 20, "too few hits");
    for (bytes, what) in [(std::fs::read(&q).unwrap(), "plain"), (std::fs::read(&qgz).unwrap(), "gzip")] {
        for t in &targets {
            assert_eq!(run_stdin(&["search", "-", t, "--evalue", "10"], &bytes), reference, "{} stdin query vs {}", what, t);
        }
    }

    // translated search: back-translated members of three families
    let nt = d.join("nt.fasta");
    let mut nts = String::new();
    for (i, r) in [recs[1], recs[7], recs[13]].iter().enumerate() {
        nts.push_str(&format!(">contig{}\nACG{}TAG\n", i, back_translate(&seq(r))));
    }
    std::fs::write(&nt, nts).unwrap();
    let tref = run(&["translatesearch", nt.to_str().unwrap(), &targets[0], "--evalue", "10"]).0;
    assert!(tref.lines().count() > 5, "translatesearch found too little: {}", tref);
    for t in &targets[1..] {
        assert_eq!(run(&["translatesearch", nt.to_str().unwrap(), t, "--evalue", "10"]).0, tref, "translatesearch vs {}", t);
    }

    // taxonomy: one genus per family under one root
    let nodes = d.join("nodes.dmp");
    let names = d.join("names.dmp");
    let seqmap = d.join("seqmap.tsv");
    let mut nd = String::from("1\t|\t1\t|\tno rank\t|\n");
    let mut nm = String::from("1\t|\troot\t|\t\t|\tscientific name\t|\n");
    for f in 0..6 {
        nd.push_str(&format!("{}\t|\t1\t|\tgenus\t|\n", 100 + f));
        nm.push_str(&format!("{}\t|\tGenus{}\t|\t\t|\tscientific name\t|\n", 100 + f, f));
    }
    let mut sm = String::new();
    for f in 0..6 {
        for m in 0..5 {
            sm.push_str(&format!("fam{}_m{}\t{}\n", f, m, 100 + f));
        }
    }
    std::fs::write(&nodes, nd).unwrap();
    std::fs::write(&names, nm).unwrap();
    std::fs::write(&seqmap, sm).unwrap();
    let tax = |t: &str| {
        run(&[
            "taxonomy", q.to_str().unwrap(), t, "--nodes", nodes.to_str().unwrap(), "--names", names.to_str().unwrap(),
            "--seqmap", seqmap.to_str().unwrap(), "--evalue", "1e-3",
        ])
        .0
    };
    let xref = tax(&targets[0]);
    assert!(xref.contains("\tgenus\tGenus"), "nothing classified to a genus: {}", xref);
    for t in &targets[1..] {
        assert_eq!(tax(t), xref, "taxonomy vs {}", t);
    }
}

/// A `.stidx` file is accepted wherever sequences are read (targets, queries,
/// makedb, cluster, convert2fasta), and a search never adopts an index built
/// with another seed or masking than its own: profile search (default k = 5)
/// gives the same hits whether or not a k = 6 index exists.
#[test]
fn seed_index_files_work_as_inputs_everywhere() {
    let d = tmpdir("stidx_in");
    let fa = d.join("t.fasta");
    let text = synthetic_fasta(0xABBA, 6, 5);
    std::fs::write(&fa, &text).unwrap();
    let db = d.join("tDB");
    run(&["makedb", fa.to_str().unwrap(), db.to_str().unwrap()]);
    let plain_db = d.join("plainDB");
    run(&["makedb", fa.to_str().unwrap(), plain_db.to_str().unwrap()]);
    run(&["createindex", db.to_str().unwrap()]);
    let stidx = format!("{}.stidx", db.to_str().unwrap());
    let f = fa.to_str().unwrap();
    let q = d.join("q.fasta");
    std::fs::write(&q, every_nth_record(&text, 4)).unwrap();
    let qs = q.to_str().unwrap();

    // targets: easysearch, align, search, profilesearch
    let first = format!(">{}", text.split('>').nth(1).unwrap());
    let q1 = d.join("q1.fasta");
    std::fs::write(&q1, first).unwrap();
    let msa = d.join("fam.a3m");
    run(&["result2msa", q1.to_str().unwrap(), f, "--out", msa.to_str().unwrap()]);
    for cmd in [
        vec!["easysearch", qs, "{t}", "--evalue", "10"],
        vec!["align", q1.to_str().unwrap(), "{t}"],
        vec!["search", qs, "{t}", "--evalue", "10", "--extra-cols"],
        vec!["profilesearch", msa.to_str().unwrap(), "{t}", "--evalue", "10"],
        vec!["result2msa", q1.to_str().unwrap(), "{t}"],
        vec!["cluster", "{t}"],
    ] {
        let with = |t: &str| -> Vec<String> { cmd.iter().map(|a| a.replace("{t}", t)).collect() };
        let reference = run(&with(f).iter().map(|s| s.as_str()).collect::<Vec<_>>()).0;
        assert!(reference.lines().count() > 3, "{:?} produced too little", cmd);
        for t in [plain_db.to_str().unwrap(), db.to_str().unwrap(), stidx.as_str()] {
            let a = with(t);
            let (out, _) = run(&a.iter().map(|s| s.as_str()).collect::<Vec<_>>());
            assert_eq!(out, reference, "{:?} with target {} differs from the FASTA target", cmd, t);
        }
    }
    // the profile search must not have adopted the k = 6 index silently
    let (_, err) = run(&["profilesearch", msa.to_str().unwrap(), &stidx, "--evalue", "10"]);
    assert!(err.contains("incompatible index"), "expected a warning about the k = 6 index: {}", err);

    // queries from an index; makedb, cluster and convert2fasta from an index
    let qdb = d.join("qDB");
    run(&["makedb", qs, qdb.to_str().unwrap()]);
    run(&["createindex", qdb.to_str().unwrap()]);
    let qidx = format!("{}.stidx", qdb.to_str().unwrap());
    assert_eq!(
        run(&["search", &qidx, f, "--evalue", "10"]).0,
        run(&["search", qs, f, "--evalue", "10"]).0
    );
    let db_from_idx = d.join("fromIdx");
    run(&["makedb", &stidx, db_from_idx.to_str().unwrap()]);
    assert_eq!(
        run(&["search", qs, db_from_idx.to_str().unwrap(), "--evalue", "10"]).0,
        run(&["search", qs, f, "--evalue", "10"]).0
    );
    // convert2fasta: identifiers only (an index keeps no descriptions)
    let ids_only: String = text
        .split('>')
        .filter(|r| !r.is_empty())
        .map(|r| {
            let mut lines = r.lines();
            let id = lines.next().unwrap().split_whitespace().next().unwrap().to_string();
            format!(">{}\n{}\n", id, lines.collect::<String>())
        })
        .collect();
    let (back, err) = run(&["convert2fasta", &stidx]);
    assert_eq!(back, ids_only);
    assert!(err.contains("identifiers only"), "{}", err);
}

#[test]
fn makedb_output_round_trips_headers_and_sequences() {
    let d = tmpdir("roundtrip");
    let fa = d.join("in.fasta");
    let text = synthetic_fasta(7, 2, 3);
    std::fs::write(&fa, &text).unwrap();
    let db = d.join("db");
    run(&["makedb", fa.to_str().unwrap(), db.to_str().unwrap()]);
    let back = run(&["convert2fasta", db.to_str().unwrap()]).0;
    // normalise the input to one sequence line per record
    let mut norm = String::new();
    let mut seq = String::new();
    for l in text.lines() {
        if let Some(h) = l.strip_prefix('>') {
            if !seq.is_empty() {
                norm.push_str(&seq);
                norm.push('\n');
                seq.clear();
            }
            norm.push('>');
            norm.push_str(h);
            norm.push('\n');
        } else {
            seq.push_str(l);
        }
    }
    norm.push_str(&seq);
    norm.push('\n');
    assert_eq!(back, norm);
}

#[test]
fn malformed_and_unsupported_inputs_fail_loudly() {
    let d = tmpdir("bad");
    let good = d.join("good.fasta");
    std::fs::write(&good, ">a\nMKVLLACDEFGHIK\n").unwrap();
    let g = good.to_str().unwrap();
    let cases: Vec<(&str, Vec<u8>, &str)> = vec![
        ("empty.fasta", b"".to_vec(), "empty"),
        ("blank.fasta", b"\n\n".to_vec(), "no FASTA records"),
        ("nohdr.fasta", b"MKVLLA\n".to_vec(), "before the first '>'"),
        ("digits.fasta", b">x\nMK12\n".to_vec(), "unexpected '1'"),
        ("emptyrec.fasta", b">x\n>y\nMK\n".to_vec(), "no residues"),
        ("bz.fasta.bz2", b"BZh91AY&SY".to_vec(), "bzip2"),
    ];
    for (name, bytes, needle) in cases {
        let p = d.join(name);
        std::fs::write(&p, bytes).unwrap();
        for args in [
            vec!["search", g, p.to_str().unwrap()],
            vec!["search", p.to_str().unwrap(), g],
            vec!["makedb", p.to_str().unwrap(), d.join("x").to_str().unwrap()],
        ] {
            let o = Command::new(bin()).args(&args).output().unwrap();
            let err = String::from_utf8_lossy(&o.stderr);
            assert!(!o.status.success(), "{:?} should fail", args);
            assert!(err.contains(needle), "{:?}: missing '{}' in: {}", args, needle, err);
        }
    }
    // a database prefix passed with an extension is not silently read as FASTA
    let db = d.join("db");
    run(&["makedb", g, db.to_str().unwrap()]);
    let idx = format!("{}.index", db.to_str().unwrap());
    let o = Command::new(bin()).args(["search", g, &idx]).output().unwrap();
    assert!(!o.status.success());
    // nucleotide database refused for protein search
    let nt = d.join("nt.fasta");
    std::fs::write(&nt, ">n\nACGTACGTACGTACGTACGTAAAACCCGGGTTT\n").unwrap();
    let ntdb = d.join("ntdb");
    run(&["makedb", nt.to_str().unwrap(), ntdb.to_str().unwrap()]);
    let o = Command::new(bin())
        .args(["search", g, ntdb.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(!o.status.success());
    assert!(String::from_utf8_lossy(&o.stderr).contains("nucleotide database"));
}

#[test]
fn duplicate_query_identifiers_stay_independent() {
    let d = tmpdir("dupq");
    let t = d.join("t.fasta");
    std::fs::write(&t, synthetic_fasta(99, 3, 3)).unwrap();
    // two different queries share the identifier "q"; output keeps input order
    let text = std::fs::read_to_string(&t).unwrap();
    let recs: Vec<&str> = text.split('>').filter(|s| !s.is_empty()).collect();
    let seq = |r: &str| r.lines().skip(1).collect::<String>();
    let q = d.join("q.fasta");
    std::fs::write(
        &q,
        format!(">q first\n{}\n>q second\n{}\n", seq(recs[6]), seq(recs[0])),
    )
    .unwrap();
    let out = run(&["search", q.to_str().unwrap(), t.to_str().unwrap(), "--evalue", "1e-3"]).0;
    let first_targets: Vec<&str> = out
        .lines()
        .skip(1)
        .map(|l| l.split('\t').nth(1).unwrap())
        .collect();
    // the first block (hits of the first "q") must be family 2, then family 0
    let fam2_end = first_targets.iter().position(|t| !t.starts_with("fam2")).unwrap();
    assert!(fam2_end >= 1);
    assert!(first_targets[fam2_end..].iter().all(|t| t.starts_with("fam0")),
        "hits of the two 'q' records were merged or reordered: {:?}", first_targets);
}

/// Every `step`-th record of a FASTA text (queries related to the targets).
fn every_nth_record(fasta: &str, step: usize) -> String {
    fasta
        .split('>')
        .filter(|r| !r.is_empty())
        .enumerate()
        .filter(|(i, _)| i % step == 0)
        .map(|(_, r)| format!(">{}", r))
        .collect()
}

fn run_env(args: &[&str], env: &[(&str, &str)]) -> String {
    let mut c = Command::new(bin());
    c.args(args);
    for (k, v) in env {
        c.env(k, v);
    }
    let o = c.output().expect("spawn sabertooth");
    assert!(o.status.success(), "sabertooth {:?} {:?} failed: {}", env, args, String::from_utf8_lossy(&o.stderr));
    String::from_utf8_lossy(&o.stdout).to_string()
}

/// Execution-mode parity: thread count, query batch size, memory limit, index
/// kind (persistent vs temporary) and SIMD kernel family change how the work
/// is scheduled, never what is reported.
#[test]
#[allow(clippy::type_complexity)]
fn execution_modes_give_byte_identical_output() {
    let d = tmpdir("modes");
    let fa = d.join("t.fasta");
    let text = synthetic_fasta(0xBEEF, 8, 6);
    std::fs::write(&fa, &text).unwrap();
    let q = d.join("q.fasta");
    std::fs::write(&q, every_nth_record(&text, 4)).unwrap();
    let db = d.join("tDB");
    run(&["makedb", fa.to_str().unwrap(), db.to_str().unwrap()]);
    run(&["createindex", db.to_str().unwrap()]);
    let base = ["search", q.to_str().unwrap(), db.to_str().unwrap(), "--evalue", "10", "--extra-cols"];
    let reference = run_env(&base, &[]);
    assert!(reference.lines().count() > 40, "too few hits: {}", reference.lines().count());
    let variants: Vec<(Vec<&str>, Vec<(&str, &str)>)> = vec![
        (vec!["--threads", "1"], vec![]),
        (vec!["--threads", "2", "--query-batch", "3"], vec![]),
        (vec!["--memory-limit", "64M"], vec![]),
        (vec!["--hit-block", "1024"], vec![]),
        (vec![], vec![("SABERTOOTH_SIMD", "sse2")]),
        (vec![], vec![("SABERTOOTH_SIMD", "avx2")]),
        (vec![], vec![("SABERTOOTH_SIMD", "scalar")]),
    ];
    for (extra, env) in &variants {
        let mut a: Vec<&str> = base.to_vec();
        a.extend(extra.iter());
        assert_eq!(run_env(&a, env), reference, "{:?} {:?} changed the output", extra, env);
    }
    // FASTA target (temporary index) gives the same records
    let fa_out = run_env(&["search", q.to_str().unwrap(), fa.to_str().unwrap(), "--evalue", "10", "--extra-cols"], &[]);
    assert_eq!(fa_out, reference);
}

#[test]
fn out_of_core_shards_give_byte_identical_output() {
    let d = tmpdir("ooc");
    let fa = d.join("t.fasta");
    let text = synthetic_fasta(0x00C5, 8, 6);
    std::fs::write(&fa, &text).unwrap();
    let q = d.join("q.fasta");
    std::fs::write(&q, every_nth_record(&text, 4)).unwrap();
    let db = d.join("tDB");
    run(&["makedb", fa.to_str().unwrap(), db.to_str().unwrap()]);
    run(&["createindex", db.to_str().unwrap()]);
    for extra in [&["--evalue", "10"][..], &["--evalue", "10", "--max-hits", "3", "--min-tcov", "0.2", "-s", "7.5"][..]] {
        let mut a = vec!["search", q.to_str().unwrap(), db.to_str().unwrap(), "--extra-cols"];
        a.extend_from_slice(extra);
        let reference = run_env(&a, &[]);
        assert!(reference.lines().count() > 20, "too few hits: {}", reference.lines().count());
        for shards in ["2", "3", "7"] {
            let mut b = a.clone();
            b.extend_from_slice(&["--shards", shards, "--query-batch", "5"]);
            assert_eq!(run_env(&b, &[]), reference, "{} local shards and {:?} differ", shards, extra);
            // second run reuses the persistent shards
            let (_, err) = run(&b);
            assert!(err.contains("reused"), "shards not reused: {}", err);
        }
    }
    // FASTA target: shards cut from the temporary index, under --tmp-dir
    let local = run_env(&["search", q.to_str().unwrap(), fa.to_str().unwrap(), "--evalue", "10"], &[]);
    let sharded = run_env(
        &["search", q.to_str().unwrap(), fa.to_str().unwrap(), "--evalue", "10", "--shards", "3", "--tmp-dir", d.to_str().unwrap()],
        &[],
    );
    assert_eq!(sharded, local);
    // stats JSON records the shards
    let sj = d.join("s.json");
    run(&["search", q.to_str().unwrap(), db.to_str().unwrap(), "--shards", "2", "--stats-json", sj.to_str().unwrap()]);
    let j = std::fs::read_to_string(&sj).unwrap();
    assert!(j.contains("\"shards\"") && j.contains("\"count\": 2"), "{}", j);
}

#[cfg(feature = "distributed")]
#[test]
fn distributed_cli_matches_local_cli() {
    let d = tmpdir("dist");
    let fa = d.join("t.fasta");
    let text = synthetic_fasta(0xD157, 8, 6);
    std::fs::write(&fa, &text).unwrap();
    let q = d.join("q.fasta");
    std::fs::write(&q, every_nth_record(&text, 4)).unwrap();
    let db = d.join("tDB");
    run(&["makedb", fa.to_str().unwrap(), db.to_str().unwrap()]);
    run(&["createindex", db.to_str().unwrap()]);
    for extra in [&["--evalue", "10"][..], &["--evalue", "10", "--max-hits", "5", "--min-qcov", "0.3", "-s", "7.5"][..]] {
        let mut a = vec!["search", q.to_str().unwrap(), db.to_str().unwrap(), "--extra-cols"];
        a.extend_from_slice(extra);
        let local = run_env(&a, &[]);
        assert!(local.lines().count() > 20, "too few hits: {}", local.lines().count());
        for shards in ["1", "3", "5"] {
            let mut b = a.clone();
            b.extend_from_slice(&["--distributed", "--shards", shards, "--hydra-quiet", "--hydra-port", "0", "--query-batch", "4"]);
            assert_eq!(run_env(&b, &[]), local, "distributed with {} shards and {:?} differs", shards, extra);
        }
    }
    // FASTA target: shards are cut from a temporary index and still match
    let local = run_env(&["search", q.to_str().unwrap(), fa.to_str().unwrap(), "--evalue", "10"], &[]);
    let dist = run_env(
        &["search", q.to_str().unwrap(), fa.to_str().unwrap(), "--evalue", "10", "--distributed", "--shards", "2", "--hydra-quiet",
          "--hydra-port", "0", "--tmp-dir", d.to_str().unwrap()],
        &[],
    );
    assert_eq!(dist, local);
}
