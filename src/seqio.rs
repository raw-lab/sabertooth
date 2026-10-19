//! Format-aware sequence input shared by every command.
//!
//! A sequence argument can be
//!
//! * a FASTA file (plain or **gzip**, detected from the magic bytes — streamed
//!   through a decompressor, never fully inflated in memory),
//! * `-` for FASTA on standard input,
//! * an **MMseqs2 database prefix** (`prefix` + `prefix.index` + `prefix.dbtype`,
//!   optional `prefix_h` headers) as written by `mmseqs createdb` or
//!   `sabertooth makedb`,
//! * a Sabertooth **seed index** (`.stidx`, protein): its sequences and
//!   identifiers (the index keeps identifiers only, not header descriptions).
//!
//! Detection order is: `<path>.dbtype` exists → database prefix; otherwise the
//! file's first bytes decide (index magic → seed index, gzip magic → gzip
//! FASTA, `>` → FASTA). Anything else
//! is rejected with an explanation instead of being parsed as an empty FASTA —
//! the failure mode that made `makedb → search` silently report zero targets.
//! bzip2 / xz / zstd inputs are recognised and rejected explicitly.

use crate::fasta::{FastaReader, LoadOptions, SeqDb};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

/// Detected on-disk format of a sequence argument.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InputFormat {
    /// FASTA (`gz` = gzip-compressed).
    Fasta { gz: bool },
    /// MMseqs2-style database prefix.
    MmseqsDb,
    /// Sabertooth seed index file (`.stidx`).
    SeedIndex,
}

impl std::fmt::Display for InputFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            InputFormat::Fasta { gz: false } => write!(f, "FASTA"),
            InputFormat::Fasta { gz: true } => write!(f, "gzip FASTA"),
            InputFormat::MmseqsDb => write!(f, "MMseqs2 database"),
            InputFormat::SeedIndex => write!(f, "seed index"),
        }
    }
}

/// What kind of residues a command expects.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SeqKind {
    Protein,
    Nucleotide,
    Any,
}

fn err(kind: io::ErrorKind, msg: String) -> io::Error {
    io::Error::new(kind, msg)
}

/// Is `path` an MMseqs2-style database prefix (has `<path>.dbtype`)?
pub fn is_db_prefix(path: &str) -> bool {
    Path::new(&format!("{}.dbtype", path)).is_file()
}

/// Detect the format of a sequence argument without consuming it.
pub fn detect_format(path: &str) -> io::Result<InputFormat> {
    if path == "-" {
        return Ok(InputFormat::Fasta { gz: false });
    }
    if is_db_prefix(path) {
        return Ok(InputFormat::MmseqsDb);
    }
    let p = Path::new(path);
    if !p.exists() {
        if Path::new(&format!("{}.index", path)).is_file() {
            return Err(err(
                io::ErrorKind::InvalidData,
                format!(
                    "{}: looks like a database prefix ({}.index exists) but {}.dbtype is missing",
                    path, path, path
                ),
            ));
        }
        return Err(err(
            io::ErrorKind::NotFound,
            format!("{}: no such file or database prefix", path),
        ));
    }
    if p.is_dir() {
        return Err(err(
            io::ErrorKind::InvalidInput,
            format!("{}: is a directory, expected a FASTA file or database prefix", path),
        ));
    }
    if crate::index::is_index_file(path) {
        return Ok(InputFormat::SeedIndex);
    }
    let mut f = File::open(p).map_err(|e| err(e.kind(), format!("{}: {}", path, e)))?;
    let mut magic = [0u8; 6];
    let mut n = 0;
    while n < magic.len() {
        let k = f.read(&mut magic[n..])?;
        if k == 0 {
            break;
        }
        n += k;
    }
    classify_magic(path, &magic[..n])
}

fn classify_magic(path: &str, m: &[u8]) -> io::Result<InputFormat> {
    if m.len() >= 2 && m[0] == 0x1f && m[1] == 0x8b {
        return Ok(InputFormat::Fasta { gz: true });
    }
    if m.len() >= 3 && &m[..3] == b"BZh" {
        return Err(unsupported_codec(path, "bzip2"));
    }
    if m.len() >= 6 && m[..6] == [0xfd, b'7', b'z', b'X', b'Z', 0x00] {
        return Err(unsupported_codec(path, "xz"));
    }
    if m.len() >= 4 && m[..4] == [0x28, 0xb5, 0x2f, 0xfd] {
        return Err(unsupported_codec(path, "zstd"));
    }
    // Plain text: the FASTA reader validates the content line by line (leading
    // blank lines and `;` comments are allowed, so do not insist on '>' here).
    if m.is_empty() {
        return Err(err(
            io::ErrorKind::InvalidData,
            format!("{}: input file is empty (0 bytes)", path),
        ));
    }
    Ok(InputFormat::Fasta { gz: false })
}

fn unsupported_codec(path: &str, codec: &str) -> io::Error {
    err(
        io::ErrorKind::InvalidData,
        format!(
            "{}: {}-compressed input is not supported; decompress it first or recompress with gzip",
            path, codec
        ),
    )
}

/// A boxed buffered reader (plain file, gzip stream, or stdin).
pub type DynReader = Box<dyn BufRead + Send>;

/// Open a FASTA (plain or gzip, or `-` for stdin) as a streaming reader.
pub fn open_fasta(path: &str) -> io::Result<FastaReader<DynReader>> {
    if path == "-" {
        let stdin = io::stdin();
        let mut br = BufReader::with_capacity(1 << 16, stdin);
        let gz = {
            let buf = br.fill_buf()?;
            buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b
        };
        let r: DynReader = if gz {
            Box::new(BufReader::with_capacity(
                1 << 16,
                flate2::bufread::MultiGzDecoder::new(br),
            ))
        } else {
            Box::new(br)
        };
        return Ok(FastaReader::new(r, "<stdin>"));
    }
    let fmt = detect_format(path)?;
    let gz = match fmt {
        InputFormat::Fasta { gz } => gz,
        InputFormat::MmseqsDb => {
            return Err(err(
                io::ErrorKind::InvalidInput,
                format!("{}: is a database prefix, not a FASTA file", path),
            ))
        }
        InputFormat::SeedIndex => {
            return Err(err(
                io::ErrorKind::InvalidInput,
                format!("{}: is a seed index, not a FASTA file", path),
            ))
        }
    };
    let f = File::open(path).map_err(|e| err(e.kind(), format!("{}: {}", path, e)))?;
    let r: DynReader = if gz {
        Box::new(BufReader::with_capacity(
            1 << 16,
            GzErrorContext {
                inner: flate2::read::MultiGzDecoder::new(BufReader::with_capacity(1 << 16, f)),
                path: path.to_string(),
            },
        ))
    } else {
        Box::new(BufReader::with_capacity(1 << 16, f))
    };
    Ok(FastaReader::new(r, path))
}

/// Adds the file name to gzip decoding errors (corrupt / truncated archives).
struct GzErrorContext<R: Read> {
    inner: R,
    path: String,
}

impl<R: Read> Read for GzErrorContext<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.inner.read(buf).map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("{}: gzip decompression failed: {}", self.path, e),
            )
        })
    }
}

/// Sequences loaded by [`load_sequences`] plus provenance and warnings.
pub struct Loaded {
    pub db: SeqDb,
    pub format: InputFormat,
    /// Non-fatal observations for the user (printed by the CLI).
    pub warnings: Vec<String>,
}

/// Load a sequence argument (FASTA, gzip FASTA, stdin, database prefix or
/// seed index).
///
/// `expect` validates the residue type: a protein command given a nucleotide
/// database is an error; a FASTA whose first record looks like DNA produces a
/// warning (FASTA carries no type information).
pub fn load_sequences(
    path: &str,
    opts: LoadOptions,
    aa2num: &[u8; 256],
    expect: SeqKind,
) -> io::Result<Loaded> {
    let format = detect_format(path)?;
    let mut warnings = Vec::new();
    let db = match format {
        InputFormat::MmseqsDb => {
            let (db, info) = crate::mmdb::read_db_with(Path::new(path), aa2num, opts, expect)?;
            if info.missing_headers > 0 {
                warnings.push(format!(
                    "{}: {} entries have no header in {}_h; their database key is used as the identifier",
                    path, info.missing_headers, path
                ));
            }
            db
        }
        InputFormat::SeedIndex => {
            if expect == SeqKind::Nucleotide {
                return Err(err(
                    io::ErrorKind::InvalidInput,
                    format!("{}: a seed index holds protein sequences, but nucleotide input is expected", path),
                ));
            }
            let idx = crate::index::SeedIndex::open(Path::new(path))?;
            let mut db = SeqDb::new(opts);
            for i in 0..idx.len() {
                db.push_encoded(idx.id(i), None, idx.seq(i), None);
            }
            if opts.keep_headers {
                warnings.push(format!(
                    "{}: a seed index stores identifiers only; header descriptions are not available",
                    path
                ));
            }
            db
        }
        InputFormat::Fasta { .. } => {
            let mut fr = open_fasta(path)?;
            let db = fr.read_all(opts, aa2num)?;
            if fr.gap_chars_removed > 0 {
                warnings.push(format!(
                    "{}: removed {} alignment gap characters ('-'/'.') from sequences",
                    path, fr.gap_chars_removed
                ));
            }
            if expect == SeqKind::Protein && looks_nucleotide_db(&db) {
                warnings.push(format!(
                    "{}: sequences look like nucleotides but a protein input is expected \
                     (use translatesearch for nucleotide queries)",
                    path
                ));
            }
            db
        }
    };
    Ok(Loaded {
        db,
        format,
        warnings,
    })
}

/// Heuristic over the first few records: ≥ 90% A/C/G/T/U/N letters.
pub fn looks_nucleotide_db(db: &SeqDb) -> bool {
    let n = db.len().min(8);
    if n == 0 {
        return false;
    }
    (0..n).all(|i| crate::nucl::looks_like_nucleotide(&db.seq(i)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::build_aa2num;
    use std::io::Write;

    fn tmpdir(name: &str) -> std::path::PathBuf {
        let mut d = std::env::temp_dir();
        d.push(format!("sabertooth_seqio_{}_{}", std::process::id(), name));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn gzip_fasta_is_streamed() {
        let d = tmpdir("gz");
        let p = d.join("x.fa.gz");
        let f = File::create(&p).unwrap();
        let mut gz = flate2::write::GzEncoder::new(f, flate2::Compression::default());
        gz.write_all(b">a one\nMKVLLA\n>b\nGGGG\n").unwrap();
        gz.finish().unwrap();
        let ps = p.to_string_lossy().to_string();
        assert_eq!(detect_format(&ps).unwrap(), InputFormat::Fasta { gz: true });
        let l = load_sequences(&ps, LoadOptions::FULL, &build_aa2num(), SeqKind::Protein).unwrap();
        assert_eq!(l.db.len(), 2);
        assert_eq!(l.db.header(0), "a one");
        assert_eq!(&*l.db.seq(1), b"GGGG");
    }

    #[test]
    fn corrupt_gzip_is_an_error() {
        let d = tmpdir("gzbad");
        let p = d.join("bad.fa.gz");
        std::fs::write(&p, [0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 3, 1, 2, 3, 4, 5]).unwrap();
        let e = load_sequences(
            &p.to_string_lossy(),
            LoadOptions::FULL,
            &build_aa2num(),
            SeqKind::Any,
        )
        .err()
        .expect("corrupt gzip must fail");
        assert!(e.to_string().contains("gzip"), "{}", e);
    }

    #[test]
    fn other_codecs_and_empty_files_are_rejected() {
        let d = tmpdir("codec");
        for (name, bytes, needle) in [
            ("a.bz2", b"BZh91AY".to_vec(), "bzip2"),
            ("a.zst", vec![0x28, 0xb5, 0x2f, 0xfd, 0, 0], "zstd"),
            ("empty.fa", vec![], "empty"),
        ] {
            let p = d.join(name);
            std::fs::write(&p, &bytes).unwrap();
            let e = detect_format(&p.to_string_lossy()).unwrap_err();
            assert!(e.to_string().contains(needle), "{}: {}", name, e);
        }
        let e = detect_format(&d.join("nope").to_string_lossy()).unwrap_err();
        assert!(e.to_string().contains("no such file"), "{}", e);
    }
}
