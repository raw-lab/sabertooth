//! MMseqs2 on-disk database format (validated read + streaming, atomic write).
//!
//! `mmseqs createdb` turns a FASTA into a small family of files. Verified byte
//! layout (from the real binary):
//!
//! * `<db>`         — data: each record is `payload + '\n' + '\0'`.
//! * `<db>.index`   — one `key<TAB>offset<TAB>length` line per record, where
//!   `length` **includes** the trailing `'\n'` and `'\0'`.
//! * `<db>.dbtype`  — 4-byte little-endian int. The low 16 bits are the main
//!   type (0 = amino acid, 1 = nucleotide, 2 = HMM profile, 12 =
//!   generic/header, ...). Bit 31 marks zstd-compressed entries; bits 17–30
//!   carry MMseqs2 "extended" flags (`(dbtype >> 16) & 0x7FFE`: 2 =
//!   index-needs-source, 4 = context pseudocounts, 8 = GPU-padded, 16 = set,
//!   32 = aux sequence) — layout taken from MMseqs2 `DBReader.h`
//!   (`getExtendedDbtype`) and `DBReader::isCompressed`.
//! * `<db>_h`       — header data (FASTA header without `>`), same framing.
//! * `<db>_h.index` — index for the header data.
//! * `<db>.lookup`  — `key<TAB>accession<TAB>set` (accession = first header token).
//! * `<db>.source`  — `0<TAB><source filename>`.
//!
//! ## Validation (reader)
//!
//! Everything a corrupt or unsupported database could do to produce silently
//! wrong results is checked and reported: missing `.dbtype`, unsupported main
//! types (profiles, alignment results, ...), compressed or GPU-padded data, split
//! data files (`<db>.0`, `<db>.1`, ...), malformed or duplicate index lines,
//! offset arithmetic overflow, entries running past the end of the data file,
//! and entries that are not NUL-terminated. Headers are matched to sequences
//! **by database key**, never by row position, so a header index in a different
//! order (or with gaps) still yields the right identifiers.
//!
//! ## Writing
//!
//! [`DbWriter`] streams records straight to disk (no in-memory copy of the whole
//! database) into temporary files that are renamed into place only when the
//! database is complete; `<db>.dbtype` is renamed last, so an interrupted
//! `makedb` never leaves something that is detected as a valid database.

use crate::fasta::{header_id, LoadOptions, SeqDb};
use crate::seqio::SeqKind;
use std::fs;
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

pub const DBTYPE_AMINO_ACIDS: i32 = 0;
pub const DBTYPE_NUCLEOTIDES: i32 = 1;
pub const DBTYPE_HMM_PROFILE: i32 = 2;
pub const DBTYPE_GENERIC: i32 = 12; // header databases use this

/// Extended dbtype flags (MMseqs2 `Parameters.h`), as returned in
/// [`DbType::extended`] (`(dbtype >> 16) & 0x7FFE`).
pub const DBTYPE_EXTENDED_INDEX_NEED_SRC: u32 = 2;
pub const DBTYPE_EXTENDED_CONTEXT_PSEUDO_COUNTS: u32 = 4;
pub const DBTYPE_EXTENDED_GPU: u32 = 8;
pub const DBTYPE_EXTENDED_SET: u32 = 16;
pub const DBTYPE_EXTENDED_AUX_SEQ: u32 = 32;

/// Human-readable MMseqs2 main dbtype name.
pub fn dbtype_name(main: i32) -> &'static str {
    match main {
        0 => "amino acid sequence",
        1 => "nucleotide sequence",
        2 => "HMM profile",
        3 => "profile state sequence",
        4 => "profile state profile",
        5 => "alignment result",
        6 => "clustering result",
        7 => "prefilter result",
        8 => "taxonomy result",
        9 => "MMseqs2 k-mer index",
        10 => "compressed a3m",
        11 => "MSA",
        12 => "generic",
        13 => "omitted",
        14 => "prefilter reverse result",
        15 => "offset",
        16 => "directory",
        17 => "flat file",
        18 => "sequence-taxonomy",
        19 => "stdin",
        20 => "URI",
        _ => "unknown",
    }
}

/// Decoded `.dbtype` word.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DbType {
    /// Main type (low 16 bits).
    pub main: i32,
    /// Extended flags, `(dbtype >> 16) & 0x7FFE`.
    pub extended: u32,
    /// Bit 31: entries are zstd-compressed.
    pub compressed: bool,
}

/// One entry of a `.index` file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IndexEntry {
    pub key: u32,
    pub offset: u64,
    pub length: u64,
}

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

fn with_path(e: io::Error, p: &str) -> io::Error {
    io::Error::new(e.kind(), format!("{}: {}", p, e))
}

fn write_dbtype(path: &Path, dbtype: i32) -> io::Result<()> {
    fs::write(path, dbtype.to_le_bytes())
}

/// Read and decode `<prefix>.dbtype`. A missing or short file is an error.
pub fn read_dbtype_file(path: &str) -> io::Result<DbType> {
    let b = fs::read(path).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!(
                "{}: {} (every MMseqs2 database has a 4-byte .dbtype file)",
                path, e
            ),
        )
    })?;
    if b.len() < 4 {
        return Err(invalid(format!(
            "{}: truncated dbtype file ({} bytes, expected 4)",
            path,
            b.len()
        )));
    }
    let raw = u32::from_le_bytes([b[0], b[1], b[2], b[3]]);
    Ok(DbType {
        main: (raw & 0xFFFF) as i32,
        extended: (raw >> 16) & 0x7FFE,
        compressed: raw & (1u32 << 31) != 0,
    })
}

/// Strictly parse a `.index` file: three tab-separated unsigned integers per
/// line, no duplicate keys, `offset + length` within `data_len`.
pub fn parse_index(text: &str, data_len: u64, what: &str) -> io::Result<Vec<IndexEntry>> {
    let mut out = Vec::new();
    for (ln, line) in text.lines().enumerate() {
        let lineno = ln + 1;
        if line.is_empty() {
            return Err(invalid(format!("{}:{}: empty index line", what, lineno)));
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 3 {
            return Err(invalid(format!(
                "{}:{}: expected 3 tab-separated fields (key, offset, length), found {}",
                what,
                lineno,
                fields.len()
            )));
        }
        let key: u32 = fields[0].parse().map_err(|_| {
            invalid(format!("{}:{}: bad key '{}'", what, lineno, fields[0]))
        })?;
        let offset: u64 = fields[1].parse().map_err(|_| {
            invalid(format!("{}:{}: bad offset '{}'", what, lineno, fields[1]))
        })?;
        let length: u64 = fields[2].parse().map_err(|_| {
            invalid(format!("{}:{}: bad length '{}'", what, lineno, fields[2]))
        })?;
        let end = offset.checked_add(length).ok_or_else(|| {
            invalid(format!(
                "{}:{}: offset {} + length {} overflows",
                what, lineno, offset, length
            ))
        })?;
        if end > data_len {
            return Err(invalid(format!(
                "{}:{}: entry key {} spans bytes {}..{} but the data file has only {} bytes \
                 (truncated or mismatched database files)",
                what, lineno, key, offset, end, data_len
            )));
        }
        out.push(IndexEntry { key, offset, length });
    }
    // duplicate keys
    let mut keys: Vec<u32> = out.iter().map(|e| e.key).collect();
    keys.sort_unstable();
    for w in keys.windows(2) {
        if w[0] == w[1] {
            return Err(invalid(format!("{}: duplicate key {}", what, w[0])));
        }
    }
    Ok(out)
}

/// Strip the trailing `'\0'` (required) and an optional `'\n'`/`"\r\n"` from an entry.
fn unframe<'a>(entry: &'a [u8], key: u32, what: &str) -> io::Result<&'a [u8]> {
    let mut end = entry.len();
    if end == 0 || entry[end - 1] != 0 {
        return Err(invalid(format!(
            "{}: entry key {} is not NUL-terminated (corrupt or not an MMseqs2 database)",
            what, key
        )));
    }
    end -= 1;
    if end > 0 && entry[end - 1] == b'\n' {
        end -= 1;
    }
    if end > 0 && entry[end - 1] == b'\r' {
        end -= 1;
    }
    Ok(&entry[..end])
}

/// Summary of a database read.
#[derive(Clone, Debug)]
pub struct DbReadInfo {
    pub dbtype: DbType,
    /// Entries whose key has no header record.
    pub missing_headers: usize,
    /// Database keys in record order.
    pub keys: Vec<u32>,
}

/// Read an MMseqs2 database at `<prefix>` into a `SeqDb` (everything retained).
pub fn read_db(prefix: &Path, aa2num: &[u8; 256]) -> io::Result<SeqDb> {
    read_db_with(prefix, aa2num, LoadOptions::FULL, SeqKind::Any).map(|(db, _)| db)
}

/// Read data entries in index order, calling `f(entry, payload)` for each.
fn for_each_entry<F>(data_path: &str, entries: &[IndexEntry], what: &str, mut f: F) -> io::Result<()>
where
    F: FnMut(&IndexEntry, &[u8]) -> io::Result<()>,
{
    let file = fs::File::open(data_path).map_err(|e| with_path(e, data_path))?;
    let mut r = BufReader::with_capacity(1 << 20, file);
    let mut pos: u64 = 0;
    let mut buf: Vec<u8> = Vec::new();
    for e in entries {
        if e.offset != pos {
            r.seek(SeekFrom::Start(e.offset))
                .map_err(|x| with_path(x, data_path))?;
            pos = e.offset;
        }
        let len = usize::try_from(e.length)
            .map_err(|_| invalid(format!("{}: entry key {} too large", what, e.key)))?;
        buf.resize(len, 0);
        r.read_exact(&mut buf).map_err(|x| with_path(x, data_path))?;
        pos += e.length;
        let payload = unframe(&buf, e.key, what)?;
        f(e, payload)?;
    }
    Ok(())
}

/// Stream an MMseqs2 database record by record (validated as in
/// [`read_db_with`]): `f(key, header, letters)` per entry in index order, where
/// `header` is the full header if `keep_headers` else only its identifier.
/// Memory is O(headers index + longest record).
pub fn stream_db<F>(prefix: &Path, keep_headers: bool, expect: SeqKind, mut f: F) -> io::Result<DbReadInfo>
where
    F: FnMut(u32, &[u8], &[u8]) -> io::Result<()>,
{
    let pstr = prefix.to_string_lossy().to_string();

    // Split data files (`<db>.0`, `<db>.1`, ...) take precedence over `<db>` in
    // MMseqs2 (`FileUtil::findDatafiles`); we only read single-file databases.
    if Path::new(&format!("{}.0", pstr)).exists() {
        return Err(invalid(format!(
            "{}: split database (data in {}.0, {}.1, ...) is not supported; \
             export it with `mmseqs convert2fasta {} out.fasta` and rebuild with `sabertooth makedb`",
            pstr, pstr, pstr, pstr
        )));
    }

    let dbt = read_dbtype_file(&format!("{}.dbtype", pstr))?;
    if dbt.compressed {
        return Err(invalid(format!(
            "{}: compressed MMseqs2 database (createdb --compressed 1) is not supported; \
             recreate it with `mmseqs createdb ... --compressed 0` or use the FASTA",
            pstr
        )));
    }
    if dbt.extended & DBTYPE_EXTENDED_GPU != 0 {
        return Err(invalid(format!(
            "{}: GPU-padded MMseqs2 database is not supported; recreate it without --gpu",
            pstr
        )));
    }
    match dbt.main {
        DBTYPE_AMINO_ACIDS => {
            if expect == SeqKind::Nucleotide {
                return Err(invalid(format!(
                    "{}: is an amino-acid database but a nucleotide input is expected",
                    pstr
                )));
            }
        }
        DBTYPE_NUCLEOTIDES => {
            if expect == SeqKind::Protein {
                return Err(invalid(format!(
                    "{}: is a nucleotide database but this command needs protein sequences \
                     (use translatesearch for nucleotide queries)",
                    pstr
                )));
            }
        }
        DBTYPE_HMM_PROFILE => {
            return Err(invalid(format!(
                "{}: is an MMseqs2 HMM-profile database; profile databases are not supported as \
                 sequence inputs (give profilesearch an MSA instead)",
                pstr
            )));
        }
        other => {
            return Err(invalid(format!(
                "{}: unsupported database type {} ({}); expected an amino-acid (0) or \
                 nucleotide (1) sequence database",
                pstr,
                other,
                dbtype_name(other)
            )));
        }
    }

    let data_len = fs::metadata(&pstr).map_err(|e| with_path(e, &pstr))?.len();
    let index_path = format!("{}.index", pstr);
    let index_text = fs::read_to_string(&index_path).map_err(|e| with_path(e, &index_path))?;
    let entries = parse_index(&index_text, data_len, &index_path)?;
    drop(index_text);
    if entries.is_empty() {
        return Err(invalid(format!("{}: database has no entries", pstr)));
    }
    if entries.len() as u64 > u32::MAX as u64 {
        return Err(invalid(format!("{}: too many entries", pstr)));
    }

    // Headers (optional), matched by key.
    let h_path = format!("{}_h", pstr);
    let hi_path = format!("{}_h.index", pstr);
    let (h_exists, hi_exists) = (Path::new(&h_path).is_file(), Path::new(&hi_path).is_file());
    if h_exists != hi_exists {
        return Err(invalid(format!(
            "{}: header database is incomplete ({} {} but {} {})",
            pstr,
            h_path,
            if h_exists { "exists" } else { "is missing" },
            hi_path,
            if hi_exists { "exists" } else { "is missing" }
        )));
    }
    // key -> header text (only first token unless headers requested)
    let mut header_of: Vec<(u32, Vec<u8>)> = Vec::new();
    if h_exists {
        let hlen = fs::metadata(&h_path).map_err(|e| with_path(e, &h_path))?.len();
        let htext = fs::read_to_string(&hi_path).map_err(|e| with_path(e, &hi_path))?;
        let hentries = parse_index(&htext, hlen, &hi_path)?;
        drop(htext);
        let mut sorted = hentries.clone();
        sorted.sort_by_key(|e| e.offset);
        header_of.reserve(sorted.len());
        for_each_entry(&h_path, &sorted, &h_path, |e, payload| {
            let text = if keep_headers {
                payload.to_vec()
            } else {
                header_id(payload).to_vec()
            };
            header_of.push((e.key, text));
            Ok(())
        })?;
        header_of.sort_unstable_by_key(|(k, _)| *k);
    }

    let mut missing = 0usize;
    let mut count = 0usize;
    let mut keys = Vec::with_capacity(entries.len());
    let mut letters: Vec<u8> = Vec::new();
    for_each_entry(&pstr, &entries, &pstr, |e, payload| {
        letters.clear();
        for (col, &b) in payload.iter().enumerate() {
            match b {
                b'A'..=b'Z' | b'*' => letters.push(b),
                b'a'..=b'z' => letters.push(b.to_ascii_uppercase()),
                b'\n' | b'\r' | b' ' | b'\t' => {}
                _ => {
                    return Err(invalid(format!(
                        "{}: entry key {} byte {}: unexpected byte 0x{:02x} in sequence data",
                        pstr, e.key, col, b
                    )))
                }
            }
        }
        if letters.last() == Some(&b'*') {
            letters.pop();
        }
        if letters.is_empty() {
            return Err(invalid(format!(
                "{}: entry key {} has an empty sequence",
                pstr, e.key
            )));
        }
        if letters.len() as u64 > u32::MAX as u64 {
            return Err(invalid(format!("{}: entry key {} is too long", pstr, e.key)));
        }
        let header: Option<&[u8]> = match header_of.binary_search_by_key(&e.key, |(k, _)| *k) {
            Ok(i) if !header_id(&header_of[i].1).is_empty() => Some(&header_of[i].1),
            _ => None,
        };
        match header {
            Some(h) => f(e.key, h, &letters)?,
            None => {
                if h_exists {
                    missing += 1;
                }
                let k = e.key.to_string();
                f(e.key, k.as_bytes(), &letters)?;
            }
        }
        keys.push(e.key);
        count += 1;
        Ok(())
    })?;
    if !h_exists {
        missing = count;
    }
    Ok(DbReadInfo {
        dbtype: dbt,
        missing_headers: missing,
        keys,
    })
}

/// Read an MMseqs2 database, validating type and structure. `expect` restricts
/// the accepted sequence type.
pub fn read_db_with(
    prefix: &Path,
    aa2num: &[u8; 256],
    opts: LoadOptions,
    expect: SeqKind,
) -> io::Result<(SeqDb, DbReadInfo)> {
    let mut db = SeqDb::new(opts);
    let info = stream_db(prefix, opts.keep_headers, expect, |_k, h, letters| {
        db.push_record(h, letters, aa2num);
        Ok(())
    })?;
    db.shrink_to_fit();
    Ok((db, info))
}

/// Streaming, atomic writer for an MMseqs2 sequence database.
pub struct DbWriter {
    prefix: String,
    dbtype: i32,
    data: BufWriter<fs::File>,
    index: BufWriter<fs::File>,
    hdata: BufWriter<fs::File>,
    hindex: BufWriter<fs::File>,
    lookup: BufWriter<fs::File>,
    offset: u64,
    hoffset: u64,
    key: u32,
    tmp: Vec<(PathBuf, PathBuf)>,
}

fn tmp_name(final_path: &str) -> String {
    format!("{}.tmp-{}", final_path, std::process::id())
}

impl DbWriter {
    /// Start writing a database at `prefix` (`nucleotide` selects the dbtype).
    pub fn create(prefix: &Path, nucleotide: bool) -> io::Result<DbWriter> {
        let pstr = prefix.to_string_lossy().to_string();
        let names = [
            pstr.clone(),
            format!("{}.index", pstr),
            format!("{}_h", pstr),
            format!("{}_h.index", pstr),
            format!("{}.lookup", pstr),
        ];
        let mut tmp = Vec::new();
        let mut files = Vec::new();
        for n in &names {
            let t = tmp_name(n);
            let f = fs::File::create(&t).map_err(|e| with_path(e, &t))?;
            files.push(BufWriter::with_capacity(1 << 20, f));
            tmp.push((PathBuf::from(t), PathBuf::from(n)));
        }
        let mut it = files.into_iter();
        Ok(DbWriter {
            prefix: pstr,
            dbtype: if nucleotide {
                DBTYPE_NUCLEOTIDES
            } else {
                DBTYPE_AMINO_ACIDS
            },
            data: it.next().unwrap(),
            index: it.next().unwrap(),
            hdata: it.next().unwrap(),
            hindex: it.next().unwrap(),
            lookup: it.next().unwrap(),
            offset: 0,
            hoffset: 0,
            key: 0,
            tmp,
        })
    }

    /// Append one record (`header` without `>`, `seq` residue letters).
    pub fn push(&mut self, header: &[u8], seq: &[u8]) -> io::Result<()> {
        let key = self.key;
        self.data.write_all(seq)?;
        self.data.write_all(b"\n\0")?;
        let len = seq.len() as u64 + 2;
        writeln!(self.index, "{}\t{}\t{}", key, self.offset, len)?;
        self.offset += len;

        self.hdata.write_all(header)?;
        self.hdata.write_all(b"\n\0")?;
        let hlen = header.len() as u64 + 2;
        writeln!(self.hindex, "{}\t{}\t{}", key, self.hoffset, hlen)?;
        self.hoffset += hlen;

        self.lookup.write_all(key.to_string().as_bytes())?;
        self.lookup.write_all(b"\t")?;
        self.lookup.write_all(header_id(header))?;
        self.lookup.write_all(b"\t0\n")?;
        self.key = self
            .key
            .checked_add(1)
            .ok_or_else(|| invalid("too many records for a u32-keyed database".into()))?;
        Ok(())
    }

    /// Records written so far.
    pub fn len(&self) -> usize {
        self.key as usize
    }

    pub fn is_empty(&self) -> bool {
        self.key == 0
    }

    /// Flush, fsync and rename everything into place. `.dbtype` goes last.
    pub fn finish(mut self, source: &str) -> io::Result<()> {
        for w in [
            &mut self.data,
            &mut self.index,
            &mut self.hdata,
            &mut self.hindex,
            &mut self.lookup,
        ] {
            w.flush()?;
            w.get_ref().sync_all()?;
        }
        let p = self.prefix.clone();
        // small files
        let extra = [
            (format!("{}_h.dbtype", p), DBTYPE_GENERIC.to_le_bytes().to_vec()),
            (
                format!("{}.source", p),
                format!("0\t{}\n", source).into_bytes(),
            ),
        ];
        for (name, bytes) in &extra {
            let t = tmp_name(name);
            fs::write(&t, bytes)?;
            self.tmp.push((PathBuf::from(t), PathBuf::from(name)));
        }
        // An existing database at this prefix stops being "a database" while
        // its files are swapped: remove the old dbtype first.
        let dbtype_final = format!("{}.dbtype", p);
        let _ = fs::remove_file(&dbtype_final);
        for (t, f) in &self.tmp {
            fs::rename(t, f).map_err(|e| with_path(e, &f.to_string_lossy()))?;
        }
        let t = tmp_name(&dbtype_final);
        write_dbtype(Path::new(&t), self.dbtype)?;
        fs::rename(&t, &dbtype_final)?;
        self.tmp.clear();
        Ok(())
    }
}

impl Drop for DbWriter {
    fn drop(&mut self) {
        // abandoned (error path): remove temporaries
        for (t, _) in &self.tmp {
            let _ = fs::remove_file(t);
        }
    }
}

/// Write a whole in-memory database (compatibility wrapper over [`DbWriter`]).
pub fn write_db(prefix: &Path, db: &SeqDb, nucleotide: bool, source: &str) -> io::Result<()> {
    let mut w = DbWriter::create(prefix, nucleotide)?;
    for r in db.iter() {
        w.push(r.header.as_bytes(), &r.seq())?;
    }
    w.finish(source)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::build_aa2num;
    use std::io::Cursor;

    fn tmp(name: &str) -> std::path::PathBuf {
        let mut d = std::env::temp_dir();
        d.push(format!("sabertooth_mmdb_test_{}_{}", std::process::id(), name));
        d
    }

    fn sample_db() -> SeqDb {
        let t = build_aa2num();
        SeqDb::from_reader(
            Cursor::new(">seq1 first protein\nMKVLLACDEF\n>seq2 second\nWYFPCDGHIK\n"),
            &t,
        )
        .unwrap()
    }

    #[test]
    fn write_then_read_roundtrips() {
        let t = build_aa2num();
        let db = sample_db();
        let prefix = tmp("rt");
        write_db(&prefix, &db, false, "in.fasta").unwrap();

        let back = read_db(&prefix, &t).unwrap();
        assert_eq!(back.len(), 2);
        assert_eq!(back.id(0), "seq1");
        assert_eq!(back.header(0), "seq1 first protein");
        assert_eq!(&*back.seq(0), b"MKVLLACDEF");
        assert_eq!(back.id(1), "seq2");
        assert_eq!(&*back.seq(1), b"WYFPCDGHIK");
        assert_eq!(back.num(0), db.num(0));
        // no temporaries left behind
        let dir = prefix.parent().unwrap();
        let stem = prefix.file_name().unwrap().to_string_lossy().to_string();
        for e in fs::read_dir(dir).unwrap() {
            let n = e.unwrap().file_name().to_string_lossy().to_string();
            assert!(!(n.starts_with(&stem) && n.contains(".tmp-")), "leftover {}", n);
        }
    }

    #[test]
    fn index_lengths_include_framing() {
        let t = build_aa2num();
        let db = SeqDb::from_reader(Cursor::new(">a\nMKVLLA\n"), &t).unwrap();
        let prefix = tmp("len");
        write_db(&prefix, &db, false, "x").unwrap();
        let idx = fs::read_to_string(format!("{}.index", prefix.to_string_lossy())).unwrap();
        assert_eq!(idx.lines().next().unwrap(), "0\t0\t8");
        let data = fs::read(prefix.to_string_lossy().to_string()).unwrap();
        assert_eq!(&data, b"MKVLLA\n\0");
    }

    #[test]
    fn dbtype_bytes_are_correct() {
        let t = build_aa2num();
        let db = SeqDb::from_reader(Cursor::new(">a\nMKVL\n"), &t).unwrap();
        let prefix = tmp("dbtype");
        write_db(&prefix, &db, false, "x").unwrap();
        let bytes = fs::read(format!("{}.dbtype", prefix.to_string_lossy())).unwrap();
        assert_eq!(bytes, vec![0, 0, 0, 0]);
        let hbytes = fs::read(format!("{}_h.dbtype", prefix.to_string_lossy())).unwrap();
        assert_eq!(hbytes, vec![12, 0, 0, 0]);
    }

    #[test]
    fn headers_are_matched_by_key_not_row() {
        // Reverse the header index row order and give the data index a gap in
        // keys: identifiers must still follow keys.
        let t = build_aa2num();
        let prefix = tmp("bykey");
        let p = prefix.to_string_lossy().to_string();
        fs::write(&p, b"MKV\n\0GGG\n\0").unwrap();
        fs::write(format!("{}.index", p), "5\t0\t5\n9\t5\t5\n").unwrap();
        fs::write(format!("{}.dbtype", p), 0i32.to_le_bytes()).unwrap();
        fs::write(format!("{}_h", p), b"nine desc\n\0five\n\0").unwrap();
        fs::write(format!("{}_h.index", p), "9\t0\t11\n5\t11\t6\n").unwrap();
        let back = read_db(&prefix, &t).unwrap();
        assert_eq!(back.id(0), "five");
        assert_eq!(&*back.seq(0), b"MKV");
        assert_eq!(back.id(1), "nine");
        assert_eq!(back.header(1), "nine desc");
    }

    fn write_raw(prefix: &Path, data: &[u8], index: &str, dbtype: i32) {
        let p = prefix.to_string_lossy().to_string();
        fs::write(&p, data).unwrap();
        fs::write(format!("{}.index", p), index).unwrap();
        fs::write(format!("{}.dbtype", p), dbtype.to_le_bytes()).unwrap();
    }

    #[test]
    fn rejects_corrupt_and_unsupported_databases() {
        let t = build_aa2num();
        let cases: Vec<(&str, Vec<u8>, &str, i32, &str)> = vec![
            ("past_end", b"MKV\n\0".to_vec(), "0\t0\t50\n", 0, "data file has only"),
            ("overflow", b"MKV\n\0".to_vec(), "0\t18446744073709551615\t5\n", 0, "overflows"),
            ("dup", b"MKV\n\0GG\n\0".to_vec(), "0\t0\t5\n0\t5\t4\n", 0, "duplicate key"),
            ("badline", b"MKV\n\0".to_vec(), "0\t0\n", 0, "3 tab-separated"),
            ("notnul", b"MKV\nX".to_vec(), "0\t0\t5\n", 0, "NUL-terminated"),
            // dbtype word written by `mmseqs createdb --compressed 1` (bytes 00 00 00 80)
            ("compressed", b"MKV\n\0".to_vec(), "0\t0\t5\n", i32::MIN, "compressed"),
            ("gpu", b"MKV\n\0".to_vec(), "0\t0\t5\n", 8 << 16, "GPU-padded"),
            ("profile", b"MKV\n\0".to_vec(), "0\t0\t5\n", 2, "HMM-profile"),
            ("result", b"MKV\n\0".to_vec(), "0\t0\t5\n", 5, "alignment result"),
            ("empty_seq", b"\n\0".to_vec(), "0\t0\t2\n", 0, "empty sequence"),
            ("no_entries", Vec::new(), "", 0, "no entries"),
        ];
        for (name, data, index, dbtype, needle) in cases {
            let prefix = tmp(&format!("bad_{}", name));
            write_raw(&prefix, &data, index, dbtype);
            let e = read_db(&prefix, &t).err().unwrap_or_else(|| panic!("{} accepted", name));
            assert!(e.to_string().contains(needle), "{}: {}", name, e);
        }
        // missing dbtype
        let prefix = tmp("nodbtype");
        write_raw(&prefix, b"MKV\n\0", "0\t0\t5\n", 0);
        fs::remove_file(format!("{}.dbtype", prefix.to_string_lossy())).unwrap();
        let e = read_db(&prefix, &t).unwrap_err();
        assert!(e.to_string().contains(".dbtype"), "{}", e);
        // split database
        let prefix = tmp("split");
        let p = prefix.to_string_lossy().to_string();
        fs::write(format!("{}.0", p), b"MKV\n\0").unwrap();
        fs::write(format!("{}.index", p), "0\t0\t5\n").unwrap();
        fs::write(format!("{}.dbtype", p), 0i32.to_le_bytes()).unwrap();
        let e = read_db(&prefix, &t).unwrap_err();
        assert!(e.to_string().contains("split database"), "{}", e);
    }

    #[test]
    fn nucleotide_db_rejected_for_protein_search() {
        let t = build_aa2num();
        let prefix = tmp("nucl");
        write_raw(&prefix, b"ACGT\n\0", "0\t0\t6\n", 1);
        let e = read_db_with(&prefix, &t, LoadOptions::SEARCH, SeqKind::Protein).unwrap_err();
        assert!(e.to_string().contains("nucleotide database"), "{}", e);
        assert!(read_db_with(&prefix, &t, LoadOptions::SEARCH, SeqKind::Nucleotide).is_ok());
    }
}
