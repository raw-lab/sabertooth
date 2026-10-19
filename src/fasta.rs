//! Compact sequence storage and validated, streaming FASTA parsing.
//!
//! ## Storage
//!
//! [`SeqDb`] keeps every sequence in **one contiguous arena** of numeric residue
//! codes (MMseqs2 order, `X` = 20) with a `u64` offset table, plus a string arena
//! of identifiers. The original ASCII letters and full FASTA headers are only
//! retained when a command needs them ([`LoadOptions`]); search never keeps a
//! second ASCII copy of the database. Compared with the previous
//! `Vec<Record{header, id, seq, num}>` layout (four heap allocations and two full
//! copies of every residue) this removes ~100 bytes of per-record overhead and
//! halves residue storage.
//!
//! ## Parsing
//!
//! [`FastaReader`] streams records one at a time from any `BufRead` and rejects
//! malformed input with actionable, line-numbered errors instead of silently
//! producing an empty database:
//!
//! * sequence data before the first `>` header,
//! * empty headers and records with no residues,
//! * characters that are not residue letters (digits, punctuation, binary data),
//! * sequences longer than `u32::MAX` residues (the seed index stores `u32`
//!   positions).
//!
//! Tolerated and normalised: lowercase letters (uppercased), Windows line endings,
//! whitespace inside sequence lines, `;` comment lines, alignment gap characters
//! `-`/`.` (removed and counted), and a terminal `*` stop symbol (removed); an
//! internal `*` is kept as a letter and encoded as `X`.

use crate::alphabet::{aa_to_num, num_to_aa};
use std::borrow::Cow;
use std::fs::File;
use std::io::{self, BufRead, Read, Write};
use std::path::Path;

/// What to retain when loading sequences.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LoadOptions {
    /// Keep the full FASTA header (identifier + description).
    pub keep_headers: bool,
    /// Keep the original residue letters (needed for nucleotide translation and
    /// for writing sequences back out; search only needs numeric codes).
    pub keep_letters: bool,
}

impl LoadOptions {
    /// Minimal footprint for search: numeric residues + identifiers only.
    pub const SEARCH: LoadOptions = LoadOptions {
        keep_headers: false,
        keep_letters: false,
    };
    /// Everything (format conversion, `createdb`, `makedb`, `convert2fasta`).
    pub const FULL: LoadOptions = LoadOptions {
        keep_headers: true,
        keep_letters: true,
    };
    /// Letters but no headers (translated search: the nucleotide letters are
    /// needed for six-frame translation).
    pub const LETTERS: LoadOptions = LoadOptions {
        keep_headers: false,
        keep_letters: true,
    };
}

/// Append-only string arena: one contiguous byte buffer plus `u64` offsets.
#[derive(Clone, Debug)]
pub struct StrArena {
    bytes: Vec<u8>,
    offsets: Vec<u64>, // len = n + 1
}

impl Default for StrArena {
    fn default() -> Self {
        StrArena {
            bytes: Vec::new(),
            offsets: vec![0],
        }
    }
}

impl StrArena {
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a string. Invalid UTF-8 is replaced (lossy) so `get` is infallible.
    pub fn push(&mut self, s: &[u8]) {
        match std::str::from_utf8(s) {
            Ok(_) => self.bytes.extend_from_slice(s),
            Err(_) => self
                .bytes
                .extend_from_slice(String::from_utf8_lossy(s).as_bytes()),
        }
        self.offsets.push(self.bytes.len() as u64);
    }

    pub fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn get(&self, i: usize) -> &str {
        let a = self.offsets[i] as usize;
        let b = self.offsets[i + 1] as usize;
        // SAFETY-free: `push` only ever stores valid UTF-8.
        std::str::from_utf8(&self.bytes[a..b]).unwrap_or("")
    }

    /// Heap bytes held (for memory accounting).
    pub fn heap_bytes(&self) -> usize {
        self.bytes.capacity() + self.offsets.capacity() * 8
    }

    fn shrink_to_fit(&mut self) {
        self.bytes.shrink_to_fit();
        self.offsets.shrink_to_fit();
    }
}

/// A borrowed view of one record.
#[derive(Clone, Copy)]
pub struct RecordRef<'a> {
    /// Ordinal of the record in its database (stable input order).
    pub index: usize,
    /// Identifier (first whitespace-delimited header token).
    pub id: &'a str,
    /// Full header if retained, otherwise the identifier.
    pub header: &'a str,
    /// Numeric residue codes.
    pub num: &'a [u8],
    letters: Option<&'a [u8]>,
}

impl<'a> RecordRef<'a> {
    pub fn len(&self) -> usize {
        self.num.len()
    }
    pub fn is_empty(&self) -> bool {
        self.num.is_empty()
    }
    /// Residue letters: the original letters when retained, otherwise decoded
    /// from the numeric codes (ambiguity codes then read back as `X`).
    pub fn seq(&self) -> Cow<'a, [u8]> {
        match self.letters {
            Some(l) => Cow::Borrowed(l),
            None => Cow::Owned(self.num.iter().map(|&r| num_to_aa(r)).collect()),
        }
    }
}

/// A collection of sequences in a compact arena plus aggregate statistics.
#[derive(Clone, Debug)]
pub struct SeqDb {
    residues: Vec<u8>,
    offsets: Vec<u64>, // len = n + 1
    ids: StrArena,
    headers: Option<StrArena>,
    letters: Option<Vec<u8>>,
    /// Total number of residues across all records (search-space size).
    pub total_residues: usize,
    max_len: usize,
}

impl Default for SeqDb {
    fn default() -> Self {
        SeqDb::new(LoadOptions::FULL)
    }
}

impl SeqDb {
    /// An empty database that retains what `opts` asks for.
    pub fn new(opts: LoadOptions) -> SeqDb {
        SeqDb {
            residues: Vec::new(),
            offsets: vec![0],
            ids: StrArena::new(),
            headers: if opts.keep_headers {
                Some(StrArena::new())
            } else {
                None
            },
            letters: if opts.keep_letters { Some(Vec::new()) } else { None },
            total_residues: 0,
            max_len: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.offsets.len() - 1
    }
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Append a record from its raw header (without `>`) and residue letters.
    /// Letters are encoded with `aa2num`.
    pub fn push_record(&mut self, header: &[u8], letters: &[u8], aa2num: &[u8; 256]) {
        let id = header_id(header);
        self.ids.push(id);
        if let Some(h) = self.headers.as_mut() {
            h.push(header);
        }
        self.residues
            .extend(letters.iter().map(|&b| aa_to_num(b, aa2num)));
        if let Some(l) = self.letters.as_mut() {
            l.extend_from_slice(letters);
        }
        self.finish_push(letters.len());
    }

    /// Append an already-encoded record (used by format readers and shards).
    pub fn push_encoded(&mut self, id: &str, header: Option<&str>, num: &[u8], letters: Option<&[u8]>) {
        self.ids.push(id.as_bytes());
        if let Some(h) = self.headers.as_mut() {
            h.push(header.unwrap_or(id).as_bytes());
        }
        self.residues.extend_from_slice(num);
        if let Some(l) = self.letters.as_mut() {
            match letters {
                Some(x) => l.extend_from_slice(x),
                None => l.extend(num.iter().map(|&r| num_to_aa(r))),
            }
        }
        self.finish_push(num.len());
    }

    fn finish_push(&mut self, n: usize) {
        self.offsets.push(self.residues.len() as u64);
        self.total_residues += n;
        if n > self.max_len {
            self.max_len = n;
        }
    }

    /// Numeric residues of record `i`.
    #[inline]
    pub fn num(&self, i: usize) -> &[u8] {
        let a = self.offsets[i] as usize;
        let b = self.offsets[i + 1] as usize;
        &self.residues[a..b]
    }

    /// Length of record `i`.
    #[inline]
    pub fn seq_len(&self, i: usize) -> usize {
        (self.offsets[i + 1] - self.offsets[i]) as usize
    }

    /// Identifier of record `i`.
    #[inline]
    pub fn id(&self, i: usize) -> &str {
        self.ids.get(i)
    }

    /// Full header of record `i` (falls back to the identifier when headers were
    /// not retained).
    pub fn header(&self, i: usize) -> &str {
        match &self.headers {
            Some(h) => h.get(i),
            None => self.ids.get(i),
        }
    }

    /// Whether full headers are retained.
    pub fn has_headers(&self) -> bool {
        self.headers.is_some()
    }

    /// Whether original letters are retained.
    pub fn has_letters(&self) -> bool {
        self.letters.is_some()
    }

    /// Residue letters of record `i` (original if retained, else decoded).
    pub fn seq(&self, i: usize) -> Cow<'_, [u8]> {
        self.record(i).seq()
    }

    /// Borrowed view of record `i`.
    pub fn record(&self, i: usize) -> RecordRef<'_> {
        let a = self.offsets[i] as usize;
        let b = self.offsets[i + 1] as usize;
        RecordRef {
            index: i,
            id: self.ids.get(i),
            header: self.header(i),
            num: &self.residues[a..b],
            letters: self.letters.as_ref().map(|l| &l[a..b]),
        }
    }

    /// Iterate records in input order.
    pub fn iter(&self) -> impl Iterator<Item = RecordRef<'_>> + '_ {
        (0..self.len()).map(move |i| self.record(i))
    }

    /// Longest record length.
    pub fn max_len(&self) -> usize {
        self.max_len
    }

    /// The whole residue arena (records are `offsets()[i]..offsets()[i+1]`).
    pub fn residues(&self) -> &[u8] {
        &self.residues
    }

    /// Record offsets into [`residues`](Self::residues) (length `len() + 1`).
    pub fn offsets(&self) -> &[u64] {
        &self.offsets
    }

    /// Heap bytes held by this database (for memory accounting).
    pub fn heap_bytes(&self) -> usize {
        self.residues.capacity()
            + self.offsets.capacity() * 8
            + self.ids.heap_bytes()
            + self.headers.as_ref().map(|h| h.heap_bytes()).unwrap_or(0)
            + self.letters.as_ref().map(|l| l.capacity()).unwrap_or(0)
    }

    /// Release spare capacity after loading.
    pub fn shrink_to_fit(&mut self) {
        self.residues.shrink_to_fit();
        self.offsets.shrink_to_fit();
        self.ids.shrink_to_fit();
        if let Some(h) = self.headers.as_mut() {
            h.shrink_to_fit();
        }
        if let Some(l) = self.letters.as_mut() {
            l.shrink_to_fit();
        }
    }

    /// Copy records `range` into a new database (distributed shards).
    pub fn subset(&self, range: std::ops::Range<usize>) -> SeqDb {
        let mut out = SeqDb::new(LoadOptions {
            keep_headers: self.headers.is_some(),
            keep_letters: self.letters.is_some(),
        });
        for i in range {
            let r = self.record(i);
            out.push_encoded(r.id, Some(r.header), r.num, r.letters);
        }
        out
    }

    /// Parse FASTA from any buffered reader, retaining everything
    /// ([`LoadOptions::FULL`]). Malformed input is an error.
    pub fn from_reader<R: BufRead>(reader: R, aa2num: &[u8; 256]) -> io::Result<SeqDb> {
        let mut fr = FastaReader::new(reader, "<input>");
        let db = fr.read_all(LoadOptions::FULL, aa2num)?;
        Ok(db)
    }

    /// Read a FASTA (plain or gzip) file. Kept for API compatibility; new code
    /// should use [`crate::seqio::load_sequences`], which also accepts databases.
    pub fn read_fasta<P: AsRef<Path>>(path: P, aa2num: &[u8; 256]) -> io::Result<SeqDb> {
        let p = path.as_ref().to_string_lossy().to_string();
        let mut fr = crate::seqio::open_fasta(&p)?;
        fr.read_all(LoadOptions::FULL, aa2num)
    }

    /// Write the database back out as FASTA (60 columns per line).
    pub fn write_fasta<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let mut out = io::BufWriter::new(File::create(path)?);
        for i in 0..self.len() {
            let r = self.record(i);
            writeln!(out, ">{}", r.header)?;
            let s = r.seq();
            for chunk in s.chunks(60) {
                out.write_all(chunk)?;
                out.write_all(b"\n")?;
            }
        }
        out.flush()
    }
}

/// First whitespace-delimited token of a header.
pub fn header_id(header: &[u8]) -> &[u8] {
    let start = header
        .iter()
        .position(|b| !b.is_ascii_whitespace())
        .unwrap_or(header.len());
    let rest = &header[start..];
    let end = rest
        .iter()
        .position(|b| b.is_ascii_whitespace())
        .unwrap_or(rest.len());
    &rest[..end]
}

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

/// Streaming, validating FASTA reader.
pub struct FastaReader<R: BufRead> {
    r: R,
    source: String,
    line: Vec<u8>,
    lineno: u64,
    /// Header of the next record (already read) and its line number.
    pending: Option<(Vec<u8>, u64)>,
    started: bool,
    eof: bool,
    records: u64,
    /// Count of alignment-gap characters (`-`, `.`) removed from sequences.
    pub gap_chars_removed: u64,
    /// Count of non-standard letters (B, Z, J, O, U, `*`, ...) encoded as `X`.
    pub ambiguous_letters: u64,
}

impl<R: BufRead> FastaReader<R> {
    pub fn new(r: R, source: &str) -> Self {
        FastaReader {
            r,
            source: source.to_string(),
            line: Vec::with_capacity(256),
            lineno: 0,
            pending: None,
            started: false,
            eof: false,
            records: 0,
            gap_chars_removed: 0,
            ambiguous_letters: 0,
        }
    }

    /// Records returned so far.
    pub fn records_read(&self) -> u64 {
        self.records
    }

    fn read_line(&mut self) -> io::Result<bool> {
        self.line.clear();
        let n = self.r.read_until(b'\n', &mut self.line).map_err(|e| {
            io::Error::new(e.kind(), format!("{}: read error: {}", self.source, e))
        })?;
        if n == 0 {
            return Ok(false);
        }
        self.lineno += 1;
        while matches!(self.line.last(), Some(b'\n') | Some(b'\r')) {
            self.line.pop();
        }
        Ok(true)
    }

    /// Read the next record. `header` receives the header text without `>`,
    /// `seq` the uppercased residue letters. Returns `Ok(false)` at end of input.
    pub fn next_record(&mut self, header: &mut Vec<u8>, seq: &mut Vec<u8>) -> io::Result<bool> {
        header.clear();
        seq.clear();
        // Locate the header of this record.
        let header_line = match self.pending.take() {
            Some((h, ln)) => {
                header.extend_from_slice(&h);
                ln
            }
            None => {
                if self.eof {
                    return Ok(false);
                }
                loop {
                    if !self.read_line()? {
                        self.eof = true;
                        return Ok(false);
                    }
                    let t = trim_ascii(&self.line);
                    if t.is_empty() || t[0] == b';' {
                        continue;
                    }
                    if t[0] == b'>' {
                        header.extend_from_slice(&t[1..]);
                        break self.lineno;
                    }
                    if !self.started {
                        return Err(invalid(format!(
                            "{}:{}: malformed FASTA: sequence data before the first '>' header \
                             (is this a FASTA file? MMseqs2 database prefixes must be given \
                             without an extension and need a matching .dbtype file)",
                            self.source, self.lineno
                        )));
                    }
                    return Err(invalid(format!(
                        "{}:{}: internal FASTA reader state error (sequence line outside a record)",
                        self.source, self.lineno
                    )));
                }
            }
        };
        self.started = true;
        let trimmed = trim_ascii(header).to_vec();
        header.clear();
        header.extend_from_slice(&trimmed);
        if header_id(header).is_empty() {
            return Err(invalid(format!(
                "{}:{}: empty FASTA header: every record needs an identifier after '>'",
                self.source, header_line
            )));
        }

        // Sequence lines until the next header or EOF.
        loop {
            if !self.read_line()? {
                self.eof = true;
                break;
            }
            let lineno = self.lineno;
            let line = std::mem::take(&mut self.line);
            let t = trim_ascii(&line);
            if t.first() == Some(&b'>') {
                self.pending = Some((t[1..].to_vec(), lineno));
                self.line = line;
                break;
            }
            if t.first() == Some(&b';') {
                self.line = line;
                continue;
            }
            for (col, &b) in t.iter().enumerate() {
                match b {
                    b'A'..=b'Z' => {
                        if matches!(b, b'B' | b'J' | b'O' | b'U' | b'X' | b'Z') {
                            self.ambiguous_letters += 1;
                        }
                        seq.push(b)
                    }
                    b'a'..=b'z' => {
                        let u = b.to_ascii_uppercase();
                        if matches!(u, b'B' | b'J' | b'O' | b'U' | b'X' | b'Z') {
                            self.ambiguous_letters += 1;
                        }
                        seq.push(u)
                    }
                    b'*' => {
                        self.ambiguous_letters += 1;
                        seq.push(b'*')
                    }
                    b'-' | b'.' => self.gap_chars_removed += 1,
                    b' ' | b'\t' | b'\r' => {}
                    _ => {
                        let shown = if b.is_ascii_graphic() {
                            format!("'{}'", b as char)
                        } else {
                            format!("byte 0x{:02x}", b)
                        };
                        return Err(invalid(format!(
                            "{}:{}:{}: unexpected {} in the sequence of record '{}' \
                             (expected residue letters; is the file binary, compressed with an \
                             unsupported codec, or not FASTA?)",
                            self.source,
                            lineno,
                            col + 1,
                            shown,
                            String::from_utf8_lossy(header_id(header))
                        )));
                    }
                }
            }
            self.line = line;
        }
        // A terminal stop symbol is not a residue.
        if seq.last() == Some(&b'*') {
            seq.pop();
            self.ambiguous_letters = self.ambiguous_letters.saturating_sub(1);
        }
        if seq.is_empty() {
            return Err(invalid(format!(
                "{}:{}: record '{}' has no residues (empty sequences cannot be searched; \
                 remove the record or fix the file)",
                self.source,
                header_line,
                String::from_utf8_lossy(header_id(header))
            )));
        }
        if seq.len() as u64 > u32::MAX as u64 {
            return Err(invalid(format!(
                "{}:{}: record '{}' has {} residues; the maximum supported length is {}",
                self.source,
                header_line,
                String::from_utf8_lossy(header_id(header)),
                seq.len(),
                u32::MAX
            )));
        }
        self.records += 1;
        Ok(true)
    }

    /// Read every remaining record into a [`SeqDb`]. An input with no records
    /// is an error.
    pub fn read_all(&mut self, opts: LoadOptions, aa2num: &[u8; 256]) -> io::Result<SeqDb> {
        let mut db = SeqDb::new(opts);
        let mut h = Vec::new();
        let mut s = Vec::new();
        while self.next_record(&mut h, &mut s)? {
            if db.len() as u64 >= u32::MAX as u64 {
                return Err(invalid(format!(
                    "{}: more than {} records are not supported",
                    self.source,
                    u32::MAX
                )));
            }
            db.push_record(&h, &s, aa2num);
        }
        if db.is_empty() {
            return Err(invalid(format!(
                "{}: no FASTA records found (the input is empty)",
                self.source
            )));
        }
        db.shrink_to_fit();
        Ok(db)
    }
}

fn trim_ascii(s: &[u8]) -> &[u8] {
    let start = s
        .iter()
        .position(|b| !b.is_ascii_whitespace())
        .unwrap_or(s.len());
    let end = s
        .iter()
        .rposition(|b| !b.is_ascii_whitespace())
        .map(|p| p + 1)
        .unwrap_or(start);
    if start >= end {
        &s[0..0]
    } else {
        &s[start..end]
    }
}

/// Read an entire file to a `String` (small helper for MSA / matrix inputs).
pub fn read_to_string<P: AsRef<Path>>(path: P) -> io::Result<String> {
    let p = path.as_ref();
    let mut s = String::new();
    File::open(p)
        .map_err(|e| io::Error::new(e.kind(), format!("{}: {}", p.display(), e)))?
        .read_to_string(&mut s)
        .map_err(|e| io::Error::new(e.kind(), format!("{}: {}", p.display(), e)))?;
    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::build_aa2num;
    use std::io::Cursor;

    fn parse(data: &str) -> io::Result<SeqDb> {
        SeqDb::from_reader(Cursor::new(data), &build_aa2num())
    }

    #[test]
    fn parses_multi_record_fasta() {
        let db = parse(">sp|P1 test one\nMKV\nLLA\n>sp|P2 test two\nGGGG\n").unwrap();
        assert_eq!(db.len(), 2);
        assert_eq!(db.id(0), "sp|P1");
        assert_eq!(db.header(0), "sp|P1 test one");
        assert_eq!(&*db.seq(0), b"MKVLLA");
        assert_eq!(db.num(0).len(), 6);
        assert_eq!(db.seq_len(1), 4);
        assert_eq!(db.total_residues, 10);
        assert_eq!(db.max_len(), 6);
    }

    #[test]
    fn rejects_empty_input() {
        let e = parse("").unwrap_err();
        assert!(e.to_string().contains("no FASTA records"), "{}", e);
        let e = parse("\n\n  \n").unwrap_err();
        assert!(e.to_string().contains("no FASTA records"), "{}", e);
    }

    #[test]
    fn rejects_sequence_before_header() {
        let e = parse("MKVLLA\n>a\nMKV\n").unwrap_err();
        assert!(e.to_string().contains("before the first '>'"), "{}", e);
        assert!(e.to_string().contains(":1:"), "line number missing: {}", e);
    }

    #[test]
    fn rejects_empty_record_and_header() {
        let e = parse(">a\n>b\nMKV\n").unwrap_err();
        assert!(e.to_string().contains("'a' has no residues"), "{}", e);
        let e = parse(">\nMKV\n").unwrap_err();
        assert!(e.to_string().contains("empty FASTA header"), "{}", e);
        let e = parse(">a\nMKV\n>c\n").unwrap_err();
        assert!(e.to_string().contains("'c' has no residues"), "{}", e);
    }

    #[test]
    fn rejects_invalid_characters_with_position() {
        let e = parse(">a\nMKV1LA\n").unwrap_err();
        let m = e.to_string();
        assert!(m.contains(":2:4:") && m.contains("'1'"), "{}", m);
        let e = parse(">a\nMK\x01V\n").unwrap_err();
        assert!(e.to_string().contains("0x01"), "{}", e);
    }

    #[test]
    fn tolerates_crlf_gaps_case_stops_and_comments() {
        let db = parse(";comment\r\n>a desc\r\nmk-v.l\r\nLA*\r\n>b\nM*K\n").unwrap();
        assert_eq!(&*db.seq(0), b"MKVLLA");
        assert_eq!(db.header(0), "a desc");
        // internal stop is kept as a letter and encodes to X
        assert_eq!(&*db.seq(1), b"M*K");
        assert_eq!(db.num(1)[1], crate::alphabet::ANY);
    }

    #[test]
    fn search_options_drop_letters_and_headers() {
        let t = build_aa2num();
        let mut fr = FastaReader::new(Cursor::new(">x long description\nMKVB\n"), "t");
        let db = fr.read_all(LoadOptions::SEARCH, &t).unwrap();
        assert!(!db.has_headers() && !db.has_letters());
        assert_eq!(db.header(0), "x");
        // decoded letters: ambiguous B reads back as X
        assert_eq!(&*db.seq(0), b"MKVX");
        assert_eq!(fr.ambiguous_letters, 1);
    }

    #[test]
    fn subset_copies_records() {
        let db = parse(">a\nMK\n>b\nVLL\n>c\nA\n").unwrap();
        let s = db.subset(1..3);
        assert_eq!(s.len(), 2);
        assert_eq!(s.id(0), "b");
        assert_eq!(&*s.seq(1), b"A");
        assert_eq!(s.total_residues, 4);
    }
}
