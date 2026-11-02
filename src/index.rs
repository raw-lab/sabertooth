//! Persistent, versioned, memory-mapped seed index with a compact sequence
//! arena.
//!
//! v0.2 rebuilt a `HashMap<u64, Vec<(u32, u32)>>` — one heap vector per distinct
//! k-mer — on every search. This module replaces it with a single file
//! (`<db>.stidx`) holding everything search needs, laid out as contiguous
//! little-endian arrays that are memory-mapped (file-backed pages: the kernel
//! can evict them under memory pressure and re-read on demand):
//!
//! | section            | type            | content                                   |
//! |--------------------|-----------------|-------------------------------------------|
//! | residues           | `u8 × N`        | numeric residues of all targets           |
//! | seq_offsets        | `u64 × (n+1)`   | record `i` = `residues[o[i]..o[i+1]]`     |
//! | id_bytes/offsets   | `u8`, `u64`     | identifiers (read only for output)        |
//! | postings           | `u64 × P`       | `(target_id << 32) | position`, grouped by k-mer, sorted by (target, position) |
//! | kmer_table         | `[u32;2] × K`   | per distinct k-mer: suffix code, offset within its bucket |
//! | bucket_kmer_start  | `u64 × (B+1)`   | bucket = first `p` residues of the k-mer  |
//! | bucket_post_start  | `u64 × (B+1)`   |                                           |
//!
//! A lookup is `bucket = code / 20^(w−p)`, then a binary search of the bucket's
//! sorted suffixes (≤ 400 entries for k = 6, p = 4), giving a contiguous
//! posting slice — no per-k-mer allocation, 8 bytes of lookup structure per
//! distinct k-mer. Databases ≥ 16·20^(w−1) residues use `p = w − 1` (one
//! bucket per `(w−1)`-prefix, no binary search; 51 MB of bucket tables for
//! k = 6). A presence bitmap over the first `min(w, 6)` residues (≤ 8 MiB) is
//! built lazily on first use so most absent k-mers cost one bit test.
//!
//! ## Integrity and compatibility
//!
//! The 4 KiB header records magic, format version, alphabet, seed pattern,
//! masking parameters, sequence/residue/posting counts, a source fingerprint
//! (size + mtime of the FASTA or of the MMseqs2 data/index/dbtype files) and a
//! section table; a header checksum and a footer marker (`SBTIDXOK` + file
//! length) are written last. Builds go to `<out>.tmp-<pid>` and are renamed into
//! place only when complete, so an interrupted build never leaves a file that
//! opens. Opening checks header checksum, footer, section bounds and table
//! end-values — O(1) work, no rescan of the database; `verify` additionally
//! checks every section checksum. Staleness is detected from the fingerprint
//! (`stat` only).
//!
//! ## Bounded-memory construction
//!
//! Sequences are streamed into the output file (never held twice), then the
//! postings are produced by partitions of buckets sized to the memory budget:
//! one counting pass, then per partition a parallel counting-sort scatter of
//! that partition's occurrences and a per-bucket sort. Peak anonymous memory is
//! the bucket tables plus one partition buffer (16 bytes per occurrence).

use crate::fasta::SeqDb;
use crate::mask::{low_complexity_mask_into, MaskParams};
use crate::prefilter::SeedPattern;
use memmap2::{Mmap, MmapOptions};
use rayon::prelude::*;
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

pub const MAGIC: [u8; 8] = *b"SBTIDX\0\x01";
pub const FOOTER_MAGIC: [u8; 8] = *b"SBTIDXOK";
pub const FORMAT_VERSION: u32 = 2;
/// Oldest format this build still reads (v1: 8-byte postings, sparse layout).
pub const MIN_FORMAT_VERSION: u32 = 1;
/// Alphabet id 1 = MMseqs2 amino-acid order `ACDEFGHIKLMNPQRSTVWY` + `X`.
pub const ALPHABET_AA21: u32 = 1;
const HEADER_SIZE: u64 = 4096;
const FOOTER_SIZE: u64 = 16;
const ALIGN: u64 = 64;
const A: u64 = crate::alphabet::PROFILE_AA_SIZE as u64;

/// Division by a power of 20 known only at run time, compiled to a
/// multiply-shift per exponent instead of a hardware `div` (64-bit division
/// was a measurable share of seed lookups and index builds).
#[derive(Clone, Copy, Debug)]
pub struct Pow20 {
    e: u32,
    v: u64,
}

impl Pow20 {
    pub fn new(e: u32) -> Pow20 {
        Pow20 { e, v: A.pow(e) }
    }
    /// The divisor, `20^e`.
    #[inline(always)]
    pub fn value(self) -> u64 {
        self.v
    }
    #[inline(always)]
    pub fn quot(self, x: u64) -> u64 {
        match self.e {
            0 => x,
            1 => x / 20,
            2 => x / 400,
            3 => x / 8_000,
            4 => x / 160_000,
            5 => x / 3_200_000,
            6 => x / 64_000_000,
            7 => x / 1_280_000_000,
            _ => x / self.v,
        }
    }
    #[inline(always)]
    pub fn modulo(self, x: u64) -> u64 {
        x - self.quot(x) * self.v
    }
}

/// Index file name for a database prefix or FASTA path.
pub fn default_index_path(source: &str) -> PathBuf {
    PathBuf::from(format!("{}.stidx", source))
}

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

// ---------------------------------------------------------------------------
// hashing (integrity checks, not cryptographic)
// ---------------------------------------------------------------------------

/// Streaming 64-bit hash over arbitrary byte chunks (word-wise multiply-mix).
#[derive(Clone)]
pub struct Hasher64 {
    state: u64,
    buf: [u8; 8],
    nbuf: usize,
    len: u64,
}

impl Default for Hasher64 {
    fn default() -> Self {
        Hasher64 {
            state: 0x9E37_79B9_7F4A_7C15,
            buf: [0; 8],
            nbuf: 0,
            len: 0,
        }
    }
}

impl Hasher64 {
    #[inline]
    fn mix(&mut self, w: u64) {
        let x = (self.state ^ w).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        self.state = x ^ (x >> 31);
    }
    pub fn update(&mut self, mut data: &[u8]) {
        self.len += data.len() as u64;
        if self.nbuf > 0 {
            while self.nbuf < 8 && !data.is_empty() {
                self.buf[self.nbuf] = data[0];
                self.nbuf += 1;
                data = &data[1..];
            }
            if self.nbuf == 8 {
                let w = u64::from_le_bytes(self.buf);
                self.mix(w);
                self.nbuf = 0;
            }
        }
        let mut chunks = data.chunks_exact(8);
        for c in &mut chunks {
            let w = u64::from_le_bytes([c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]]);
            self.mix(w);
        }
        for &b in chunks.remainder() {
            self.buf[self.nbuf] = b;
            self.nbuf += 1;
        }
    }
    pub fn finish(&self) -> u64 {
        let mut h = self.clone();
        if h.nbuf > 0 {
            let mut b = [0u8; 8];
            b[..h.nbuf].copy_from_slice(&h.buf[..h.nbuf]);
            h.mix(u64::from_le_bytes(b));
        }
        let l = h.len;
        h.mix(l);
        let mut z = h.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

/// A writer that counts bytes and hashes them.
struct HashingWriter<W: Write> {
    inner: W,
    hasher: Hasher64,
    written: u64,
}

impl<W: Write> HashingWriter<W> {
    fn new(inner: W) -> Self {
        HashingWriter {
            inner,
            hasher: Hasher64::default(),
            written: 0,
        }
    }
}

impl<W: Write> Write for HashingWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let n = self.inner.write(buf)?;
        self.hasher.update(&buf[..n]);
        self.written += n as u64;
        Ok(n)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

// ---------------------------------------------------------------------------
// metadata
// ---------------------------------------------------------------------------

/// Section identifiers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u32)]
pub enum SectionKind {
    Residues = 1,
    SeqOffsets = 2,
    IdBytes = 3,
    IdOffsets = 4,
    Postings = 5,
    KmerTable = 6,
    BucketKmerStart = 7,
    BucketPostStart = 8,
    /// direct layout: posting offset of every possible k-mer code (4 or 8
    /// bytes), or — grouped offsets — its 16-bit offset relative to the start
    /// of its group of [`OFFSET_GROUP`] codes
    KmerOffsets = 9,
    /// grouped direct offsets: posting offset of the first code of every group
    KmerGroupBase = 10,
}

/// Codes per group of the grouped direct offset table.
pub const OFFSET_GROUP: usize = 64;

/// The automatic layout choice takes the direct table whenever the sparse
/// tables would exceed this size: they no longer fit any CPU cache, so every
/// sparse lookup costs several dependent memory accesses where the direct
/// table needs one (1.7× faster searches at 17 M residues, `bench/README.md`).
pub const AUTO_DIRECT_MIN_SPARSE_BYTES: u64 = 64 << 20;

const SECTION_KINDS: [SectionKind; 10] = [
    SectionKind::Residues,
    SectionKind::SeqOffsets,
    SectionKind::IdBytes,
    SectionKind::IdOffsets,
    SectionKind::Postings,
    SectionKind::KmerTable,
    SectionKind::BucketKmerStart,
    SectionKind::BucketPostStart,
    SectionKind::KmerOffsets,
    SectionKind::KmerGroupBase,
];

/// One section table entry.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Section {
    pub offset: u64,
    pub len: u64,
    pub hash: u64,
}

/// Fingerprint of the files an index was built from (`stat` data only).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SourceFingerprint {
    /// 0 = none (transient), 1 = FASTA file, 2 = MMseqs2 database prefix.
    pub kind: u32,
    /// (size, mtime in ns) of up to three files (FASTA: the file; DB: data,
    /// `.index`, `.dbtype`).
    pub files: [(u64, i64); 3],
    pub path: String,
}

fn stat_of(p: &str) -> io::Result<(u64, i64)> {
    let md = fs::metadata(p).map_err(|e| io::Error::new(e.kind(), format!("{}: {}", p, e)))?;
    let mtime = md
        .modified()
        .ok()
        .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
        .map(|d| d.as_nanos() as i64)
        .unwrap_or(0);
    Ok((md.len(), mtime))
}

impl SourceFingerprint {
    /// Fingerprint a FASTA path or MMseqs2 database prefix.
    pub fn of(source: &str) -> io::Result<SourceFingerprint> {
        if crate::seqio::is_db_prefix(source) {
            Ok(SourceFingerprint {
                kind: 2,
                files: [
                    stat_of(source)?,
                    stat_of(&format!("{}.index", source))?,
                    stat_of(&format!("{}.dbtype", source))?,
                ],
                path: source.to_string(),
            })
        } else {
            Ok(SourceFingerprint {
                kind: 1,
                files: [stat_of(source)?, (0, 0), (0, 0)],
                path: source.to_string(),
            })
        }
    }

    /// No source (transient indexes, distributed shards: their staleness is
    /// tracked by their owner).
    pub fn transient() -> SourceFingerprint {
        SourceFingerprint::default()
    }
}

/// Everything recorded in an index header.
#[derive(Clone, Debug)]
pub struct IndexMeta {
    pub version: u32,
    pub alphabet: u32,
    pub pattern: SeedPattern,
    /// Residues per bucket prefix (`p`).
    pub prefix_len: u32,
    pub masked: bool,
    pub mask: MaskParams,
    pub num_seqs: u64,
    pub total_residues: u64,
    pub max_seq_len: u64,
    pub num_postings: u64,
    pub num_kmers: u64,
    pub num_buckets: u64,
    pub masked_residues: u64,
    pub source: SourceFingerprint,
    pub created_unix: u64,
    pub tool_version: String,
    pub sections: Vec<(SectionKind, Section)>,
    /// k-mer lookup layout (format ≥ 2; v1 files are sparse)
    pub layout: IndexLayout,
    /// bytes per posting (1..=8) and bits of the position field: a posting is
    /// the little-endian integer `(target << pos_bits) | position`
    /// (v1 files: 8 bytes, 32 position bits)
    pub posting_width: u32,
    pub posting_pos_bits: u32,
    /// bytes per entry of the direct offset table: 4 or 8, or 2 for grouped
    /// offsets (16-bit offsets relative to group bases of 4 bytes, or 8 when
    /// there are 2^32 postings or more); 0 = sparse
    pub offset_width: u32,
}

/// How k-mer codes are resolved to posting ranges.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IndexLayout {
    /// Bucket tables + a sorted (suffix, offset) entry per distinct k-mer:
    /// 8 bytes per distinct k-mer (small and medium databases).
    Sparse,
    /// One offset per possible k-mer code (`20^w + 1` entries): a lookup is
    /// one contiguous read. Used when it is no larger than the sparse tables
    /// (dense k-mer spaces, i.e. large databases).
    Direct,
}

/// Layout requested for a build.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LayoutChoice {
    /// The smaller of the two (ties: direct).
    Auto,
    Sparse,
    Direct,
}

impl IndexLayout {
    pub fn name(self) -> &'static str {
        match self {
            IndexLayout::Sparse => "sparse",
            IndexLayout::Direct => "direct",
        }
    }
}

impl IndexMeta {
    fn section(&self, k: SectionKind) -> Section {
        self.sections
            .iter()
            .find(|(kk, _)| *kk == k)
            .map(|(_, s)| *s)
            .unwrap_or_default()
    }
}

/// Bucket prefix length for a seed weight and database size. The default is
/// `min(w, 4)`, raised so the suffix (`w − p` residues) fits `u32`
/// (`20^7 < 2^32`). For large databases the prefix grows to `w − 1` (at most
/// 5): every bucket is then exactly one `(w−1)`-residue run, so
/// [`SeedIndex::prefix_ranges`] needs no binary search (seed lookups ≈1.4×
/// faster at 68 M residues). The bucket tables then cost `16 · 20^(w−1)`
/// bytes (51 MB for k = 6), so the layout is used only from
/// `threshold · 20^(w−1)` residues on (≈51 M residues for k = 6 at the
/// default threshold, where the tables are ~5% of the index).
fn prefix_len_for_size(weight: usize, total_residues: u64, threshold: u64) -> u32 {
    let mut p = weight.min(4);
    let big = weight.saturating_sub(1).min(5);
    if big > p && total_residues >= threshold.saturating_mul(A.pow(big as u32)) {
        p = big;
    }
    p.max(weight.saturating_sub(7)) as u32
}

/// Default residues per `(w−1)`-prefix bucket for the large layout.
pub const LARGE_PREFIX_MIN_RESIDUES_PER_BUCKET: u64 = 16;

struct ByteWriter {
    b: Vec<u8>,
}

impl ByteWriter {
    fn u32(&mut self, v: u32) {
        self.b.extend_from_slice(&v.to_le_bytes());
    }
    fn u64(&mut self, v: u64) {
        self.b.extend_from_slice(&v.to_le_bytes());
    }
    fn i64(&mut self, v: i64) {
        self.b.extend_from_slice(&v.to_le_bytes());
    }
    fn f64(&mut self, v: f64) {
        self.b.extend_from_slice(&v.to_le_bytes());
    }
    fn fixed(&mut self, s: &[u8], n: usize) {
        let k = s.len().min(n);
        self.b.extend_from_slice(&s[..k]);
        self.b.extend(std::iter::repeat_n(0u8, n - k));
    }
}

struct ByteReader<'a> {
    b: &'a [u8],
    p: usize,
}

impl<'a> ByteReader<'a> {
    fn take(&mut self, n: usize) -> &'a [u8] {
        let s = &self.b[self.p..self.p + n];
        self.p += n;
        s
    }
    fn u32(&mut self) -> u32 {
        let s = self.take(4);
        u32::from_le_bytes([s[0], s[1], s[2], s[3]])
    }
    fn u64(&mut self) -> u64 {
        let s = self.take(8);
        u64::from_le_bytes([s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]])
    }
    fn i64(&mut self) -> i64 {
        self.u64() as i64
    }
    fn f64(&mut self) -> f64 {
        f64::from_bits(self.u64())
    }
    fn fixed_str(&mut self, n: usize) -> String {
        let s = self.take(n);
        let end = s.iter().position(|&b| b == 0).unwrap_or(n);
        String::from_utf8_lossy(&s[..end]).to_string()
    }
}

const MAX_SECTIONS: usize = 16;

fn encode_header(m: &IndexMeta) -> Vec<u8> {
    let mut w = ByteWriter { b: Vec::with_capacity(HEADER_SIZE as usize) };
    w.b.extend_from_slice(&MAGIC);
    w.u32(m.version);
    w.u32(HEADER_SIZE as u32);
    w.u32(m.alphabet);
    w.u32(m.prefix_len);
    w.fixed(m.pattern.mask_string().as_bytes(), 64);
    w.u32(m.masked as u32);
    w.u32(m.mask.window as u32);
    w.f64(m.mask.locut);
    w.f64(m.mask.hicut);
    w.u64(m.num_seqs);
    w.u64(m.total_residues);
    w.u64(m.max_seq_len);
    w.u64(m.num_postings);
    w.u64(m.num_kmers);
    w.u64(m.num_buckets);
    w.u64(m.masked_residues);
    w.u32(m.source.kind);
    w.u32(0);
    for (len, mt) in m.source.files {
        w.u64(len);
        w.i64(mt);
    }
    w.fixed(m.source.path.as_bytes(), 1024);
    w.u64(m.created_unix);
    w.fixed(m.tool_version.as_bytes(), 32);
    w.u32(m.sections.len() as u32);
    w.u32(0);
    for i in 0..MAX_SECTIONS {
        match m.sections.get(i) {
            Some((k, s)) => {
                w.u32(*k as u32);
                w.u32(0);
                w.u64(s.offset);
                w.u64(s.len);
                w.u64(s.hash);
            }
            None => w.fixed(&[], 32),
        }
    }
    // format 2
    w.u32(match m.layout {
        IndexLayout::Sparse => 0,
        IndexLayout::Direct => 1,
    });
    w.u32(m.posting_width);
    w.u32(m.posting_pos_bits);
    w.u32(m.offset_width);
    let used = w.b.len();
    assert!(used + 8 <= HEADER_SIZE as usize, "header overflow");
    w.b.resize(HEADER_SIZE as usize - 8, 0);
    let mut h = Hasher64::default();
    h.update(&w.b);
    let sum = h.finish();
    w.u64(sum);
    w.b
}

fn decode_header(b: &[u8], path: &str) -> io::Result<IndexMeta> {
    if b.len() < HEADER_SIZE as usize {
        return Err(invalid(format!("{}: truncated index header", path)));
    }
    if b[..8] != MAGIC {
        return Err(invalid(format!("{}: not a Sabertooth index (bad magic)", path)));
    }
    let mut h = Hasher64::default();
    h.update(&b[..HEADER_SIZE as usize - 8]);
    let stored = u64::from_le_bytes(b[HEADER_SIZE as usize - 8..HEADER_SIZE as usize].try_into().unwrap());
    if h.finish() != stored {
        return Err(invalid(format!("{}: index header checksum mismatch (corrupt file)", path)));
    }
    let mut r = ByteReader { b, p: 8 };
    let version = r.u32();
    if !(MIN_FORMAT_VERSION..=FORMAT_VERSION).contains(&version) {
        return Err(invalid(format!(
            "{}: index format version {} is not supported by this build (reads {}..={}); rebuild with `sabertooth createindex`",
            path, version, MIN_FORMAT_VERSION, FORMAT_VERSION
        )));
    }
    let _hs = r.u32();
    let alphabet = r.u32();
    let prefix_len = r.u32();
    let mask_str = r.fixed_str(64);
    let pattern = SeedPattern::from_mask(&mask_str)
        .map_err(|e| invalid(format!("{}: bad seed pattern in header: {}", path, e)))?;
    let masked = r.u32() != 0;
    let window = r.u32() as usize;
    let locut = r.f64();
    let hicut = r.f64();
    let num_seqs = r.u64();
    let total_residues = r.u64();
    let max_seq_len = r.u64();
    let num_postings = r.u64();
    let num_kmers = r.u64();
    let num_buckets = r.u64();
    let masked_residues = r.u64();
    let kind = r.u32();
    let _ = r.u32();
    let mut files = [(0u64, 0i64); 3];
    for f in files.iter_mut() {
        let len = r.u64();
        let mt = r.i64();
        *f = (len, mt);
    }
    let spath = r.fixed_str(1024);
    let created_unix = r.u64();
    let tool_version = r.fixed_str(32);
    let ns = r.u32() as usize;
    let _ = r.u32();
    if ns > MAX_SECTIONS {
        return Err(invalid(format!("{}: bad section count", path)));
    }
    let mut sections = Vec::new();
    for i in 0..MAX_SECTIONS {
        let k = r.u32();
        let _ = r.u32();
        let offset = r.u64();
        let len = r.u64();
        let hash = r.u64();
        if i < ns {
            let kind = SECTION_KINDS
                .iter()
                .find(|x| **x as u32 == k)
                .copied()
                .ok_or_else(|| invalid(format!("{}: unknown section kind {}", path, k)))?;
            sections.push((kind, Section { offset, len, hash }));
        }
    }
    let (layout, posting_width, posting_pos_bits, offset_width) = if version >= 2 {
        let l = r.u32();
        let pw = r.u32();
        let pb = r.u32();
        let ow = r.u32();
        let layout = match l {
            0 => IndexLayout::Sparse,
            1 => IndexLayout::Direct,
            _ => return Err(invalid(format!("{}: unknown index layout {}", path, l))),
        };
        if !(1..=8).contains(&pw) || pb > 32 || pb > pw * 8 || (layout == IndexLayout::Direct && ow != 2 && ow != 4 && ow != 8) {
            return Err(invalid(format!("{}: bad posting/offset encoding in header (corrupt index)", path)));
        }
        (layout, pw, pb, ow)
    } else {
        (IndexLayout::Sparse, 8, 32, 0)
    };
    Ok(IndexMeta {
        version,
        alphabet,
        pattern,
        prefix_len,
        masked,
        mask: MaskParams { window, locut, hicut },
        num_seqs,
        total_residues,
        max_seq_len,
        num_postings,
        num_kmers,
        num_buckets,
        masked_residues,
        source: SourceFingerprint {
            kind,
            files,
            path: spath,
        },
        created_unix,
        tool_version,
        sections,
        layout,
        posting_width,
        posting_pos_bits,
        offset_width,
    })
}

// ---------------------------------------------------------------------------
// build
// ---------------------------------------------------------------------------

/// Index construction options.
#[derive(Clone, Debug)]
pub struct BuildOptions {
    pub pattern: SeedPattern,
    /// Low-complexity masking of target seeds (`None` = off).
    pub mask: Option<MaskParams>,
    /// Memory budget (bytes) for the posting partition buffer.
    pub memory_budget: usize,
    /// Directory for temporary files.
    pub tmp_dir: PathBuf,
    /// Residues per `(w−1)`-prefix bucket at which the large-database bucket
    /// layout is used (see `prefix_len_for_size`; 0 forces it).
    pub large_prefix_threshold: u64,
    /// k-mer lookup layout (default: see [`AUTO_DIRECT_MIN_SPARSE_BYTES`]).
    pub layout: LayoutChoice,
    /// Direct layout: store 16-bit offsets in groups of [`OFFSET_GROUP`] codes
    /// (half the table) when every group's postings fit; `false` forces full
    /// 4/8-byte offsets (tests).
    pub direct_grouped: bool,
    /// Force a posting width in bytes (tests; default: the smallest that
    /// holds every (target, position)).
    pub posting_width: Option<u32>,
}

impl BuildOptions {
    pub fn new(pattern: SeedPattern) -> Self {
        BuildOptions {
            pattern,
            mask: None,
            memory_budget: 1 << 30,
            tmp_dir: std::env::temp_dir(),
            large_prefix_threshold: LARGE_PREFIX_MIN_RESIDUES_PER_BUCKET,
            layout: LayoutChoice::Auto,
            direct_grouped: true,
            posting_width: None,
        }
    }
}

/// What a build did (for reports and benchmarks).
#[derive(Clone, Debug, Default)]
pub struct BuildStats {
    pub num_seqs: u64,
    pub total_residues: u64,
    pub num_postings: u64,
    pub num_kmers: u64,
    pub masked_residues: u64,
    pub partitions: usize,
    pub partition_buffer_bytes: usize,
    pub bucket_table_bytes: usize,
    /// "sparse" or "direct"
    pub layout: &'static str,
    pub posting_width: u32,
    /// bytes of the k-mer lookup structures (sparse tables or direct offsets)
    pub lookup_bytes: u64,
    pub file_bytes: u64,
    pub seconds_sequences: f64,
    pub seconds_count: f64,
    pub seconds_postings: f64,
    pub seconds_finish: f64,
}

/// A record source for the builder.
pub trait RecordSource {
    /// Call `f(id, numeric_residues)` for every record in order.
    #[allow(clippy::type_complexity)]
    fn for_each(&mut self, f: &mut dyn FnMut(&[u8], &[u8]) -> io::Result<()>) -> io::Result<()>;
}

/// Records from an in-memory database.
pub struct DbSource<'a>(pub &'a SeqDb);

impl<'a> RecordSource for DbSource<'a> {
    fn for_each(&mut self, f: &mut dyn FnMut(&[u8], &[u8]) -> io::Result<()>) -> io::Result<()> {
        for i in 0..self.0.len() {
            f(self.0.id(i).as_bytes(), self.0.num(i))?;
        }
        Ok(())
    }
}

/// Records streamed from an MMseqs2 database prefix — never fully loaded.
pub struct MmseqsDbSource {
    pub prefix: String,
    pub aa2num: [u8; 256],
}

impl RecordSource for MmseqsDbSource {
    fn for_each(&mut self, f: &mut dyn FnMut(&[u8], &[u8]) -> io::Result<()>) -> io::Result<()> {
        let mut num = Vec::new();
        let t = self.aa2num;
        crate::mmdb::stream_db(
            Path::new(&self.prefix),
            false,
            crate::seqio::SeqKind::Protein,
            |_k, h, letters| {
                num.clear();
                num.extend(letters.iter().map(|&b| t[b as usize]));
                f(crate::fasta::header_id(h), &num)
            },
        )
        .map(|_| ())
    }
}

/// Re-index the records of an existing index (e.g. with another seed pattern).
pub struct IndexSource<'a>(pub &'a SeedIndex);

impl<'a> RecordSource for IndexSource<'a> {
    fn for_each(&mut self, f: &mut dyn FnMut(&[u8], &[u8]) -> io::Result<()>) -> io::Result<()> {
        for i in 0..self.0.len() {
            f(self.0.id(i).as_bytes(), self.0.seq(i))?;
        }
        Ok(())
    }
}

/// Does `path` start with the index magic?
pub fn is_index_file(path: &str) -> bool {
    let mut b = [0u8; 8];
    match File::open(path) {
        Ok(mut f) => f.read_exact(&mut b).is_ok() && b == MAGIC,
        Err(_) => false,
    }
}

/// Records of an opened seed index (owned).
pub struct OwnedIndexSource(pub SeedIndex);

impl RecordSource for OwnedIndexSource {
    fn for_each(&mut self, f: &mut dyn FnMut(&[u8], &[u8]) -> io::Result<()>) -> io::Result<()> {
        IndexSource(&self.0).for_each(f)
    }
}

/// A record source for any supported sequence input (FASTA, gzip FASTA, MMseqs2
/// database prefix, seed index), streamed.
pub fn source_for(path: &str) -> io::Result<Box<dyn RecordSource>> {
    let aa2num = crate::alphabet::build_aa2num();
    match crate::seqio::detect_format(path)? {
        crate::seqio::InputFormat::SeedIndex => Ok(Box::new(OwnedIndexSource(SeedIndex::open(Path::new(path))?))),
        crate::seqio::InputFormat::MmseqsDb => Ok(Box::new(MmseqsDbSource {
            prefix: path.to_string(),
            aa2num,
        })),
        crate::seqio::InputFormat::Fasta { .. } => Ok(Box::new(FastaSource {
            path: path.to_string(),
            aa2num,
        })),
    }
}

/// Records streamed from a FASTA (plain/gzip) — never fully loaded.
pub struct FastaSource {
    pub path: String,
    pub aa2num: [u8; 256],
}

impl RecordSource for FastaSource {
    fn for_each(&mut self, f: &mut dyn FnMut(&[u8], &[u8]) -> io::Result<()>) -> io::Result<()> {
        let mut fr = crate::seqio::open_fasta(&self.path)?;
        let (mut h, mut s) = (Vec::new(), Vec::new());
        let mut num = Vec::new();
        let mut any = false;
        while fr.next_record(&mut h, &mut s)? {
            any = true;
            num.clear();
            num.extend(s.iter().map(|&b| self.aa2num[b as usize]));
            f(crate::fasta::header_id(&h), &num)?;
        }
        if !any {
            return Err(invalid(format!("{}: no FASTA records found (the input is empty)", self.path)));
        }
        Ok(())
    }
}

fn align_up(x: u64) -> u64 {
    x.div_ceil(ALIGN) * ALIGN
}

fn pad_to<W: Write>(w: &mut W, cur: u64, target: u64) -> io::Result<()> {
    let n = (target - cur) as usize;
    if n > 0 {
        w.write_all(&vec![0u8; n])?;
    }
    Ok(())
}

/// Temporary file removed on drop.
struct TempFile {
    path: PathBuf,
}

impl TempFile {
    fn new(dir: &Path, tag: &str) -> io::Result<(TempFile, File)> {
        static COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let n = COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let path = dir.join(format!(
            ".sabertooth-{}-{}-{}.tmp",
            std::process::id(),
            n,
            tag
        ));
        let f = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(true)
            .open(&path)
            .map_err(|e| io::Error::new(e.kind(), format!("{}: {}", path.display(), e)))?;
        Ok((TempFile { path }, f))
    }
}

impl Drop for TempFile {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.path);
    }
}

/// Encode the care positions of a window into a base-20 code; `None` if a
/// care position is not a real residue or is masked.
#[inline]
pub fn encode_window(seq: &[u8], start: usize, pattern: &SeedPattern, masked: Option<&[bool]>) -> Option<u64> {
    let mut code = 0u64;
    for &off in pattern.offsets() {
        let r = seq[start + off];
        if r as u64 >= A {
            return None;
        }
        if let Some(m) = masked {
            if m[start + off] {
                return None;
            }
        }
        code = code * A + r as u64;
    }
    Some(code)
}

/// Build an index from `source` into `out` (atomically).
pub fn build_index(
    source: &mut dyn RecordSource,
    fingerprint: SourceFingerprint,
    out: &Path,
    opts: &BuildOptions,
) -> io::Result<BuildStats> {
    if !cfg!(target_endian = "little") {
        return Err(invalid("the index format requires a little-endian platform".into()));
    }
    let w = opts.pattern.weight();
    crate::prefilter::check_k(w).map_err(invalid)?;
    let span = opts.pattern.span();
    let mut stats = BuildStats::default();

    let out_tmp = PathBuf::from(format!("{}.tmp-{}", out.display(), std::process::id()));
    let tmp_guard = TempFileGuard { path: out_tmp.clone() };
    let mut file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(true)
        .open(&out_tmp)
        .map_err(|e| io::Error::new(e.kind(), format!("{}: {}", out_tmp.display(), e)))?;

    // ---- stage A: stream sequences ----------------------------------------
    let t0 = Instant::now();
    let tdir = &opts.tmp_dir;
    let (_t_off, f_off) = TempFile::new(tdir, "seqoff")?;
    let (_t_ids, f_ids) = TempFile::new(tdir, "ids")?;
    let (_t_idoff, f_idoff) = TempFile::new(tdir, "idoff")?;
    let (_t_mask, f_mask) = TempFile::new(tdir, "mask")?;
    let mut w_res = HashingWriter::new(BufWriter::with_capacity(1 << 20, &mut file));
    w_res.inner.write_all(&vec![0u8; HEADER_SIZE as usize])?; // header placeholder (not hashed)
    w_res.hasher = Hasher64::default();
    w_res.written = 0;
    let mut w_off = HashingWriter::new(BufWriter::with_capacity(1 << 20, &f_off));
    let mut w_ids = HashingWriter::new(BufWriter::with_capacity(1 << 20, &f_ids));
    let mut w_idoff = HashingWriter::new(BufWriter::with_capacity(1 << 20, &f_idoff));
    let mut w_mask = BufWriter::with_capacity(1 << 20, &f_mask);
    let mut n_seqs: u64 = 0;
    let mut total: u64 = 0;
    let mut max_len: u64 = 0;
    let mut id_total: u64 = 0;
    let mut masked_res: u64 = 0;
    let mut maskbuf: Vec<bool> = Vec::new();
    let mut bitacc: u64 = 0;
    let mut nbits: u32 = 0;
    w_off.write_all(&0u64.to_le_bytes())?;
    w_idoff.write_all(&0u64.to_le_bytes())?;
    source.for_each(&mut |id: &[u8], num: &[u8]| -> io::Result<()> {
        if n_seqs >= u32::MAX as u64 {
            return Err(invalid(format!("more than {} sequences are not supported", u32::MAX)));
        }
        if num.len() as u64 > u32::MAX as u64 {
            return Err(invalid("sequence longer than u32::MAX residues".into()));
        }
        w_res.write_all(num)?;
        total += num.len() as u64;
        max_len = max_len.max(num.len() as u64);
        w_off.write_all(&total.to_le_bytes())?;
        let idv: &[u8] = if std::str::from_utf8(id).is_ok() { id } else { b"?" };
        w_ids.write_all(idv)?;
        id_total += idv.len() as u64;
        w_idoff.write_all(&id_total.to_le_bytes())?;
        // seed mask bits (1 = excluded)
        match opts.mask {
            Some(mp) => {
                masked_res += low_complexity_mask_into(num, mp, &mut maskbuf) as u64;
            }
            None => {
                maskbuf.clear();
                maskbuf.resize(num.len(), false);
            }
        }
        for &b in &maskbuf {
            bitacc |= (b as u64) << nbits;
            nbits += 1;
            if nbits == 64 {
                w_mask.write_all(&bitacc.to_le_bytes())?;
                bitacc = 0;
                nbits = 0;
            }
        }
        n_seqs += 1;
        Ok(())
    })?;
    if nbits > 0 {
        w_mask.write_all(&bitacc.to_le_bytes())?;
    }
    w_mask.flush()?;
    drop(w_mask);
    if n_seqs == 0 {
        return Err(invalid("cannot index an empty database".into()));
    }
    let res_hash = w_res.hasher.finish();
    let res_len = w_res.written;
    w_res.flush()?;
    drop(w_res);
    w_off.flush()?;
    w_ids.flush()?;
    w_idoff.flush()?;
    let (off_hash, off_len) = (w_off.hasher.finish(), w_off.written);
    let (ids_hash, ids_len) = (w_ids.hasher.finish(), w_ids.written);
    let (idoff_hash, idoff_len) = (w_idoff.hasher.finish(), w_idoff.written);
    drop((w_off, w_ids, w_idoff));

    let mut sections: Vec<(SectionKind, Section)> = Vec::new();
    sections.push((
        SectionKind::Residues,
        Section {
            offset: HEADER_SIZE,
            len: res_len,
            hash: res_hash,
        },
    ));
    // append side files
    let mut cur = HEADER_SIZE + res_len;
    file.seek(SeekFrom::Start(cur))?;
    for (kind, mut f, len, hash) in [
        (SectionKind::SeqOffsets, &f_off, off_len, off_hash),
        (SectionKind::IdBytes, &f_ids, ids_len, ids_hash),
        (SectionKind::IdOffsets, &f_idoff, idoff_len, idoff_hash),
    ] {
        let start = align_up(cur);
        pad_to(&mut file, cur, start)?;
        f.seek(SeekFrom::Start(0))?;
        let copied = io::copy(&mut f.take(len), &mut file)?;
        if copied != len {
            return Err(invalid("short copy while assembling index".into()));
        }
        sections.push((kind, Section { offset: start, len, hash }));
        cur = start + len;
    }
    file.flush()?;
    stats.seconds_sequences = t0.elapsed().as_secs_f64();

    // ---- map sequences + mask for the seeding passes -----------------------
    // SAFETY: the mapped temporary files are private to this build and are not
    // modified while mapped.
    let map = unsafe { MmapOptions::new().map(&file)? };
    let res = &map[HEADER_SIZE as usize..(HEADER_SIZE + res_len) as usize];
    let soff = sections[1].1;
    let seq_offsets: &[u64] = cast_u64(&map[soff.offset as usize..(soff.offset + soff.len) as usize])?;
    let mask_map = if f_mask.metadata()?.len() > 0 {
        Some(unsafe { MmapOptions::new().map(&f_mask)? })
    } else {
        None
    };
    let mask_words: &[u64] = match &mask_map {
        Some(m) => cast_u64(&m[..])?,
        None => &[],
    };
    let n = n_seqs as usize;

    // the bucket layout depends on the database size (known now)
    let p = prefix_len_for_size(w, total, opts.large_prefix_threshold) as usize;
    let suffix_space = Pow20::new((w - p) as u32);
    let num_buckets = A.pow(p as u32) as usize;

    // ---- stage B: count occurrences per bucket ------------------------------
    let t1 = Instant::now();
    let threads = rayon::current_num_threads().max(1);
    let chunks = seq_chunks(seq_offsets, threads * 4);
    let window_ok = |gpos: u64, start: usize, s: &[u8]| -> Option<u64> {
        // gpos = global residue position of window start
        let mut code = 0u64;
        for &off in opts.pattern.offsets() {
            let r = s[start + off];
            if r as u64 >= A {
                return None;
            }
            let g = gpos + off as u64;
            if !mask_words.is_empty() && (mask_words[(g >> 6) as usize] >> (g & 63)) & 1 == 1 {
                return None;
            }
            code = code * A + r as u64;
        }
        Some(code)
    };
    let bucket_counts: Vec<u64> = chunks
        .par_iter()
        .map(|&(a, b)| {
            let mut c = vec![0u64; num_buckets];
            for i in a..b {
                let (s0, s1) = (seq_offsets[i] as usize, seq_offsets[i + 1] as usize);
                let s = &res[s0..s1];
                if s.len() < span {
                    continue;
                }
                for st in 0..=(s.len() - span) {
                    if let Some(code) = window_ok(s0 as u64 + st as u64, st, s) {
                        c[suffix_space.quot(code) as usize] += 1;
                    }
                }
            }
            c
        })
        .reduce(
            || vec![0u64; num_buckets],
            |mut x, y| {
                for (a, b) in x.iter_mut().zip(y.iter()) {
                    *a += *b;
                }
                x
            },
        );
    let num_postings: u64 = bucket_counts.iter().sum();
    stats.seconds_count = t1.elapsed().as_secs_f64();

    // ---- posting encoding: (target << pos_bits) | position, `pw` bytes -----
    let bits_for = |x: u64| 64 - x.leading_zeros();
    let pos_bits = bits_for(max_len.saturating_sub(1));
    let need_bits = bits_for(n_seqs - 1) + pos_bits;
    let pw: u32 = match opts.posting_width {
        Some(wd) => {
            if !(1..=8).contains(&wd) || wd * 8 < need_bits {
                return Err(invalid(format!(
                    "posting width {} bytes cannot hold {} bits (targets and positions)",
                    wd, need_bits
                )));
            }
            wd
        }
        None => need_bits.div_ceil(8).clamp(1, 8),
    };
    let encode_posting = |p: u64| -> [u8; 8] { (((p >> 32) << pos_bits) | (p & 0xFFFF_FFFF)).to_le_bytes() };

    // ---- stage C: partitions of buckets -> postings + k-mer table ----------
    let t2 = Instant::now();
    let post_off = align_up(cur);
    pad_to(&mut file, cur, post_off)?;
    let (_t_kmer, f_kmer) = TempFile::new(tdir, "kmer")?;
    let mut w_kmer = HashingWriter::new(BufWriter::with_capacity(1 << 20, &f_kmer));
    let mut w_post = HashingWriter::new(BufWriter::with_capacity(1 << 20, &mut file));
    let mut bucket_kmer_start = vec![0u64; num_buckets + 1];
    let mut bucket_post_start = vec![0u64; num_buckets + 1];
    let per_entry = 16usize; // (u32 suffix, u64 posting) with padding
    let budget_entries = (opts.memory_budget / per_entry).max(1 << 16) as u64;
    let mut kmers: u64 = 0;
    let mut posts: u64 = 0;
    let mut b_lo = 0usize;
    let mut max_part = 0usize;
    while b_lo < num_buckets {
        let mut b_hi = b_lo;
        let mut part: u64 = 0;
        while b_hi < num_buckets && (part + bucket_counts[b_hi] <= budget_entries || b_hi == b_lo) {
            part += bucket_counts[b_hi];
            b_hi += 1;
        }
        stats.partitions += 1;
        let nb = b_hi - b_lo;
        let part_n = part as usize;
        max_part = max_part.max(part_n);
        let mut buf: Vec<(u32, u64)> = vec![(0, 0); part_n];
        if part_n > 0 {
            // per-chunk counts within the partition
            let chunk_counts: Vec<Vec<u32>> = chunks
                .par_iter()
                .map(|&(a, b)| {
                    let mut c = vec![0u32; nb];
                    for i in a..b {
                        let (s0, s1) = (seq_offsets[i] as usize, seq_offsets[i + 1] as usize);
                        let s = &res[s0..s1];
                        if s.len() < span {
                            continue;
                        }
                        for st in 0..=(s.len() - span) {
                            if let Some(code) = window_ok(s0 as u64 + st as u64, st, s) {
                                let bk = suffix_space.quot(code) as usize;
                                if bk >= b_lo && bk < b_hi {
                                    c[bk - b_lo] += 1;
                                }
                            }
                        }
                    }
                    c
                })
                .collect();
            // write cursors: bucket-major, chunk-minor (keeps (seq,pos) order)
            let mut starts: Vec<Vec<usize>> = vec![vec![0usize; nb]; chunks.len()];
            let mut acc = 0usize;
            for bk in 0..nb {
                for (ci, cc) in chunk_counts.iter().enumerate() {
                    starts[ci][bk] = acc;
                    acc += cc[bk] as usize;
                }
            }
            debug_assert_eq!(acc, part_n);
            // scatter into disjoint slots
            let buf_ptr = SyncPtr(buf.as_mut_ptr());
            chunks.par_iter().zip(starts.into_par_iter()).for_each(|(&(a, b), mut cursor)| {
                let bp = buf_ptr;
                for i in a..b {
                    let (s0, s1) = (seq_offsets[i] as usize, seq_offsets[i + 1] as usize);
                    let s = &res[s0..s1];
                    if s.len() < span {
                        continue;
                    }
                    for st in 0..=(s.len() - span) {
                        if let Some(code) = window_ok(s0 as u64 + st as u64, st, s) {
                            let bk = suffix_space.quot(code) as usize;
                            if bk >= b_lo && bk < b_hi {
                                let slot = cursor[bk - b_lo];
                                cursor[bk - b_lo] += 1;
                                // SAFETY: every (chunk, bucket) owns a disjoint
                                // slot range computed by the prefix sums above.
                                unsafe {
                                    *bp.0.add(slot) = (suffix_space.modulo(code) as u32, ((i as u64) << 32) | st as u64);
                                }
                            }
                        }
                    }
                }
            });
        }
        // per-bucket sort by suffix (posting order preserved: postings unique
        // and already increasing within a bucket)
        let mut bstart = Vec::with_capacity(nb + 1);
        let mut acc = 0usize;
        for bk in b_lo..b_hi {
            bstart.push(acc);
            acc += bucket_counts[bk] as usize;
        }
        bstart.push(acc);
        {
            let mut rest: &mut [(u32, u64)] = &mut buf;
            let mut slices: Vec<&mut [(u32, u64)]> = Vec::with_capacity(nb);
            for k in 0..nb {
                let len = bstart[k + 1] - bstart[k];
                let (a, b) = rest.split_at_mut(len);
                slices.push(a);
                rest = b;
            }
            slices.par_iter_mut().for_each(|s| s.sort_unstable());
        }
        // emit
        for (k, bk) in (b_lo..b_hi).enumerate() {
            bucket_kmer_start[bk] = kmers;
            bucket_post_start[bk] = posts;
            let sl = &buf[bstart[k]..bstart[k + 1]];
            let mut j = 0;
            while j < sl.len() {
                let suf = sl[j].0;
                let rel: u32 = u32::try_from(j)
                    .map_err(|_| invalid("bucket with more than 2^32 occurrences".into()))?;
                w_kmer.write_all(&suf.to_le_bytes())?;
                w_kmer.write_all(&rel.to_le_bytes())?;
                kmers += 1;
                while j < sl.len() && sl[j].0 == suf {
                    w_post.write_all(&encode_posting(sl[j].1)[..pw as usize])?;
                    j += 1;
                }
            }
            posts += sl.len() as u64;
        }
        drop(buf);
        b_lo = b_hi;
    }
    bucket_kmer_start[num_buckets] = kmers;
    bucket_post_start[num_buckets] = posts;
    debug_assert_eq!(posts, num_postings);
    let post_hash = w_post.hasher.finish();
    let post_len = w_post.written;
    w_post.flush()?;
    drop(w_post);
    w_kmer.flush()?;
    let (kmer_hash, kmer_len) = (w_kmer.hasher.finish(), w_kmer.written);
    drop(w_kmer);
    drop(map);
    drop(mask_map);
    sections.push((
        SectionKind::Postings,
        Section {
            offset: post_off,
            len: post_len,
            hash: post_hash,
        },
    ));
    stats.seconds_postings = t2.elapsed().as_secs_f64();

    // ---- stage D: lookup tables, header, footer -----------------------------
    let t3 = Instant::now();
    cur = post_off + post_len;
    file.seek(SeekFrom::Start(cur))?;
    // decoders read 8 bytes per posting: keep 8 bytes of slack after the section
    let codes = A.pow(w as u32);
    let ow: u64 = if num_postings < u32::MAX as u64 { 4 } else { 8 };
    let direct_ok = w <= 6;
    let ss = suffix_space.value();
    // offsets of every code 0..=codes (the last is the sentinel), streamed in
    // code order from the (bucket, suffix, offset) entries
    let dense_offsets = |f: &mut dyn FnMut(u64, u64) -> io::Result<()>| -> io::Result<()> {
        let mut fk = io::BufReader::with_capacity(1 << 20, &f_kmer);
        fk.seek(SeekFrom::Start(0))?;
        let mut next: u64 = 0;
        let mut e = [0u8; 8];
        for bk in 0..num_buckets {
            let base_post = bucket_post_start[bk];
            for _ in bucket_kmer_start[bk]..bucket_kmer_start[bk + 1] {
                fk.read_exact(&mut e)?;
                let suf = u32::from_le_bytes([e[0], e[1], e[2], e[3]]) as u64;
                let rel = u32::from_le_bytes([e[4], e[5], e[6], e[7]]) as u64;
                let code = bk as u64 * ss + suf;
                while next <= code {
                    f(next, base_post + rel)?;
                    next += 1;
                }
            }
        }
        while next <= codes {
            f(next, num_postings)?;
            next += 1;
        }
        Ok(())
    };
    let ngroups = (codes + 1).div_ceil(OFFSET_GROUP as u64);
    let grouped_bytes = 2 * (codes + 1) + ow * ngroups;
    let sparse_bytes = kmer_len + 16 * (num_buckets as u64 + 1);
    let auto_direct = |direct_bytes: u64| direct_bytes <= sparse_bytes || sparse_bytes >= AUTO_DIRECT_MIN_SPARSE_BYTES;
    let may_direct = direct_ok
        && match opts.layout {
            LayoutChoice::Sparse => false,
            LayoutChoice::Direct => true,
            LayoutChoice::Auto => auto_direct(grouped_bytes),
        };
    // grouped offsets need every group's postings (but its last code's) to
    // fit 16 bits
    let grouped = may_direct
        && opts.direct_grouped
        && {
            let mut ok = true;
            let mut gbase = 0u64;
            dense_offsets(&mut |c, off| {
                if (c as usize).is_multiple_of(OFFSET_GROUP) {
                    gbase = off;
                }
                ok &= off - gbase <= u16::MAX as u64;
                Ok(())
            })?;
            ok
        };
    let direct_bytes = if grouped { grouped_bytes } else { (codes + 1) * ow };
    let layout = match opts.layout {
        LayoutChoice::Sparse => IndexLayout::Sparse,
        LayoutChoice::Direct => {
            if !direct_ok {
                return Err(invalid(format!("the direct layout supports seed weights up to 6 (got {})", w)));
            }
            IndexLayout::Direct
        }
        LayoutChoice::Auto => {
            if may_direct && auto_direct(direct_bytes) {
                IndexLayout::Direct
            } else {
                IndexLayout::Sparse
            }
        }
    };
    let kstart = align_up(cur + 8);
    pad_to(&mut file, cur, kstart)?;
    if layout == IndexLayout::Sparse {
        {
            let mut fk = &f_kmer;
            fk.seek(SeekFrom::Start(0))?;
            let copied = io::copy(&mut fk.take(kmer_len), &mut file)?;
            if copied != kmer_len {
                return Err(invalid("short copy of k-mer table".into()));
            }
        }
        sections.push((
            SectionKind::KmerTable,
            Section {
                offset: kstart,
                len: kmer_len,
                hash: kmer_hash,
            },
        ));
        cur = kstart + kmer_len;
        for (kind, table) in [
            (SectionKind::BucketKmerStart, &bucket_kmer_start),
            (SectionKind::BucketPostStart, &bucket_post_start),
        ] {
            let start = align_up(cur);
            pad_to(&mut file, cur, start)?;
            let mut hw = HashingWriter::new(BufWriter::with_capacity(1 << 20, &mut file));
            for v in table.iter() {
                hw.write_all(&v.to_le_bytes())?;
            }
            hw.flush()?;
            let (h, l) = (hw.hasher.finish(), hw.written);
            drop(hw);
            sections.push((kind, Section { offset: start, len: l, hash: h }));
            cur = start + l;
        }
        stats.lookup_bytes = kmer_len + 16 * (num_buckets as u64 + 1);
    } else {
        // offsets[c] = postings of all codes < c, for every code and a final
        // sentinel: full-width entries, or 16-bit offsets relative to the
        // first code of each group plus a table of group bases
        let mut hw = HashingWriter::new(BufWriter::with_capacity(1 << 20, &mut file));
        let mut chunk: Vec<u8> = Vec::with_capacity(1 << 20);
        let mut bases: Vec<u64> = Vec::with_capacity(if grouped { ngroups as usize } else { 0 });
        dense_offsets(&mut |c, off| {
            if grouped {
                if (c as usize).is_multiple_of(OFFSET_GROUP) {
                    bases.push(off);
                }
                let rel = off - *bases.last().unwrap();
                chunk.extend_from_slice(&(rel as u16).to_le_bytes());
            } else if ow == 4 {
                chunk.extend_from_slice(&(off as u32).to_le_bytes());
            } else {
                chunk.extend_from_slice(&off.to_le_bytes());
            }
            if chunk.len() >= (1 << 20) - 8 {
                hw.write_all(&chunk)?;
                chunk.clear();
            }
            Ok(())
        })?;
        hw.write_all(&chunk)?;
        hw.flush()?;
        let (h, l) = (hw.hasher.finish(), hw.written);
        drop(hw);
        sections.push((SectionKind::KmerOffsets, Section { offset: kstart, len: l, hash: h }));
        cur = kstart + l;
        stats.lookup_bytes = l;
        if grouped {
            debug_assert_eq!(bases.len() as u64, ngroups);
            let start = align_up(cur);
            pad_to(&mut file, cur, start)?;
            let mut hw = HashingWriter::new(BufWriter::with_capacity(1 << 20, &mut file));
            for &b in &bases {
                if ow == 4 {
                    hw.write_all(&(b as u32).to_le_bytes())?;
                } else {
                    hw.write_all(&b.to_le_bytes())?;
                }
            }
            hw.flush()?;
            let (h, l) = (hw.hasher.finish(), hw.written);
            drop(hw);
            sections.push((SectionKind::KmerGroupBase, Section { offset: start, len: l, hash: h }));
            cur = start + l;
            stats.lookup_bytes += l;
        }
        debug_assert_eq!(stats.lookup_bytes, direct_bytes);
    }
    stats.bucket_table_bytes = (num_buckets + 1) * 16 + num_buckets * 8;
    stats.layout = layout.name();
    stats.posting_width = pw;
    // footer
    let file_len = cur + FOOTER_SIZE;
    file.write_all(&FOOTER_MAGIC)?;
    file.write_all(&file_len.to_le_bytes())?;
    let meta = IndexMeta {
        version: FORMAT_VERSION,
        alphabet: ALPHABET_AA21,
        pattern: opts.pattern.clone(),
        prefix_len: p as u32,
        masked: opts.mask.is_some(),
        mask: opts.mask.unwrap_or_default(),
        num_seqs: n as u64,
        total_residues: total,
        max_seq_len: max_len,
        num_postings,
        num_kmers: kmers,
        num_buckets: num_buckets as u64,
        masked_residues: masked_res,
        source: fingerprint,
        created_unix: std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_secs())
            .unwrap_or(0),
        tool_version: crate::VERSION.to_string(),
        sections,
        layout,
        posting_width: pw,
        posting_pos_bits: pos_bits,
        offset_width: match layout {
            IndexLayout::Direct if grouped => 2,
            IndexLayout::Direct => ow as u32,
            IndexLayout::Sparse => 0,
        },
    };
    file.seek(SeekFrom::Start(0))?;
    file.write_all(&encode_header(&meta))?;
    file.flush()?;
    file.sync_all()?;
    drop(file);
    fs::rename(&out_tmp, out).map_err(|e| io::Error::new(e.kind(), format!("{}: {}", out.display(), e)))?;
    std::mem::forget(tmp_guard);
    stats.seconds_finish = t3.elapsed().as_secs_f64();
    stats.num_seqs = n as u64;
    stats.total_residues = total;
    stats.num_postings = num_postings;
    stats.num_kmers = kmers;
    stats.masked_residues = masked_res;
    stats.partition_buffer_bytes = max_part * per_entry;
    stats.file_bytes = file_len;
    Ok(stats)
}

/// Removes the in-progress output on error paths.
struct TempFileGuard {
    path: PathBuf,
}

impl Drop for TempFileGuard {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.path);
    }
}

#[derive(Clone, Copy)]
struct SyncPtr(*mut (u32, u64));
// SAFETY: used only for disjoint writes from rayon tasks (see scatter above).
unsafe impl Send for SyncPtr {}
unsafe impl Sync for SyncPtr {}

/// Split `n` sequences into up to `parts` contiguous chunks of similar residue
/// counts.
fn seq_chunks(offsets: &[u64], parts: usize) -> Vec<(usize, usize)> {
    let n = offsets.len() - 1;
    let total = offsets[n].max(1);
    let parts = parts.max(1).min(n.max(1));
    let mut out = Vec::with_capacity(parts);
    let mut a = 0usize;
    for k in 1..=parts {
        let target = total * k as u64 / parts as u64;
        let mut b = a;
        while b < n && (offsets[b] < target || b == a) {
            b += 1;
        }
        if k == parts {
            b = n;
        }
        if b > a {
            out.push((a, b));
        }
        a = b;
        if a >= n {
            break;
        }
    }
    if out.is_empty() {
        out.push((0, n));
    }
    out
}

fn cast_u64(b: &[u8]) -> io::Result<&[u64]> {
    // SAFETY: u64 has no invalid bit patterns; alignment is checked below.
    let (pre, mid, post) = unsafe { b.align_to::<u64>() };
    if !pre.is_empty() || !post.is_empty() {
        return Err(invalid("misaligned u64 section".into()));
    }
    Ok(mid)
}

fn cast_pairs(b: &[u8]) -> io::Result<&[[u32; 2]]> {
    // SAFETY: [u32; 2] has no invalid bit patterns; alignment checked below.
    let (pre, mid, post) = unsafe { b.align_to::<[u32; 2]>() };
    if !pre.is_empty() || !post.is_empty() {
        return Err(invalid("misaligned k-mer table".into()));
    }
    Ok(mid)
}

// ---------------------------------------------------------------------------
// open / lookup
// ---------------------------------------------------------------------------

/// An opened (memory-mapped) seed index.
pub struct SeedIndex {
    map: Mmap,
    meta: IndexMeta,
    path: PathBuf,
    suffix_space: Pow20,
    res: (usize, usize),
    seq_off: (usize, usize),
    id_bytes: (usize, usize),
    id_off: (usize, usize),
    postings: (usize, usize),
    /// sparse layout: k-mer table and bucket tables (empty for direct)
    kmer: (usize, usize),
    bk: (usize, usize),
    bp: (usize, usize),
    /// direct layout: offset table (empty for sparse) and its entry width
    /// (2 = grouped 16-bit offsets relative to the group bases `kgrp`)
    koff: (usize, usize),
    koff_w: usize,
    kgrp: (usize, usize),
    kgrp_w: usize,
    direct: bool,
    /// posting decoding: `width` bytes, `(target << pos_bits) | position`
    post_w: usize,
    post_mask: u64,
    pos_bits: u32,
    pos_mask: u64,
    _transient: Option<TempFileGuard>,
    /// Lazily built presence bitmap (see [`SeedIndex::presence`]).
    presence: std::sync::OnceLock<Presence>,
}

/// Bitmap of the k-mer codes present in an index, over the first
/// `min(w, PRESENCE_MAX_RESIDUES)` residues of a code: exact presence for seed
/// weights ≤ 6, a necessary condition (prefix present) beyond. At most
/// `20^6` bits = 8 MiB, built once in `O(distinct k-mers)`.
pub struct Presence {
    bits: Vec<u64>,
    /// code → bit index divisor (`20^(w − residues)`)
    div: Pow20,
}

/// Residues covered by the presence bitmap.
pub const PRESENCE_MAX_RESIDUES: u32 = 6;

impl Presence {
    /// Whether `code` may be present (exact for weights ≤ 6).
    #[inline]
    pub fn may_contain(&self, code: u64) -> bool {
        // b < 20^residues ≤ 64·bits.len() for every valid code
        let b = self.div.quot(code);
        (self.bits[(b >> 6) as usize] >> (b & 63)) & 1 == 1
    }

    /// Bytes held.
    pub fn bytes(&self) -> usize {
        self.bits.len() * 8
    }
}

impl std::fmt::Debug for SeedIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "SeedIndex({}, {} seqs)", self.path.display(), self.meta.num_seqs)
    }
}

/// Why an existing index cannot be used for a search.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Unusable {
    Stale(String),
    Incompatible(String),
}

impl std::fmt::Display for Unusable {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Unusable::Stale(s) => write!(f, "stale index: {}", s),
            Unusable::Incompatible(s) => write!(f, "incompatible index: {}", s),
        }
    }
}

impl SeedIndex {
    /// Open and validate an index file (O(1) checks; see module docs).
    pub fn open(path: &Path) -> io::Result<SeedIndex> {
        let p = path.display().to_string();
        let file = File::open(path).map_err(|e| io::Error::new(e.kind(), format!("{}: {}", p, e)))?;
        let flen = file.metadata()?.len();
        if flen < HEADER_SIZE + FOOTER_SIZE {
            return Err(invalid(format!("{}: too small to be an index (incomplete build?)", p)));
        }
        // SAFETY: index files are immutable after their atomic rename; callers
        // must not modify an index while it is in use.
        let map = unsafe { MmapOptions::new().map(&file)? };
        // Ask for transparent huge pages: when the kernel reads the file into
        // large folios (Linux, file THP), random seed lookups take far fewer
        // TLB misses. Advisory only; harmless where unsupported.
        #[cfg(target_os = "linux")]
        {
            let _ = map.advise(memmap2::Advice::HugePage);
        }
        let meta = decode_header(&map[..HEADER_SIZE as usize], &p)?;
        let foot = &map[(flen - FOOTER_SIZE) as usize..];
        if foot[..8] != FOOTER_MAGIC
            || u64::from_le_bytes(foot[8..16].try_into().unwrap()) != flen
        {
            return Err(invalid(format!("{}: missing completion footer (incomplete or truncated index)", p)));
        }
        if meta.alphabet != ALPHABET_AA21 {
            return Err(invalid(format!("{}: unsupported alphabet id {}", p, meta.alphabet)));
        }
        let get = |k: SectionKind| -> io::Result<(usize, usize)> {
            let s = meta.section(k);
            let end = s
                .offset
                .checked_add(s.len)
                .ok_or_else(|| invalid(format!("{}: section {:?} overflows", p, k)))?;
            if s.offset < HEADER_SIZE || end > flen - FOOTER_SIZE {
                return Err(invalid(format!("{}: section {:?} out of bounds", p, k)));
            }
            Ok((s.offset as usize, end as usize))
        };
        let res = get(SectionKind::Residues)?;
        let seq_off = get(SectionKind::SeqOffsets)?;
        let id_bytes = get(SectionKind::IdBytes)?;
        let id_off = get(SectionKind::IdOffsets)?;
        let postings = get(SectionKind::Postings)?;
        let direct = meta.layout == IndexLayout::Direct;
        let (kmer, bk, bp, koff) = if direct {
            ((0, 0), (0, 0), (0, 0), get(SectionKind::KmerOffsets)?)
        } else {
            (
                get(SectionKind::KmerTable)?,
                get(SectionKind::BucketKmerStart)?,
                get(SectionKind::BucketPostStart)?,
                (0, 0),
            )
        };
        // alignment / size validation once; accessors then build slices directly
        for r in [seq_off, id_off, bk, bp] {
            cast_u64(&map[r.0..r.1]).map_err(|e| invalid(format!("{}: {}", p, e)))?;
        }
        cast_pairs(&map[kmer.0..kmer.1]).map_err(|e| invalid(format!("{}: {}", p, e)))?;
        let koff_w = meta.offset_width as usize;
        let (kgrp, kgrp_w) = if direct && koff_w == 2 {
            (get(SectionKind::KmerGroupBase)?, if meta.num_postings < u32::MAX as u64 { 4 } else { 8 })
        } else {
            ((0, 0), 1)
        };
        if direct && (koff.0 % koff_w != 0 || kgrp.0 % kgrp_w != 0) {
            return Err(invalid(format!("{}: misaligned k-mer offset table", p)));
        }
        let post_w = meta.posting_width as usize;
        // 8-byte unaligned posting reads need 8 bytes after the section's start
        // of the last entry (writers leave the slack)
        if post_w < 8 && meta.num_postings > 0 && postings.1 + 8 > (flen - FOOTER_SIZE) as usize {
            return Err(invalid(format!("{}: posting section lacks read slack (corrupt index)", p)));
        }
        let idx = SeedIndex {
            suffix_space: Pow20::new(meta.pattern.weight() as u32 - meta.prefix_len),
            post_w,
            post_mask: if post_w >= 8 { u64::MAX } else { (1u64 << (8 * post_w)) - 1 },
            pos_bits: meta.posting_pos_bits,
            pos_mask: (1u64 << meta.posting_pos_bits) - 1,
            map,
            meta,
            path: path.to_path_buf(),
            res,
            seq_off,
            id_bytes,
            id_off,
            postings,
            kmer,
            bk,
            bp,
            koff,
            koff_w,
            kgrp,
            kgrp_w,
            direct,
            _transient: None,
            presence: std::sync::OnceLock::new(),
        };
        let m = &idx.meta;
        let n = m.num_seqs as usize;
        let w = m.pattern.weight() as u32;
        let mut checks: Vec<(bool, &str)> = vec![
            ((idx.res.1 - idx.res.0) as u64 == m.total_residues, "residue count"),
            (idx.seq_offsets().len() == n + 1, "sequence offsets"),
            (idx.id_offsets().len() == n + 1, "identifier offsets"),
            (
                (idx.postings.1 - idx.postings.0) as u64 == m.num_postings * post_w as u64,
                "posting count",
            ),
            (*idx.seq_offsets().last().unwrap() == m.total_residues, "sequence offsets"),
            (
                *idx.id_offsets().last().unwrap() == (idx.id_bytes.1 - idx.id_bytes.0) as u64,
                "identifier offsets",
            ),
        ];
        if direct {
            let entries = (koff.1 - koff.0) / koff_w.max(1);
            let groups = (kgrp.1 - kgrp.0) / kgrp_w;
            let sized = w <= 6
                && entries as u64 == A.pow(w) + 1
                && (koff_w != 2 || groups == entries.div_ceil(OFFSET_GROUP));
            checks.push((sized, "k-mer offset table"));
            // the sentinel is read only from a table of the right size
            checks.push((sized && idx.kmer_offset(entries - 1) == m.num_postings, "k-mer offset table"));
        } else {
            checks.extend([
                (idx.kmer_table().len() as u64 == m.num_kmers, "k-mer count"),
                (idx.bucket_kmer_start().len() as u64 == m.num_buckets + 1, "bucket table"),
                (idx.bucket_post_start().len() as u64 == m.num_buckets + 1, "bucket table"),
                (m.num_buckets == A.pow(m.prefix_len), "bucket count"),
                (*idx.bucket_kmer_start().last().unwrap() == m.num_kmers, "bucket table"),
                (*idx.bucket_post_start().last().unwrap() == m.num_postings, "bucket table"),
            ]);
        }
        for (ok, what) in checks {
            if !ok {
                return Err(invalid(format!("{}: inconsistent {} (corrupt index)", p, what)));
            }
        }
        Ok(idx)
    }

    /// Build a transient index (temporary file, removed when dropped) from an
    /// in-memory database.
    pub fn build_transient(db: &SeqDb, opts: &BuildOptions) -> io::Result<(SeedIndex, BuildStats)> {
        static COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let k = COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let path = opts
            .tmp_dir
            .join(format!(".sabertooth-{}-{}-transient.stidx", std::process::id(), k));
        let stats = build_index(&mut DbSource(db), SourceFingerprint::transient(), &path, opts)?;
        let mut idx = SeedIndex::open(&path)?;
        // remove now on Unix (the mapping stays valid); elsewhere on drop
        #[cfg(unix)]
        {
            let _ = fs::remove_file(&path);
        }
        idx._transient = Some(TempFileGuard { path });
        Ok((idx, stats))
    }

    /// Build a transient index from any record source.
    pub fn build_transient_from(source: &mut dyn RecordSource, opts: &BuildOptions) -> io::Result<(SeedIndex, BuildStats)> {
        static COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let k = COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let path = opts
            .tmp_dir
            .join(format!(".sabertooth-{}-{}-transient-src.stidx", std::process::id(), k));
        let stats = build_index(source, SourceFingerprint::transient(), &path, opts)?;
        let mut idx = SeedIndex::open(&path)?;
        #[cfg(unix)]
        {
            let _ = fs::remove_file(&path);
        }
        idx._transient = Some(TempFileGuard { path });
        Ok((idx, stats))
    }

    /// Whether this index lives only in a removed temporary file.
    pub fn is_transient(&self) -> bool {
        self._transient.is_some()
    }

    /// Check that this index can serve a search with `pattern` / `masked`, and
    /// (if `source` is given) that it is not stale relative to it.
    pub fn usable_for(&self, pattern: &SeedPattern, masked: Option<MaskParams>, source: Option<&str>) -> Result<(), Unusable> {
        if self.meta.pattern.mask_string() != pattern.mask_string() {
            return Err(Unusable::Incompatible(format!(
                "built with seed {} but the search uses {}",
                self.meta.pattern.mask_string(),
                pattern.mask_string()
            )));
        }
        let want = masked.is_some();
        if self.meta.masked != want || (want && self.meta.mask != masked.unwrap()) {
            return Err(Unusable::Incompatible(format!(
                "built with seed masking {} but the search requests {}",
                if self.meta.masked { "on" } else { "off" },
                if want { "on" } else { "off" }
            )));
        }
        if let Some(src) = source {
            match SourceFingerprint::of(src) {
                Ok(fp) => {
                    if fp.files != self.meta.source.files || fp.kind != self.meta.source.kind {
                        return Err(Unusable::Stale(format!(
                            "{} changed after the index was built (size/mtime differ); rebuild with `sabertooth createindex {}`",
                            src, src
                        )));
                    }
                }
                Err(e) => return Err(Unusable::Stale(e.to_string())),
            }
        }
        Ok(())
    }

    /// Verify every section checksum (full read of the file).
    pub fn verify(&self) -> io::Result<()> {
        for (k, s) in &self.meta.sections {
            let mut h = Hasher64::default();
            h.update(&self.map[s.offset as usize..(s.offset + s.len) as usize]);
            if h.finish() != s.hash {
                return Err(invalid(format!("{}: section {:?} checksum mismatch", self.path.display(), k)));
            }
        }
        Ok(())
    }

    pub fn meta(&self) -> &IndexMeta {
        &self.meta
    }
    pub fn path(&self) -> &Path {
        &self.path
    }
    pub fn pattern(&self) -> &SeedPattern {
        &self.meta.pattern
    }
    /// Effective k (seed weight).
    pub fn k(&self) -> usize {
        self.meta.pattern.weight()
    }
    pub fn len(&self) -> usize {
        self.meta.num_seqs as usize
    }
    pub fn is_empty(&self) -> bool {
        self.meta.num_seqs == 0
    }
    pub fn total_residues(&self) -> usize {
        self.meta.total_residues as usize
    }
    pub fn distinct_kmers(&self) -> usize {
        self.meta.num_kmers as usize
    }
    /// Bytes of the mapped file.
    pub fn mapped_bytes(&self) -> usize {
        self.map.len()
    }

    #[inline]
    fn residues(&self) -> &[u8] {
        &self.map[self.res.0..self.res.1]
    }
    #[inline]
    fn u64s(&self, r: (usize, usize)) -> &[u64] {
        // SAFETY: alignment and length were validated in `open`; the map lives
        // as long as `self` and is never written.
        unsafe { std::slice::from_raw_parts(self.map.as_ptr().add(r.0) as *const u64, (r.1 - r.0) / 8) }
    }
    #[inline]
    fn seq_offsets(&self) -> &[u64] {
        self.u64s(self.seq_off)
    }
    #[inline]
    fn id_offsets(&self) -> &[u64] {
        self.u64s(self.id_off)
    }
    /// Entry `c` of the direct offset table.
    #[inline]
    fn kmer_offset(&self, c: usize) -> u64 {
        let base = self.map.as_ptr();
        // SAFETY: c < table entries (callers pass valid codes / the sentinel);
        // the tables were validated in `open` and are aligned to their widths
        unsafe {
            match self.koff_w {
                2 => {
                    let g = c / OFFSET_GROUP;
                    let gb = if self.kgrp_w == 4 {
                        *(base.add(self.kgrp.0 + g * 4) as *const u32) as u64
                    } else {
                        *(base.add(self.kgrp.0 + g * 8) as *const u64)
                    };
                    gb + *(base.add(self.koff.0 + c * 2) as *const u16) as u64
                }
                4 => *(base.add(self.koff.0 + c * 4) as *const u32) as u64,
                _ => *(base.add(self.koff.0 + c * 8) as *const u64),
            }
        }
    }
    #[inline]
    fn kmer_table(&self) -> &[[u32; 2]] {
        // SAFETY: as for `u64s` (validated in `open`).
        unsafe {
            std::slice::from_raw_parts(
                self.map.as_ptr().add(self.kmer.0) as *const [u32; 2],
                (self.kmer.1 - self.kmer.0) / 8,
            )
        }
    }
    #[inline]
    fn bucket_kmer_start(&self) -> &[u64] {
        self.u64s(self.bk)
    }
    #[inline]
    fn bucket_post_start(&self) -> &[u64] {
        self.u64s(self.bp)
    }

    /// Numeric residues of target `i`.
    #[inline]
    pub fn seq(&self, i: usize) -> &[u8] {
        let o = self.seq_offsets();
        &self.residues()[o[i] as usize..o[i + 1] as usize]
    }

    /// Hint the CPU to fetch the residues of target `i` around position `pos`
    /// (prefilter look-ahead; no effect on results).
    #[inline]
    pub fn prefetch_residues(&self, i: usize, pos: usize) {
        #[cfg(target_arch = "x86_64")]
        {
            let o = self.seq_offsets();
            if i + 1 < o.len() {
                let start = o[i] as usize + pos;
                let res = self.residues();
                if start < res.len() {
                    // SAFETY: prefetch has no architectural effect; the address is in bounds
                    unsafe {
                        std::arch::x86_64::_mm_prefetch::<{ std::arch::x86_64::_MM_HINT_T0 }>(res.as_ptr().add(start) as *const i8);
                    }
                }
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            let _ = (i, pos);
        }
    }

    /// Hint the CPU to fetch the offset-table entry of target `i`.
    #[inline]
    pub fn prefetch_offset(&self, i: usize) {
        #[cfg(target_arch = "x86_64")]
        {
            let o = self.seq_offsets();
            if i < o.len() {
                // SAFETY: prefetch has no architectural effect; the address is in bounds
                unsafe {
                    std::arch::x86_64::_mm_prefetch::<{ std::arch::x86_64::_MM_HINT_T0 }>(o.as_ptr().add(i) as *const i8);
                }
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            let _ = i;
        }
    }

    /// Length of target `i`.
    #[inline]
    pub fn seq_len(&self, i: usize) -> usize {
        let o = self.seq_offsets();
        (o[i + 1] - o[i]) as usize
    }

    /// Identifier of target `i`.
    pub fn id(&self, i: usize) -> &str {
        let o = self.id_offsets();
        let b = &self.map[self.id_bytes.0 + o[i] as usize..self.id_bytes.0 + o[i + 1] as usize];
        std::str::from_utf8(b).unwrap_or("?")
    }

    /// Lengths of every `step`-th target (statistics sampling).
    pub fn length_sample(&self, max: usize) -> Vec<usize> {
        let n = self.len();
        let step = (n / max.max(1)).max(1);
        (0..n).step_by(step).map(|i| self.seq_len(i)).collect()
    }

    /// The k-mer presence bitmap (built on first use; `O(distinct k-mers)`,
    /// or one pass over the offset table for the direct layout).
    pub fn presence(&self) -> &Presence {
        self.presence.get_or_init(|| {
            let w = self.meta.pattern.weight() as u32;
            let r = w.min(PRESENCE_MAX_RESIDUES);
            let nbits = A.pow(r);
            let div = Pow20::new(w - r);
            let mut bits = vec![0u64; nbits.div_ceil(64) as usize];
            if self.direct {
                let codes = A.pow(w) as usize;
                let mut prev = self.kmer_offset(0);
                for c in 0..codes {
                    let next = self.kmer_offset(c + 1);
                    if next > prev {
                        let i = div.quot(c as u64);
                        bits[(i >> 6) as usize] |= 1 << (i & 63);
                    }
                    prev = next;
                }
            } else {
                let bks = self.bucket_kmer_start();
                let tab = self.kmer_table();
                for b in 0..bks.len().saturating_sub(1) {
                    let base = b as u64 * self.suffix_space.value();
                    for e in &tab[bks[b] as usize..bks[b + 1] as usize] {
                        let code = base + e[0] as u64;
                        let i = div.quot(code);
                        bits[(i >> 6) as usize] |= 1 << (i & 63);
                    }
                }
            }
            Presence { bits, div }
        })
    }

    /// Whether lookups use the direct offset table.
    pub fn is_direct(&self) -> bool {
        self.direct
    }

    /// Posting range of `code` in the direct layout (`[a, b)`, empty if the
    /// k-mer is absent). Requires [`is_direct`](Self::is_direct).
    #[inline]
    pub fn direct_range(&self, code: u64) -> (usize, usize) {
        debug_assert!(self.direct);
        let c = code as usize;
        let base = self.map.as_ptr();
        // SAFETY: `code < 20^w` (callers pass k-mer codes of this index's
        // weight), so entries c and c + 1 lie in the validated offset table
        unsafe {
            match self.koff_w {
                2 => (self.kmer_offset(c) as usize, self.kmer_offset(c + 1) as usize),
                4 => {
                    let p = base.add(self.koff.0 + c * 4) as *const u32;
                    (*p as usize, *p.add(1) as usize)
                }
                _ => {
                    let p = base.add(self.koff.0 + c * 8) as *const u64;
                    (*p as usize, *p.add(1) as usize)
                }
            }
        }
    }

    /// Prefetch the direct-table entries of `code` (no effect on results).
    #[inline]
    pub fn prefetch_code(&self, code: u64) {
        #[cfg(target_arch = "x86_64")]
        {
            if self.direct {
                use std::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
                let c = code as usize;
                // SAFETY: prefetch has no architectural effect; in bounds
                unsafe {
                    _mm_prefetch::<{ _MM_HINT_T0 }>(self.map.as_ptr().add(self.koff.0 + c * self.koff_w) as *const i8);
                    if self.koff_w == 2 {
                        _mm_prefetch::<{ _MM_HINT_T0 }>(self.map.as_ptr().add(self.kgrp.0 + c / OFFSET_GROUP * self.kgrp_w) as *const i8);
                    }
                }
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            let _ = code;
        }
    }

    /// Bytes per stored posting.
    pub fn posting_width(&self) -> usize {
        self.post_w
    }

    /// Prefetch the lookup-table entries a [`prefix_ranges`](Self::prefix_ranges)
    /// call for `prefix` will read first (no effect on results).
    #[inline]
    pub fn prefetch_prefix_bucket(&self, prefix: u64) {
        #[cfg(target_arch = "x86_64")]
        {
            use std::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
            if self.direct {
                let c = (prefix * A) as usize;
                let base = self.map.as_ptr();
                // SAFETY: prefetch has no architectural effect; the 21 entries
                // of this prefix lie inside the offset table
                unsafe {
                    let p = base.add(self.koff.0 + c * self.koff_w);
                    _mm_prefetch::<{ _MM_HINT_T0 }>(p as *const i8);
                    _mm_prefetch::<{ _MM_HINT_T0 }>(p.add(20 * self.koff_w) as *const i8);
                }
                return;
            }
            let b = self.suffix_space.quot(prefix * A) as usize;
            let bks = self.bucket_kmer_start();
            let bps = self.bucket_post_start();
            if b < bks.len() && b < bps.len() {
                // SAFETY: prefetch has no architectural effect; addresses in bounds
                unsafe {
                    _mm_prefetch::<{ _MM_HINT_T0 }>(bks.as_ptr().add(b) as *const i8);
                    _mm_prefetch::<{ _MM_HINT_T0 }>(bps.as_ptr().add(b) as *const i8);
                }
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            let _ = prefix;
        }
    }

    /// Prefetch the first k-mer-table entries of `prefix`'s bucket (reads the
    /// bucket start, ideally already prefetched; no-op for the direct layout).
    #[inline]
    pub fn prefetch_prefix_entries(&self, prefix: u64) {
        #[cfg(target_arch = "x86_64")]
        {
            if self.direct {
                return;
            }
            let b = self.suffix_space.quot(prefix * A) as usize;
            let bks = self.bucket_kmer_start();
            if b + 1 < bks.len() {
                let ks = bks[b] as usize;
                let tab = self.kmer_table();
                if ks < tab.len() {
                    // SAFETY: prefetch has no architectural effect; address in bounds
                    unsafe {
                        use std::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
                        _mm_prefetch::<{ _MM_HINT_T0 }>(tab.as_ptr().add(ks) as *const i8);
                    }
                }
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            let _ = prefix;
        }
    }

    /// Posting ranges of the 20 k-mers `prefix·20 + r` (`prefix` = the first
    /// `w − 1` residues of a code, i.e. `code / 20`): for every present k-mer
    /// sets bit `r` of the returned mask and `ranges[r] = [start, end)` in the
    /// global posting array. One lookup serves all 20 codes — the seed
    /// generator emits k-mers grouped by this prefix.
    #[inline]
    pub fn prefix_ranges(&self, prefix: u64, ranges: &mut [(usize, usize); 20]) -> u32 {
        let code0 = prefix * A;
        let mut mask = 0u32;
        if self.direct {
            let c = code0 as usize;
            let mut prev = self.kmer_offset(c);
            for (r, slot) in ranges.iter_mut().enumerate() {
                let next = self.kmer_offset(c + r + 1);
                if next > prev {
                    *slot = (prev as usize, next as usize);
                    mask |= 1 << r;
                }
                prev = next;
            }
            return mask;
        }
        let bks = self.bucket_kmer_start();
        let bps = self.bucket_post_start();
        if self.suffix_space.value() == 1 {
            // every k-mer is its own bucket (w ≤ 4)
            for (r, slot) in ranges.iter_mut().enumerate() {
                let b = (code0 + r as u64) as usize;
                if b + 1 < bks.len() && bks[b] != bks[b + 1] {
                    *slot = (bps[b] as usize, bps[b + 1] as usize);
                    mask |= 1 << r;
                }
            }
            return mask;
        }
        // suffix_space is a power of 20 ≥ 20: the 20 codes share a bucket and
        // occupy suffixes [s0, s0 + 20)
        let b = self.suffix_space.quot(code0) as usize;
        let s0 = self.suffix_space.modulo(code0) as u32;
        let (ks, ke) = (bks[b] as usize, bks[b + 1] as usize);
        if ks == ke {
            return 0;
        }
        let tab = &self.kmer_table()[ks..ke];
        let mut e = tab.partition_point(|x| x[0] < s0);
        let base = bps[b] as usize;
        let bend = bps[b + 1] as usize;
        while e < tab.len() && tab[e][0] < s0 + A as u32 {
            let r = (tab[e][0] - s0) as usize;
            let st = base + tab[e][1] as usize;
            let en = if e + 1 < tab.len() { base + tab[e + 1][1] as usize } else { bend };
            ranges[r] = (st, en);
            mask |= 1 << r;
            e += 1;
        }
        mask
    }

    /// Global posting-array range `[a, b)` of a k-mer code (empty if absent).
    #[inline]
    pub fn posting_range(&self, code: u64) -> (usize, usize) {
        if self.direct {
            let c = code as usize;
            let (a, b) = (self.kmer_offset(c), self.kmer_offset(c + 1));
            return if b > a { (a as usize, b as usize) } else { (0, 0) };
        }
        let b = self.suffix_space.quot(code) as usize;
        let suf = self.suffix_space.modulo(code) as u32;
        let bks = self.bucket_kmer_start();
        let (ks, ke) = (bks[b] as usize, bks[b + 1] as usize);
        if ks == ke {
            return (0, 0);
        }
        let tab = &self.kmer_table()[ks..ke];
        match tab.binary_search_by_key(&suf, |e| e[0]) {
            Ok(e) => {
                let bps = self.bucket_post_start();
                let base = bps[b] as usize;
                let s = base + tab[e][1] as usize;
                let end = if e + 1 < tab.len() {
                    base + tab[e + 1][1] as usize
                } else {
                    bps[b + 1] as usize
                };
                (s, end)
            }
            Err(_) => (0, 0),
        }
    }

    /// Posting `k` of the global posting array as `(target, position)`.
    #[inline]
    pub fn posting(&self, k: usize) -> (u64, u64) {
        debug_assert!((k as u64) < self.meta.num_postings);
        // SAFETY: k < num_postings (callers use ranges from the lookup tables),
        // so the 8-byte read starts inside the posting section; writers leave
        // 8 bytes of slack after it (checked in `open`)
        let v = unsafe {
            (self.map.as_ptr().add(self.postings.0 + k * self.post_w) as *const u64).read_unaligned()
        } & self.post_mask;
        (v >> self.pos_bits, v & self.pos_mask)
    }

    /// Postings `[a, b)` decoded as `(target << 32) | position` values.
    pub fn postings_in(&self, a: usize, b: usize) -> Vec<u64> {
        (a..b)
            .map(|k| {
                let (t, p) = self.posting(k);
                (t << 32) | p
            })
            .collect()
    }

    /// Posting list of a k-mer code as `(target << 32) | position` values,
    /// sorted by (target, position).
    pub fn postings(&self, code: u64) -> Vec<u64> {
        let (a, b) = self.posting_range(code);
        self.postings_in(a, b)
    }

    /// `(offset, length)` of target `i` in [`all_residues`](Self::all_residues).
    #[inline]
    pub fn seq_span(&self, i: usize) -> (usize, usize) {
        let o = self.seq_offsets();
        (o[i] as usize, (o[i + 1] - o[i]) as usize)
    }

    /// The concatenated residues of every target.
    #[inline]
    pub fn all_residues(&self) -> &[u8] {
        self.residues()
    }

    /// Hint the CPU to fetch `all_residues()[a..b]` (no effect on results).
    #[inline]
    pub fn prefetch_residue_range(&self, a: usize, b: usize) {
        #[cfg(target_arch = "x86_64")]
        {
            let res = self.residues();
            let b = b.min(res.len());
            let mut p = a & !63;
            while p < b {
                // SAFETY: prefetch has no architectural effect; p < res.len()
                unsafe {
                    std::arch::x86_64::_mm_prefetch::<{ std::arch::x86_64::_MM_HINT_T0 }>(res.as_ptr().add(p) as *const i8);
                }
                p += 64;
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            let _ = (a, b);
        }
    }

    /// Hint the CPU to fetch posting `k` (no effect on results).
    #[inline]
    pub fn prefetch_posting(&self, k: usize) {
        #[cfg(target_arch = "x86_64")]
        {
            if (k as u64) < self.meta.num_postings {
                // SAFETY: prefetch has no architectural effect; k is in bounds
                unsafe {
                    std::arch::x86_64::_mm_prefetch::<{ std::arch::x86_64::_MM_HINT_T0 }>(
                        self.map.as_ptr().add(self.postings.0 + k * self.post_w) as *const i8,
                    );
                }
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            let _ = k;
        }
    }

    /// Length of the longest target.
    #[inline]
    pub fn max_seq_len(&self) -> usize {
        self.meta.max_seq_len as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::build_aa2num;
    use std::collections::HashMap;
    use std::io::Cursor;

    fn db_from(fasta: &str) -> SeqDb {
        SeqDb::from_reader(Cursor::new(fasta), &build_aa2num()).unwrap()
    }

    fn tmp(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("sabertooth_index_test_{}", std::process::id()));
        fs::create_dir_all(&d).unwrap();
        d.join(name)
    }

    /// reference: naive hash-map index (the v0.2 structure)
    fn naive(db: &SeqDb, pat: &SeedPattern, mask: Option<MaskParams>) -> HashMap<u64, Vec<u64>> {
        let mut m: HashMap<u64, Vec<u64>> = HashMap::new();
        for i in 0..db.len() {
            let s = db.num(i);
            let mk = mask.map(|p| crate::mask::low_complexity_mask(s, p));
            if s.len() < pat.span() {
                continue;
            }
            for st in 0..=(s.len() - pat.span()) {
                if let Some(c) = encode_window(s, st, pat, mk.as_deref()) {
                    m.entry(c).or_default().push(((i as u64) << 32) | st as u64);
                }
            }
        }
        m
    }

    fn random_db(seed: u64, n: usize) -> SeqDb {
        let mut rng = crate::evalue::SplitMix64::new(seed);
        let aa = b"ACDEFGHIKLMNPQRSTVWYX";
        let mut s = String::new();
        for i in 0..n {
            let len = 1 + rng.below(300);
            s.push_str(&format!(">s{} d\n", i));
            for _ in 0..len {
                // skewed composition + occasional X and repeats
                let r = if rng.below(10) == 0 {
                    20
                } else {
                    let alpha = 8 + rng.below(12);
                    rng.below(alpha)
                };
                s.push(aa[r] as char);
            }
            if i % 7 == 0 {
                s.push_str(&"Q".repeat(25));
            }
            s.push('\n');
        }
        db_from(&s)
    }

    #[test]
    fn index_matches_naive_hashmap_for_patterns_masks_and_partitions() {
        let db = random_db(3, 400);
        for mask_s in ["111111", "11011", "1101011", "11", "1111111"] {
            let pat = SeedPattern::from_mask(mask_s).unwrap();
            for mask in [None, Some(MaskParams::default())] {
                for budget in [1usize << 30, 1 << 12] {
                    let mut o = BuildOptions::new(pat.clone());
                    o.mask = mask;
                    o.memory_budget = budget;
                    let (idx, st) = SeedIndex::build_transient(&db, &o).unwrap();
                    if budget < 1 << 20 {
                        assert!(st.partitions >= 1);
                    }
                    let want = naive(&db, &pat, mask);
                    let total: usize = want.values().map(|v| v.len()).sum();
                    assert_eq!(idx.meta().num_postings as usize, total);
                    assert_eq!(idx.distinct_kmers(), want.len(), "pattern {}", mask_s);
                    for (code, list) in &want {
                        assert_eq!(idx.postings(*code), &list[..], "code {} pattern {}", code, mask_s);
                    }
                    // absent codes
                    let mut rng = crate::evalue::SplitMix64::new(9);
                    let space = 20u64.pow(pat.weight() as u32);
                    for _ in 0..2000 {
                        let c = rng.next_u64() % space;
                        assert_eq!(idx.postings(c).len(), want.get(&c).map(|v| v.len()).unwrap_or(0));
                    }
                    // sequences and ids round-trip
                    for i in (0..db.len()).step_by(37) {
                        assert_eq!(idx.seq(i), db.num(i));
                        assert_eq!(idx.id(i), db.id(i));
                    }
                }
            }
        }
    }

    #[test]
    fn layouts_and_posting_widths_agree() {
        // every lookup path (sparse / direct) and posting encoding (minimal
        // width, forced wider widths, the v1-style 8 bytes) returns the same
        // postings, prefix ranges and presence answers
        let db = random_db(5, 300);
        for mask_s in ["111111", "11011", "1111", "11"] {
            let pat = SeedPattern::from_mask(mask_s).unwrap();
            let want = naive(&db, &pat, None);
            let w = pat.weight() as u32;
            for (layout, grouped) in [
                (LayoutChoice::Sparse, true),
                (LayoutChoice::Direct, true),
                (LayoutChoice::Direct, false),
                (LayoutChoice::Auto, true),
            ] {
                for width in [None, Some(4), Some(8)] {
                    let mut o = BuildOptions::new(pat.clone());
                    o.layout = layout;
                    o.direct_grouped = grouped;
                    o.posting_width = width;
                    let (idx, st) = SeedIndex::build_transient(&db, &o).unwrap();
                    let tag = format!("{} {:?} grouped {} {:?}", mask_s, layout, grouped, width);
                    match layout {
                        LayoutChoice::Sparse => assert!(!idx.is_direct(), "{}", tag),
                        LayoutChoice::Direct => {
                            assert!(idx.is_direct(), "{}", tag);
                            assert_eq!(idx.meta().offset_width, if grouped { 2 } else { 4 }, "{}", tag);
                        }
                        LayoutChoice::Auto => {}
                    }
                    if let Some(wd) = width {
                        assert_eq!(idx.posting_width(), wd as usize, "{}", tag);
                    }
                    assert_eq!(st.posting_width as usize, idx.posting_width());
                    idx.verify().unwrap();
                    assert_eq!(idx.distinct_kmers(), want.len(), "{}", tag);
                    for (code, list) in &want {
                        assert_eq!(idx.postings(*code), &list[..], "{} code {}", tag, code);
                    }
                    let mut rng = crate::evalue::SplitMix64::new(77);
                    let space = 20u64.pow(w);
                    let mut ranges = [(0usize, 0usize); 20];
                    let presence = idx.presence();
                    for _ in 0..3000 {
                        let c = rng.next_u64() % space;
                        let present = want.contains_key(&c);
                        assert_eq!(idx.postings(c).len(), want.get(&c).map(|v| v.len()).unwrap_or(0), "{}", tag);
                        assert_eq!(presence.may_contain(c), present, "{} presence {}", tag, c);
                        let mask = idx.prefix_ranges(c / 20, &mut ranges);
                        assert_eq!((mask >> (c % 20)) & 1 == 1, present, "{} prefix {}", tag, c);
                        if present {
                            assert_eq!(ranges[(c % 20) as usize], idx.posting_range(c), "{}", tag);
                        }
                    }
                    for i in (0..db.len()).step_by(29) {
                        assert_eq!(idx.seq(i), db.num(i));
                        assert_eq!(idx.id(i), db.id(i));
                    }
                }
            }
        }
        // a group whose postings overflow 16 bits falls back to full offsets
        let t = crate::alphabet::build_aa2num();
        let big = SeqDb::from_reader(io::Cursor::new(format!(">poly\n{}\n>mix\nMKVLLACDEFGHIKLMNPQRSTVWY\n", "A".repeat(70_000))), &t).unwrap();
        let mut o = BuildOptions::new(SeedPattern::contiguous(2));
        o.mask = None;
        o.layout = LayoutChoice::Direct;
        let (idx, _) = SeedIndex::build_transient(&big, &o).unwrap();
        assert_eq!(idx.meta().offset_width, 4, "overflowing group must not be grouped");
        idx.verify().unwrap();
        let want = naive(&big, &SeedPattern::contiguous(2), None);
        for (code, list) in &want {
            assert_eq!(idx.postings(*code), &list[..], "code {}", code);
        }
        // impossible requests fail loudly
        let mut o = BuildOptions::new(SeedPattern::from_mask("1111111").unwrap());
        o.layout = LayoutChoice::Direct;
        assert!(SeedIndex::build_transient(&db, &o).is_err(), "direct layout for weight 7");
        let mut o = BuildOptions::new(SeedPattern::contiguous(6));
        o.posting_width = Some(1);
        assert!(SeedIndex::build_transient(&db, &o).is_err(), "1-byte postings cannot hold 300 targets");
    }

    #[test]
    fn format_v1_indexes_remain_readable() {
        // tests/data/v1_fixture_k2.stidx was written by sabertooth 0.3.0
        // (format 1: sparse tables, 8-byte postings) from v1_fixture.fa
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        let v1 = SeedIndex::open(&dir.join("v1_fixture_k2.stidx")).unwrap();
        assert_eq!(v1.meta().version, 1);
        assert_eq!(v1.posting_width(), 8);
        assert!(!v1.is_direct());
        v1.verify().unwrap();
        let text = std::fs::read_to_string(dir.join("v1_fixture.fa")).unwrap();
        let db = db_from(&text);
        let pat = SeedPattern::contiguous(2);
        let (v2, _) = SeedIndex::build_transient(&db, &BuildOptions::new(pat.clone())).unwrap();
        assert_eq!(v2.meta().version, FORMAT_VERSION);
        assert_eq!(v1.len(), v2.len());
        assert_eq!(v1.meta().num_postings, v2.meta().num_postings);
        let want = naive(&db, &pat, None);
        for c in 0..400u64 {
            assert_eq!(v1.postings(c), v2.postings(c), "code {}", c);
            assert_eq!(v1.postings(c).len(), want.get(&c).map(|v| v.len()).unwrap_or(0));
        }
        for i in 0..db.len() {
            assert_eq!(v1.seq(i), v2.seq(i));
            assert_eq!(v1.id(i), v2.id(i));
        }
    }

    #[test]
    fn prefix_ranges_agree_with_single_lookups() {
        let db = random_db(13, 300);
        for mask_s in ["11", "111", "1111", "11111", "111111", "1101011", "11111111"] {
            let pat = SeedPattern::from_mask(mask_s).unwrap();
            let (idx, _) = SeedIndex::build_transient(&db, &BuildOptions::new(pat.clone())).unwrap();
            let w = pat.weight() as u32;
            let want = naive(&db, &pat, None);
            // every present prefix plus random (mostly absent) ones
            let mut prefixes: Vec<u64> = want.keys().map(|c| c / 20).collect();
            let mut rng = crate::evalue::SplitMix64::new(21);
            for _ in 0..500 {
                prefixes.push(rng.next_u64() % 20u64.pow(w - 1));
            }
            prefixes.sort_unstable();
            prefixes.dedup();
            let mut ranges = [(0usize, 0usize); 20];
            let presence = idx.presence();
            for &pre in &prefixes {
                let mask = idx.prefix_ranges(pre, &mut ranges);
                for r in 0..20u64 {
                    let code = pre * 20 + r;
                    let single = idx.posting_range(code);
                    let present = (mask >> r) & 1 == 1;
                    assert_eq!(present, single.1 > single.0, "pattern {} code {}", mask_s, code);
                    // the bitmap never rejects a present code, and is exact up to weight 6
                    if present {
                        assert!(presence.may_contain(code));
                    } else if w <= PRESENCE_MAX_RESIDUES {
                        assert!(!presence.may_contain(code), "pattern {} code {}", mask_s, code);
                    }
                    if present {
                        assert_eq!(ranges[r as usize], single, "pattern {} code {}", mask_s, code);
                        assert_eq!(idx.postings_in(single.0, single.1), &want[&code][..]);
                    } else {
                        assert!(!want.contains_key(&code));
                    }
                }
            }
        }
    }

    #[test]
    fn small_budget_forces_multiple_partitions_with_identical_output() {
        let db = random_db(5, 1500);
        let pat = SeedPattern::contiguous(5);
        let mut big = BuildOptions::new(pat.clone());
        big.memory_budget = 1 << 30;
        let mut small = BuildOptions::new(pat);
        small.memory_budget = 1; // floor of 65,536 entries per partition
        let a = tmp("big.stidx");
        let b = tmp("small.stidx");
        let sa = build_index(&mut DbSource(&db), SourceFingerprint::transient(), &a, &big).unwrap();
        let sb = build_index(&mut DbSource(&db), SourceFingerprint::transient(), &b, &small).unwrap();
        assert_eq!(sa.partitions, 1);
        assert!(sb.partitions > 1, "expected several partitions, got {}", sb.partitions);
        let (ia, ib) = (SeedIndex::open(&a).unwrap(), SeedIndex::open(&b).unwrap());
        for k in [SectionKind::Postings, SectionKind::KmerTable, SectionKind::Residues] {
            assert_eq!(ia.meta().section(k).hash, ib.meta().section(k).hash, "{:?}", k);
        }
        ia.verify().unwrap();
        ib.verify().unwrap();
    }

    #[test]
    fn corrupt_truncated_and_stale_indexes_are_detected() {
        let db = random_db(7, 50);
        let pat = SeedPattern::contiguous(6);
        let fa = tmp("src.fasta");
        db.write_fasta(&fa).unwrap();
        let src = fa.to_string_lossy().to_string();
        let out = tmp("src.stidx");
        let mut fsrc = FastaSource {
            path: src.clone(),
            aa2num: build_aa2num(),
        };
        build_index(&mut fsrc, SourceFingerprint::of(&src).unwrap(), &out, &BuildOptions::new(pat.clone())).unwrap();
        let idx = SeedIndex::open(&out).unwrap();
        assert!(idx.usable_for(&pat, None, Some(&src)).is_ok());
        // incompatible requests
        assert!(matches!(idx.usable_for(&SeedPattern::contiguous(5), None, None), Err(Unusable::Incompatible(_))));
        assert!(matches!(idx.usable_for(&pat, Some(MaskParams::default()), None), Err(Unusable::Incompatible(_))));
        // stale: touch the source
        std::thread::sleep(std::time::Duration::from_millis(20));
        let mut text = fs::read_to_string(&fa).unwrap();
        text.push_str(">extra\nMKV\n");
        fs::write(&fa, text).unwrap();
        assert!(matches!(idx.usable_for(&pat, None, Some(&src)), Err(Unusable::Stale(_))));
        drop(idx);
        // truncated (footer missing)
        let bytes = fs::read(&out).unwrap();
        let tr = tmp("trunc.stidx");
        fs::write(&tr, &bytes[..bytes.len() - 9]).unwrap();
        let e = SeedIndex::open(&tr).unwrap_err();
        assert!(e.to_string().contains("footer") || e.to_string().contains("incomplete"), "{}", e);
        // header corruption
        let mut bad = bytes.clone();
        bad[100] ^= 0xFF;
        let bp = tmp("bad.stidx");
        fs::write(&bp, &bad).unwrap();
        assert!(SeedIndex::open(&bp).unwrap_err().to_string().contains("checksum"));
        // payload corruption is caught by verify()
        let mut bad2 = bytes.clone();
        let post = decode_header(&bytes[..4096], "x").unwrap().section(SectionKind::Postings);
        bad2[post.offset as usize + 3] ^= 0x10;
        let bp2 = tmp("bad2.stidx");
        fs::write(&bp2, &bad2).unwrap();
        let i2 = SeedIndex::open(&bp2).unwrap();
        assert!(i2.verify().is_err());
        // not an index at all
        let np = tmp("not.stidx");
        fs::write(&np, vec![7u8; 8192]).unwrap();
        assert!(SeedIndex::open(&np).unwrap_err().to_string().contains("magic"));
    }

    #[test]
    fn no_partial_file_after_failed_build() {
        struct Failing;
        impl RecordSource for Failing {
            fn for_each(&mut self, f: &mut dyn FnMut(&[u8], &[u8]) -> io::Result<()>) -> io::Result<()> {
                f(b"a", &[0, 1, 2, 3, 4, 5, 6, 7])?;
                Err(io::Error::other("simulated read failure"))
            }
        }
        let out = tmp("fail.stidx");
        let _ = fs::remove_file(&out);
        let e = build_index(&mut Failing, SourceFingerprint::transient(), &out, &BuildOptions::new(SeedPattern::contiguous(4)));
        assert!(e.is_err());
        assert!(!out.exists(), "final path must not exist after a failed build");
        let dir = out.parent().unwrap();
        for ent in fs::read_dir(dir).unwrap() {
            let n = ent.unwrap().file_name().to_string_lossy().to_string();
            assert!(!n.starts_with("fail.stidx.tmp"), "leftover {}", n);
        }
    }

    #[test]
    fn prefix_lengths_cover_all_weights() {
        for w in 2..=12usize {
            for total in [0u64, 1 << 20, 1 << 26, 1 << 36] {
                let p = prefix_len_for_size(w, total, LARGE_PREFIX_MIN_RESIDUES_PER_BUCKET) as usize;
                assert!(p <= w);
                assert!(p <= 5 || w > 9, "w={} p={}", w, p);
                assert!(20u64.pow((w - p) as u32) <= u32::MAX as u64, "w={}", w);
            }
        }
        assert_eq!(prefix_len_for_size(6, 1_000_000, LARGE_PREFIX_MIN_RESIDUES_PER_BUCKET), 4);
        assert_eq!(prefix_len_for_size(6, 68_000_000, LARGE_PREFIX_MIN_RESIDUES_PER_BUCKET), 5);
        assert_eq!(prefix_len_for_size(6, 10, 0), 5);
    }

    #[test]
    fn large_prefix_layout_gives_identical_lookups() {
        // force the (w−1)-prefix layout with a zero threshold
        let db = random_db(17, 600);
        for w in [5usize, 6, 7] {
            let pat = SeedPattern::contiguous(w);
            let mut o = BuildOptions::new(pat.clone());
            o.large_prefix_threshold = 0;
            let (idx, _) = SeedIndex::build_transient(&db, &o).unwrap();
            assert_eq!(idx.meta().prefix_len as usize, (w - 1).min(5));
            let want = naive(&db, &pat, None);
            assert_eq!(idx.distinct_kmers(), want.len());
            for (code, list) in &want {
                assert_eq!(idx.postings(*code), &list[..]);
            }
            let mut ranges = [(0usize, 0usize); 20];
            // every present prefix plus a deterministic sample of others
            let mut prefixes: Vec<u64> = want.keys().map(|c| c / 20).collect();
            let space = 20u64.pow(w as u32 - 1);
            prefixes.extend((0..4000u64).map(|i| (i * 2_654_435_761) % space));
            prefixes.sort_unstable();
            prefixes.dedup();
            for &pre in &prefixes {
                let mask = idx.prefix_ranges(pre, &mut ranges);
                for r in 0..20u64 {
                    let c = pre * 20 + r;
                    assert_eq!((mask >> r) & 1 == 1, want.contains_key(&c));
                }
            }
        }
    }
}
