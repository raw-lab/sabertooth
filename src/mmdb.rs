//! MMseqs2 on-disk database format (read + write).
//!
//! `mmseqs createdb` turns a FASTA into a small family of files. Verified byte
//! layout (from the real binary):
//!
//! * `<db>`         — data: each record is `payload + '\n' + '\0'`.
//! * `<db>.index`   — one `key<TAB>offset<TAB>length` line per record, where
//!                    `length` **includes** the trailing `'\n'` and `'\0'`.
//! * `<db>.dbtype`  — 4-byte little-endian int: 0 = amino acid, 1 = nucleotide,
//!                    2 = profile, 12 = generic/header.
//! * `<db>_h`       — header data (FASTA header without `>`), same `\n\0` framing.
//! * `<db>_h.index` — index for the header data.
//! * `<db>.lookup`  — `key<TAB>accession<TAB>set` (accession = first header token).
//! * `<db>.source`  — `0<TAB><source filename>`.
//!
//! We write all of these so a Sabertooth DB is readable by `mmseqs`, and we read
//! the same layout so an `mmseqs createdb` output is readable by Sabertooth.

use crate::fasta::{Record, SeqDb};
use crate::matrix::SubstitutionMatrix;
use crate::profile::Profile;
use std::fs;
use std::io::{self, Write};
use std::path::Path;

pub const DBTYPE_AMINO_ACIDS: i32 = 0;
pub const DBTYPE_NUCLEOTIDES: i32 = 1;
/// Profile database (`mmseqs msa2profile` output). Entries are 25-byte columns.
pub const DBTYPE_PROFILE: i32 = 2;
pub const DBTYPE_GENERIC: i32 = 12; // header databases use this

/// One entry of a `.index` file.
struct IndexEntry {
    key: u32,
    offset: u64,
    length: u64,
}

fn write_dbtype(path: &Path, dbtype: i32) -> io::Result<()> {
    fs::write(path, dbtype.to_le_bytes())
}

/// Read a 4-byte little-endian dbtype; missing file defaults to amino acids.
fn read_dbtype(path: &Path) -> i32 {
    match fs::read(path) {
        Ok(b) if b.len() >= 4 => i32::from_le_bytes([b[0], b[1], b[2], b[3]]),
        _ => DBTYPE_AMINO_ACIDS,
    }
}

fn parse_index(text: &str) -> Vec<IndexEntry> {
    let mut out = Vec::new();
    for line in text.lines() {
        let mut it = line.split('\t');
        let key = it.next().and_then(|s| s.trim().parse().ok());
        let offset = it.next().and_then(|s| s.trim().parse().ok());
        let length = it.next().and_then(|s| s.trim().parse().ok());
        if let (Some(key), Some(offset), Some(length)) = (key, offset, length) {
            out.push(IndexEntry { key, offset, length });
        }
    }
    out
}

/// Append `payload + '\n' + '\0'` to `data`, recording the index entry.
fn push_entry(data: &mut Vec<u8>, index: &mut Vec<IndexEntry>, key: u32, payload: &[u8]) {
    let offset = data.len() as u64;
    data.extend_from_slice(payload);
    data.push(b'\n');
    data.push(b'\0');
    let length = payload.len() as u64 + 2;
    index.push(IndexEntry { key, offset, length });
}

fn write_index(path: &Path, index: &[IndexEntry]) -> io::Result<()> {
    let mut f = io::BufWriter::new(fs::File::create(path)?);
    for e in index {
        writeln!(f, "{}\t{}\t{}", e.key, e.offset, e.length)?;
    }
    f.flush()
}

/// Write a sequence database to `<prefix>` in MMseqs2 format. `nucleotide`
/// selects the data dbtype. `source` is recorded in `<prefix>.source`.
pub fn write_db(prefix: &Path, db: &SeqDb, nucleotide: bool, source: &str) -> io::Result<()> {
    let mut data = Vec::new();
    let mut dindex = Vec::new();
    let mut hdata = Vec::new();
    let mut hindex = Vec::new();
    let mut lookup = String::new();

    for (i, rec) in db.records.iter().enumerate() {
        let key = i as u32;
        push_entry(&mut data, &mut dindex, key, &rec.seq);
        push_entry(&mut hdata, &mut hindex, key, rec.header.as_bytes());
        lookup.push_str(&format!("{}\t{}\t0\n", key, rec.id));
    }

    let pstr = prefix.to_string_lossy().to_string();
    fs::write(&pstr, &data)?;
    write_index(Path::new(&format!("{}.index", pstr)), &dindex)?;
    write_dbtype(
        Path::new(&format!("{}.dbtype", pstr)),
        if nucleotide {
            DBTYPE_NUCLEOTIDES
        } else {
            DBTYPE_AMINO_ACIDS
        },
    )?;
    fs::write(format!("{}_h", pstr), &hdata)?;
    write_index(Path::new(&format!("{}_h.index", pstr)), &hindex)?;
    write_dbtype(Path::new(&format!("{}_h.dbtype", pstr)), DBTYPE_GENERIC)?;
    fs::write(format!("{}.lookup", pstr), lookup)?;
    fs::write(format!("{}.source", pstr), "0\t".to_string() + source + "\n")?;
    Ok(())
}

/// Strip the trailing `'\n''\0'` (or a bare `'\0'`) framing from a data slice.
fn unframe(entry: &[u8]) -> &[u8] {
    let mut end = entry.len();
    if end > 0 && entry[end - 1] == 0 {
        end -= 1;
    }
    if end > 0 && entry[end - 1] == b'\n' {
        end -= 1;
    }
    &entry[..end]
}

/// Read an MMseqs2 database at `<prefix>` back into a `SeqDb`. Recognises the
/// data dbtype (amino acid vs nucleotide) to map residues; unknown letters
/// become `X`. Header text is taken from `<prefix>_h` when present, otherwise
/// the numeric key is used as the id.
pub fn read_db(prefix: &Path, aa2num: &[u8; 256]) -> io::Result<SeqDb> {
    let pstr = prefix.to_string_lossy().to_string();
    let data = fs::read(&pstr)?;
    let dindex = parse_index(&fs::read_to_string(format!("{}.index", pstr))?);
    let dbtype = read_dbtype(Path::new(&format!("{}.dbtype", pstr)));
    if dbtype == DBTYPE_PROFILE {
        // Profile payloads are binary and legitimately contain 0x0A, so the
        // text `\n\0` unframing below would silently corrupt them.
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "this is a profile DB (dbtype 2) — use read_profile_db, not read_db",
        ));
    }

    // headers are optional
    let hdata = fs::read(format!("{}_h", pstr)).ok();
    let hindex = fs::read_to_string(format!("{}_h.index", pstr))
        .ok()
        .map(|s| parse_index(&s));

    let mut records = Vec::with_capacity(dindex.len());
    let mut total = 0usize;
    for (i, e) in dindex.iter().enumerate() {
        let start = e.offset as usize;
        let end = start + e.length as usize;
        if end > data.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("index entry {} runs past end of data file", e.key),
            ));
        }
        let seq = unframe(&data[start..end]).to_vec();

        // matching header, if available
        let header = match (&hdata, &hindex) {
            (Some(hd), Some(hi)) if i < hi.len() => {
                let hs = hi[i].offset as usize;
                let he = hs + hi[i].length as usize;
                if he <= hd.len() {
                    String::from_utf8_lossy(unframe(&hd[hs..he])).to_string()
                } else {
                    e.key.to_string()
                }
            }
            _ => e.key.to_string(),
        };
        let id = header
            .split_whitespace()
            .next()
            .unwrap_or(&header)
            .to_string();
        let num: Vec<u8> = seq.iter().map(|&b| aa2num[b as usize]).collect();
        total += num.len();
        records.push(Record {
            header,
            id,
            seq,
            num,
        });
    }
    Ok(SeqDb {
        records,
        total_residues: total,
    })
}

/// Read an MMseqs2 **profile** database (`dbtype == 2`, e.g. the output of
/// `mmseqs msa2profile`) into `(name, Profile)` pairs, in index order.
///
/// Profile entries are 25-byte columns terminated by a single `NUL`. Unlike text
/// entries they are **not** `\n\0`-framed: the binary payload contains `0x0A`
/// bytes as ordinary score/residue values, so only the trailing `NUL` is
/// stripped. Decoding of each column is done by [`Profile::from_pssm_columns`].
pub fn read_profile_db(
    prefix: &Path,
    mat: &SubstitutionMatrix,
) -> io::Result<Vec<(String, Profile)>> {
    let pstr = prefix.to_string_lossy().to_string();
    let dbtype = read_dbtype(Path::new(&format!("{}.dbtype", pstr)));
    if dbtype != DBTYPE_PROFILE {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("expected a profile DB (dbtype {}), found dbtype {}", DBTYPE_PROFILE, dbtype),
        ));
    }
    let data = fs::read(&pstr)?;
    let dindex = parse_index(&fs::read_to_string(format!("{}.index", pstr))?);
    let hdata = fs::read(format!("{}_h", pstr)).ok();
    let hindex = fs::read_to_string(format!("{}_h.index", pstr))
        .ok()
        .map(|s| parse_index(&s));

    let mut out = Vec::with_capacity(dindex.len());
    for (i, e) in dindex.iter().enumerate() {
        let start = e.offset as usize;
        let end = start + e.length as usize;
        if end > data.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("profile entry {} runs past end of data file", e.key),
            ));
        }
        // strip exactly one trailing NUL if present; never strip '\n'
        let mut payload = &data[start..end];
        if payload.last() == Some(&0) {
            payload = &payload[..payload.len() - 1];
        }
        let profile = Profile::from_pssm_columns(payload, mat).map_err(|msg| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("profile entry {}: {}", e.key, msg),
            )
        })?;
        let name = match (&hdata, &hindex) {
            (Some(hd), Some(hi)) if i < hi.len() => {
                let hs = hi[i].offset as usize;
                let he = hs + hi[i].length as usize;
                if he <= hd.len() {
                    String::from_utf8_lossy(unframe(&hd[hs..he]))
                        .split_whitespace()
                        .next()
                        .unwrap_or("")
                        .to_string()
                } else {
                    e.key.to_string()
                }
            }
            _ => e.key.to_string(),
        };
        out.push((name, profile));
    }
    Ok(out)
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

    #[test]
    fn write_then_read_roundtrips() {
        let t = build_aa2num();
        let fasta = ">seq1 first protein\nMKVLLACDEF\n>seq2 second\nWYFPCDGHIK\n";
        let db = SeqDb::from_reader(Cursor::new(fasta), &t).unwrap();
        let prefix = tmp("rt");
        write_db(&prefix, &db, false, "in.fasta").unwrap();

        let back = read_db(&prefix, &t).unwrap();
        assert_eq!(back.records.len(), 2);
        assert_eq!(back.records[0].id, "seq1");
        assert_eq!(back.records[0].header, "seq1 first protein");
        assert_eq!(back.records[0].seq, b"MKVLLACDEF");
        assert_eq!(back.records[1].id, "seq2");
        assert_eq!(back.records[1].seq, b"WYFPCDGHIK");
        assert_eq!(back.records[0].num, db.records[0].num);
    }

    #[test]
    fn index_lengths_include_framing() {
        // an entry's index length must be payload + 2 (\n and \0)
        let t = build_aa2num();
        let db = SeqDb::from_reader(Cursor::new(">a\nMKVLLA\n"), &t).unwrap();
        let prefix = tmp("len");
        write_db(&prefix, &db, false, "x").unwrap();
        let idx = fs::read_to_string(format!("{}.index", prefix.to_string_lossy())).unwrap();
        // "MKVLLA" is 6 -> length 8
        assert_eq!(idx.lines().next().unwrap(), "0\t0\t8");
        let data = fs::read(prefix.to_string_lossy().to_string()).unwrap();
        assert_eq!(&data, b"MKVLLA\n\0");
    }

    /// Write a minimal profile DB (dbtype 2) with bare-NUL entry framing.
    fn write_profile_db_fixture(prefix: &std::path::Path, entries: &[Vec<u8>]) {
        let pstr = prefix.to_string_lossy().to_string();
        let mut data = Vec::new();
        let mut index = String::new();
        for (key, payload) in entries.iter().enumerate() {
            let off = data.len();
            data.extend_from_slice(payload);
            data.push(0); // single NUL, no '\n'
            index.push_str(&format!("{}\t{}\t{}\n", key, off, payload.len() + 1));
        }
        fs::write(&pstr, &data).unwrap();
        fs::write(format!("{}.index", pstr), index).unwrap();
        fs::write(format!("{}.dbtype", pstr), DBTYPE_PROFILE.to_le_bytes()).unwrap();
    }

    fn fixture_column(peak: usize, q: u8) -> Vec<u8> {
        let mut c = vec![0u8; crate::profile::PROFILE_READIN_SIZE];
        c[peak] = 25u8; // a positive i8 score
        c[crate::profile::PROFILE_QUERY_OFFSET] = q;
        c[crate::profile::PROFILE_CONSENSUS_OFFSET] = q;
        c[crate::profile::PROFILE_NEFF_OFFSET] = 65; // Neff 2.0
        c
    }

    #[test]
    fn read_profile_db_decodes_entries() {
        let m = SubstitutionMatrix::blosum62();
        let prefix = tmp("profdb");
        let mut e0 = fixture_column(5, 5);
        e0.extend(fixture_column(0, 0)); // 2 columns
        let e1 = fixture_column(19, 19); // 1 column
        write_profile_db_fixture(&prefix, &[e0, e1]);

        let profiles = read_profile_db(&prefix, &m).unwrap();
        assert_eq!(profiles.len(), 2);
        assert_eq!(profiles[0].1.query_len, 2);
        assert_eq!(profiles[1].1.query_len, 1);
        assert_eq!(profiles[0].1.pssm[5], 25);
        assert_eq!(profiles[1].1.consensus[0], 19);
        assert!((profiles[0].1.neff[0] - 2.0).abs() < 1e-5);
        assert_eq!(profiles[0].1.set_size, 0); // documented: not recoverable
    }

    #[test]
    fn read_db_refuses_a_profile_db() {
        let m = SubstitutionMatrix::blosum62();
        let prefix = tmp("guard_prof");
        write_profile_db_fixture(&prefix, &[fixture_column(3, 3)]);
        let err = read_db(&prefix, &m.aa2num).unwrap_err();
        assert!(
            err.to_string().contains("profile DB"),
            "read_db must refuse a dbtype-2 DB, got: {}",
            err
        );
    }

    #[test]
    fn read_profile_db_refuses_a_sequence_db() {
        let t = build_aa2num();
        let m = SubstitutionMatrix::blosum62();
        let db = SeqDb::from_reader(Cursor::new(">a\nMKVLLA\n"), &t).unwrap();
        let prefix = tmp("guard_seq");
        write_db(&prefix, &db, false, "x").unwrap();
        let err = read_profile_db(&prefix, &m).unwrap_err();
        assert!(err.to_string().contains("expected a profile DB"), "got: {}", err);
    }

    #[test]
    fn dbtype_bytes_are_correct() {
        let t = build_aa2num();
        let db = SeqDb::from_reader(Cursor::new(">a\nMKVL\n"), &t).unwrap();
        let prefix = tmp("dbtype");
        write_db(&prefix, &db, false, "x").unwrap();
        let bytes = fs::read(format!("{}.dbtype", prefix.to_string_lossy())).unwrap();
        assert_eq!(bytes, vec![0, 0, 0, 0]); // amino acid
        let hbytes = fs::read(format!("{}_h.dbtype", prefix.to_string_lossy())).unwrap();
        assert_eq!(hbytes, vec![12, 0, 0, 0]); // header/generic
    }
}
