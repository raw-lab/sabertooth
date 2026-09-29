//! FASTA parsing and a compact in-memory sequence database.
//!
//! Sequences are stored twice: the raw uppercase residue bytes (for output) and
//! a numeric encoding (internal residue indices) used by the aligner and the
//! k-mer prefilter. This mirrors the split MMseqs2 keeps between its sequence
//! DB and the numeric `Sequence` objects fed to the compute kernels.

use crate::alphabet::aa_to_num;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

/// A single sequence record.
#[derive(Clone, Debug)]
pub struct Record {
    /// FASTA header without the leading `>` (identifier + description).
    pub header: String,
    /// The identifier only (first whitespace-delimited token of the header).
    pub id: String,
    /// Uppercased residue bytes as read from disk.
    pub seq: Vec<u8>,
    /// Numeric residue encoding (internal indices, `X` for anything unknown).
    pub num: Vec<u8>,
}

impl Record {
    pub fn len(&self) -> usize {
        self.num.len()
    }
    pub fn is_empty(&self) -> bool {
        self.num.is_empty()
    }
}

/// A collection of sequence records plus aggregate statistics.
#[derive(Clone, Debug)]
pub struct SeqDb {
    pub records: Vec<Record>,
    /// Total number of residues across all records (search-space size).
    pub total_residues: usize,
}

impl SeqDb {
    pub fn len(&self) -> usize {
        self.records.len()
    }
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Read a FASTA file into a sequence database, encoding residues with the
    /// supplied ASCII->index table.
    pub fn read_fasta<P: AsRef<Path>>(path: P, aa2num: &[u8; 256]) -> io::Result<SeqDb> {
        let file = File::open(&path)?;
        Self::from_reader(BufReader::new(file), aa2num)
    }

    /// Parse FASTA from any buffered reader.
    pub fn from_reader<R: BufRead>(reader: R, aa2num: &[u8; 256]) -> io::Result<SeqDb> {
        let mut records: Vec<Record> = Vec::new();
        let mut header: Option<String> = None;
        let mut seq: Vec<u8> = Vec::new();

        let flush = |header: &mut Option<String>, seq: &mut Vec<u8>, records: &mut Vec<Record>| {
            if let Some(h) = header.take() {
                let id = h.split_whitespace().next().unwrap_or("").to_string();
                let num: Vec<u8> = seq.iter().map(|&b| aa_to_num(b, aa2num)).collect();
                records.push(Record {
                    header: h,
                    id,
                    seq: std::mem::take(seq),
                    num,
                });
            } else {
                seq.clear();
            }
        };

        for line in reader.lines() {
            let line = line?;
            let line = line.trim_end();
            if line.is_empty() {
                continue;
            }
            if let Some(stripped) = line.strip_prefix('>') {
                flush(&mut header, &mut seq, &mut records);
                header = Some(stripped.to_string());
            } else {
                for &b in line.as_bytes() {
                    if b.is_ascii_alphabetic() {
                        seq.push(b.to_ascii_uppercase());
                    }
                }
            }
        }
        flush(&mut header, &mut seq, &mut records);

        let total_residues = records.iter().map(|r| r.num.len()).sum();
        Ok(SeqDb {
            records,
            total_residues,
        })
    }

    /// Write the database back out as FASTA (used by `createdb` round-tripping).
    pub fn write_fasta<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let mut out = File::create(path)?;
        for rec in &self.records {
            writeln!(out, ">{}", rec.header)?;
            for chunk in rec.seq.chunks(60) {
                out.write_all(chunk)?;
                out.write_all(b"\n")?;
            }
        }
        Ok(())
    }
}

/// Read an entire file to a `String` (small helper for MSA / matrix inputs).
pub fn read_to_string<P: AsRef<Path>>(path: P) -> io::Result<String> {
    let mut s = String::new();
    File::open(path)?.read_to_string(&mut s)?;
    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::build_aa2num;
    use std::io::Cursor;

    #[test]
    fn parses_multi_record_fasta() {
        let t = build_aa2num();
        let data = ">sp|P1 test one\nMKV\nLLA\n>sp|P2 test two\nGGGG\n";
        let db = SeqDb::from_reader(Cursor::new(data), &t).unwrap();
        assert_eq!(db.len(), 2);
        assert_eq!(db.records[0].id, "sp|P1");
        assert_eq!(db.records[0].seq, b"MKVLLA");
        assert_eq!(db.records[0].num.len(), 6);
        assert_eq!(db.records[1].len(), 4);
        assert_eq!(db.total_residues, 10);
    }
}
