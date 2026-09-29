//! # Sabertooth
//!
//! A pure-Rust reimplementation of the profile/PSSM sensitivity core of
//! [MMseqs2](https://github.com/soedinglab/MMseqs2). It provides:
//!
//! - a faithful port of MMseqs2's substitution-matrix reconstruction
//!   ([`matrix`]) and profile/PSSM construction ([`profile`]),
//! - a k-mer prefilter with similar-k-mer branch-and-bound expansion
//!   ([`prefilter`]),
//! - Smith–Waterman–Gotoh local alignment for sequences and profiles
//!   ([`align`]),
//! - Karlin–Altschul E-value statistics ([`evalue`]), and
//! - a search pipeline tying them together ([`search`]).
//!
//! The design goal is remote-homology sensitivity via PSSM search, suitable as
//! the search backend for a Rust reimplementation of geNomad-style profile
//! scanning.

pub mod align;
pub mod alphabet;
pub mod banner;
pub mod cluster;
#[cfg(feature = "distributed")]
pub mod distributed;
pub mod evalue;
pub mod fasta;
pub mod matrix;
pub mod mmdb;
pub mod msa;
pub mod nucl;
pub mod prefilter;
pub mod profile;
pub mod search;
pub mod taxonomy;
pub mod simd;

/// Crate version string (from Cargo).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
