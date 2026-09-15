//! Amino-acid alphabet handling.
//!
//! Sabertooth uses the exact same internal residue ordering as MMseqs2 so that
//! substitution scores, background frequencies and profile columns line up
//! byte-for-byte with the upstream tool:
//!
//! ```text
//! index:  0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19  20
//! residue A C D E F G H I K L M  N  P  Q  R  S  T  V  W  Y   X
//! ```
//!
//! `X` (index 20) is the "any / unknown" residue. There are 20 real amino
//! acids (`PROFILE_AA_SIZE`), matching MMseqs2's `Sequence::PROFILE_AA_SIZE`.

/// Number of real amino acids used for profiles/PSSMs (mirrors
/// `Sequence::PROFILE_AA_SIZE` in MMseqs2).
pub const PROFILE_AA_SIZE: usize = 20;

/// Alphabet size including the `X` any-residue symbol.
pub const ALPHABET_SIZE: usize = 21;

/// Internal residue index reserved for the `X` / unknown residue.
pub const ANY: u8 = 20;

/// Internal residue index used to mark a gap column in an alignment/MSA.
///
/// MMseqs2 uses the sentinel value `20` for `MultipleAlignment::GAP` in the
/// numeric MSA representation, distinct from the printable `-`. Because our
/// numeric residues only run 0..=20, we use `0xFF` as an out-of-band gap marker
/// in numeric MSA rows and treat any value `>= PROFILE_AA_SIZE` as "not a
/// countable amino acid" (gap or X) when building profiles — exactly matching
/// the MMseqs2 `aa_pos < PROFILE_AA_SIZE` guard.
pub const GAP: u8 = 0xFF;

/// Residue letters in Sabertooth/MMseqs2 internal order.
pub const NUM2AA: [u8; ALPHABET_SIZE] = *b"ACDEFGHIKLMNPQRSTVWYX";

/// Build the 256-entry lookup table mapping an ASCII byte to an internal
/// residue index. Ambiguity codes (`B`, `Z`, `J`, `U`, `O`, `*`) and anything
/// unrecognised collapse to `X` (`ANY`), matching MMseqs2's letter mapping.
pub fn build_aa2num() -> [u8; 256] {
    let mut table = [ANY; 256];
    for (idx, &letter) in NUM2AA.iter().enumerate() {
        table[letter as usize] = idx as u8;
        // lowercase maps to the same residue (a3m lowercase = insert state, but
        // the residue identity is preserved when we uppercase during parsing)
        table[letter.to_ascii_lowercase() as usize] = idx as u8;
    }
    table
}

/// Convenience: translate an ASCII residue byte to an internal index.
#[inline]
pub fn aa_to_num(byte: u8, table: &[u8; 256]) -> u8 {
    table[byte as usize]
}

/// Convenience: translate an internal index back to its residue letter.
#[inline]
pub fn num_to_aa(num: u8) -> u8 {
    if (num as usize) < ALPHABET_SIZE {
        NUM2AA[num as usize]
    } else {
        b'-'
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_letters() {
        let t = build_aa2num();
        for (idx, &letter) in NUM2AA.iter().enumerate() {
            assert_eq!(aa_to_num(letter, &t) as usize, idx);
            assert_eq!(num_to_aa(idx as u8), letter);
        }
    }

    #[test]
    fn ambiguity_collapses_to_x() {
        let t = build_aa2num();
        for &b in b"BZJUO*?." {
            assert_eq!(aa_to_num(b, &t), ANY);
        }
    }
}
