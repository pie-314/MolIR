pub mod smiles;

pub use smiles::{parse_smiles, MolecularGraph};
use serde::{Deserialize, Serialize};

/// Number of 64-bit words required to represent a 2048-bit fingerprint.
pub const FINGERPRINT_WORDS: usize = 32;

/// A 2048-bit molecular fingerprint aligned to 64-byte cache line boundaries.
#[repr(C, align(64))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct MolecularFingerprint {
    pub words: [u64; FINGERPRINT_WORDS],
}

impl Default for MolecularFingerprint {
    #[inline]
    fn default() -> Self {
        Self {
            words: [0u64; FINGERPRINT_WORDS],
        }
    }
}

impl MolecularFingerprint {
    /// Creates a new fingerprint from an array of 32 64-bit words.
    #[inline]
    pub const fn from_words(words: [u64; FINGERPRINT_WORDS]) -> Self {
        Self { words }
    }

    /// Parses a chemical SMILES string and computes its 2048-bit ECFP4 fingerprint.
    pub fn from_smiles(smiles: &str) -> anyhow::Result<Self> {
        let mut graph = parse_smiles(smiles)?;
        Ok(graph.to_ecfp4())
    }

    /// Creates a fingerprint with all bits cleared (all zeroes).
    #[inline]
    pub const fn zeros() -> Self {
        Self {
            words: [0u64; FINGERPRINT_WORDS],
        }
    }

    /// Creates a fingerprint with all bits set (all ones).
    #[inline]
    pub const fn ones() -> Self {
        Self {
            words: [u64::MAX; FINGERPRINT_WORDS],
        }
    }

    /// Sets the bit at the given index (0..2048).
    #[inline]
    pub fn set_bit(&mut self, bit_index: usize) {
        if bit_index < 2048 {
            let word_idx = bit_index / 64;
            let bit_offset = bit_index % 64;
            self.words[word_idx] |= 1u64 << bit_offset;
        }
    }

    /// Checks whether the bit at the given index is set.
    #[inline]
    pub fn get_bit(&self, bit_index: usize) -> bool {
        if bit_index < 2048 {
            let word_idx = bit_index / 64;
            let bit_offset = bit_index % 64;
            (self.words[word_idx] & (1u64 << bit_offset)) != 0
        } else {
            false
        }
    }

    /// Computes the total number of set bits (population count).
    #[inline]
    pub fn popcount(&self) -> u32 {
        let mut count = 0u32;
        for &w in &self.words {
            count += w.count_ones();
        }
        count
    }

    /// Calculates the Tanimoto (Jaccard) similarity score against another fingerprint.
    ///
    /// Tanimoto(A, B) = popcount(A & B) / (popcount(A) + popcount(B) - popcount(A & B))
    #[inline]
    pub fn tanimoto(&self, other: &Self) -> f32 {
        let mut intersection = 0u32;
        let mut union_count = 0u32;

        for i in 0..FINGERPRINT_WORDS {
            intersection += (self.words[i] & other.words[i]).count_ones();
            union_count += (self.words[i] | other.words[i]).count_ones();
        }

        if union_count == 0 {
            1.0
        } else {
            intersection as f32 / union_count as f32
        }
    }

    /// Calculates Tanimoto similarity when the other fingerprint's popcount is already known.
    #[inline]
    pub fn tanimoto_with_popcounts(&self, self_popcount: u32, other: &Self, other_popcount: u32) -> f32 {
        let mut intersection = 0u32;
        for i in 0..FINGERPRINT_WORDS {
            intersection += (self.words[i] & other.words[i]).count_ones();
        }

        let union_count = self_popcount + other_popcount - intersection;
        if union_count == 0 {
            1.0
        } else {
            intersection as f32 / union_count as f32
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fingerprint_bit_operations() {
        let mut fp = MolecularFingerprint::zeros();
        assert_eq!(fp.popcount(), 0);

        fp.set_bit(0);
        fp.set_bit(63);
        fp.set_bit(64);
        fp.set_bit(2047);

        assert!(fp.get_bit(0));
        assert!(fp.get_bit(63));
        assert!(fp.get_bit(64));
        assert!(fp.get_bit(2047));
        assert!(!fp.get_bit(1));
        assert_eq!(fp.popcount(), 4);
    }

    #[test]
    fn test_tanimoto_identity_and_orthogonality() {
        let mut fp1 = MolecularFingerprint::zeros();
        fp1.set_bit(10);
        fp1.set_bit(20);

        let mut fp2 = MolecularFingerprint::zeros();
        fp2.set_bit(10);
        fp2.set_bit(20);

        let mut fp3 = MolecularFingerprint::zeros();
        fp3.set_bit(30);
        fp3.set_bit(40);

        assert!((fp1.tanimoto(&fp2) - 1.0).abs() < 1e-6);
        assert!((fp1.tanimoto(&fp3) - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_from_smiles_convenience() {
        let fp = MolecularFingerprint::from_smiles("c1ccccc1").expect("Benzene parse");
        assert!(fp.popcount() > 0);
    }
}
