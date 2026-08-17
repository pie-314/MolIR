use std::fs::File;
use std::path::Path;
use memmap2::Mmap;
use serde::{Deserialize, Serialize};
use crate::fingerprint::MolecularFingerprint;

/// Cache-aligned binary record representing a single molecule in the index.
#[repr(C, align(64))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct FingerprintRecord {
    pub cid: u32,
    pub popcount: u16,
    pub _reserved: [u8; 10],
    pub fingerprint: MolecularFingerprint,
}

impl FingerprintRecord {
    pub fn new(cid: u32, fingerprint: MolecularFingerprint) -> Self {
        let popcount = fingerprint.popcount() as u16;
        Self {
            cid,
            popcount,
            _reserved: [0u8; 10],
            fingerprint,
        }
    }
}

/// Metadata manifest describing a packaged molecular index dataset.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct DatasetManifest {
    pub dataset_name: String,
    pub version: String,
    pub fingerprint_type: String,
    pub radius: u32,
    pub bit_length: u32,
    pub record_count: u64,
    pub endianness: String,
}

impl Default for DatasetManifest {
    fn default() -> Self {
        Self {
            dataset_name: "molir-dataset".to_string(),
            version: "0.1.0".to_string(),
            fingerprint_type: "ECFP4".to_string(),
            radius: 2,
            bit_length: 2048,
            record_count: 0,
            endianness: "little".to_string(),
        }
    }
}

/// Zero-copy memory-mapped store for fast sequential and parallel search.
pub struct MmapFingerprintStore {
    _file: File,
    mmap: Mmap,
    record_count: usize,
}

impl MmapFingerprintStore {
    /// Opens and memory-maps a binary fingerprint store file.
    pub fn open<P: AsRef<Path>>(path: P) -> anyhow::Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        
        let record_size = std::mem::size_of::<FingerprintRecord>();
        if mmap.len() % record_size != 0 {
            anyhow::bail!(
                "Invalid binary store size: {} bytes is not a multiple of record size ({} bytes)",
                mmap.len(),
                record_size
            );
        }

        let record_count = mmap.len() / record_size;
        Ok(Self {
            _file: file,
            mmap,
            record_count,
        })
    }

    /// Returns a typed slice over the memory-mapped records with zero copying.
    pub fn as_slice(&self) -> &[FingerprintRecord] {
        let ptr = self.mmap.as_ptr() as *const FingerprintRecord;
        unsafe { std::slice::from_raw_parts(ptr, self.record_count) }
    }

    /// Number of molecular records loaded in the memory map.
    #[inline]
    pub fn len(&self) -> usize {
        self.record_count
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.record_count == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_record_size_and_alignment() {
        assert_eq!(std::mem::align_of::<FingerprintRecord>(), 64);
        assert_eq!(std::mem::size_of::<FingerprintRecord>() % 64, 0);
    }
}
