use crate::fingerprint::MolecularFingerprint;
use crate::storage::FingerprintRecord;
use serde::{Deserialize, Serialize};

/// Supported SIMD execution targets.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SimdBackend {
    Scalar,
    Avx2,
    Avx512,
}

impl SimdBackend {
    /// Detects the highest performance instruction set supported by the host CPU.
    pub fn detect() -> Self {
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx512f") && is_x86_feature_detected!("avx512vpopcntdq") {
                return SimdBackend::Avx512;
            }
            if is_x86_feature_detected!("avx2") {
                return SimdBackend::Avx2;
            }
        }
        SimdBackend::Scalar
    }
}

/// Baseline scalar linear scan over records.
#[inline]
pub fn scan_scalar(
    query: &MolecularFingerprint,
    query_popcount: u32,
    records: &[FingerprintRecord],
    threshold: f32,
) -> Vec<(u32, f32)> {
    let mut results = Vec::new();
    let min_popcount = (query_popcount as f32 * threshold).ceil() as u16;
    let max_popcount = if threshold > 0.0 {
        (query_popcount as f32 / threshold).floor() as u16
    } else {
        u16::MAX
    };

    for record in records {
        // Fast popcount bound filtering
        if record.popcount < min_popcount || record.popcount > max_popcount {
            continue;
        }

        let score = query.tanimoto_with_popcounts(
            query_popcount,
            &record.fingerprint,
            record.popcount as u32,
        );

        if score >= threshold {
            results.push((record.cid, score));
        }
    }

    results
}

/// Dispatches fingerprint scan to the best available SIMD backend.
#[inline]
pub fn scan_records(
    query: &MolecularFingerprint,
    query_popcount: u32,
    records: &[FingerprintRecord],
    threshold: f32,
    backend: SimdBackend,
) -> Vec<(u32, f32)> {
    match backend {
        SimdBackend::Scalar => scan_scalar(query, query_popcount, records, threshold),
        SimdBackend::Avx2 => {
            // For now fallback to scalar; AVX2 Harley-Seal kernel will be implemented in Phase 4
            scan_scalar(query, query_popcount, records, threshold)
        }
        SimdBackend::Avx512 => {
            // For now fallback to scalar; AVX-512 VPOPCNT kernel will be implemented in Phase 4
            scan_scalar(query, query_popcount, records, threshold)
        }
    }
}
