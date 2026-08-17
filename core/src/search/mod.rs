use std::cmp::Ordering;
use std::collections::BinaryHeap;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::fingerprint::MolecularFingerprint;
use crate::simd::{scan_records, SimdBackend};
use crate::storage::FingerprintRecord;

/// Candidate match item containing molecule CID and computed similarity score.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SearchHit {
    pub cid: u32,
    pub score: f32,
}

impl Eq for SearchHit {}

impl PartialOrd for SearchHit {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SearchHit {
    #[inline]
    fn cmp(&self, other: &Self) -> Ordering {
        // Min-heap ordering by score (lowest score at the top for easy eviction)
        self.score
            .partial_cmp(&other.score)
            .unwrap_or(Ordering::Equal)
            .reverse()
    }
}

/// Bounded Min-Heap accumulator for retrieving Top-K highest scoring hits.
#[derive(Debug, Clone)]
pub struct TopKAccumulator {
    capacity: usize,
    heap: BinaryHeap<SearchHit>,
}

impl TopKAccumulator {
    pub fn new(capacity: usize) -> Self {
        Self {
            capacity,
            heap: BinaryHeap::with_capacity(capacity + 1),
        }
    }

    #[inline]
    pub fn push(&mut self, hit: SearchHit) {
        if self.heap.len() < self.capacity {
            self.heap.push(hit);
        } else if let Some(min_hit) = self.heap.peek() {
            if hit.score > min_hit.score {
                self.heap.pop();
                self.heap.push(hit);
            }
        }
    }

    pub fn into_sorted_vec(self) -> Vec<SearchHit> {
        let mut vec = self.heap.into_vec();
        // Sort descending by score
        vec.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(Ordering::Equal));
        vec
    }
}

/// Query configuration parameters for similarity search.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SearchQuery {
    pub fingerprint: MolecularFingerprint,
    pub threshold: f32,
    pub top_k: usize,
    pub simd_backend: Option<SimdBackend>,
}

impl SearchQuery {
    pub fn new(fingerprint: MolecularFingerprint, threshold: f32, top_k: usize) -> Self {
        Self {
            fingerprint,
            threshold,
            top_k,
            simd_backend: None,
        }
    }
}

/// Executes parallel similarity search using Rayon work-stealing parallelism.
pub fn search_parallel(
    records: &[FingerprintRecord],
    query: &SearchQuery,
    chunk_size: usize,
) -> Vec<SearchHit> {
    let query_popcount = query.fingerprint.popcount();
    let backend = query.simd_backend.unwrap_or_else(SimdBackend::detect);

    let local_results: Vec<Vec<(u32, f32)>> = records
        .par_chunks(chunk_size.max(1024))
        .map(|chunk| {
            scan_records(
                &query.fingerprint,
                query_popcount,
                chunk,
                query.threshold,
                backend,
            )
        })
        .collect();

    let mut accumulator = TopKAccumulator::new(query.top_k);
    for chunk_res in local_results {
        for (cid, score) in chunk_res {
            accumulator.push(SearchHit { cid, score });
        }
    }

    accumulator.into_sorted_vec()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_top_k_accumulator() {
        let mut acc = TopKAccumulator::new(3);
        acc.push(SearchHit { cid: 1, score: 0.5 });
        acc.push(SearchHit { cid: 2, score: 0.9 });
        acc.push(SearchHit { cid: 3, score: 0.2 });
        acc.push(SearchHit { cid: 4, score: 0.8 });

        let results = acc.into_sorted_vec();
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].cid, 2); // 0.9
        assert_eq!(results[1].cid, 4); // 0.8
        assert_eq!(results[2].cid, 1); // 0.5
    }
}
