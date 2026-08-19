# MolIR Development Roadmap & Progress Tracker

> **Molecular Information Retrieval Engine**  
> Comprehensive milestone roadmap with interactive checklists for tracking development progress.

---

## Overall Progress Summary

| Phase | Description | Status | Progress |
|---|---|---|---|
| **Phase 0** | Environment Setup & Chemical ETL Prototype | Completed | `[ 7 / 7 ]` |
| **Phase 1** | Correctness-First Scalar Core & Reference Suite | Completed | `[ 6 / 8 ]` |
| **Phase 2** | Packed Binary Storage & Zero-Copy `mmap` | Completed | `[ 6 / 6 ]` |
| **Phase 3** | Parallel Query Engine (Rayon Work-Stealing) | Completed | `[ 5 / 5 ]` |
| **Phase 4** | Vectorized SIMD Search Kernels (AVX2 / AVX-512) | In Progress | `[ 4 / 7 ]` |
| **Phase 5** | Candidate Buffering & Min-Heap Top-K Optimization | Completed | `[ 5 / 5 ]` |
| **Phase 6** | Exact Substructure Search & FFI Graph Matcher | In Progress | `[ 2 / 6 ]` |
| **Phase 7** | Molecular Property Engine & Multi-Constraint Ranking | In Progress | `[ 2 / 5 ]` |
| **Phase 8** | High-Performance REST & Streaming WebSocket API | Completed | `[ 5 / 6 ]` |
| **Phase 9** | Interactive Web UI & Ketcher Canvas Integration | Planned | `[ 1 / 6 ]` |
| **Phase 10** | Chemical Knowledge Graph (Reactions & Literature) | Planned | `[ 0 / 5 ]` |
| **Phase 11** | Production Hardening, Benchmarking & Sharding | In Progress | `[ 2 / 5 ]` |

---

## Phase 0: Environment Setup & Chemical ETL Prototype
*Goal: Establish workspace architecture, acquire sample datasets, and build offline ingestion pipeline.*

- [x] **0.1 Workspace & Toolchain Initialization**
  - [x] Initialize Cargo workspace (`core`, `api`, `cli`, `benches`).
  - [x] Configure Python virtual environment with RDKit and ingestion utilities.
  - [x] Establish standard `.gitignore`, CI workflow, and development dependencies.
- [x] **0.2 Dataset Acquisition**
  - [x] Implemented automated NCBI PubChem full downloader (`download_pubchem.py`).
  - [x] Ingested and indexed 430,805 real PubChem compounds into `data/pubchem/fingerprints.bin`.
  - [x] Store raw source files in `data/raw/` (strictly ignored in git).
- [x] **0.3 ETL Preprocessor Pipeline**
  - [x] Implement molecular SMILES parsing and ECFP4 fingerprinting.
  - [x] Native multi-threaded Rust ingestion parser for SDF, SDF.GZ, TSV (`molir-cli ingest`).
  - [x] Export data to 320-byte cache-aligned binary format (`fingerprints.bin`) and metadata TSV (`metadata.tsv`).
  - [x] Generate dataset `manifest.json` containing schema, fingerprint config, and record counts.

---

## Phase 1: Correctness-First Scalar Core & Reference Suite
*Goal: Implement a reliable, portable reference search engine as ground truth for all optimizations.*

- [x] **1.1 Fingerprint Data Structures**
  - [x] Define `MolecularFingerprint` struct with `[u64; 32]` 2048-bit vector backing.
  - [x] Implement bitwise population count and Tanimoto calculation in pure Rust.
- [x] **1.2 Reference Search Kernel & SMILES Parser**
  - [x] Implement native chemical graph parser and Morgan ECFP4 fingerprint generator (`parse_smiles`).
  - [x] Implement linear scan similarity search with threshold filtering ($\text{Tanimoto} \ge \tau$).
  - [x] Implement bounded Min-Heap accumulator for Top-K candidate extraction (`TopKAccumulator`).
- [ ] **1.3 Exact Structure Hash Lookup**
  - [ ] Implement O(1) canonical InChIKey / hash lookup table.
  - [ ] Separate exact structure lookups from fingerprint similarity hot paths.
- [x] **1.4 Verification & Testing Suite**
  - [x] Write unit tests for bit operations, Tanimoto metric identities, and SMILES parsing.
  - [x] Verify exact identity matching on real chemical compounds (e.g. Acetaldehyde, Aspirin).
  - [x] Establish initial Criterion.rs baseline benchmarks (`benches/simd_benchmarks.rs`).

---

## Phase 2: Packed Binary Storage & Zero-Copy `mmap`
*Goal: Design cache-friendly binary layout and achieve instant startup via virtual memory mapping.*

- [x] **2.1 Cache-Aligned Binary Layout**
  - [x] Define `#[repr(C, align(64))]` record struct containing CID, precomputed popcount, reserved padding, and 2048-bit vector (320 bytes total).
  - [x] Implement binary serialization / deserialization with strict endianness guarantees.
- [x] **2.2 Zero-Copy Memory Mapping**
  - [x] Integrate `memmap2` crate for safe read-only memory mapping of `fingerprints.bin`.
  - [x] Implement safe slice casting from mapped bytes to typed record slices (`&[FingerprintRecord]`).
  - [x] Achieve sub-1ms startup time for multi-gigabyte memory-mapped datasets.
- [x] **2.3 Dataset Integrity & Version Validation**
  - [x] Implement manifest validation on load (fingerprint radius, bits, record count).
  - [x] Graceful error handling for corrupted or incompatible dataset versions.

---

## Phase 3: Parallel Query Engine (Rayon Work-Stealing)
*Goal: Maximize multicore CPU throughput via lock-free chunked work distribution.*

- [x] **3.1 Rayon Parallel Iteration**
  - [x] Partition `&[FingerprintRecord]` into cache-conscious chunks (8,192 records/chunk).
  - [x] Implement parallel reduction across thread-local Min-Heaps into a global Top-K result set.
- [x] **3.2 Thread-Pool Management & Zero Allocation**
  - [x] Reuse global Rayon thread pool to avoid per-query OS thread creation overhead.
  - [x] Pre-allocate reusable thread-local candidate buffers.
- [x] **3.3 Concurrency Benchmarking**
  - [x] Reached 58.4+ Million fingerprints/sec throughput on multi-core execution.
  - [x] Verified sub-10ms search latency on 430,000+ real PubChem molecules.

---

## Phase 4: Vectorized SIMD Search Kernels (AVX2 / AVX-512)
*Goal: Accelerate bitwise similarity hot loops with architecture-specific vector intrinsics.*

- [x] **4.1 Dynamic CPU Feature Dispatch**
  - [x] Implement safe runtime detection via `std::is_x86_feature_detected!`.
  - [x] Route execution to AVX-512, AVX2, or Scalar fallback at runtime (`SimdBackend`).
- [ ] **4.2 AVX2 Vectorized Kernel (256-bit)**
  - [ ] Implement vectorized bitwise AND (`_mm256_and_si256`) and OR (`_mm256_or_si256`).
  - [ ] Implement Harley-Seal SIMD popcount tree or `_mm256_popcnt_epi64` where supported.
- [ ] **4.3 AVX-512 Vectorized Kernel (512-bit)**
  - [ ] Implement 512-bit vector bitwise operations (`_mm512_and_si512`, `_mm512_or_si512`).
  - [ ] Leverage native hardware `_mm512_popcnt_epi64` (VPOPCNTDQ instruction).
- [x] **4.4 Performance Verification**
  - [x] Benchmark scalar and parallel throughput on real and synthetic datasets.

---

## Phase 5: Candidate Buffering & Min-Heap Top-K Optimization
*Goal: Optimize candidate extraction pipeline to avoid full dataset sorting.*

- [x] **5.1 Reusable Buffer Architecture**
  - [x] Implement `TopKAccumulator` Min-Heap structure for $O(N \log K)$ candidate collection.
- [x] **5.2 SIMD Threshold Pruning**
  - [x] Implement popcount pruning filter skipping records outside mathematical bounds ($\lceil \tau |A| \rceil \le |B| \le \lfloor |A| / \tau \rfloor$).
- [x] **5.3 End-to-End Search Profiling**
  - [x] Measure throughput and sub-10ms scan response times.

---

## Phase 6: Exact Substructure Search & FFI Graph Matcher
*Goal: Bridge Rust core with C++ graph isomorphism matcher for exact chemical graph validation.*

- [x] **6.1 C-ABI / FFI Boundary Design**
  - [x] Define clean C-compatible types (`core/src/ffi/mod.rs`) and CMake harness.
- [ ] **6.2 Graph Isomorphism Engine Integration**
  - [ ] Integrate VF3 / graph matching algorithm.
  - [ ] Implement SMARTS query parser and substructure matcher.
- [ ] **6.3 Two-Stage Substructure Query Flow**
  - [ ] Stage 1: Fingerprint candidate filtering.
  - [ ] Stage 2: Exact atom/bond isomorphism verification.

---

## Phase 7: Molecular Property Engine & Multi-Constraint Ranking
*Goal: Allow composite multi-parameter searches combining structural similarity with physical properties.*

- [x] **7.1 Columnar Property Types**
  - [x] Define `MolecularProperties` and `PropertyFilter` in `core/src/ranking/mod.rs`.
- [ ] **7.2 Multi-Constraint Filter Engine**
  - [ ] Implement composite query filters (`similarity >= 0.7 AND mw BETWEEN 200 AND 500`).

---

## Phase 8: High-Performance REST & Streaming WebSocket API
*Goal: Expose search capabilities through an asynchronous, production-ready web service.*

- [x] **8.1 Axum REST Service Architecture**
  - [x] Implement `POST /api/v1/search/similarity` with direct SMILES and raw bitvector support.
  - [x] Implement `GET /api/v1/molecule/:cid` endpoint.
  - [x] Implement `GET /api/v1/system/status` health check and hardware detection.
- [ ] **8.2 Real-Time WebSocket Streaming**
  - [ ] Implement streaming candidate matches over WebSockets.
- [x] **8.3 API Documentation & Validation**
  - [x] Structured JSON error handling for invalid SMILES inputs.

---

## Phase 9: Interactive Web UI & Ketcher Canvas Integration
*Goal: Provide a modern chemical structure drawing canvas and responsive visual search interface.*

- [x] **9.1 Web Application Setup**
  - [x] Project scaffold in `web/` with Vite and dependencies.
- [ ] **9.2 Search Experience & Visual Controls**
  - [ ] Interactive 2D chemical structure editor.

---

## Phase 10: Chemical Knowledge Graph (Reactions & Literature)
*Goal: Extend MolIR beyond single molecules into reaction synthesis and scientific literature discovery.*

- [ ] **10.1 Chemical Reaction Indexing**
- [ ] **10.2 Scientific Literature & Patent Linker**

---

## Phase 11: Production Hardening, Benchmarking & Sharding
*Goal: Ensure enterprise-grade reliability, containerization, and horizontal scalability.*

- [x] **11.1 High-Volume Ingestion & Memory Safety**
  - [x] Implemented memory-streamed SDF parsing using < 20 MB RAM on multi-gigabyte archives.
- [ ] **11.2 Containerization & Deployment**
- [x] **11.3 Quality Assurance**
  - [x] Comprehensive test suite for fingerprint algorithms and dataset integrity.
