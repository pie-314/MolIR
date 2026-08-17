# MolIR Development Roadmap & Progress Tracker

> **Molecular Information Retrieval Engine**  
> Comprehensive milestone roadmap with interactive checklists for tracking development progress.

---

## Overall Progress Summary

| Phase | Description | Status | Progress |
|---|---|---|---|
| **Phase 0** | Environment Setup & Chemical ETL Prototype | In Progress | `[ 4 / 7 ]` |
| **Phase 1** | Correctness-First Scalar Core & Reference Suite | Planned | `[ 0 / 8 ]` |
| **Phase 2** | Packed Binary Storage & Zero-Copy `mmap` | Planned | `[ 0 / 6 ]` |
| **Phase 3** | Parallel Query Engine (Rayon Work-Stealing) | Planned | `[ 0 / 5 ]` |
| **Phase 4** | Vectorized SIMD Search Kernels (AVX2 / AVX-512) | Planned | `[ 0 / 7 ]` |
| **Phase 5** | Candidate Buffering & Min-Heap Top-K Optimization | Planned | `[ 0 / 5 ]` |
| **Phase 6** | Exact Substructure Search & FFI Graph Matcher | Planned | `[ 0 / 6 ]` |
| **Phase 7** | Molecular Property Engine & Multi-Constraint Ranking | Planned | `[ 0 / 5 ]` |
| **Phase 8** | High-Performance REST & Streaming WebSocket API | Planned | `[ 0 / 6 ]` |
| **Phase 9** | Interactive Web UI & Ketcher Canvas Integration | Planned | `[ 0 / 6 ]` |
| **Phase 10** | Chemical Knowledge Graph (Reactions & Literature) | Planned | `[ 0 / 5 ]` |
| **Phase 11** | Production Hardening, Benchmarking & Sharding | Planned | `[ 0 / 5 ]` |

---

## Phase 0: Environment Setup & Chemical ETL Prototype
*Goal: Establish workspace architecture, acquire sample datasets, and build offline ingestion pipeline.*

- [x] **0.1 Workspace & Toolchain Initialization**
  - [x] Initialize Cargo workspace (`core`, `api`, `cli`, `benches`).
  - [x] Configure Python virtual environment with RDKit and NumPy.
  - [x] Establish standard `.gitignore`, CI workflow, and development dependencies.
- [ ] **0.2 Dataset Acquisition**
  - [ ] Download ChEMBL / PubChem sample subsets (100k, 1M, and 10M compounds).
  - [ ] Store raw source files in `.data/raw/` (ignored in git).
- [x] **0.3 ETL Preprocessor Pipeline**
  - [ ] Implement molecular sanitizer and canonicalizer (`canonical_smiles`, `inchikey`).
  - [ ] Generate 2048-bit ECFP4 / Morgan fingerprints (radius = 2).
  - [ ] Calculate baseline physicochemical properties (MW, LogP, TPSA, HBD, HBA).
  - [x] Export data to binary packed test format (`fingerprints.bin`) and metadata database (`metadata.db`).
  - [x] Generate dataset `manifest.json` containing schema, fingerprint config, and record counts.

---

## Phase 1: Correctness-First Scalar Core & Reference Suite
*Goal: Implement a reliable, portable reference search engine as ground truth for all optimizations.*

- [ ] **1.1 Fingerprint Data Structures**
  - [ ] Define `MolecularFingerprint` struct with `[u64; 32]` bitvector backing.
  - [ ] Implement bitwise population count and Jaccard / Tanimoto calculation in pure Rust.
- [ ] **1.2 Reference Scalar Search Kernel**
  - [ ] Implement linear scan similarity search with threshold filtering ($\text{Tanimoto} \ge \tau$).
  - [ ] Implement bounded Min-Heap accumulator for Top-K candidate extraction.
- [ ] **1.3 Exact Structure Hash Lookup**
  - [ ] Implement O(1) canonical InChIKey / hash lookup table.
  - [ ] Separate exact structure lookups from fingerprint similarity hot paths.
- [ ] **1.4 Verification & Testing Suite**
  - [ ] Write unit tests comparing Tanimoto calculations against RDKit ground truth.
  - [ ] Implement property-based testing (`proptest` / `quickcheck`) for mathematical invariants.
  - [ ] Establish initial Criterion.rs baseline benchmarks (`benches/scalar_baseline.rs`).

---

## Phase 2: Packed Binary Storage & Zero-Copy `mmap`
*Goal: Design cache-friendly binary layout and achieve instant startup via virtual memory mapping.*

- [ ] **2.1 Cache-Aligned Binary Layout**
  - [ ] Define `#[repr(C, align(64))]` record struct containing CID, precomputed popcount, and 2048-bit vector.
  - [ ] Implement binary serialization / deserialization with strict endianness guarantees.
- [ ] **2.2 Zero-Copy Memory Mapping**
  - [ ] Integrate `memmap2` crate for safe read-only memory mapping of `fingerprints.bin`.
  - [ ] Implement safe slice casting from mapped bytes to typed record slices (`&[FingerprintRecord]`).
  - [ ] Measure startup time and memory footprint (target < 100ms startup for 10M molecules).
- [ ] **2.3 Dataset Integrity & Version Validation**
  - [ ] Implement manifest validation on load (fingerprint radius, bits, record count, checksum).
  - [ ] Graceful error handling for corrupted or incompatible dataset versions.

---

## Phase 3: Parallel Query Engine (Rayon Work-Stealing)
*Goal: Maximize multicore CPU throughput via lock-free chunked work distribution.*

- [ ] **3.1 Rayon Parallel Iteration**
  - [ ] Partition `&[FingerprintRecord]` into cache-conscious chunks (e.g., 8,192 records/chunk).
  - [ ] Implement parallel reduction across thread-local Min-Heaps into a global Top-K result set.
- [ ] **3.2 Thread-Pool Management & Zero Allocation**
  - [ ] Reuse global Rayon thread pool to avoid per-query OS thread creation overhead.
  - [ ] Pre-allocate reusable thread-local candidate buffers.
- [ ] **3.3 Concurrency Benchmarking**
  - [ ] Benchmark scaling efficiency from 1 to 32+ cores.
  - [ ] Profile thread contention, cache invalidation, and NUMA node locality.

---

## Phase 4: Vectorized SIMD Search Kernels (AVX2 / AVX-512)
*Goal: Accelerate bitwise similarity hot loops with architecture-specific vector intrinsics.*

- [ ] **4.1 Dynamic CPU Feature Dispatch**
  - [ ] Implement safe runtime detection via `std::is_x86_feature_detected!`.
  - [ ] Route execution to AVX-512, AVX2, or Scalar fallback at runtime.
- [ ] **4.2 AVX2 Vectorized Kernel (256-bit)**
  - [ ] Implement vectorized bitwise AND (`_mm256_and_si256`) and OR (`_mm256_or_si256`).
  - [ ] Implement Harley-Seal SIMD popcount tree or `_mm256_popcnt_epi64` where supported.
- [ ] **4.3 AVX-512 Vectorized Kernel (512-bit)**
  - [ ] Implement 512-bit vector bitwise operations (`_mm512_and_si512`, `_mm512_or_si512`).
  - [ ] Leverage native hardware `_mm512_popcnt_epi64` (VPOPCNTDQ instruction).
- [ ] **4.4 Differential Fuzzing & Performance Verification**
  - [ ] Run automated differential fuzz tests comparing SIMD outputs against scalar reference across $10^7$ random fingerprints.
  - [ ] Benchmark speedup multipliers (Target: 4x-8x over scalar on AVX2, 8x-15x on AVX-512).

---

## Phase 5: Candidate Buffering & Min-Heap Top-K Optimization
*Goal: Optimize the candidate extraction pipeline to avoid full dataset sorting.*

- [ ] **5.1 Reusable Buffer Architecture**
  - [ ] Implement thread-local `CandidateBuffer` pooling to eliminate heap allocations per query.
  - [ ] Use `SmallVec` or fixed-capacity circular heaps for Top-K candidate collection.
- [ ] **5.2 SIMD Threshold Pruning**
  - [ ] Implement early exit pruning using upper-bound popcount bounds ($\text{Tanimoto}_{\max} = \frac{\text{popcount}(A)}{\text{popcount}(B)}$).
  - [ ] Skip vector evaluations for candidates whose precomputed popcount violates threshold bounds.
- [ ] **5.3 End-to-End Search Profiling**
  - [ ] Profile L1/L2/L3 cache misses with `perf` / `valgrind --tool=cachegrind`.
  - [ ] Measure throughput (Target: > 100M fingerprints/sec on modern desktop hardware).

---

## Phase 6: Exact Substructure Search & FFI Graph Matcher
*Goal: Bridge Rust core with C++ graph isomorphism matcher for exact chemical graph validation.*

- [ ] **6.1 C-ABI / FFI Boundary Design**
  - [ ] Define clean C-compatible header (`molir_ffi.h`) with POD data structs.
  - [ ] Create `cxx` or `bindgen` Rust FFI wrapper crate (`molir-ffi`).
- [ ] **6.2 Graph Isomorphism Engine Integration**
  - [ ] Integrate VF3 / RDKit C++ Subgraph Isomorphism engine.
  - [ ] Implement SMARTS query parser and substructure matcher.
- [ ] **6.3 Two-Stage Substructure Query Flow**
  - [ ] Stage 1: Generate substructure candidate subset using pattern/subgraph fingerprint pre-filters.
  - [ ] Stage 2: Pass candidate structures to FFI graph matcher for exact atom/bond isomorphism verification.
  - [ ] Benchmark false positive elimination and speedup over unindexed graph matching.

---

## Phase 7: Molecular Property Engine & Multi-Constraint Ranking
*Goal: Allow composite multi-parameter searches combining structural similarity with physical properties.*

- [ ] **7.1 Columnar Property Store**
  - [ ] Store aligned arrays for MW, LogP, TPSA, HBD, HBA, and Rotatable Bonds.
  - [ ] Memory-map property columns for high-speed SIMD range filtering.
- [ ] **7.2 Multi-Constraint Filter Engine**
  - [ ] Implement boolean query expressions (e.g., `similarity >= 0.75 AND mw BETWEEN 200 AND 500 AND logp <= 3.5`).
  - [ ] Vectorize range filter checks with SIMD comparisons (`_mm256_cmp_ps`).
- [ ] **7.3 Composite Multi-Objective Scoring**
  - [ ] Support custom scoring functions combining Tanimoto distance with property penalty functions.

---

## Phase 8: High-Performance REST & Streaming WebSocket API
*Goal: Expose search capabilities through an asynchronous, production-ready web service.*

- [ ] **8.1 Axum REST Service Architecture**
  - [ ] Implement `POST /api/v1/search/similarity` (JSON query payload, Top-K, threshold).
  - [ ] Implement `POST /api/v1/search/exact` (SMILES / InChIKey exact match).
  - [ ] Implement `POST /api/v1/search/substructure` (SMARTS / MOL block substructure match).
  - [ ] Implement `GET /api/v1/molecule/:cid` (Hydrate chemical structure, 2D coordinates, properties).
  - [ ] Implement `GET /api/v1/system/status` (Healthcheck, dataset metrics, CPU vector extensions).
- [ ] **8.2 Real-Time WebSocket Streaming**
  - [ ] Implement `WS /api/v1/search/stream` for streaming candidate matches as Rayon chunks resolve.
- [ ] **8.3 API Documentation & Validation**
  - [ ] Generate OpenAPI 3.0 / Swagger schema.
  - [ ] Add request validation, rate limiting, and structured JSON error responses.

---

## Phase 9: Interactive Web UI & Ketcher Canvas Integration
*Goal: Provide a modern chemical structure drawing canvas and responsive visual search interface.*

- [ ] **9.1 Web Application Setup**
  - [ ] Set up modern web application frontend.
  - [ ] Integrate Ketcher 2D chemical structure editor canvas.
- [ ] **9.2 Search Experience & Visual Controls**
  - [ ] Support interactive structure drawing, SMILES copy-pasting, and MOL file upload.
  - [ ] Dynamic similarity threshold slider with live result count preview.
  - [ ] Tabbed search mode selector (Similarity, Exact, Substructure, Property Filters).
- [ ] **9.3 Chemical Result Explorer**
  - [ ] Render 2D SVG molecular diagrams for retrieved Top-K hits.
  - [ ] Interactive property inspection drawer (MW, LogP, TPSA, Lipinski Rule of 5 status).
  - [ ] Export result sets to CSV, SDF, and SMILES lists.

---

## Phase 10: Chemical Knowledge Graph (Reactions & Literature)
*Goal: Extend MolIR beyond single molecules into reaction synthesis and scientific literature discovery.*

- [ ] **10.1 Chemical Reaction Indexing**
  - [ ] Parse Reaction SMILES / SMIRKS datasets (e.g., USPTO reaction dataset).
  - [ ] Index reactant, product, catalyst, and solvent roles.
  - [ ] Implement reaction similarity and transformation search.
- [ ] **10.2 Scientific Literature & Patent Linker**
  - [ ] Ingest PubMed / Patent metadata linking compound CIDs to DOIs and patent numbers.
  - [ ] Query interface connecting retrieved molecules to relevant research papers.

---

## Phase 11: Production Hardening, Benchmarking & Sharding
*Goal: Ensure enterprise-grade reliability, containerization, and horizontal scalability.*

- [ ] **11.1 Sharding & Distributed Partitions**
  - [ ] Partition database into discrete shards for multi-node horizontal scaling.
  - [ ] Implement scatter-gather query coordinator across cluster nodes.
- [ ] **11.2 Containerization & Deployment**
  - [ ] Build multi-stage `Dockerfile` with optimized native CPU target flags (`target-cpu=native`).
  - [ ] Provide Docker Compose setup with preloaded sample dataset and web UI.
- [ ] **11.3 Continuous Benchmarking & Quality Assurance**
  - [ ] Implement automated CI performance regression tracking with Criterion.
  - [ ] Conduct memory leak audits (Valgrind, AddressSanitizer, Miri).
