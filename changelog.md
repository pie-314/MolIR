# Changelog

All notable changes to the **MolIR** project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Planned
- **Core Engine**:
  - `MolecularFingerprint` 2048-bit bitpacked representation.
  - Scalar reference Tanimoto similarity and Min-Heap Top-K accumulator.
  - Binary packed file reader and zero-copy `memmap2` integration.
  - Rayon chunked work-stealing parallel iterator.
  - AVX2 and AVX-512 SIMD popcount vector kernels with dynamic CPU dispatch.
  - C-ABI boundary and VF3 / RDKit C++ graph isomorphism matcher.
- **ETL Pipeline**:
  - RDKit-based dataset preprocessing, canonicalization, and binary packing scripts.
  - Manifest generation and dataset versioning support.
- **API & Frontend**:
  - Axum HTTP REST and WebSocket streaming server.
  - Web frontend with integrated Ketcher 2D molecular drawing canvas.
- **Benchmarking**:
  - Criterion.rs suite for SIMD throughput, cache efficiency, and latency profiling.

---

## [0.1.0-alpha] - 2026-08-17

### Added
- **System Architecture Blueprint**: Documented complete architecture in [`architecture.md`](./architecture.md), including:
  - Vectorized SIMD two-stage retrieval pipeline.
  - Memory-mapped cache-aligned binary storage format.
  - Rayon work-stealing parallel execution model.
  - C-ABI graph isomorphism verification boundary.
  - REST and WebSocket streaming API specifications.
- **Interactive Development Roadmap**: Created [`roadmap.md`](./roadmap.md) featuring granular progress checklists across Phases 0 through 11.
- **Project Documentation**: Created comprehensive [`README.md`](./README.md) with quick start guide, system diagrams, and repository layout.
- **Project Governance**: Added [`CONTRIBUTING.md`](./CONTRIBUTING.md) and open source [`LICENSE`](./LICENSE).
- **Repository Setup**: Initialized project structure and standard `.gitignore`.

