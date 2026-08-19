# Changelog

All notable changes to the MolIR project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Added
- **Native Chemical SMILES Parser & ECFP4 Generator**:
  - Added pure Rust chemical graph parser supporting organic and aromatic atoms, charges, rings, and branching.
  - Added radius-2 circular Morgan fingerprinting with canonical FNV-1a neighborhood hashing.
  - Exposed `MolecularFingerprint::from_smiles(&str)` in `molir-core`.
- **SMILES Query Support in API & CLI**:
  - Updated `POST /api/v1/search/similarity` in `molir-api` to accept raw chemical SMILES queries directly in JSON payloads.
  - Added `--smiles` / `-s` flag to `molir-cli scan` for instant terminal-based molecular similarity search.
- **High-Performance Dataset Ingestion**:
  - Added native multi-threaded `molir-cli ingest` capable of streaming and indexing SDF, SDF.GZ, and TSV files at over 50,000 compounds/sec.
  - Added `download_pubchem.py` for automated chunked downloads and streaming ingestion of NCBI PubChem datasets.
  - Standardized `FingerprintRecord` binary layout to 320 bytes (64-byte cache-aligned header + 256-byte vector).
- **Verified Dataset Ingestion**:
  - Successfully indexed 430,805 real PubChem compounds into `data/pubchem/fingerprints.bin` (131.47 MB).
  - Verified exact 100% identity and analogue discovery for real compounds (e.g. acetaldehyde, indole pyrroles) in sub-10ms latency.

---

## [0.1.0-alpha] - 2026-08-17

### Added
- **Core Engine Architecture**:
  - 2048-bit cache-aligned `MolecularFingerprint` struct.
  - Zero-copy virtual memory mapping via `memmap2`.
  - Multi-threaded search with Rayon work-stealing (`search_parallel`).
  - Min-Heap `TopKAccumulator` candidate ranking.
  - Dynamic CPU SIMD architecture detection (`AVX-512`, `AVX2`, `Scalar`).
- **REST API Server**:
  - Axum web server with health check, status, similarity search, and molecule hydration endpoints.
- **Command Line Tool**:
  - CLI commands for hardware inspection, synthetic benchmarking, dataset ingestion, and similarity scans.
- **Documentation & Research**:
  - Comprehensive architectural specification in `architecture.md`.
  - Technical mathematical monograph in `MolIR_Technical_Monograph.tex`.
  - Interactive project roadmap in `roadmap.md`.
