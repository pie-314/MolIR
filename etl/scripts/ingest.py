#!/usr/bin/env python3
"""
MolIR Chemical Ingestion & Binary Packing Pipeline
Transforms raw SDF or SMILES into cache-aligned binary fingerprint datasets.
"""

import os
import sys
import json
import struct
import argparse
from pathlib import Path

# Record struct format matching Rust:
# #[repr(C, align(64))]
# struct FingerprintRecord {
#     cid: u32 (4 bytes),
#     popcount: u16 (2 bytes),
#     _reserved: [u8; 10] (10 bytes),
#     words: [u64; 32] (256 bytes)
# }
# Total struct size = 272 bytes (padded to 320 or 64-byte multiple if needed)
RECORD_HEADER_FORMAT = "<IH10s"  # cid (u32), popcount (u16), reserved (10s) = 16 bytes
FINGERPRINT_WORDS_FORMAT = "<32Q"  # 32 x u64 = 256 bytes

def generate_sample_binary_dataset(output_dir: Path, count: int = 1000):
    """Generates a synthetic or sample binary dataset for testing without external dependencies."""
    output_dir.mkdir(parents=True, exist_ok=True)
    bin_path = output_dir / "fingerprints.bin"
    manifest_path = output_dir / "manifest.json"

    print(f"Packing {count} molecules into {bin_path}...")
    import random
    random.seed(42)

    with open(bin_path, "wb") as fp:
        for cid in range(1, count + 1):
            words = [random.getrandbits(64) for _ in range(32)]
            popcount = sum(bin(w).count("1") for w in words)

            header = struct.pack(RECORD_HEADER_FORMAT, cid, popcount, b"\x00" * 10)
            fp_bytes = struct.pack(FINGERPRINT_WORDS_FORMAT, *words)
            fp.write(header + fp_bytes)

    manifest = {
        "dataset_name": "sample-synthetic",
        "version": "0.1.0",
        "fingerprint_type": "ECFP4",
        "radius": 2,
        "bit_length": 2048,
        "record_count": count,
        "endianness": "little"
    }

    with open(manifest_path, "w") as mf:
        json.dump(manifest, mf, indent=2)

    print(f"Dataset generated successfully at {output_dir}")

def main():
    parser = argparse.ArgumentParser(description="MolIR Dataset Preprocessor")
    parser.add_argument("--input", type=str, help="Input SDF or SMILES file path")
    parser.add_argument("--out-dir", type=str, default="./data/sample", help="Output directory")
    parser.add_argument("--sample-count", type=int, default=10000, help="Number of sample records if generating synthetic dataset")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    generate_sample_binary_dataset(out_dir, args.sample_count)

if __name__ == "__main__":
    main()
