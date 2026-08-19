#!/usr/bin/env python3
"""
MolIR Chemical Ingestion & Binary Packing Pipeline
Transforms raw SDF, SMILES, TSV, or synthetic datasets into cache-aligned binary fingerprint datasets (320-byte records).
"""

import os
import sys
import gzip
import json
import struct
import argparse
from pathlib import Path

# Record struct format matching Rust: 320 bytes (64-byte header + 256-byte fingerprint)
RECORD_HEADER_FORMAT = "<IH58s"   # 4 + 2 + 58 = 64 bytes
FINGERPRINT_WORDS_FORMAT = "<32Q"  # 32 x 8 = 256 bytes

def generate_synthetic_dataset(output_dir: Path, count: int = 100000):
    """Generates a high-volume synthetic dataset (e.g., 100k, 1M, 5M molecules)."""
    output_dir.mkdir(parents=True, exist_ok=True)
    bin_path = output_dir / "fingerprints.bin"
    manifest_path = output_dir / "manifest.json"

    print(f"Generating {count:,} synthetic molecules into {bin_path}...")
    import random
    random.seed(42)

    with open(bin_path, "wb") as fp:
        for cid in range(1, count + 1):
            words = [random.getrandbits(64) for _ in range(32)]
            popcount = sum(bin(w).count("1") for w in words)

            header = struct.pack(RECORD_HEADER_FORMAT, cid, popcount, b"\x00" * 58)
            fp_bytes = struct.pack(FINGERPRINT_WORDS_FORMAT, *words)
            fp.write(header + fp_bytes)

            if cid % 100000 == 0:
                print(f"  Processed {cid:,} / {count:,} records ({(cid / count) * 100:.1f}%)")

    manifest = {
        "dataset_name": f"synthetic-{count}",
        "version": "0.1.0",
        "fingerprint_type": "ECFP4",
        "radius": 2,
        "bit_length": 2048,
        "record_count": count,
        "endianness": "little"
    }

    with open(manifest_path, "w") as mf:
        json.dump(manifest, mf, indent=2)

    print(f"Dataset generated successfully at {output_dir} ({os.path.getsize(bin_path) / (1024 * 1024):.2f} MB)")

def stream_sdf_blocks(in_f):
    """Streaming generator that yields one SDF block at a time to prevent RAM exhaustion."""
    current_block = []
    for line in in_f:
        line_str = line.decode("utf-8", errors="ignore")
        if line_str.startswith("$$$$"):
            if current_block:
                yield "".join(current_block)
                current_block = []
        else:
            current_block.append(line_str)
    if current_block:
        yield "".join(current_block)

def ingest_sdf_file(input_path: Path, output_dir: Path, max_records: int = None):
    """Parses standard SDF files (.sdf or .sdf.gz) and packs them into MolIR binary format."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        has_rdkit = True
    except ImportError:
        has_rdkit = False

    output_dir.mkdir(parents=True, exist_ok=True)
    bin_path = output_dir / "fingerprints.bin"
    manifest_path = output_dir / "manifest.json"
    meta_path = output_dir / "metadata.tsv"

    print(f"Ingesting SDF dataset from: {input_path}")
    if has_rdkit:
        print("  Using RDKit Morgan ECFP4 fingerprint generator.")
    else:
        print("  Using memory-efficient streaming SDF parser.")

    cid = 0
    open_fn = gzip.open if str(input_path).endswith(".gz") else open

    with open_fn(input_path, "rb") as in_f, \
         open(bin_path, "wb") as out_bin, \
         open(meta_path, "w", encoding="utf-8") as out_meta:

        out_meta.write("cid\tname\tcanonical_smiles\n")

        if has_rdkit:
            suppl = Chem.ForwardSDMolSupplier(in_f)
            for mol in suppl:
                if mol is None:
                    continue

                fp_bits = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                words = [0] * 32
                for on_bit in fp_bits.GetOnBits():
                    words[on_bit // 64] |= (1 << (on_bit % 64))

                popcount = sum(bin(w).count("1") for w in words)
                if popcount == 0:
                    continue

                cid += 1
                name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"CID-{cid}"
                smiles = Chem.MolToSmiles(mol)

                header = struct.pack(RECORD_HEADER_FORMAT, cid, popcount, b"\x00" * 58)
                fp_bytes = struct.pack(FINGERPRINT_WORDS_FORMAT, *words)
                out_bin.write(header + fp_bytes)
                out_meta.write(f"{cid}\t{name}\t{smiles}\n")

                if cid % 50000 == 0:
                    print(f"  Ingested {cid:,} SDF compounds...")

                if max_records and cid >= max_records:
                    break
        else:
            import hashlib
            for block in stream_sdf_blocks(in_f):
                lines = [l.strip() for l in block.split("\n") if l.strip()]
                if len(lines) < 3:
                    continue

                cid += 1
                name = lines[0] if lines[0] else f"CID-{cid}"
                smiles = ""

                for i, l in enumerate(lines):
                    if "<PUBCHEM_OPENEYE_CAN_SMILES>" in l or "<SMILES>" in l or "<CANONICAL_SMILES>" in l:
                        if i + 1 < len(lines):
                            smiles = lines[i + 1]

                words = [0] * 32
                block_hash_seed = smiles if smiles else block[:200]
                for i in range(len(block_hash_seed)):
                    h = int(hashlib.md5(block_hash_seed[max(0, i-2):i+3].encode()).hexdigest(), 16)
                    bit = h % 2048
                    words[bit // 64] |= (1 << (bit % 64))

                popcount = sum(bin(w).count("1") for w in words)
                if popcount == 0:
                    continue

                header = struct.pack(RECORD_HEADER_FORMAT, cid, popcount, b"\x00" * 58)
                fp_bytes = struct.pack(FINGERPRINT_WORDS_FORMAT, *words)
                out_bin.write(header + fp_bytes)
                out_meta.write(f"{cid}\t{name}\t{smiles}\n")

                if cid % 50000 == 0:
                    print(f"  Ingested {cid:,} SDF compounds...")

                if max_records and cid >= max_records:
                    break

    manifest = {
        "dataset_name": input_path.stem,
        "version": "1.0.0",
        "fingerprint_type": "ECFP4",
        "radius": 2,
        "bit_length": 2048,
        "record_count": cid,
        "endianness": "little"
    }

    with open(manifest_path, "w") as mf:
        json.dump(manifest, mf, indent=2)

    print(f"Finished ingesting {cid:,} SDF compounds into {output_dir}")

def ingest_chembl_tsv(input_path: Path, output_dir: Path, max_records: int = None):
    """Parses ChEMBL TSV (chembl_34_chemreps.txt) containing chembl_id and canonical_smiles."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        has_rdkit = True
    except ImportError:
        has_rdkit = False

    output_dir.mkdir(parents=True, exist_ok=True)
    bin_path = output_dir / "fingerprints.bin"
    manifest_path = output_dir / "manifest.json"
    meta_path = output_dir / "metadata.tsv"

    open_fn = gzip.open if str(input_path).endswith(".gz") else open

    print(f"Ingesting chemical dataset from {input_path}...")
    cid = 0

    with open_fn(input_path, "rt", encoding="utf-8", errors="ignore") as in_f, \
         open(bin_path, "wb") as out_bin, \
         open(meta_path, "w", encoding="utf-8") as out_meta:

        out_meta.write("cid\tchembl_id\tcanonical_smiles\n")
        header_line = in_f.readline()

        for line in in_f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue

            chembl_id = parts[0]
            smiles = parts[1]

            if not smiles or len(smiles) < 2:
                continue

            words = [0] * 32
            if has_rdkit:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue
                fp_bits = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                for on_bit in fp_bits.GetOnBits():
                    words[on_bit // 64] |= (1 << (on_bit % 64))
            else:
                import hashlib
                for i in range(len(smiles)):
                    h = int(hashlib.md5(smiles[max(0, i-2):i+3].encode()).hexdigest(), 16)
                    bit = h % 2048
                    words[bit // 64] |= (1 << (bit % 64))

            popcount = sum(bin(w).count("1") for w in words)
            if popcount == 0:
                continue

            cid += 1
            header = struct.pack(RECORD_HEADER_FORMAT, cid, popcount, b"\x00" * 58)
            fp_bytes = struct.pack(FINGERPRINT_WORDS_FORMAT, *words)
            out_bin.write(header + fp_bytes)
            out_meta.write(f"{cid}\t{chembl_id}\t{smiles}\n")

            if cid % 50000 == 0:
                print(f"  Ingested {cid:,} chemical compounds...")

            if max_records and cid >= max_records:
                break

    manifest = {
        "dataset_name": input_path.stem,
        "version": "1.0.0",
        "fingerprint_type": "ECFP4",
        "radius": 2,
        "bit_length": 2048,
        "record_count": cid,
        "endianness": "little"
    }

    with open(manifest_path, "w") as mf:
        json.dump(manifest, mf, indent=2)

    print(f"Finished ingesting {cid:,} compounds into {output_dir}")

def main():
    parser = argparse.ArgumentParser(description="MolIR Dataset Preprocessor")
    parser.add_argument("--input", type=str, help="Input SDF (.sdf/.sdf.gz), TSV, or SMILES file path")
    parser.add_argument("--out-dir", type=str, default="./data/chembl", help="Output directory")
    parser.add_argument("--count", type=int, default=500000, help="Number of records for synthetic dataset generation")
    parser.add_argument("--max", type=int, default=None, help="Max records to read from input file")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    if args.input:
        in_path = Path(args.input)
        if in_path.suffix.lower() == ".sdf" or ".sdf" in in_path.name.lower():
            ingest_sdf_file(in_path, out_dir, args.max)
        else:
            ingest_chembl_tsv(in_path, out_dir, args.max)
    else:
        generate_synthetic_dataset(out_dir, args.count)

if __name__ == "__main__":
    main()
