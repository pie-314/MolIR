#!/usr/bin/env python3
"""
PubChem Full Dataset Downloader & Multi-Chunk Ingestor
Downloads SDF chunks from NCBI FTP and streams them into fingerprints.bin using native fast Rust ingestion.
"""

import os
import sys
import subprocess
import urllib.request
import re
import argparse
from pathlib import Path

BASE_FTP_URL = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/"

def get_sdf_file_list():
    """Fetches the list of all Compound_*.sdf.gz files from NCBI FTP directory listing."""
    print(f"Fetching directory index from {BASE_FTP_URL}...")
    req = urllib.request.Request(BASE_FTP_URL, headers={"User-Agent": "MolIR-Downloader/1.0"})
    with urllib.request.urlopen(req) as resp:
        html = resp.read().decode("utf-8")

    filenames = re.findall(r'href="(Compound_\d+_\d+\.sdf\.gz)"', html)
    filenames = sorted(list(set(filenames)))
    print(f"Found {len(filenames)} total PubChem SDF chunks (~500,000 compounds each).")
    return filenames

def download_and_ingest_chunks(chunks_to_process, out_dir: Path, raw_dir: Path, keep_gz: bool = False):
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)

    cli_path = Path("./target/release/molir-cli")
    if not cli_path.exists():
        print("Building release molir-cli binary...")
        subprocess.run(["cargo", "build", "--release", "-p", "molir-cli"], check=True)

    print("============================================================")
    print(f"PubChem Ingestion: Processing {len(chunks_to_process)} chunks")
    print(f"  Output Directory: {out_dir}")
    print(f"  Raw Directory:    {raw_dir}")
    print("============================================================")

    for idx, filename in enumerate(chunks_to_process, 1):
        url = BASE_FTP_URL + filename
        gz_path = raw_dir / filename

        print(f"\n[{idx}/{len(chunks_to_process)}] Downloading {filename}...")
        urllib.request.urlretrieve(url, gz_path)
        file_size_mb = os.path.getsize(gz_path) / (1024 * 1024)
        print(f"  Downloaded {filename} ({file_size_mb:.1f} MB)")

        print(f"  Ingesting {filename} with native Rust parallel parser...")
        cmd = [
            str(cli_path), "ingest",
            "--input", str(gz_path),
            "--out-dir", str(out_dir)
        ]
        subprocess.run(cmd, check=True)

        if not keep_gz:
            try:
                os.remove(gz_path)
            except OSError:
                pass

def main():
    parser = argparse.ArgumentParser(description="PubChem Full Downloader & Ingestor")
    parser.add_argument("--out-dir", type=str, default="./data/pubchem", help="Output directory for unified fingerprints.bin")
    parser.add_argument("--raw-dir", type=str, default="./data/raw/pubchem", help="Temporary directory for downloaded .sdf.gz files")
    parser.add_argument("--chunks", type=int, default=None, help="Number of chunks to download (default: all)")
    parser.add_argument("--keep-gz", action="store_true", help="Keep downloaded .sdf.gz files on disk")
    args = parser.parse_args()

    file_list = get_sdf_file_list()
    if args.chunks:
        file_list = file_list[:args.chunks]

    download_and_ingest_chunks(file_list, Path(args.out_dir), Path(args.raw_dir), args.keep_gz)

if __name__ == "__main__":
    main()
