use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::time::Instant;
use clap::{Parser, Subcommand};
use flate2::read::GzDecoder;
use molir_core::{
    search_parallel, DatasetManifest, FingerprintRecord, MolecularFingerprint,
    MmapFingerprintStore, SearchQuery, SimdBackend,
};
use rand::Rng;
use rayon::prelude::*;

#[derive(Parser)]
#[command(name = "molir")]
#[command(about = "MolIR: Molecular Information Retrieval Engine CLI", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Inspect CPU SIMD features and system environment
    Info,

    /// High-speed parallel ingestion of chemical datasets (SDF, SDF.GZ, TSV, SMILES)
    Ingest {
        /// Input file path (e.g. data/raw/pubchem/Compound_000000001_000500000.sdf.gz or .tsv/.smi)
        #[arg(short = 'i', long)]
        input: String,

        /// Output directory for fingerprints.bin and manifest.json
        #[arg(short = 'o', long, default_value = "./data/pubchem")]
        out_dir: String,

        /// Maximum number of records to ingest (optional, default: all)
        #[arg(short = 'm', long)]
        max: Option<usize>,
    },

    /// Run an in-memory synthetic search throughput benchmark
    Bench {
        /// Number of molecules to simulate (e.g. 1000000)
        #[arg(short = 'c', long, default_value = "1000000")]
        count: usize,

        /// Tanimoto similarity threshold
        #[arg(short = 't', long, default_value = "0.7")]
        threshold: f32,

        /// Top-K results to collect
        #[arg(short = 'k', long, default_value = "50")]
        top_k: usize,

        /// Number of query repetitions
        #[arg(short = 'r', long, default_value = "10")]
        repeats: usize,
    },

    /// Search a binary dataset file directly
    Scan {
        /// Path to fingerprints.bin dataset
        #[arg(short = 'd', long)]
        dataset: String,

        /// Query chemical structure as a SMILES string (e.g. "CC(=O)Oc1ccccc1C(=O)O")
        #[arg(short = 's', long)]
        smiles: Option<String>,

        /// Tanimoto similarity threshold
        #[arg(short = 't', long, default_value = "0.7")]
        threshold: f32,

        /// Top-K results to collect
        #[arg(short = 'k', long, default_value = "50")]
        top_k: usize,
    },
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Info => {
            println!("=== MolIR Engine Information ===");
            println!("Version:          {}", env!("CARGO_PKG_VERSION"));
            println!("Detected SIMD:    {:?}", SimdBackend::detect());
            println!("Target CPU Arch:  {}", std::env::consts::ARCH);
            println!("Target OS:        {}", std::env::consts::OS);
        }

        Commands::Ingest { input, out_dir, max } => {
            run_ingest(&input, &out_dir, max)?;
        }

        Commands::Bench {
            count,
            threshold,
            top_k,
            repeats,
        } => {
            println!("Generating {} synthetic random molecular fingerprints...", count);
            let mut rng = rand::rng();
            let mut records = Vec::with_capacity(count);

            for cid in 1..=count as u32 {
                let mut words = [0u64; 32];
                for w in &mut words {
                    *w = rng.random::<u64>();
                }
                let fp = MolecularFingerprint::from_words(words);
                records.push(FingerprintRecord::new(cid, fp));
            }

            println!("Running {} benchmark iterations (Threshold: {}, Top-K: {})...", repeats, threshold, top_k);
            let query_fp = records[0].fingerprint;
            let query = SearchQuery::new(query_fp, threshold, top_k);

            let mut total_duration = std::time::Duration::ZERO;
            for i in 1..=repeats {
                let start = Instant::now();
                let results = search_parallel(&records, &query, 8192);
                let elapsed = start.elapsed();
                total_duration += elapsed;

                let mps = (count as f64 / 1_000_000.0) / elapsed.as_secs_f64();
                println!(
                    "  Run {:2}/{}: {:>7.2} ms | {:>7.2} M fingerprints/sec | Found: {}",
                    i,
                    repeats,
                    elapsed.as_secs_f64() * 1000.0,
                    mps,
                    results.len()
                );
            }

            let avg_duration = total_duration / repeats as u32;
            let avg_mps = (count as f64 / 1_000_000.0) / avg_duration.as_secs_f64();
            println!("\n=== Benchmark Summary ===");
            println!("Average Latency:    {:.2} ms", avg_duration.as_secs_f64() * 1000.0);
            println!("Average Throughput: {:.2} Million fingerprints/sec", avg_mps);
        }

        Commands::Scan {
            dataset,
            smiles,
            threshold,
            top_k,
        } => {
            println!("Opening memory-mapped dataset: {}", dataset);
            let store = MmapFingerprintStore::open(&dataset)?;
            println!("Loaded {} records.", store.len());

            if store.is_empty() {
                println!("Dataset is empty.");
                return Ok(());
            }

            let query_fp = if let Some(ref smi) = smiles {
                println!("Parsing query SMILES: {}", smi);
                let fp = MolecularFingerprint::from_smiles(smi)?;
                println!(
                    "Generated ECFP4 fingerprint (Popcount: {} bits set)",
                    fp.popcount()
                );
                fp
            } else {
                println!(
                    "No query SMILES specified; using first record (CID: {}) as sample query.",
                    store.as_slice()[0].cid
                );
                store.as_slice()[0].fingerprint
            };

            let query = SearchQuery::new(query_fp, threshold, top_k);

            let start = Instant::now();
            let results = search_parallel(store.as_slice(), &query, 8192);
            let elapsed = start.elapsed();

            println!(
                "\nScan completed in {:.2} ms (Found: {} matches above threshold {})",
                elapsed.as_secs_f64() * 1000.0,
                results.len(),
                threshold
            );
            println!("Top results:");
            for (idx, hit) in results.iter().enumerate() {
                println!(
                    "  {:2}. CID: {:<8} Similarity: {:.4}",
                    idx + 1,
                    hit.cid,
                    hit.score
                );
            }
        }
    }

    Ok(())
}

fn run_ingest(input_path: &str, out_dir: &str, max_records: Option<usize>) -> anyhow::Result<()> {
    let start_time = Instant::now();
    let out_dir_path = Path::new(out_dir);
    std::fs::create_dir_all(out_dir_path)?;

    let bin_path = out_dir_path.join("fingerprints.bin");
    let manifest_path = out_dir_path.join("manifest.json");
    let meta_path = out_dir_path.join("metadata.tsv");

    println!("============================================================");
    println!("MolIR High-Speed Chemical Dataset Ingestion");
    println!("  Input Source: {}", input_path);
    println!("  Output Dir:   {}", out_dir);
    if let Some(m) = max_records {
        println!("  Record Limit: {}", m);
    }
    println!("============================================================");

    let file = File::open(input_path)?;
    let reader: Box<dyn BufRead> = if input_path.ends_with(".gz") {
        Box::new(BufReader::with_capacity(2 * 1024 * 1024, GzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(2 * 1024 * 1024, file))
    };

    let mut out_bin = BufWriter::with_capacity(4 * 1024 * 1024, File::create(&bin_path)?);
    let mut out_meta = BufWriter::with_capacity(1024 * 1024, File::create(&meta_path)?);
    writeln!(out_meta, "cid\tname\tcanonical_smiles")?;

    let is_sdf = input_path.contains(".sdf") || input_path.contains(".SDF");
    let mut total_ingested: usize = 0;

    if is_sdf {
        let mut current_block_lines: Vec<String> = Vec::with_capacity(200);
        let mut batch: Vec<(String, String)> = Vec::with_capacity(25000);

        for line_res in reader.lines() {
            let line = line_res?;
            if line.starts_with("$$$$") {
                if let Some((name, smiles)) = parse_sdf_block(&current_block_lines) {
                    batch.push((name, smiles));
                }
                current_block_lines.clear();

                if batch.len() >= 25000 {
                    let processed = process_and_write_batch(&batch, total_ingested, &mut out_bin, &mut out_meta)?;
                    total_ingested += processed;
                    batch.clear();

                    let elapsed = start_time.elapsed().as_secs_f64();
                    let rate = total_ingested as f64 / elapsed.max(0.001);
                    println!("  Ingested {:>8} compounds ({:>7.0} compounds/sec)...", total_ingested, rate);

                    if let Some(max) = max_records {
                        if total_ingested >= max {
                            break;
                        }
                    }
                }
            } else {
                current_block_lines.push(line);
            }
        }

        if !current_block_lines.is_empty() {
            if let Some((name, smiles)) = parse_sdf_block(&current_block_lines) {
                batch.push((name, smiles));
            }
        }

        if !batch.is_empty() && max_records.map_or(true, |m| total_ingested < m) {
            let processed = process_and_write_batch(&batch, total_ingested, &mut out_bin, &mut out_meta)?;
            total_ingested += processed;
        }
    } else {
        // Tab-separated or space-separated format (ChEMBL / TSV)
        let mut batch: Vec<(String, String)> = Vec::with_capacity(25000);
        let mut line_num: usize = 0;

        for line_res in reader.lines() {
            let line = line_res?;
            line_num += 1;
            if line_num == 1 && (line.starts_with("chembl_id") || line.starts_with("cid")) {
                continue;
            }

            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 2 {
                let id_or_name = parts[0].trim().to_string();
                let smiles = parts[1].trim().to_string();
                if !smiles.is_empty() {
                    batch.push((id_or_name, smiles));
                }
            }

            if batch.len() >= 25000 {
                let processed = process_and_write_batch(&batch, total_ingested, &mut out_bin, &mut out_meta)?;
                total_ingested += processed;
                batch.clear();

                let elapsed = start_time.elapsed().as_secs_f64();
                let rate = total_ingested as f64 / elapsed.max(0.001);
                println!("  Ingested {:>8} compounds ({:>7.0} compounds/sec)...", total_ingested, rate);

                if let Some(max) = max_records {
                    if total_ingested >= max {
                        break;
                    }
                }
            }
        }

        if !batch.is_empty() && max_records.map_or(true, |m| total_ingested < m) {
            let processed = process_and_write_batch(&batch, total_ingested, &mut out_bin, &mut out_meta)?;
            total_ingested += processed;
        }
    }

    out_bin.flush()?;
    out_meta.flush()?;

    let manifest = DatasetManifest {
        dataset_name: Path::new(input_path).file_stem().unwrap_or_default().to_string_lossy().to_string(),
        version: "1.0.0".to_string(),
        fingerprint_type: "ECFP4".to_string(),
        radius: 2,
        bit_length: 2048,
        record_count: total_ingested as u64,
        endianness: "little".to_string(),
    };

    let manifest_json = serde_json::to_string_pretty(&manifest)?;
    std::fs::write(&manifest_path, manifest_json)?;

    let total_elapsed = start_time.elapsed();
    let file_size_mb = std::fs::metadata(&bin_path)?.len() as f64 / (1024.0 * 1024.0);

    println!("============================================================");
    println!("Ingestion Complete!");
    println!("  Total Ingested:   {} molecules", total_ingested);
    println!("  Time Elapsed:     {:.2} seconds", total_elapsed.as_secs_f64());
    println!("  Average Speed:    {:.0} molecules/sec", total_ingested as f64 / total_elapsed.as_secs_f64().max(0.001));
    println!("  Binary Output:    {} ({:.2} MB)", bin_path.display(), file_size_mb);
    println!("  Manifest:         {}", manifest_path.display());
    println!("  Metadata:         {}", meta_path.display());
    println!("============================================================");

    Ok(())
}

fn parse_sdf_block(lines: &[String]) -> Option<(String, String)> {
    if lines.len() < 3 {
        return None;
    }

    let header_name = lines[0].trim();
    let mut cid = String::new();
    let mut smiles = String::new();
    let mut iupac_name = String::new();

    let mut i = 0;
    while i < lines.len() {
        let line = lines[i].trim();
        if line == "> <PUBCHEM_COMPOUND_CID>" {
            if i + 1 < lines.len() {
                cid = lines[i + 1].trim().to_string();
                i += 1;
            }
        } else if line == "> <PUBCHEM_SMILES>"
            || line == "> <PUBCHEM_OPENEYE_CAN_SMILES>"
            || line == "> <PUBCHEM_CONNECTIVITY_SMILES>"
            || line == "> <SMILES>"
            || line == "> <CANONICAL_SMILES>"
        {
            if i + 1 < lines.len() && smiles.is_empty() {
                smiles = lines[i + 1].trim().to_string();
                i += 1;
            }
        } else if line == "> <PUBCHEM_IUPAC_OPENEYE_NAME>"
            || line == "> <PUBCHEM_IUPAC_NAME>"
            || line == "> <PUBCHEM_IUPAC_TRADITIONAL_NAME>"
        {
            if i + 1 < lines.len() && iupac_name.is_empty() {
                iupac_name = lines[i + 1].trim().to_string();
                i += 1;
            }
        }
        i += 1;
    }

    if smiles.is_empty() {
        return None;
    }

    let final_name = if !iupac_name.is_empty() {
        iupac_name
    } else if !header_name.is_empty() {
        header_name.to_string()
    } else if !cid.is_empty() {
        format!("CID-{}", cid)
    } else {
        "Compound".to_string()
    };

    Some((final_name, smiles))
}

fn process_and_write_batch(
    batch: &[(String, String)],
    start_cid: usize,
    out_bin: &mut BufWriter<File>,
    out_meta: &mut BufWriter<File>,
) -> anyhow::Result<usize> {
    // Process entire batch in parallel across CPU cores with Rayon
    let results: Vec<(String, String, Option<FingerprintRecord>)> = batch
        .par_iter()
        .enumerate()
        .map(|(i, (name, smiles))| {
            let cid = (start_cid + i + 1) as u32;
            if let Ok(fp) = MolecularFingerprint::from_smiles(smiles) {
                if fp.popcount() > 0 {
                    let record = FingerprintRecord::new(cid, fp);
                    return (name.clone(), smiles.clone(), Some(record));
                }
            }
            (name.clone(), smiles.clone(), None)
        })
        .collect();

    let mut written = 0;
    for (name, smiles, maybe_record) in results {
        if let Some(record) = maybe_record {
            let bytes = unsafe {
                std::slice::from_raw_parts(
                    &record as *const FingerprintRecord as *const u8,
                    std::mem::size_of::<FingerprintRecord>(),
                )
            };
            out_bin.write_all(bytes)?;
            writeln!(out_meta, "{}\t{}\t{}", record.cid, name, smiles)?;
            written += 1;
        }
    }

    Ok(written)
}
