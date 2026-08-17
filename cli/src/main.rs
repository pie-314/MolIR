use std::time::Instant;
use clap::{Parser, Subcommand};
use molir_core::{
    search_parallel, FingerprintRecord, MolecularFingerprint, MmapFingerprintStore, SearchQuery,
    SimdBackend,
};
use rand::Rng;

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

            let query_fp = store.as_slice()[0].fingerprint;
            let query = SearchQuery::new(query_fp, threshold, top_k);

            let start = Instant::now();
            let results = search_parallel(store.as_slice(), &query, 8192);
            let elapsed = start.elapsed();

            println!("\nScan completed in {:.2} ms", elapsed.as_secs_f64() * 1000.0);
            println!("Top results:");
            for (idx, hit) in results.iter().enumerate() {
                println!("  {:2}. CID: {:<8} Similarity: {:.4}", idx + 1, hit.cid, hit.score);
            }
        }
    }

    Ok(())
}
