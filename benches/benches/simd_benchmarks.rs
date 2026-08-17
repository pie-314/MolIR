use criterion::{black_box, criterion_group, criterion_main, Criterion};
use molir_core::{
    search_parallel, FingerprintRecord, MolecularFingerprint, SearchQuery, SimdBackend,
};
use rand::Rng;

fn generate_dataset(size: usize) -> Vec<FingerprintRecord> {
    let mut rng = rand::rng();
    let mut records = Vec::with_capacity(size);

    for cid in 1..=size as u32 {
        let mut words = [0u64; 32];
        for w in &mut words {
            *w = rng.random::<u64>();
        }
        let fp = MolecularFingerprint::from_words(words);
        records.push(FingerprintRecord::new(cid, fp));
    }
    records
}

fn bench_fingerprint_scan(c: &mut Criterion) {
    let dataset = generate_dataset(100_000);
    let query_fp = dataset[0].fingerprint;
    let query = SearchQuery::new(query_fp, 0.7, 50);

    c.bench_function("scan_100k_parallel", |b| {
        b.iter(|| {
            search_parallel(
                black_box(&dataset),
                black_box(&query),
                black_box(8192),
            )
        })
    });
}

criterion_group!(benches, bench_fingerprint_scan);
criterion_main!(benches);
