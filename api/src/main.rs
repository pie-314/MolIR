use axum::{
    extract::{Path, State},
    http::StatusCode,
    response::IntoResponse,
    routing::{get, post},
    Json, Router,
};
use clap::Parser;
use molir_core::{
    search_parallel, FingerprintRecord, MmapFingerprintStore, MolecularFingerprint, SearchHit,
    SearchQuery, SimdBackend,
};
use serde::{Deserialize, Serialize};
use std::sync::Arc;
use tower_http::cors::CorsLayer;
use tower_http::trace::TraceLayer;
use tracing_subscriber::{layer::SubscriberExt, util::SubscriberInitExt};

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "MolIR High-Performance Chemical Search API Server"
)]
struct Args {
    #[arg(short, long, default_value = "8080")]
    port: u16,

    #[arg(short, long)]
    dataset: Option<String>,
}

struct AppState {
    store: Option<MmapFingerprintStore>,
    sample_records: Vec<FingerprintRecord>,
}

#[derive(Serialize)]
struct SystemStatusResponse {
    status: &'static str,
    version: &'static str,
    simd_backend: SimdBackend,
    loaded_molecules: usize,
}

#[derive(Deserialize)]
struct SimilaritySearchRequest {
    words: Option<[u64; 32]>,
    threshold: Option<f32>,
    top_k: Option<usize>,
}

#[derive(Serialize)]
struct SimilaritySearchResponse {
    total_scanned: usize,
    matched_count: usize,
    results: Vec<SearchHit>,
}

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    tracing_subscriber::registry()
        .with(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| "molir_api=info,tower_http=info".into()),
        )
        .with(tracing_subscriber::fmt::layer())
        .init();

    let args = Args::parse();
    tracing::info!("Initializing MolIR Search Engine API...");

    let store = if let Some(path) = &args.dataset {
        tracing::info!("Loading dataset from: {}", path);
        Some(MmapFingerprintStore::open(path)?)
    } else {
        tracing::warn!("No dataset path specified, running in memory-mock mode");
        None
    };

    let app_state = Arc::new(AppState {
        store,
        sample_records: Vec::new(),
    });

    let app = Router::new()
        .route("/api/v1/system/status", get(system_status))
        .route("/api/v1/search/similarity", post(search_similarity))
        .route("/api/v1/molecule/{cid}", get(get_molecule_info))
        .layer(CorsLayer::permissive())
        .layer(TraceLayer::new_for_http())
        .with_state(app_state);

    let listener = tokio::net::TcpListener::bind(format!("0.0.0.0:{}", args.port)).await?;
    tracing::info!(
        "MolIR API server listening on http://localhost:{}",
        args.port
    );
    axum::serve(listener, app).await?;

    Ok(())
}

async fn system_status(State(state): State<Arc<AppState>>) -> Json<SystemStatusResponse> {
    let loaded = state
        .store
        .as_ref()
        .map(|s| s.len())
        .unwrap_or(state.sample_records.len());

    Json(SystemStatusResponse {
        status: "healthy",
        version: env!("CARGO_PKG_VERSION"),
        simd_backend: SimdBackend::detect(),
        loaded_molecules: loaded,
    })
}

async fn search_similarity(
    State(state): State<Arc<AppState>>,
    Json(payload): Json<SimilaritySearchRequest>,
) -> Json<SimilaritySearchResponse> {
    let fp = if let Some(words) = payload.words {
        MolecularFingerprint::from_words(words)
    } else {
        MolecularFingerprint::zeros()
    };

    let threshold = payload.threshold.unwrap_or(0.7);
    let top_k = payload.top_k.unwrap_or(50);
    let query = SearchQuery::new(fp, threshold, top_k);

    let records: &[FingerprintRecord] = if let Some(store) = &state.store {
        store.as_slice()
    } else {
        &state.sample_records
    };

    let total_scanned = records.len();
    let results = search_parallel(records, &query, 8192);
    let matched_count = results.len();

    Json(SimilaritySearchResponse {
        total_scanned,
        matched_count,
        results,
    })
}

async fn get_molecule_info(Path(cid): Path<u32>) -> impl IntoResponse {
    (
        StatusCode::OK,
        Json(serde_json::json!({
            "cid": cid,
            "name": format!("Compound-{}", cid),
            "status": "placeholder_hydrated"
        })),
    )
}
