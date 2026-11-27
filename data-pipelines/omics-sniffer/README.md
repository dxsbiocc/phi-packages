# Omics Sniffer Data Pipeline

Curated cross-repository public omics metadata index and Parquet-native incremental auto-sync pipeline backing `mcp-servers/omics-sniffer` (`/omics-sniffer/mcp` on `phi-omics-gateway`).

**Strict Code / Data Separation**:
- Only source code (`schema_v1.sql`, `pipeline/`, `scanners/`, `scripts/`, `deploy/`) lives in git.
- Mutable runtime directories (`data/`, `serving/`, `logs/`) live outside the git checkout via `PHI_OMICS_SNIFFER_STATE_DIR` (e.g. `/var/lib/phi-omics-sniffer`), and `serving/` (`releases/v_*` + atomic relative symlink `current`) is mounted read-only into `phi-omics-gateway` via `PHI_OMICS_SNIFFER_DATA_DIR=/var/lib/phi-omics-sniffer/serving`.

## Quick start

```bash
pip install -e '.[test]'
pytest

# Check incremental sync schedule against external state dir
PHI_OMICS_SNIFFER_STATE_DIR=/var/lib/phi-omics-sniffer \
  python3 pipeline/incremental_parquet_sync.py --status

# Run incremental Parquet sync (Tier-1 daily GEO / ArrayExpress / PRIDE + atomic symlink swap)
PHI_OMICS_SNIFFER_STATE_DIR=/var/lib/phi-omics-sniffer \
  bash scripts/run_sync_pipeline.sh

# Bootstrap initial Parquet snapshot from an existing SQLite omics_v1.db (one-time)
PHI_OMICS_SNIFFER_STATE_DIR=/var/lib/phi-omics-sniffer \
  python3 pipeline/export_parquet.py --full
```

## Co-deployment with `phi-site` & `mcp-servers`

1. Place initial Parquet snapshot under `/var/lib/phi-omics-sniffer/serving/releases/v_<timestamp>` with relative symlink `/var/lib/phi-omics-sniffer/serving/current -> releases/v_<timestamp>`.
2. Start the unified MCP gateway (`mcp-servers/docker-compose.yml`) with `PHI_OMICS_SNIFFER_DATA_DIR=/var/lib/phi-omics-sniffer/serving`:
   - Because `serving/` is mounted as a parent directory and `current` is a relative symlink (`releases/v_...`), atomic swaps performed by `incremental_parquet_sync.py` become visible to `phi-omics-gateway` immediately with **zero container restarts**.
3. Enable 6-hour automatic incremental updates either via:
   - **Docker Compose**: `docker compose --profile sync up -d` in `mcp-servers/docker-compose.yml`, or
   - **systemd timer**: copy `deploy/systemd/phi-omics-sniffer-sync.{service,timer}` to `/etc/systemd/system/` and run `systemctl enable --now phi-omics-sniffer-sync.timer`.

## Architecture

- **`schema_v1.sql`**: Canonical tables (`datasets`, `samples`, `dataset_artifacts`, `dataset_members`, `publications`, `dataset_publications`) + `sync_cursors`.
- **`pipeline/paths.py`**: Resolves `CODE_ROOT` (`schema_v1.sql`, `sync_policy.json`) separately from `PHI_OMICS_SNIFFER_STATE_DIR` (`data/`, `serving/`, `logs/`).
- **`pipeline/incremental_parquet_sync.py`**: Serverless Parquet-native incremental updater (no 15+ GB main SQLite needed in production). Fetches deltas into a temporary SQLite, runs two-stage DuckDB L1 cross-repository deduplication (`dataset_members` accession & DOI map), rewrites only touched `source=*` partitions, hardlinks unchanged partitions, validates row counts against regression, and atomically swaps `serving/current`.
- **`pipeline/export_parquet.py` & `pipeline/merge_staging_to_main.py`**: Initial/batch staging SQLite merger and source-partitioned ZSTD Parquet exporter.
- **`scanners/`**: 37 repository scanners covering GEO, SRA, ArrayExpress, PRIDE, MassIVE, MetaboLights, Metabolomics Workbench, CELLxGENE, Single Cell Portal, ENCODE, Roadmap/4DN, ChIP-Atlas, GTEx, GDC/TCGA, ICGC, cBioPortal, CPTAC, PDC, ICPC, DepMap, Drug Screens, GenomeCRISPR, BioGRID ORCS, HPA, ProteomicsDB, PeptideAtlas, iProX/jPOST, HuBMAP, HTAN, ImmPort, GWAS Catalog, FinnGen, Synapse, CNSA, and STOmics.
