from __future__ import annotations

import importlib
import json
import os
import sqlite3
import subprocess
import sys
from pathlib import Path

import duckdb
import pytest

CODE_DIR = Path(__file__).resolve().parents[1]
SCHEMA_SQL = (CODE_DIR / "schema_v1.sql").read_text(encoding="utf-8")


def _seed_sqlite(db_path: Path) -> None:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.executescript(SCHEMA_SQL)
    ts = "2026-11-27T10:00:00+0800"
    conn.execute(
        """
        INSERT INTO datasets(
            canonical_id, title, title_norm, organism, tissue, disease_tags,
            assay, platforms, n_samples, access, license, confidence,
            content_hash, first_seen_at, updated_at, probe_status, last_checked
        ) VALUES
        ('geo:GSE100001', 'Single-cell liver HCC atlas', 'single cell liver hcc atlas',
         'Homo sapiens', 'Liver', '["HCC","liver cancer"]', 'scRNA-seq', '["GPL24676"]',
         24, 'public', 'CC0', 'verified_live', 'h1', ?, ?, 'ok', ?),
        ('pride:PXD000123', 'HCC phosphoproteomics cohort', 'hcc phosphoproteomics cohort',
         'Homo sapiens', 'Liver', '["HCC"]', 'proteomics', '["Orbitrap"]',
         40, 'public', 'CC-BY-4.0', 'metadata_only', 'h2', ?, ?, 'ok', ?)
        """,
        (ts, ts, ts, ts, ts, ts),
    )
    conn.execute(
        """
        INSERT INTO dataset_members(dataset_id, source, accession, role, external_url) VALUES
        ('geo:GSE100001', 'geo', 'GSE100001', 'series', 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100001'),
        ('pride:PXD000123', 'pride', 'PXD000123', 'project', 'https://www.ebi.ac.uk/pride/archive/projects/PXD000123')
        """
    )
    conn.execute(
        """
        INSERT INTO publications(pmid, doi, title, authors, journal, year, abstract, url, fetched_at) VALUES
        ('38000001', '10.1038/s41586-026-00001-1', 'Atlas of HCC', 'Li X et al.', 'Nature', 2026,
         'Single-cell liver cancer study.', 'https://pubmed.ncbi.nlm.nih.gov/38000001/', ?)
        """,
        (ts,),
    )
    conn.execute(
        """
        INSERT INTO dataset_publications(dataset_id, pmid, relation, source) VALUES
        ('geo:GSE100001', '38000001', 'primary', 'geo_link')
        """
    )
    conn.execute(
        """
        INSERT INTO samples(
            sample_id, dataset_id, title, organism, tissue, disease,
            condition, sex, age, attributes, run_ids, content_hash, updated_at
        ) VALUES
        ('geo:GSM1001', 'geo:GSE100001', 'Tumor_1', 'Homo sapiens', 'Liver', 'HCC', 'tumor', 'male', '58', '{"stage":"II"}', '["SRR1001"]', 'sh1', ?),
        ('pride:S1', 'pride:PXD000123', 'Proteome_1', 'Homo sapiens', 'Liver', 'HCC', 'tumor', 'female', '61', '{}', '[]', 'sh2', ?)
        """,
        (ts, ts),
    )
    conn.execute(
        """
        INSERT INTO dataset_artifacts(
            artifact_id, dataset_id, sample_id, kind, role, title, filename,
            format, uri, size_bytes, access, license, source, derived_from,
            probe_status, last_checked, content_hash
        ) VALUES
        ('geo:GSE100001:matrix', 'geo:GSE100001', NULL, 'processed', 'counts',
         'Count matrix', 'GSE100001_counts.tsv.gz', 'tsv',
         'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100001/suppl/GSE100001_counts.tsv.gz',
         1048576, 'public', 'CC0', 'geo', NULL, 'ok', ?, 'ah1')
        """,
        (ts,),
    )
    conn.commit()
    conn.close()


def test_code_and_data_strictly_separated_export_and_incremental_merge(tmp_path: Path) -> None:
    state_dir = tmp_path / "external_state"
    env = {
        **os.environ,
        "PHI_OMICS_SNIFFER_STATE_DIR": str(state_dir),
    }

    main_db = state_dir / "data" / "omics_v1.db"
    _seed_sqlite(main_db)

    # 1. Run export_parquet.py into external_state/serving (git code tree stays untouched)
    res_export = subprocess.run(
        [sys.executable, "-B", str(CODE_DIR / "pipeline" / "export_parquet.py"), "--full"],
        cwd=str(CODE_DIR),
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )
    assert "Published atomic Parquet release" in res_export.stdout
    assert not (CODE_DIR / "data").exists()
    assert not (CODE_DIR / "serving").exists()

    current = state_dir / "serving" / "current"
    assert current.is_symlink()
    manifest_1 = json.loads((current / "_manifest.json").read_text(encoding="utf-8"))
    assert manifest_1["datasets_count"] == 2
    assert manifest_1["samples_sources"] == 2

    # 2. Build a delta SQLite with:
    #    (a) an ArrayExpress mirror of GSE100001 (should L1-deduplicate into geo:GSE100001!)
    #    (b) a brand-new GEO series geo:GSE100002
    delta_db = tmp_path / "delta.db"
    dconn = sqlite3.connect(delta_db)
    dconn.executescript(SCHEMA_SQL)
    ts2 = "2026-11-27T16:00:00+0800"
    dconn.execute(
        """
        INSERT INTO datasets(
            canonical_id, title, title_norm, organism, tissue, disease_tags,
            assay, platforms, n_samples, access, license, confidence,
            content_hash, first_seen_at, updated_at, probe_status, last_checked
        ) VALUES
        ('ae:E-GEOD-100001', 'Mirror of GSE100001', 'mirror of gse100001',
         'Homo sapiens', 'Liver', NULL, 'scRNA-seq', NULL, 24, 'public', NULL, 'metadata_only', 'dh1', ?, ?, 'ok', ?),
        ('geo:GSE100002', 'Spatial transcriptomics of cholangiocarcinoma', 'spatial transcriptomics of cholangiocarcinoma',
         'Homo sapiens', 'Bile duct', '["CCA"]', 'Visium', '["GPL24676"]', 12, 'public', 'CC0', 'metadata_only', 'dh2', ?, ?, 'ok', ?)
        """,
        (ts2, ts2, ts2, ts2, ts2, ts2),
    )
    dconn.execute(
        """
        INSERT INTO dataset_members(dataset_id, source, accession, role, external_url) VALUES
        ('ae:E-GEOD-100001', 'arrayexpress', 'E-GEOD-100001', 'entry', 'https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-GEOD-100001'),
        ('ae:E-GEOD-100001', 'geo', 'GSE100001', 'series_link', 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100001'),
        ('geo:GSE100002', 'geo', 'GSE100002', 'series', 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100002')
        """
    )
    dconn.execute(
        """
        INSERT INTO samples(
            sample_id, dataset_id, title, organism, tissue, disease,
            condition, sex, age, attributes, run_ids, content_hash, updated_at
        ) VALUES
        ('geo:GSM2001', 'geo:GSE100002', 'CCA_Spot_1', 'Homo sapiens', 'Bile duct', 'CCA', 'tumor', NULL, NULL, '{}', '[]', 'sh3', ?)
        """,
        (ts2,),
    )
    dconn.commit()
    dconn.close()

    # Delete main SQLite to prove incremental_parquet_sync.py is 100% Parquet-native!
    main_db.unlink()

    res_merge = subprocess.run(
        [
            sys.executable,
            "-B",
            str(CODE_DIR / "pipeline" / "incremental_parquet_sync.py"),
            "--delta-sqlite",
            str(delta_db),
            "--source",
            "geo",
        ],
        cwd=str(CODE_DIR),
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )
    assert "Published atomic incremental release" in res_merge.stdout

    manifest_2 = json.loads((current / "_manifest.json").read_text(encoding="utf-8"))
    # 2 original + 1 new (ae:E-GEOD-100001 deduplicated into geo:GSE100001!) = 3 total datasets
    assert manifest_2["datasets_count"] == 3
    assert manifest_2["l1_deduped"] == 1
    assert manifest_2["rewritten_sample_partitions"] == ["geo"]

    con = duckdb.connect(":memory:")
    cids = [
        r[0]
        for r in con.execute(
            f"SELECT canonical_id FROM read_parquet('{current / 'datasets.parquet'}') ORDER BY canonical_id"
        ).fetchall()
    ]
    assert cids == ["geo:GSE100001", "geo:GSE100002", "pride:PXD000123"]
    con.close()


def test_status_cli_works_without_existing_data(tmp_path: Path) -> None:
    env = {**os.environ, "PHI_OMICS_SNIFFER_STATE_DIR": str(tmp_path / "empty_state")}
    res = subprocess.run(
        [sys.executable, "-B", str(CODE_DIR / "pipeline" / "incremental_parquet_sync.py"), "--status"],
        cwd=str(CODE_DIR),
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )
    assert "geo" in res.stdout
    assert "arrayexpress" in res.stdout
    assert "pride" in res.stdout
