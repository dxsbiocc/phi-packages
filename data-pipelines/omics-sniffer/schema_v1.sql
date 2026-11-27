-- omics-sniffer schema v1
-- GEO / SRA / ArrayExpress curation database
--
-- IMPORTANT: Application code MUST run: PRAGMA foreign_keys = ON;
-- SQLite does not enforce FK constraints unless this pragma is enabled per connection.
--
-- Enum-like columns are TEXT (e.g. access, confidence, kind, probe_status).
-- JSON-valued columns are stored as TEXT; validate in application layer.
-- Timestamps are ISO-8601 TEXT (UTC+8 / local app convention).

------------------------------------------------------------------------
-- Unified business tables
------------------------------------------------------------------------

CREATE TABLE datasets (
    canonical_id   TEXT PRIMARY KEY,
    title          TEXT,
    title_norm     TEXT,
    organism       TEXT,
    tissue         TEXT,
    disease_tags   TEXT,   -- JSON array, e.g. ["HCC","liver cancer"]
    assay          TEXT,
    platforms      TEXT,   -- JSON array of platform IDs/names
    n_samples      INTEGER,
    access         TEXT,   -- public | controlled | unknown
    license        TEXT,
    confidence     TEXT,   -- verified_live | metadata_only | paper_claimed
    content_hash   TEXT,
    first_seen_at  TEXT,
    updated_at     TEXT,
    probe_status   TEXT,
    last_checked   TEXT
);

CREATE INDEX idx_datasets_organism_assay ON datasets(organism, assay);
CREATE INDEX idx_datasets_title_norm ON datasets(title_norm);

CREATE TABLE samples (
    sample_id      TEXT PRIMARY KEY,
    dataset_id     TEXT NOT NULL REFERENCES datasets(canonical_id),
    title          TEXT,
    organism       TEXT,
    tissue         TEXT,
    disease        TEXT,
    condition      TEXT,
    sex            TEXT,
    age            TEXT,
    attributes     TEXT,   -- JSON object of raw key-value attributes
    run_ids        TEXT,   -- JSON array of SRA Run accessions
    content_hash   TEXT,
    updated_at     TEXT
);

CREATE INDEX idx_samples_dataset_id ON samples(dataset_id);

CREATE TABLE dataset_artifacts (
    artifact_id    TEXT PRIMARY KEY,
    dataset_id     TEXT NOT NULL REFERENCES datasets(canonical_id),
    sample_id      TEXT REFERENCES samples(sample_id),  -- NULL = dataset-level
    kind           TEXT,   -- raw | processed | supplementary | paper_supplement
    role           TEXT,   -- counts | tpm | fpkm | clinical_table | phenotype | ...
    title          TEXT,
    filename       TEXT,
    format         TEXT,
    uri            TEXT,
    size_bytes     INTEGER,
    access         TEXT,   -- may differ from dataset.access
    license        TEXT,
    source         TEXT,   -- geo | sra | arrayexpress | paper | ...
    derived_from   TEXT,
    probe_status   TEXT,
    last_checked   TEXT,
    content_hash   TEXT
);

CREATE INDEX idx_artifacts_dataset_kind_access
    ON dataset_artifacts(dataset_id, kind, access);
CREATE INDEX idx_artifacts_dataset_id ON dataset_artifacts(dataset_id);

CREATE TABLE dataset_members (
    id             INTEGER PRIMARY KEY AUTOINCREMENT,
    dataset_id     TEXT NOT NULL REFERENCES datasets(canonical_id),
    source         TEXT NOT NULL,   -- geo | sra | arrayexpress | ...
    accession      TEXT NOT NULL,
    role           TEXT,            -- series | bioproject | experiment | entry | ...
    external_url   TEXT,
    UNIQUE(source, accession)
);

CREATE INDEX idx_members_dataset_id ON dataset_members(dataset_id);
CREATE INDEX idx_members_accession ON dataset_members(accession);

------------------------------------------------------------------------
-- Publications
-- PK is pmid TEXT. For DOI-only rows (no PubMed ID), store a synthetic
-- primary key of the form 'DOI:<doi>', e.g. 'DOI:10.1038/s41586-020-0000'.
-- Real PMIDs are numeric strings without prefix. Application must enforce
-- uniqueness of the bare DOI when using the synthetic form.
------------------------------------------------------------------------

CREATE TABLE publications (
    pmid           TEXT PRIMARY KEY,  -- real PMID, or 'DOI:<doi>' if no PMID
    doi            TEXT UNIQUE,
    title          TEXT,
    authors        TEXT,
    journal        TEXT,
    year           INTEGER,
    abstract       TEXT,
    url            TEXT,
    fetched_at     TEXT
);

CREATE TABLE dataset_publications (
    dataset_id     TEXT NOT NULL REFERENCES datasets(canonical_id),
    pmid           TEXT NOT NULL REFERENCES publications(pmid),
    relation       TEXT,   -- primary | secondary | cited_by
    source         TEXT,   -- geo_link | paper_claim | manual | ...
    PRIMARY KEY (dataset_id, pmid)
);

------------------------------------------------------------------------
-- Source mirrors (raw pull snapshots; not exposed via MCP)
------------------------------------------------------------------------

CREATE TABLE src_geo_series (
    accession              TEXT PRIMARY KEY,  -- GSE...
    title                  TEXT,
    summary                TEXT,
    organism               TEXT,
    platform               TEXT,
    sample_count           INTEGER,
    supplementary_files    TEXT,   -- JSON
    raw_xml                TEXT,
    fetched_at             TEXT,
    updatedate             TEXT
);

CREATE INDEX idx_src_geo_fetched_at ON src_geo_series(fetched_at);

CREATE TABLE src_sra_experiment (
    accession              TEXT PRIMARY KEY,  -- SRX...
    bioproject             TEXT,
    biosample              TEXT,
    title                  TEXT,
    platform               TEXT,
    spots                  INTEGER,
    bases                  INTEGER,
    run_accessions         TEXT,   -- JSON array
    ftp_paths              TEXT,   -- JSON array
    raw_xml                TEXT,
    fetched_at             TEXT
);

CREATE INDEX idx_src_sra_fetched_at ON src_sra_experiment(fetched_at);

CREATE TABLE src_arrayexpress_entry (
    accession              TEXT PRIMARY KEY,  -- E-MTAB-...
    title                  TEXT,
    organism               TEXT,
    assay                  TEXT,
    sample_count           INTEGER,
    related_geo            TEXT,
    related_sra            TEXT,
    files                  TEXT,   -- JSON
    raw_json               TEXT,
    fetched_at             TEXT
);

CREATE INDEX idx_src_ae_fetched_at ON src_arrayexpress_entry(fetched_at);

------------------------------------------------------------------------
-- Incremental sync cursors
------------------------------------------------------------------------

CREATE TABLE sync_cursors (
    source         TEXT PRIMARY KEY,  -- geo | sra | arrayexpress | pubmed | ...
    cursor_value   TEXT,
    cursor_kind    TEXT,              -- pdat | updatedate | timestamp | ...
    last_run_at    TEXT,
    last_error     TEXT
);
