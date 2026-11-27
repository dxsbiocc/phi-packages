#!/usr/bin/env python3
"""GDSC / CTRP / PRISM drug-sensitivity catalogs -> omics_drug_screen_staging.db.

Metadata-only: cell lines, compounds, screen descriptors, and file pointers to
dose-response matrices. Does NOT download large fitted/raw matrices
(GDSC fitted CSVs ~40-50MB; PRISM LFC/matrices 100MB+).

Sources:
  GDSC  – Sanger FTP current_release (Cell_Lines_Details.xlsx, screened_compounds CSV)
  PRISM – figshare Repurposing Public 24Q2 (cell + treatment meta CSVs)
  CTRP  – catalog pointers only (NCI CTD2 FTP/HTTPS currently unreachable from this box)

Overlap with omics_depmap_staging.db: record DepMap model IDs / CCLE names when present.
Never writes omics_v1.db.
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)
import csv, io, json
from stage_common import *
from http_util import http_get

try:
    import openpyxl
except ImportError:
    openpyxl = None

DB = ROOT / "data" / "omics_drug_screen_staging.db"
SRC = "drug_screen_scan"
GDSC_FTP = "https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release"
PRISM_ARTICLE = 25917643  # Repurposing Public 24Q2
DEPMAP_DB = ROOT / "data" / "omics_depmap_staging.db"
log = get_logger(SRC)
thr = Throttle(0.5)

DDL = """
CREATE TABLE IF NOT EXISTS src_gdsc_cell_line(
  cosmic_id TEXT PRIMARY KEY, sample_name TEXT, ccle_name TEXT, gdsc_name TEXT,
  tcga_desc TEXT, tissue TEXT, tissue_sub TEXT, microsatellite TEXT,
  is_liver INTEGER, liver_tags TEXT, raw_json TEXT, fetched_at TEXT);
CREATE TABLE IF NOT EXISTS src_gdsc_compound(
  drug_id TEXT PRIMARY KEY, drug_name TEXT, synonyms TEXT, targets TEXT,
  target_pathway TEXT, pubchem TEXT, raw_json TEXT, fetched_at TEXT);
CREATE TABLE IF NOT EXISTS src_prism_cell_line(
  depmap_id TEXT PRIMARY KEY, ccle_name TEXT, cell_line_name TEXT, lineage TEXT,
  primary_disease TEXT, subtype TEXT, is_liver INTEGER, liver_tags TEXT,
  raw_json TEXT, fetched_at TEXT);
CREATE TABLE IF NOT EXISTS src_prism_compound(
  broad_id TEXT PRIMARY KEY, name TEXT, dose TEXT, screen_id TEXT,
  moa TEXT, target TEXT, disease_area TEXT, indication TEXT,
  raw_json TEXT, fetched_at TEXT);
CREATE TABLE IF NOT EXISTS src_ctrp_catalog(
  key TEXT PRIMARY KEY, title TEXT, uri TEXT, note TEXT, fetched_at TEXT);
CREATE TABLE IF NOT EXISTS src_drug_file(
  file_id TEXT PRIMARY KEY, source TEXT, release_name TEXT, name TEXT,
  size_bytes INTEGER, url TEXT, role TEXT, downloaded INTEGER, fetched_at TEXT);
"""

GDSC_FILES = [
    ("Cell_Lines_Details.xlsx", "cell_line_annotation", True),
    ("screened_compounds_rel_8.4.csv", "compound_annotation", True),
    ("GDSC1_fitted_dose_response_24Jul22.csv", "fitted_dose_response", False),
    ("GDSC1_fitted_dose_response_24Jul22.xlsx", "fitted_dose_response", False),
    ("GDSC2_fitted_dose_response_24Jul22.csv", "fitted_dose_response", False),
    ("GDSC2_fitted_dose_response_24Jul22.xlsx", "fitted_dose_response", False),
    ("ANOVA_results_GDSC1_24Jul22.xlsx", "anova_markers", False),
    ("ANOVA_results_GDSC2_24Jul22.xlsx", "anova_markers", False),
    ("GDSC_Fitted_Data_Description.pdf", "docs", False),
    ("GDSC1_public_raw_data_24Jul22.csv", "raw_viability", False),
    ("GDSC2_public_raw_data_24Jul22.csv", "raw_viability", False),
]

CTRP_POINTERS = [
    ("CTRPv2_ExpandedDataset_zip",
     "ftp://anonymous:guest@caftpd.nci.nih.gov/pub/OCG-DCC/CTD2/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/CTRPv2.0_2015_ctd2_ExpandedDataset.zip",
     "NCI CTD2 FTP zip (meta+curves); FTP/HTTPS unreachable from this box as of scan"),
    ("CTRPv2_NCI_web_dir",
     "https://ctd2-data.nci.nih.gov/Public/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/",
     "Directory now serves SPA HTML; direct v20.meta.*.txt URLs return HTML not files"),
    ("CTRPv2_Orcestra_PSet",
     "https://zenodo.org/record/7826870/files/PSet_CTRPv2.rds?download=1",
     "Orcestra/Zenodo PharmacoSet RDS (R binary) — pointer only, not parsed"),
    ("CTRPv2_Broad_portal",
     "https://portals.broadinstitute.org/ctrp/",
     "Cancer Therapeutics Response Portal UI"),
]

def liver_tags_blob(*parts):
    blob = " ".join(str(p) for p in parts if p)
    b = blob.lower()
    t = []
    if re.search(r"\bhcc\b|hepatocellular", b): t.append("HCC")
    if "cholangi" in b or "bile duct" in b or "biliary" in b: t.append("CHOL")
    if "hepatoblastoma" in b: t.append("Hepatoblastoma")
    if re.search(r"\bliver\b", b) or "hepat" in b: t.append("liver")
    if is_liver(blob) and not t: t.append("liver_related")
    return sorted(set(t))

def load_depmap_ids():
    """DepMap ModelID -> basic fields for overlap annotation."""
    m = {}
    if not DEPMAP_DB.exists():
        return m
    try:
        dc = sqlite3.connect(DEPMAP_DB, timeout=30)
        for mid, raw in dc.execute("SELECT model_id, raw_json FROM src_depmap_model"):
            try:
                j = json.loads(raw)
            except Exception:
                j = {}
            m[mid] = j
            ccle = j.get("CCLEName") or j.get("CellLineName")
            if ccle:
                m.setdefault(f"ccle:{ccle.upper()}", mid)
        dc.close()
    except Exception as e:
        log.warning("depmap overlap skip: %s", e)
    return m

def register_file(c, source, release, name, size, url, role, downloaded):
    fid = f"{source}:{release}:{name}"
    c.execute("INSERT OR REPLACE INTO src_drug_file VALUES(?,?,?,?,?,?,?,?,?)",
              (fid, source, release, name, size, url, role, int(downloaded), now_iso()))
    ds = f"DRUG:{source}:{release}"
    add_artifact(c, dict(
        artifact_id=f"DRUGFILE:{fid}", dataset_id=ds, kind="processed",
        role=role, title=name, filename=name,
        format=name.rsplit(".", 1)[-1], uri=url, size_bytes=size,
        access="public", source=source.lower()))

def ingest_gdsc(c, depmap):
    release = "GDSC_current_release"
    ds = f"DRUG:GDSC:{release}"
    upsert_dataset(c, dict(
        canonical_id=ds, title="GDSC drug sensitivity (Sanger current_release)",
        organism="Homo sapiens", assay="drug sensitivity (GDSC1/GDSC2)",
        platforms=["GDSC"], access="public", license="GDSC data access terms",
        confidence="metadata_only", probe_status="ok"))
    add_member(c, ds, "gdsc", release, "release", GDSC_FTP + "/")

    # register all file pointers; download only small annotation files
    for name, role, do_fetch in GDSC_FILES:
        url = f"{GDSC_FTP}/{name}"
        size = None
        try:
            hr = http_get(url, timeout=(10, 30), stream=True)
            size = int(hr.headers.get("content-length") or 0) or None
            hr.close()
        except Exception as e:
            log.warning("HEAD-ish %s: %s", name, e)
        register_file(c, "GDSC", release, name, size, url, role, False)

    # Cell lines
    if openpyxl is None:
        raise RuntimeError("openpyxl required for Cell_Lines_Details.xlsx")
    xurl = f"{GDSC_FTP}/Cell_Lines_Details.xlsx"
    log.info("fetching %s", xurl)
    thr()
    thr(); xb = http_get(xurl, timeout=(10, 120)).content
    register_file(c, "GDSC", release, "Cell_Lines_Details.xlsx", len(xb), xurl, "cell_line_annotation", True)
    wb = openpyxl.load_workbook(io.BytesIO(xb), read_only=True, data_only=True)
    # sheet name historically "Cell line details"
    ws = wb[wb.sheetnames[0]]
    rows = list(ws.iter_rows(values_only=True))
    headers = [str(h).strip() if h is not None else f"col{i}" for i, h in enumerate(rows[0])]
    # normalize header keys
    def hget(row, *cands):
        for cand in cands:
            for i, h in enumerate(headers):
                if h.lower().replace("\n", " ") == cand.lower() or cand.lower() in h.lower():
                    return row[i] if i < len(row) else None
        return None
    n_cell = 0; n_liv = 0
    for row in rows[1:]:
        if not row or all(v is None for v in row):
            continue
        cosmic = hget(row, "cosmic identifier", "cosmic_id", "cosmic")
        if cosmic is None:
            continue
        cosmic = str(int(cosmic)) if isinstance(cosmic, float) else str(cosmic).strip()
        sample = hget(row, "sample name", "cell line name", "cell_line")
        ccle = hget(row, "ccle name", "ccle")
        gdsc_n = hget(row, "gdsc name", "gdsc")
        tcga = hget(row, "tcga classification", "tcga desc", "cancer type", "tcga")
        tissue = hget(row, "tissue", "gdsc tissue descriptor 1", "site")
        tissue_sub = hget(row, "tissue sub-type", "gdsc tissue descriptor 2", "histology")
        msi = hget(row, "microsatellite", "msi")
        tags = liver_tags_blob(sample, ccle, tcga, tissue, tissue_sub)
        # also match LIHC / LIHC-like TCGA
        if tcga and str(tcga).upper() in ("LIHC", "CHOL", "LIHC/CHOL"):
            extra = ["liver"]
            tu = str(tcga).upper()
            if "LIHC" in tu: extra.append("HCC")
            if "CHOL" in tu: extra.append("CHOL")
            tags = sorted(set(tags + extra))
        is_liv = int(bool(tags))
        n_cell += 1; n_liv += is_liv
        raw = {headers[i]: row[i] for i in range(min(len(headers), len(row)))}
        c.execute("INSERT OR REPLACE INTO src_gdsc_cell_line VALUES(?,?,?,?,?,?,?,?,?,?,?,?)",
                  (cosmic, sample, ccle, gdsc_n, str(tcga) if tcga is not None else None,
                   str(tissue) if tissue is not None else None,
                   str(tissue_sub) if tissue_sub is not None else None,
                   str(msi) if msi is not None else None,
                   is_liv, json.dumps(tags), json.dumps(raw, default=str)[:8000], now_iso()))
        depmap_id = None
        if ccle:
            depmap_id = depmap.get(f"ccle:{str(ccle).upper()}")
        attrs = {"cosmic_id": cosmic, "ccle_name": ccle, "gdsc_name": gdsc_n,
                 "tcga": tcga, "is_liver": is_liv, "depmap_model_id": depmap_id,
                 "overlap_depmap": bool(depmap_id)}
        upsert_sample(c, dict(
            sample_id=f"GDSC:COSMIC:{cosmic}", dataset_id=ds, title=sample or cosmic,
            organism="Homo sapiens", tissue=str(tissue) if tissue else None,
            disease=str(tcga) if tcga else None, condition=str(tissue_sub) if tissue_sub else None,
            attributes=attrs))
    wb.close()
    c.execute("UPDATE datasets SET n_samples=? WHERE canonical_id=?", (n_cell, ds))

    # Compounds
    curl = f"{GDSC_FTP}/screened_compounds_rel_8.4.csv"
    log.info("fetching %s", curl)
    thr(); ctxt = http_get(curl, timeout=(10, 60)).text
    register_file(c, "GDSC", release, "screened_compounds_rel_8.4.csv", len(ctxt.encode()), curl, "compound_annotation", True)
    reader = csv.DictReader(io.StringIO(ctxt))
    n_cmp = 0
    for row in reader:
        did = str(row.get("DRUG_ID") or row.get("drug_id") or "").strip()
        if not did:
            continue
        n_cmp += 1
        c.execute("INSERT OR REPLACE INTO src_gdsc_compound VALUES(?,?,?,?,?,?,?,?)",
                  (did, row.get("DRUG_NAME") or row.get("drug_name"),
                   row.get("SYNONYMS") or row.get("synonyms"),
                   row.get("TARGET") or row.get("targets") or row.get("PUTATIVE_TARGET"),
                   row.get("TARGET_PATHWAY") or row.get("target_pathway") or row.get("PATHWAY_NAME"),
                   row.get("PUBCHEM") or row.get("pubchem"),
                   json.dumps(row)[:4000], now_iso()))
    log.info("GDSC cell_lines=%d liver=%d compounds=%d", n_cell, n_liv, n_cmp)
    return dict(cells=n_cell, liver=n_liv, compounds=n_cmp, release=release)

def ingest_prism(c, depmap):
    art = fetch(f"https://api.figshare.com/v2/articles/{PRISM_ARTICLE}", log, thr).json()
    release = "PRISM_Repurposing_24Q2"
    ds = f"DRUG:PRISM:{release}"
    upsert_dataset(c, dict(
        canonical_id=ds, title=art.get("title") or "PRISM Repurposing Public 24Q2",
        organism="Homo sapiens", assay="PRISM drug repurposing screen",
        platforms=["PRISM"], access="public", license="CC BY 4.0",
        confidence="metadata_only", probe_status="ok"))
    add_member(c, ds, "figshare", str(PRISM_ARTICLE), "release", art.get("url_public_html"))
    add_member(c, ds, "prism", release, "screen", "https://depmap.org/repurposing/")

    cell_url = treat_url = None
    for f in art.get("files") or []:
        name, size, url = f["name"], f["size"], f["download_url"]
        role = "other"
        download = False
        if "Cell_Line_Meta" in name:
            role, download, cell_url = "cell_line_annotation", True, url
        elif "Treatment_Meta" in name or "Compound_List" in name:
            role = "compound_annotation"
            if "Treatment_Meta" in name:
                download, treat_url = True, url
        elif "LFC" in name or "matrix" in name.lower() or "LMFI" in name:
            role = "dose_response_matrix"
        elif "README" in name.upper() or "Readme" in name:
            role = "docs"
        register_file(c, "PRISM", release, name, size, url, role, False)

    # cell lines (small)
    log.info("fetching PRISM cell meta")
    thr(); ctxt = http_get(cell_url, timeout=(10, 90)).text
    register_file(c, "PRISM", release, "Repurposing_Public_24Q2_Cell_Line_Meta_Data.csv",
                  len(ctxt.encode()), cell_url, "cell_line_annotation", True)
    n_cell = 0; n_liv = 0
    for row in csv.DictReader(io.StringIO(ctxt)):
        # column names vary slightly across releases
        mid = row.get("depmap_id") or row.get("DepMap_ID") or row.get("ModelID") or row.get("depmap_ID")
        if not mid:
            # sometimes first column
            mid = next(iter(row.values()), None)
        if not mid or not str(mid).startswith("ACH-"):
            # try find ACH-
            for v in row.values():
                if v and str(v).startswith("ACH-"):
                    mid = v; break
        if not mid:
            continue
        mid = str(mid).strip()
        ccle = row.get("ccle_name") or row.get("CCLE_Name") or row.get("ccle")
        name = row.get("cell_line_name") or row.get("cell_line") or row.get("stripped_cell_line_name") or ccle
        lineage = row.get("lineage") or row.get("OncotreeLineage") or row.get("lineage_1")
        disease = row.get("primary_disease") or row.get("OncotreePrimaryDisease") or row.get("disease")
        subtype = row.get("subtype") or row.get("OncotreeSubtype") or row.get("lineage_2")
        tags = liver_tags_blob(name, ccle, lineage, disease, subtype)
        is_liv = int(bool(tags))
        n_cell += 1; n_liv += is_liv
        c.execute("INSERT OR REPLACE INTO src_prism_cell_line VALUES(?,?,?,?,?,?,?,?,?,?)",
                  (mid, ccle, name, lineage, disease, subtype, is_liv, json.dumps(tags),
                   json.dumps(row)[:6000], now_iso()))
        overlap = mid in depmap
        upsert_sample(c, dict(
            sample_id=f"PRISM:{mid}", dataset_id=ds, title=name or mid,
            organism="Homo sapiens", tissue=lineage, disease=disease, condition=subtype,
            attributes={"depmap_model_id": mid, "ccle_name": ccle, "is_liver": is_liv,
                        "overlap_depmap": overlap, "source": "PRISM_24Q2"}))
    c.execute("UPDATE datasets SET n_samples=? WHERE canonical_id=?", (n_cell, ds))

    # treatments / compounds (meta only; may be ~1.5MB)
    log.info("fetching PRISM treatment meta")
    thr(); ttxt = http_get(treat_url, timeout=(10, 180)).text
    register_file(c, "PRISM", release, "Repurposing_Public_24Q2_Treatment_Meta_Data.csv",
                  len(ttxt.encode()), treat_url, "compound_annotation", True)
    n_cmp = 0
    seen = set()
    for row in csv.DictReader(io.StringIO(ttxt)):
        broad = row.get("broad_id") or row.get("BROAD_ID") or row.get("broad_cpd_id") or row.get("drug_id")
        name = row.get("name") or row.get("drug_name") or row.get("compound_name") or row.get("pert_iname")
        key = str(broad or name or "").strip()
        if not key or key in seen:
            # still count unique compound names
            if not key:
                continue
        seen.add(key)
        n_cmp += 1
        c.execute("INSERT OR REPLACE INTO src_prism_compound VALUES(?,?,?,?,?,?,?,?,?,?)",
                  (key, name, row.get("dose") or row.get("dose_unit"),
                   row.get("screen_id") or row.get("screen"),
                   row.get("moa") or row.get("MOA"),
                   row.get("target") or row.get("Target"),
                   row.get("disease_area") or row.get("disease.area"),
                   row.get("indication"),
                   json.dumps(row)[:4000], now_iso()))
    log.info("PRISM cell_lines=%d liver=%d compounds_rows=%d", n_cell, n_liv, n_cmp)
    return dict(cells=n_cell, liver=n_liv, compounds=n_cmp, release=release, article=PRISM_ARTICLE)

def ingest_ctrp(c):
    release = "CTRPv2_2015"
    ds = f"DRUG:CTRP:{release}"
    upsert_dataset(c, dict(
        canonical_id=ds, title="CTRPv2 (CTD2) drug sensitivity — catalog pointers",
        organism="Homo sapiens", assay="drug sensitivity (CTRP)",
        platforms=["CTRP", "CTD2"], access="public",
        confidence="metadata_only", probe_status="blocked_download"))
    add_member(c, ds, "ctrp", release, "screen", "https://portals.broadinstitute.org/ctrp/")
    for key, uri, note in CTRP_POINTERS:
        c.execute("INSERT OR REPLACE INTO src_ctrp_catalog VALUES(?,?,?,?,?)",
                  (key, key, uri, note, now_iso()))
        add_artifact(c, dict(
            artifact_id=f"CTRP:{key}", dataset_id=ds, kind="processed",
            role="catalog_pointer", title=key, format="url",
            uri=uri, access="public", source="ctrp"))
        register_file(c, "CTRP", release, key, None, uri, "catalog_pointer", False)
    log.info("CTRP pointers only (NCI FTP/HTTPS blocked)")
    return dict(blocker="NCI CTD2 FTP/HTTPS unreachable; Orcestra RDS pointer recorded", release=release)

def main():
    c = connect(DB, DDL)
    cur = load_cursor(c, SRC)
    if cur.get("complete"):
        log.info("already complete; exit"); return

    depmap = load_depmap_ids()
    log.info("depmap overlap index keys=%d", len(depmap))

    if not cur.get("gdsc_done"):
        cur["gdsc"] = ingest_gdsc(c, depmap)
        cur["gdsc_done"] = True
        save_cursor(c, SRC, cur)
    if not cur.get("prism_done"):
        cur["prism"] = ingest_prism(c, depmap)
        cur["prism_done"] = True
        save_cursor(c, SRC, cur)
    if not cur.get("ctrp_done"):
        cur["ctrp"] = ingest_ctrp(c)
        cur["ctrp_done"] = True
        save_cursor(c, SRC, cur)

    cur["complete"] = True
    cur["note"] = ("Large dose-response matrices pointer-only. "
                   "PRISM absent from omics_depmap_staging 24Q4 — recorded here from figshare 24Q2. "
                   "CTRP meta files not fetched (NCI host blocked).")
    save_cursor(c, SRC, cur)
    log.info("DONE gdsc=%s prism=%s ctrp=%s", cur.get("gdsc"), cur.get("prism"), cur.get("ctrp"))

if __name__ == "__main__":
    main()
