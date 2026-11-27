#!/usr/bin/env python3
"""STOmicsDB (CNGB) dataset catalog -> omics_stomics_staging.db (never omics_v1.db).

Public API: GET https://db.cngb.org/stomics/api/datasets/
Requires JSON-string query params: filters=[] & pagination={page,per_page}.
Stores dataset cards + file_id / portal / FTP pointers only (no raw matrices).
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)
import json
from stage_common import *

DB = ROOT / "data" / "omics_stomics_staging.db"
SRC = "stomics_scan"
API = "https://db.cngb.org/stomics/api/datasets/"
WEB = "https://db.cngb.org/stomics"
FTP_PUB = "https://ftp.cngb.org/pub/stomics/"
FTP_SCIRAID = "https://ftp.cngb.org/pub/SciRAID/stomics/"
log = get_logger(SRC)
thr = Throttle(0.4)

DDL = """
CREATE TABLE IF NOT EXISTS src_stomics_dataset(
  dataset_id TEXT PRIMARY KEY, title TEXT, summary TEXT, status TEXT,
  species TEXT, tissues TEXT, organ_parts TEXT, disease TEXT,
  platforms TEXT, technologies TEXT, sample_number INTEGER, section_number INTEGER,
  cells INTEGER, genes INTEGER, spots INTEGER, submission_date TEXT,
  file_ids TEXT, accessions TEXT, pmids TEXT, portal_url TEXT,
  is_liver INTEGER, liver_tags TEXT, raw_json TEXT, fetched_at TEXT);
"""

def _names(xs):
    out = []
    for x in xs or []:
        if isinstance(x, dict):
            out.append(x.get("name") or x.get("scientific_name") or x.get("common_names", [None])[0] or "")
        else:
            out.append(str(x))
    return [n for n in out if n]

def liver_tags(d):
    blob = " ".join([
        d.get("title") or "", d.get("summary") or "", d.get("overall_design") or "",
        " ".join(_names(d.get("tissues"))), " ".join(_names(d.get("organ_parts"))),
        " ".join(_names(d.get("disease"))), " ".join(_names(d.get("cell_types"))),
    ])
    b = blob.lower()
    t = []
    if re.search(r"\bhcc\b|hepatocellular", b): t.append("HCC")
    if "cholangi" in b or "bile duct" in b or "biliary" in b: t.append("CHOL")
    if "masld" in b or "nafld" in b or "nash" in b: t.append("MASLD")
    if "hepatocyte" in b: t.append("hepatocyte")
    if "hepatoblastoma" in b: t.append("Hepatoblastoma")
    if "gallbladder" in b: t.append("gallbladder")
    if re.search(r"\bliver\b", b) or "hepat" in b: t.append("liver")
    if is_liver(blob) and not t: t.append("liver_related")
    return sorted(set(t))

def fetch_page(page, per_page=200):
    params = {
        "filters": json.dumps([]),
        "pagination": json.dumps({"page": page, "per_page": per_page}),
    }
    return fetch(API, log, thr, params=params, headers={"Accept": "application/json"}).json()

def main():
    c = connect(DB, DDL)
    cur = load_cursor(c, SRC)
    if cur.get("complete"):
        log.info("already complete; exit"); return
    done = set(cur.get("done", []))

    first = fetch_page(1, 200)
    meta = (first.get("meta") or {}).get("pagination") or {}
    total = int(meta.get("total") or 0)
    per_page = 200
    total_pages = max(1, (total + per_page - 1) // per_page)
    log.info("total=%d pages=%d done=%d", total, total_pages, len(done))
    cur.update(total=total, total_pages=total_pages)
    save_cursor(c, SRC, cur)

    for page in range(1, total_pages + 1):
        j = first if page == 1 else fetch_page(page, per_page)
        for d in j.get("data") or []:
            did = d.get("dataset_id")
            if not did or did in done:
                continue
            tags = liver_tags(d)
            species = _names(d.get("species"))
            tissues = _names(d.get("tissues"))
            organs = _names(d.get("organ_parts"))
            diseases = _names(d.get("disease"))
            platforms = _names(d.get("platforms"))
            techs = _names(d.get("stomics_technologies"))
            org = species[0] if len(species) == 1 else (", ".join(species) if species else None)
            tissue = tissues[0] if len(tissues) == 1 else (", ".join(tissues) if tissues else None)
            portal = f"{WEB}/datasets/{did}"
            ds = f"STOMICS:{did}"
            upsert_dataset(c, dict(
                canonical_id=ds, title=d.get("title"), organism=org, tissue=tissue,
                disease_tags=tags or diseases, assay="spatial transcriptomics (STOmicsDB)",
                platforms=platforms or techs, n_samples=d.get("sample_number"),
                access="public" if (d.get("status") or "").lower() in ("public", "released", "") else "unknown",
                confidence="metadata_only", probe_status="ok"))
            add_member(c, ds, "stomics", did, "dataset", portal)
            for acc in d.get("accessions") or []:
                if isinstance(acc, dict) and acc.get("value"):
                    add_member(c, ds, "cnsa" if str(acc["value"]).startswith("CNP") else "external",
                               acc["value"], acc.get("name") or "accession", acc.get("url"))
            add_artifact(c, dict(
                artifact_id=f"STOMICS:{did}:portal", dataset_id=ds, kind="processed",
                role="dataset_portal", title=d.get("title"), format="html",
                uri=portal, access="public", source="stomics"))
            add_artifact(c, dict(
                artifact_id=f"STOMICS:{did}:ftp_pub", dataset_id=ds, kind="raw",
                role="ftp_pointer", title="STOmics public FTP root (pointer)",
                format="html", uri=FTP_PUB, access="public", source="stomics"))
            for fid in d.get("file_ids") or []:
                add_artifact(c, dict(
                    artifact_id=f"STOMICS:{did}:file:{fid}", dataset_id=ds, kind="processed",
                    role="stomics_file_id", title=fid, filename=fid, format="id",
                    uri=f"{portal}#file={fid}", access="public", source="stomics"))
            n = d.get("sample_number") or 0
            if n:
                upsert_sample(c, dict(
                    sample_id=f"STOMICS:{did}:cohort", dataset_id=ds,
                    title=f"{did} sample cohort", organism=org, tissue=tissue,
                    disease=";".join(diseases) if diseases else None,
                    attributes={"sample_number": n, "section_number": d.get("section_number"),
                                "cells": d.get("cells"), "genes": d.get("genes"),
                                "spots": d.get("spots"), "organ_parts": organs,
                                "technologies": techs, "is_liver": int(bool(tags))}))
            c.execute(
                "INSERT OR REPLACE INTO src_stomics_dataset VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (did, d.get("title"), (d.get("summary") or "")[:5000], d.get("status"),
                 json.dumps(species), json.dumps(tissues), json.dumps(organs), json.dumps(diseases),
                 json.dumps(platforms), json.dumps(techs), d.get("sample_number"), d.get("section_number"),
                 d.get("cells"), d.get("genes"), d.get("spots"), str(d.get("submission_date") or ""),
                 json.dumps(d.get("file_ids") or []), json.dumps(d.get("accessions") or [])[:8000],
                 json.dumps(d.get("pmids") or []), portal,
                 int(bool(tags)), json.dumps(tags), json.dumps(d, ensure_ascii=False)[:30000], now_iso()))
            done.add(did)
        cur["done"] = sorted(done)
        cur["last_page"] = page
        save_cursor(c, SRC, cur)
        log.info("page %d/%d stored=%d", page, total_pages, len(done))

    cur["complete"] = True
    cur["ftp_note"] = f"raw FTP trees: {FTP_PUB} and {FTP_SCIRAID} (pointers only)"
    save_cursor(c, SRC, cur)
    nliv = c.execute("SELECT COUNT(*) FROM src_stomics_dataset WHERE is_liver=1").fetchone()[0]
    log.info("DONE datasets=%d liver=%d", len(done), nliv)

if __name__ == "__main__":
    main()
