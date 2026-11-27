#!/usr/bin/env python3
"""Broad Single Cell Portal study catalog -> omics_scp_staging.db (never omics_v1.db).

Uses public search API only. GET /studies/{id} returns 401 without auth, so
file/cluster pointers are limited to study page URLs from search hits.
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)
import json
from stage_common import *

DB = ROOT / "data" / "omics_scp_staging.db"
SRC = "scp_scan"
API = "https://singlecell.broadinstitute.org/single_cell/api/v1"
WEB = "https://singlecell.broadinstitute.org"
log = get_logger(SRC)
thr = Throttle(0.35)

DDL = """
CREATE TABLE IF NOT EXISTS src_scp_study(
  accession TEXT PRIMARY KEY, name TEXT, description TEXT, study_url TEXT,
  cell_count INTEGER, gene_count INTEGER, public INTEGER, study_source TEXT,
  species TEXT, organs TEXT, diseases TEXT, sex TEXT, protocols TEXT,
  is_liver INTEGER, liver_tags TEXT, raw_json TEXT, fetched_at TEXT);
"""

def tags_for(study):
    meta = study.get("metadata") or {}
    organs = meta.get("organ") or []
    diseases = meta.get("disease") or []
    blob = " ".join([
        study.get("accession") or "", study.get("name") or "",
        study.get("description") or "", " ".join(organs), " ".join(diseases),
    ])
    t = []
    b = blob.lower()
    if any(k in b for k in ("hepatocellular", " hcc", "hcc ", "liver cancer")):
        t.append("HCC")
    if any(k in b for k in ("cholangi", "bile duct", "biliary")):
        t.append("CHOL")
    if "hepatocyte" in b:
        t.append("hepatocyte")
    if "hepatoblastoma" in b:
        t.append("Hepatoblastoma")
    if any(o.lower() == "liver" for o in organs) or re.search(r"\bliver\b", b):
        t.append("liver")
    if is_liver(blob) and not t:
        t.append("liver_related")
    return sorted(set(t)), organs, diseases

def search_page(page):
    # API caps page size at 10 regardless of size=
    return fetch(
        f"{API}/search", log, thr,
        params={"type": "study", "page": page, "size": 10},
        headers={"Accept": "application/json"},
    ).json()

def main():
    c = connect(DB, DDL)
    cur = load_cursor(c, SRC)
    if cur.get("complete"):
        log.info("already complete; exit"); return
    done = set(cur.get("done", []))
    first = search_page(1)
    total_pages = int(first.get("total_pages") or 1)
    total_studies = int(first.get("total_studies") or 0)
    log.info("total_studies=%d total_pages=%d done=%d", total_studies, total_pages, len(done))
    cur["total_studies"] = total_studies
    save_cursor(c, SRC, cur)

    pages = list(range(1, total_pages + 1))
    # resume: skip fully-done pages by checking accessions; still iterate
    for page in pages:
        j = first if page == 1 else search_page(page)
        studies = j.get("studies") or []
        for s in studies:
            acc = s.get("accession")
            if not acc or acc in done:
                continue
            lt, organs, diseases = tags_for(s)
            meta = s.get("metadata") or {}
            species = meta.get("species") or []
            sex = meta.get("sex") or []
            protocols = meta.get("library_preparation_protocol") or []
            org = species[0] if len(species) == 1 else (", ".join(species) if species else None)
            tissue = organs[0] if len(organs) == 1 else (", ".join(organs) if organs else None)
            path = s.get("study_url") or f"/single_cell/study/{acc}"
            if not path.startswith("http"):
                url = WEB + path
            else:
                url = path
            ds = f"SCP:{acc}"
            upsert_dataset(c, dict(
                canonical_id=ds, title=s.get("name"), organism=org, tissue=tissue,
                disease_tags=lt or diseases, assay="single-cell (Broad SCP)",
                platforms=protocols, n_samples=s.get("cell_count"),
                access="public" if s.get("public") else "unknown",
                confidence="metadata_only", probe_status="ok"))
            add_member(c, ds, "scp", acc, "study", url)
            # pointer only — detail API is 401 anonymously
            add_artifact(c, dict(
                artifact_id=f"SCP:{acc}:portal", dataset_id=ds, kind="processed",
                role="study_portal", title=s.get("name"), format="html",
                uri=url, access="public", source="scp"))
            add_artifact(c, dict(
                artifact_id=f"SCP:{acc}:api_search", dataset_id=ds, kind="processed",
                role="api_search_record", title="SCP search API record (pointer)",
                format="json", uri=f"{API}/search?type=study&terms={acc}",
                access="public", source="scp"))
            # one pseudo-sample summarizing study cell cohort
            if s.get("cell_count"):
                upsert_sample(c, dict(
                    sample_id=f"SCP:{acc}:cells", dataset_id=ds,
                    title=f"{acc} cell cohort", organism=org, tissue=tissue,
                    disease=";".join(diseases) if diseases else None,
                    attributes={"cell_count": s.get("cell_count"), "gene_count": s.get("gene_count"),
                                "sex": sex, "protocols": protocols}))
            c.execute(
                "INSERT OR REPLACE INTO src_scp_study VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (acc, s.get("name"), (s.get("description") or "")[:5000], url,
                 s.get("cell_count"), s.get("gene_count"), int(bool(s.get("public"))),
                 s.get("study_source"), json.dumps(species), json.dumps(organs),
                 json.dumps(diseases), json.dumps(sex), json.dumps(protocols),
                 int(bool(lt)), json.dumps(lt), json.dumps(s)[:20000], now_iso()))
            done.add(acc)
        cur["done"] = sorted(done)
        cur["last_page"] = page
        save_cursor(c, SRC, cur)
        if page % 10 == 0 or page == total_pages:
            log.info("page %d/%d studies_stored=%d", page, total_pages, len(done))

    cur["complete"] = True
    cur["note"] = "GET /studies/{accession} returns 401 without auth; file/cluster metadata limited to portal URLs"
    save_cursor(c, SRC, cur)
    log.info("DONE studies=%d", len(done))

if __name__ == "__main__":
    main()
