#!/usr/bin/env python3
"""ImmPort Shared Data study catalog -> omics_immport_staging.db (never omics_v1.db).

Uses the public search API (no token):
  GET https://www.immport.org/data/query/api/search/study

Study-detail endpoints under /data/query/api/study/{SDY*} return 401 without
an ImmPort API bearer token, so we store catalog cards + portal / download-browser
pointers only (no biosample-level rows, no raw file pulls).
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)
import json
from stage_common import *

DB = ROOT / "data" / "omics_immport_staging.db"
SRC = "immport_scan"
API = "https://www.immport.org/data/query/api/search/study"
DETAIL_API = "https://www.immport.org/data/query/api/study"  # auth required
WEB = "https://www.immport.org/shared/study"
DOWNLOAD_BROWSER = "https://www.immport.org/shared/study/all/download"
PAGE_SIZE = 200
SOURCE_FIELDS = ",".join([
    "study_accession", "brief_title", "brief_description",
    "condition_or_disease", "assay_method_count", "shared_subject_count",
    "research_focus", "species", "program_name", "actual_enrollment",
    "pubmed_id", "study_pi", "doi", "clinicaltrials_link", "biosample_type",
    "clinical_trial", "latest_data_release_version", "initial_data_release_date",
    "gender_included", "age_range",
])
log = get_logger(SRC)
thr = Throttle(0.35)

DDL = """
CREATE TABLE IF NOT EXISTS src_immport_study(
  study_accession TEXT PRIMARY KEY,
  brief_title TEXT,
  brief_description TEXT,
  condition_or_disease TEXT,
  research_focus TEXT,
  program_name TEXT,
  species TEXT,
  biosample_type TEXT,
  assay_method_count TEXT,
  actual_enrollment INTEGER,
  shared_subject_count TEXT,
  clinical_trial TEXT,
  clinicaltrials_link TEXT,
  pubmed_id TEXT,
  study_pi TEXT,
  doi TEXT,
  age_range TEXT,
  gender_included TEXT,
  latest_data_release_version TEXT,
  initial_data_release_date TEXT,
  portal_url TEXT,
  is_liver INTEGER,
  liver_tags TEXT,
  raw_json TEXT,
  fetched_at TEXT
);
CREATE INDEX IF NOT EXISTS idx_src_immport_liver ON src_immport_study(is_liver);
"""


def _as_list(v):
    if v is None:
        return []
    if isinstance(v, list):
        return [str(x) for x in v if x is not None and str(x).strip()]
    return [str(v)] if str(v).strip() else []


def liver_tags(src: dict) -> list[str]:
    blob = " ".join([
        src.get("brief_title") or "",
        src.get("brief_description") or "",
        " ".join(_as_list(src.get("condition_or_disease"))),
        " ".join(_as_list(src.get("research_focus"))),
        " ".join(_as_list(src.get("biosample_type"))),
        " ".join(_as_list(src.get("program_name"))),
        src.get("study_accession") or "",
    ])
    b = blob.lower()
    t: list[str] = []
    if re.search(r"\bhcc\b|hepatocellular", b) or "liver cancer" in b:
        t.append("HCC")
    if "cholangi" in b or "bile duct" in b or re.search(r"\bbiliary\b", b):
        t.append("CHOL")
    if "masld" in b or "nafld" in b or "nash" in b:
        t.append("MASLD")
    if "cirrhosis" in b or "cirrhotic" in b:
        t.append("cirrhosis")
    if "hepatitis" in b or re.search(r"\bhbv\b|\bhcv\b", b):
        t.append("hepatitis")
    if "liver transplantation" in b or "liver transplant" in b or (
        "transplant" in b and "liver" in b
    ):
        t.append("transplant-liver")
    if re.search(r"\bliver\b", b) or "hepatic" in b or "hepat" in b:
        t.append("liver")
    if is_liver(blob) and not t:
        t.append("liver_related")
    return sorted(set(t))


def search_page(from_record: int, page_size: int = PAGE_SIZE) -> dict:
    return fetch(
        API, log, thr,
        params={
            "pageSize": page_size,
            "fromRecord": from_record,
            "sourceFields": SOURCE_FIELDS,
            "sortField": "study_accession",
            "sortFieldDirection": "asc",
        },
        headers={"Accept": "application/json"},
    ).json()


def assay_summary(assay_method_count) -> str | None:
    items = _as_list(assay_method_count)
    if not items:
        return None
    names = []
    for it in items:
        names.append(it.split("(")[0].strip() if "(" in it else it)
    return "; ".join(names[:12])


def store_hit(c, h, done: set):
    src = h.get("_source") or {}
    acc = src.get("study_accession") or h.get("_id")
    if not acc or acc in done:
        return False
    tags = liver_tags(src)
    species = _as_list(src.get("species"))
    conditions = _as_list(src.get("condition_or_disease"))
    focuses = _as_list(src.get("research_focus"))
    biosamples = _as_list(src.get("biosample_type"))
    programs = _as_list(src.get("program_name"))
    assays = _as_list(src.get("assay_method_count"))
    org = species[0] if len(species) == 1 else (", ".join(species) if species else None)
    liverish = any(
        t in ("liver", "transplant-liver", "HCC", "CHOL", "cirrhosis", "MASLD", "hepatitis")
        for t in tags
    )
    if liverish and ("liver" in tags or "transplant-liver" in tags):
        tissue = "liver"
    elif biosamples:
        tissue = biosamples[0] if len(biosamples) == 1 else ", ".join(biosamples[:3])
    else:
        tissue = None
    portal = f"{WEB}/{acc}"
    enroll = src.get("actual_enrollment")
    try:
        enroll_i = int(enroll) if enroll is not None and str(enroll).strip() != "" else None
    except (TypeError, ValueError):
        enroll_i = None
    ds = f"IMMPORT:{acc}"
    upsert_dataset(c, dict(
        canonical_id=ds,
        title=src.get("brief_title"),
        organism=org,
        tissue=tissue,
        disease_tags=tags or conditions,
        assay=assay_summary(assays) or "ImmPort study (mixed assays)",
        platforms=assays,
        n_samples=enroll_i,
        access="public",
        confidence="metadata_only",
        probe_status="ok",
    ))
    add_member(c, ds, "immport", acc, "study", portal)
    for pmid in _as_list(src.get("pubmed_id")):
        add_member(c, ds, "pubmed", str(pmid), "publication",
                   f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/")
    doi = src.get("doi")
    if doi:
        add_member(c, ds, "doi", str(doi), "doi", f"https://doi.org/{doi}")
    for ct in _as_list(src.get("clinicaltrials_link")):
        add_artifact(c, dict(
            artifact_id=f"IMMPORT:{acc}:ctg:{abs(hash(ct)) % 10_000_000}",
            dataset_id=ds, kind="supplementary", role="clinicaltrials_link",
            title="ClinicalTrials.gov link", format="html", uri=ct,
            access="public", source="immport"))
    add_artifact(c, dict(
        artifact_id=f"IMMPORT:{acc}:portal", dataset_id=ds, kind="processed",
        role="study_portal", title=src.get("brief_title"), format="html",
        uri=portal, access="public", source="immport"))
    add_artifact(c, dict(
        artifact_id=f"IMMPORT:{acc}:search_api", dataset_id=ds, kind="processed",
        role="api_search_record", title="ImmPort public search API record (pointer)",
        format="json",
        uri=f"{API}?studyAccession={acc}&pageSize=1",
        access="public", source="immport"))
    add_artifact(c, dict(
        artifact_id=f"IMMPORT:{acc}:download_browser", dataset_id=ds, kind="raw",
        role="download_browser_pointer",
        title="ImmPort Shared Data download browser (all studies; pointer only)",
        format="html", uri=DOWNLOAD_BROWSER, access="public", source="immport"))
    add_artifact(c, dict(
        artifact_id=f"IMMPORT:{acc}:detail_api_auth", dataset_id=ds, kind="processed",
        role="auth_gated_detail_api",
        title="Study detail API (401 without ImmPort token)",
        format="json", uri=f"{DETAIL_API}/{acc}", access="controlled", source="immport"))
    if enroll_i:
        upsert_sample(c, dict(
            sample_id=f"IMMPORT:{acc}:subjects",
            dataset_id=ds,
            title=f"{acc} enrolled subject cohort",
            organism=org,
            tissue=tissue,
            disease="; ".join(conditions) if conditions else None,
            condition="; ".join(focuses) if focuses else None,
            sex=src.get("gender_included"),
            age=src.get("age_range"),
            attributes={
                "actual_enrollment": enroll_i,
                "biosample_type": biosamples,
                "assay_method_count": assays,
                "program_name": programs,
                "is_liver": int(bool(tags)),
                "liver_tags": tags,
                "note": "cohort row from search catalog; per-subject detail needs auth",
            },
        ))
    c.execute(
        "INSERT OR REPLACE INTO src_immport_study VALUES("
        "?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
        (
            acc,
            src.get("brief_title"),
            (src.get("brief_description") or "")[:8000],
            json.dumps(conditions, ensure_ascii=False),
            json.dumps(focuses, ensure_ascii=False),
            json.dumps(programs, ensure_ascii=False),
            json.dumps(species, ensure_ascii=False),
            json.dumps(biosamples, ensure_ascii=False),
            json.dumps(assays, ensure_ascii=False),
            enroll_i,
            json.dumps(_as_list(src.get("shared_subject_count")), ensure_ascii=False),
            src.get("clinical_trial"),
            json.dumps(_as_list(src.get("clinicaltrials_link")), ensure_ascii=False),
            json.dumps(_as_list(src.get("pubmed_id")), ensure_ascii=False),
            json.dumps(_as_list(src.get("study_pi")), ensure_ascii=False),
            src.get("doi"),
            src.get("age_range"),
            src.get("gender_included"),
            src.get("latest_data_release_version"),
            src.get("initial_data_release_date"),
            portal,
            int(bool(tags)),
            json.dumps(tags, ensure_ascii=False),
            json.dumps(src, ensure_ascii=False)[:40000],
            now_iso(),
        ),
    )
    done.add(acc)
    return True


def main():
    c = connect(DB, DDL)
    cur = load_cursor(c, SRC)
    if cur.get("complete"):
        log.info("already complete; exit")
        return
    done = set(cur.get("done", []))
    start_offset = int(cur.get("from_record", 0))

    probe = search_page(0, page_size=1)
    total = int(((probe.get("hits") or {}).get("total") or {}).get("value") or 0)
    log.info("total_studies=%d start_offset=%d done=%d", total, start_offset, len(done))
    cur["total"] = total
    cur["detail_api_auth"] = (
        f"{DETAIL_API}/{{SDY}} returns 401 without ImmPort bearer token; "
        "catalog-only ingest via public search API"
    )
    save_cursor(c, SRC, cur)

    for offset in range(start_offset, max(total, 1), PAGE_SIZE):
        j = search_page(offset)
        hits = ((j.get("hits") or {}).get("hits")) or []
        n_new = 0
        for h in hits:
            if store_hit(c, h, done):
                n_new += 1
        cur["done"] = sorted(done)
        cur["from_record"] = offset + PAGE_SIZE
        cur["last_offset"] = offset
        save_cursor(c, SRC, cur)
        log.info("offset=%d/%d stored=%d new=%d page_hits=%d",
                 offset, total, len(done), n_new, len(hits))
        if not hits:
            break

    cur["complete"] = True
    cur["note"] = (
        "Public catalog via search API only. "
        f"Detail/biosample/file endpoints at {DETAIL_API}/{{SDY}} need ImmPort API token. "
        f"Shared download browser pointer: {DOWNLOAD_BROWSER}"
    )
    save_cursor(c, SRC, cur)
    nliv = c.execute("SELECT COUNT(*) FROM src_immport_study WHERE is_liver=1").fetchone()[0]
    nsubj = c.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
    log.info("DONE studies=%d liver=%d cohort_samples=%d", len(done), nliv, nsubj)


if __name__ == "__main__":
    main()
