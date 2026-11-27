#!/usr/bin/env python3
"""GWAS Catalog study catalog -> omics_gwas_catalog_staging.db (never omics_v1.db).

Primary ingest: official FTP release TSVs under
  https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/
which match REST API v2 /studies cardinality for published studies (~232k) and
include FULL SUMMARY STATISTICS + SUMMARY STATS LOCATION columns.

Legacy REST /rest/api/studies is 410 Gone. REST v2 works but paginating 232k
studies at ~20–25s/page is impractical through this proxy; FTP bulk is the
supported file-download path (same Catalog). Metadata + sumstats URI pointers
only — no association dump, no summary-stat file downloads.
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)
import csv
import json
import re
from stage_common import *
from http_util import http_get as _http_get

DB = ROOT / "data" / "omics_gwas_catalog_staging.db"
SRC = "gwas_catalog_scan"
FTP_LATEST = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest"
PUB_TSV = "gwas-catalog-download-studies-v1.0.3.1.tsv"
UNPUB_TSV = "gwas-catalog-unpublished-studies-v1.0.3.1.tsv"
CACHE = ROOT / "cache_gwas"
PORTAL = "https://www.ebi.ac.uk/gwas/studies"
API_STUDY = "https://www.ebi.ac.uk/gwas/rest/api/v2/studies"
log = get_logger(SRC)
thr = Throttle(0.2)

DDL = """
CREATE TABLE IF NOT EXISTS src_gwas_study(
  accession_id TEXT PRIMARY KEY,
  disease_trait TEXT,
  study_title TEXT,
  efo_traits TEXT,
  bg_efo_traits TEXT,
  pubmed_id TEXT,
  initial_sample_size TEXT,
  replication_sample_size TEXT,
  association_count INTEGER,
  genotyping_technologies TEXT,
  platforms TEXT,
  cohort TEXT,
  full_summary_stats_available INTEGER,
  full_summary_stats TEXT,
  terms_of_license TEXT,
  catalog_status TEXT,
  portal_url TEXT,
  is_liver INTEGER,
  liver_tags TEXT,
  raw_json TEXT,
  fetched_at TEXT
);
CREATE INDEX IF NOT EXISTS idx_src_gwas_liver ON src_gwas_study(is_liver);
CREATE INDEX IF NOT EXISTS idx_src_gwas_sumstats ON src_gwas_study(full_summary_stats_available);
CREATE INDEX IF NOT EXISTS idx_src_gwas_status ON src_gwas_study(catalog_status);
"""


def liver_tags_for(disease: str, mapped: str, bg: str, title: str, acc: str) -> list[str]:
    blob = " ".join([disease or "", mapped or "", bg or "", title or "", acc or ""])
    b = blob.lower()
    t: list[str] = []
    if re.search(r"\bhcc\b|hepatocellular", b) or "liver cancer" in b:
        t.append("HCC")
    if "cholangi" in b or "bile duct" in b or re.search(r"\bbiliary\b", b) or "gallbladder" in b:
        t.append("CHOL")
    if "masld" in b or "nafld" in b or "nash" in b or "fatty liver" in b:
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


def parse_n_samples(text) -> int | None:
    if not text or str(text).strip().upper() in ("NA", "N/A", "", "NOT YET CURATED"):
        return None
    nums = [int(x.replace(",", "")) for x in re.findall(r"\d[\d,]*", str(text))]
    return sum(nums) if nums else None


def norm_ss_uri(uri: str | None) -> str | None:
    if not uri or not str(uri).strip() or str(uri).strip().upper() in ("NA", "N/A"):
        return None
    u = str(uri).strip()
    if u.startswith("ftp://ftp."):
        u = "https://" + u[len("ftp://"):]
    elif u.startswith("http://ftp."):
        u = "https://" + u[len("http://"):]
    return u


def yesish(v) -> bool:
    return str(v or "").strip().lower() in ("yes", "y", "true", "1")


def efo_from_mapped(mapped: str, uri: str) -> list[dict]:
    names = [x.strip() for x in (mapped or "").split(",") if x.strip()]
    uris = [x.strip() for x in (uri or "").split(",") if x.strip()]
    out = []
    for i, name in enumerate(names):
        efo_id = None
        if i < len(uris):
            m = re.search(r"/(EFO_[0-9]+|MONDO_[0-9]+|HP_[0-9]+|Orphanet_[0-9]+|GO_[0-9]+|CHEBI_[0-9]+)", uris[i])
            if m:
                efo_id = m.group(1)
        out.append({"efo_trait": name, "efo_id": efo_id, "uri": uris[i] if i < len(uris) else None})
    return out


def download_tsv(name: str) -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    dest = CACHE / name
    if dest.exists() and dest.stat().st_size > 1_000_000:
        log.info("cache hit %s (%d bytes)", dest.name, dest.stat().st_size)
        return dest
    url = f"{FTP_LATEST}/{name}"
    log.info("downloading %s", url)
    thr()
    for a in range(1, MAX_RETRIES + 1):
        try:
            r = _http_get(url, timeout=(15, 300))
            r.raise_for_status()
            dest.write_bytes(r.content)
            log.info("saved %s (%d bytes)", dest, len(r.content))
            return dest
        except Exception as e:
            log.warning("download retry %d/%d %s: %s", a, MAX_RETRIES, url, e)
            time.sleep(min(60, 2 ** a))
    raise RuntimeError(f"failed download {url}")


def row_get(row: dict, *keys):
    for k in keys:
        if k in row and row[k] is not None:
            return row[k]
    # case-insensitive fallback
    lower = {k.lower(): v for k, v in row.items()}
    for k in keys:
        if k.lower() in lower:
            return lower[k.lower()]
    return None


def store_row(c, row: dict, catalog_status: str) -> bool:
    acc = (row_get(row, "STUDY ACCESSION") or "").strip()
    if not acc or acc.upper() == "NA":
        return False
    disease = row_get(row, "DISEASE/TRAIT") or ""
    title = row_get(row, "STUDY") or ""
    mapped = row_get(row, "MAPPED_TRAIT") or ""
    mapped_uri = row_get(row, "MAPPED_TRAIT_URI") or ""
    bg = row_get(row, "MAPPED BACKGROUND TRAIT", "BACKGROUND TRAIT") or ""
    bg_uri = row_get(row, "MAPPED BACKGROUND TRAIT URI") or ""
    tags = liver_tags_for(disease, mapped, bg, title, acc)
    efo = efo_from_mapped(mapped, mapped_uri)
    bg_efo = efo_from_mapped(bg, bg_uri) if bg else []
    pmid = row_get(row, "PUBMED ID", "PUBMEDID")
    if pmid and str(pmid).strip().lower() in ("na", "not yet curated", ""):
        pmid = None
    geno = row_get(row, "GENOTYPING TECHNOLOGY") or ""
    platform = row_get(row, "PLATFORM [SNPS PASSING QC]") or ""
    has_ss = yesish(row_get(row, "FULL SUMMARY STATISTICS"))
    ss_uri = norm_ss_uri(row_get(row, "SUMMARY STATS LOCATION"))
    portal = f"{PORTAL}/{acc}"
    n = parse_n_samples(row_get(row, "INITIAL SAMPLE SIZE"))
    try:
        assoc = int(row_get(row, "ASSOCIATION COUNT") or 0)
    except (TypeError, ValueError):
        assoc = None
    ds = f"GWAS:{acc}"
    geno_list = [g.strip() for g in str(geno).split(",") if g.strip()] if geno else []
    assay = "GWAS"
    if geno_list:
        assay = f"GWAS ({'; '.join(geno_list[:3])})"
    upsert_dataset(c, dict(
        canonical_id=ds,
        title=disease or title or acc,
        organism="Homo sapiens",
        tissue=None,
        disease_tags=tags or ([mapped] if mapped else ([disease] if disease else [])),
        assay=assay,
        platforms=[platform] if platform else geno_list,
        n_samples=n,
        access="public",
        confidence="metadata_only",
        probe_status="ok",
    ))
    add_member(c, ds, "gwas_catalog", acc, "study", portal)
    if pmid:
        add_member(c, ds, "pubmed", str(pmid).strip(), "publication",
                   f"https://pubmed.ncbi.nlm.nih.gov/{str(pmid).strip()}/")
    add_artifact(c, dict(
        artifact_id=f"GWAS:{acc}:portal", dataset_id=ds, kind="processed",
        role="study_portal", title=disease or acc, format="html",
        uri=portal, access="public", source="gwas_catalog"))
    add_artifact(c, dict(
        artifact_id=f"GWAS:{acc}:api", dataset_id=ds, kind="processed",
        role="api_study_record", title="GWAS Catalog REST API v2 study pointer",
        format="json", uri=f"{API_STUDY}/{acc}", access="public", source="gwas_catalog"))
    if has_ss and ss_uri:
        add_artifact(c, dict(
            artifact_id=f"GWAS:{acc}:sumstats", dataset_id=ds, kind="processed",
            role="summary_statistics",
            title=f"{acc} full summary statistics (FTP pointer; not downloaded)",
            format="ftp_dir", uri=ss_uri, access="public", source="gwas_catalog"))
    raw = {k: row.get(k) for k in row}
    c.execute(
        "INSERT OR REPLACE INTO src_gwas_study VALUES("
        "?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
        (
            acc,
            disease,
            title,
            json.dumps(efo, ensure_ascii=False),
            json.dumps(bg_efo, ensure_ascii=False),
            str(pmid).strip() if pmid else None,
            row_get(row, "INITIAL SAMPLE SIZE"),
            row_get(row, "REPLICATION SAMPLE SIZE"),
            assoc,
            json.dumps(geno_list, ensure_ascii=False),
            platform,
            row_get(row, "COHORT"),
            int(has_ss),
            ss_uri,
            None,
            catalog_status,
            portal,
            int(bool(tags)),
            json.dumps(tags, ensure_ascii=False),
            json.dumps(raw, ensure_ascii=False)[:20000],
            now_iso(),
        ),
    )
    return True


def ingest_tsv(c, path: Path, catalog_status: str, cur: dict, start_offset: int = 0):
    scanned = int(cur.get("scanned", 0))
    liver = int(cur.get("liver", 0))
    errors = int(cur.get("errors", 0))
    with path.open(newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for i, row in enumerate(reader):
            if i < start_offset:
                continue
            try:
                if store_row(c, row, catalog_status):
                    scanned += 1
                    tags = liver_tags_for(
                        row_get(row, "DISEASE/TRAIT") or "",
                        row_get(row, "MAPPED_TRAIT") or "",
                        row_get(row, "MAPPED BACKGROUND TRAIT", "BACKGROUND TRAIT") or "",
                        row_get(row, "STUDY") or "",
                        row_get(row, "STUDY ACCESSION") or "",
                    )
                    if tags:
                        liver += 1
            except Exception as e:
                errors += 1
                log.warning("%s row %d error: %s", catalog_status, i, e)
            if (i + 1) % 2000 == 0:
                c.commit()
                cur.update({
                    "phase": catalog_status,
                    "offset": i + 1,
                    "scanned": scanned,
                    "liver": liver,
                    "errors": errors,
                })
                save_cursor(c, SRC, cur)
                log.info("%s offset=%d scanned=%d liver=%d errors=%d",
                         catalog_status, i + 1, scanned, liver, errors)
    c.commit()
    cur.update({"scanned": scanned, "liver": liver, "errors": errors, "offset": 0})
    save_cursor(c, SRC, cur)
    return scanned, liver, errors


def main():
    c = connect(DB, DDL)
    cur = load_cursor(c, SRC)
    if cur.get("complete"):
        log.info("already complete; exit")
        return

    cur["api_note"] = (
        "Ingest via FTP release TSVs (official Catalog downloads). "
        "REST API v2 /studies reports the same published totalElements; "
        "legacy /rest/api/studies returns 410. Associations not ingested."
    )
    cur["ftp_latest"] = FTP_LATEST
    save_cursor(c, SRC, cur)

    phase = cur.get("phase") or "download"
    if phase in ("download", "published"):
        pub_path = download_tsv(PUB_TSV)
        cur["published_tsv"] = str(pub_path)
        cur["phase"] = "published"
        cur["offset"] = int(cur.get("offset", 0)) if cur.get("phase") == "published" else 0
        save_cursor(c, SRC, cur)
        log.info("ingest published from %s offset=%s", pub_path, cur.get("offset", 0))
        ingest_tsv(c, pub_path, "published", cur, start_offset=int(cur.get("offset", 0)))
        cur["phase"] = "unpublished"
        cur["offset"] = 0
        save_cursor(c, SRC, cur)

    if cur.get("phase") == "unpublished":
        unpub_path = download_tsv(UNPUB_TSV)
        cur["unpublished_tsv"] = str(unpub_path)
        save_cursor(c, SRC, cur)
        log.info("ingest unpublished from %s offset=%s", unpub_path, cur.get("offset", 0))
        # Prefer published on conflict: only insert unpublished-only rows
        # (INSERT OR REPLACE would overwrite published with unpublished — avoid).
        # Use a temp approach: skip if accession already present.
        start = int(cur.get("offset", 0))
        scanned = int(cur.get("scanned", 0))
        liver = int(cur.get("liver", 0))
        errors = int(cur.get("errors", 0))
        with unpub_path.open(newline="", encoding="utf-8", errors="replace") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for i, row in enumerate(reader):
                if i < start:
                    continue
                acc = (row_get(row, "STUDY ACCESSION") or "").strip()
                if not acc:
                    continue
                exists = c.execute(
                    "SELECT 1 FROM src_gwas_study WHERE accession_id=?", (acc,)
                ).fetchone()
                if exists:
                    continue
                try:
                    if store_row(c, row, "unpublished"):
                        scanned += 1
                        tags = liver_tags_for(
                            row_get(row, "DISEASE/TRAIT") or "",
                            row_get(row, "MAPPED_TRAIT") or "",
                            row_get(row, "MAPPED BACKGROUND TRAIT", "BACKGROUND TRAIT") or "",
                            row_get(row, "STUDY") or "",
                            acc,
                        )
                        if tags:
                            liver += 1
                except Exception as e:
                    errors += 1
                    log.warning("unpublished row %d error: %s", i, e)
                if (i + 1) % 2000 == 0:
                    c.commit()
                    cur.update({
                        "phase": "unpublished",
                        "offset": i + 1,
                        "scanned": scanned,
                        "liver": liver,
                        "errors": errors,
                    })
                    save_cursor(c, SRC, cur)
                    log.info("unpublished offset=%d scanned=%d liver=%d errors=%d",
                             i + 1, scanned, liver, errors)
        c.commit()
        cur.update({"scanned": scanned, "liver": liver, "errors": errors})

    cur["complete"] = True
    cur["finished_at"] = now_iso()
    cur["phase"] = "done"
    save_cursor(c, SRC, cur)
    n_db = c.execute("SELECT COUNT(*) FROM src_gwas_study").fetchone()[0]
    n_liv = c.execute("SELECT COUNT(*) FROM src_gwas_study WHERE is_liver=1").fetchone()[0]
    n_ss = c.execute(
        "SELECT COUNT(*) FROM src_gwas_study WHERE full_summary_stats_available=1"
    ).fetchone()[0]
    n_pub = c.execute(
        "SELECT COUNT(*) FROM src_gwas_study WHERE catalog_status='published'"
    ).fetchone()[0]
    n_un = c.execute(
        "SELECT COUNT(*) FROM src_gwas_study WHERE catalog_status='unpublished'"
    ).fetchone()[0]
    log.info(
        "DONE src_rows=%d published=%d unpublished_only=%d liver=%d sumstats_ptr=%d errors=%s",
        n_db, n_pub, n_un, n_liv, n_ss, cur.get("errors", 0),
    )


if __name__ == "__main__":
    main()
