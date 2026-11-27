#!/usr/bin/env python3
"""GenomeCRISPR experiment catalog -> omics_genomecrispr_staging.db (never omics_v1.db).

Uses live POST https://genomecrispr.dkfz.de/api/experiments with PubMed ID range
queries (API returns at most ~20 docs per call; empty GET /api/experiments serves
the SPA HTML; GET /api/experiments/bibliography currently 502).

Stores publication + per-condition screen metadata and portal/API pointers ONLY.
Does NOT pull /api/sgrnas hit matrices (multi-MB per gene).
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)
import json
import re
from stage_common import *
import time

DB = ROOT / "data" / "omics_genomecrispr_staging.db"
SRC = "genomecrispr_scan"
API = "https://genomecrispr.dkfz.de/api/experiments"
PORTAL = "https://genomecrispr.dkfz.de"
log = get_logger(SRC)
thr = Throttle(0.35)

DDL = """
CREATE TABLE IF NOT EXISTS src_gc_publication(
  pubmed TEXT PRIMARY KEY,
  title TEXT,
  authors TEXT,
  abstract TEXT,
  n_screens INTEGER,
  n_cell_lines INTEGER,
  is_liver INTEGER,
  liver_tags TEXT,
  portal_url TEXT,
  raw_json TEXT,
  fetched_at TEXT
);
CREATE INDEX IF NOT EXISTS idx_src_gc_pub_liver ON src_gc_publication(is_liver);

CREATE TABLE IF NOT EXISTS src_gc_screen(
  screen_id TEXT PRIMARY KEY,
  pubmed TEXT,
  cell_line TEXT,
  condition_name TEXT,
  maxfc REAL,
  minfc REAL,
  title TEXT,
  is_liver INTEGER,
  liver_tags TEXT,
  portal_url TEXT,
  raw_json TEXT,
  fetched_at TEXT
);
CREATE INDEX IF NOT EXISTS idx_src_gc_screen_liver ON src_gc_screen(is_liver);
CREATE INDEX IF NOT EXISTS idx_src_gc_screen_pmid ON src_gc_screen(pubmed);
"""

LIVER_CELL_RE = re.compile(
    r"hepg2|huh-?7|heparg|hep3b|sk-?hep|plc/?prf|huh-?6|li-?7|egi-?1|"
    r"hepatocyte|hepatoma|hepatoblast|hepatocellular|"
    r"liver[- ]derived|liver cell|\bliver\b|"
    r"cholangiocarcin|\bcholangi|bile[- ]duct|biliary",
    re.I,
)


def liver_tags_for(*parts) -> list[str]:
    blob = " ".join(str(x) for x in parts if x)
    if not blob:
        return []
    b = blob.lower()
    t: list[str] = []
    if LIVER_CELL_RE.search(blob):
        t.append("liver")
    if re.search(r"\bhcc\b|hepatocellular", b) or "liver cancer" in b:
        t.append("HCC")
    if "cholangi" in b or "bile duct" in b or re.search(r"\bbiliary\b", b):
        t.append("CHOL")
    if "masld" in b or "nafld" in b or "nash" in b or "fatty liver" in b:
        t.append("MASLD")
    if "cirrhosis" in b or "cirrhotic" in b:
        t.append("cirrhosis")
    if "hepatitis" in b or re.search(r"\bhbv\b|\bhcv\b", b):
        t.append("hepatitis")
    if re.search(r"\bliver\b", b) or "hepatic" in b or "hepatoma" in b or "hepatocyte" in b:
        t.append("liver")
    if is_liver(blob) and not t:
        t.append("liver_related")
    return sorted(set(t))



def post_ex(query: dict) -> list:
    from http_util import http_post
    for attempt in range(1, MAX_RETRIES + 1):
        thr()
        try:
            r = http_post(
                API,
                data=json.dumps(query),
                timeout=TIMEOUT,
                headers={
                    "Content-Type": "application/json",
                    "Accept": "application/json",
                    "Connection": "close",
                    "User-Agent": "omics-sniffer/1.0",
                },
            )
            if r.status_code >= 500 or r.status_code == 429:
                raise RuntimeError(f"HTTP {r.status_code}")
            r.raise_for_status()
            if not r.content or r.content.strip() in (b"", b"null"):
                return []
            # SPA HTML fallback (Angular catch-all) — treat as empty/unavailable
            ct = (r.headers.get("content-type") or "").lower()
            if "html" in ct and r.text.lstrip().startswith("<!"):
                raise RuntimeError("API returned SPA HTML instead of JSON")
            j = r.json()
            if isinstance(j, list):
                return j
            if j is None:
                return []
            raise RuntimeError(f"unexpected JSON type {type(j)}")
        except Exception as e:
            log.warning("post retry %d/%d query=%s: %s", attempt, MAX_RETRIES, query, e)
            time.sleep(min(60, 2 ** attempt))
    raise RuntimeError(f"giving up POST {query}")


def expand_screens(pub: dict) -> list[dict]:
    out = []
    cond = pub.get("condition") or {}
    if not isinstance(cond, dict):
        return out
    for cell, phenotypes in cond.items():
        if not isinstance(phenotypes, dict):
            continue
        for pname, meta in phenotypes.items():
            if not isinstance(meta, dict):
                continue
            sid = meta.get("screenid")
            if not sid:
                continue
            out.append({
                "screen_id": sid,
                "pubmed": pub.get("pubmed"),
                "cell_line": cell,
                "condition_name": pname,
                "maxfc": meta.get("maxfc"),
                "minfc": meta.get("minfc"),
                "title": pub.get("title"),
            })
    return out


def store_publication(c, pub: dict) -> tuple[int, int]:
    """Returns (n_screens_stored, n_liver_screens)."""
    pmid = str(pub.get("pubmed") or "").strip()
    if not pmid:
        return 0, 0
    screens = expand_screens(pub)
    # tags from title/abstract/authors + cell lines
    cells = sorted({s["cell_line"] for s in screens if s.get("cell_line")})
    conds = sorted({s["condition_name"] for s in screens if s.get("condition_name")})
    tags = liver_tags_for(
        pub.get("title"), pub.get("abstract"), pub.get("authors"),
        " ".join(cells), " ".join(conds),
    )
    portal = f"{PORTAL}/#/experiments"  # hash SPA; PMID deep-link not stable
    ds = f"GENOMECRISPR:PMID:{pmid}"
    upsert_dataset(c, dict(
        canonical_id=ds,
        title=f"GenomeCRISPR PMID:{pmid}: {pub.get('title') or ''}"[:500],
        organism="Homo sapiens",
        tissue=None,
        disease_tags=tags or [],
        assay="CRISPR/Cas9 screen (GenomeCRISPR curated)",
        platforms=["GenomeCRISPR"],
        n_samples=len(screens) or None,
        access="public",
        confidence="metadata_only",
        probe_status="ok",
    ))
    add_member(c, ds, "pubmed", pmid, "publication",
               f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/")
    add_member(c, ds, "genomecrispr", pmid, "publication_record", PORTAL)
    add_artifact(c, dict(
        artifact_id=f"GENOMECRISPR:{pmid}:portal", dataset_id=ds,
        kind="processed", role="portal", title=pub.get("title"),
        format="html", uri=PORTAL, access="public", source="genomecrispr"))
    add_artifact(c, dict(
        artifact_id=f"GENOMECRISPR:{pmid}:api", dataset_id=ds,
        kind="processed", role="api_experiments_by_pubmed",
        title="POST /api/experiments pubmed $in pointer",
        format="json",
        uri=f"{API}?pubmed={pmid}",
        access="public", source="genomecrispr"))
    authors = pub.get("authors")
    if isinstance(authors, (list, dict)):
        authors = json.dumps(authors, ensure_ascii=False)
    elif authors is not None:
        authors = str(authors)
    # strip bulky nested condition from raw stored on pub (screens have it)
    raw_pub = {
        "pubmed": pmid,
        "title": pub.get("title"),
        "authors": authors,
        "abstract": (pub.get("abstract") or "")[:4000],
        "n_screens": len(screens),
        "cell_lines": cells,
        "conditions": conds,
        "_id": pub.get("_id"),
    }
    c.execute(
        "INSERT OR REPLACE INTO src_gc_publication VALUES(?,?,?,?,?,?,?,?,?,?,?)",
        (
            pmid,
            pub.get("title"),
            authors,
            (pub.get("abstract") or "")[:8000] or None,
            len(screens),
            len(cells),
            int(bool(tags)),
            json.dumps(tags, ensure_ascii=False),
            portal,
            json.dumps(raw_pub, ensure_ascii=False)[:20000],
            now_iso(),
        ),
    )
    n_liv = 0
    for s in screens:
        # Screen-level tags from cell/condition only (avoid paper-abstract bleed)
        stags = liver_tags_for(s.get("cell_line"), s.get("condition_name"))
        if stags:
            n_liv += 1
        sds = f"GENOMECRISPR:{s['screen_id']}"
        upsert_dataset(c, dict(
            canonical_id=sds,
            title=f"GenomeCRISPR {s['screen_id']}: {s.get('cell_line')} / {s.get('condition_name')}"[:500],
            organism="Homo sapiens",
            tissue=s.get("cell_line"),
            disease_tags=stags or [],
            assay=f"CRISPR screen — {s.get('condition_name')}",
            platforms=["GenomeCRISPR"],
            access="public",
            confidence="metadata_only",
            probe_status="ok",
        ))
        add_member(c, sds, "genomecrispr", s["screen_id"], "screen", PORTAL)
        add_member(c, sds, "pubmed", pmid, "publication",
                   f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/")
        add_artifact(c, dict(
            artifact_id=f"GENOMECRISPR:{s['screen_id']}:api_ptr",
            dataset_id=sds, kind="processed", role="api_screen_pointer",
            title="GenomeCRISPR screen via parent publication API (pointer; sgRNA matrix NOT fetched)",
            format="json", uri=f"{API}", access="public", source="genomecrispr"))
        c.execute(
            "INSERT OR REPLACE INTO src_gc_screen VALUES(?,?,?,?,?,?,?,?,?,?,?,?)",
            (
                s["screen_id"],
                pmid,
                s.get("cell_line"),
                s.get("condition_name"),
                s.get("maxfc"),
                s.get("minfc"),
                pub.get("title"),
                int(bool(stags)),
                json.dumps(stags, ensure_ascii=False),
                PORTAL,
                json.dumps(s, ensure_ascii=False)[:8000],
                now_iso(),
            ),
        )
    return len(screens), n_liv


def pmid_ranges(start: int, end: int, step: int) -> list[tuple[int, int]]:
    out = []
    lo = start
    while lo <= end:
        hi = min(lo + step - 1, end)
        out.append((lo, hi))
        lo = hi + 1
    return out


def harvest_range(c, lo: int, hi: int, seen: set, depth: int = 0) -> tuple[int, int, int]:
    """Adaptive split when API returns the 20-doc cap. Returns (pubs, screens, liver_screens)."""
    docs = post_ex({"pubmed": {"$gte": str(lo), "$lte": str(hi)}})
    if len(docs) >= 20 and lo < hi and depth < 12:
        mid = (lo + hi) // 2
        log.info("cap=20 split [%d,%d] -> [%d,%d]+[%d,%d] depth=%d", lo, hi, lo, mid, mid + 1, hi, depth)
        a = harvest_range(c, lo, mid, seen, depth + 1)
        b = harvest_range(c, mid + 1, hi, seen, depth + 1)
        return (a[0] + b[0], a[1] + b[1], a[2] + b[2])
    n_pub = n_scr = n_liv = 0
    for pub in docs:
        pmid = str(pub.get("pubmed") or "")
        if not pmid or pmid in seen:
            continue
        seen.add(pmid)
        s, lv = store_publication(c, pub)
        n_pub += 1
        n_scr += s
        n_liv += lv
    return n_pub, n_scr, n_liv


def probe_api() -> dict:
    """Record API health notes."""
    notes = {
        "api": API,
        "portal": PORTAL,
        "bibliography_get": "GET /api/experiments/bibliography -> 502 Proxy Error (probed)",
        "experiments_get": "GET /api/experiments serves SPA HTML (Angular); use POST",
        "sgrnas": "POST /api/sgrnas returns large matrices — intentionally NOT fetched",
    }
    try:
        docs = post_ex({"pubmed": {"$exists": True}})
        notes["probe_exists_n"] = len(docs)
        notes["probe_ok"] = True
    except Exception as e:
        notes["probe_ok"] = False
        notes["probe_error"] = str(e)
    return notes


def main():
    c = connect(DB, DDL)
    cur = load_cursor(c, SRC)
    if cur.get("complete"):
        log.info("already complete; exit")
        return

    notes = probe_api()
    cur["api_notes"] = notes
    save_cursor(c, SRC, cur)
    if not notes.get("probe_ok"):
        cur["complete"] = True
        cur["blocker"] = notes.get("probe_error")
        cur["finished_at"] = now_iso()
        save_cursor(c, SRC, cur, err=notes.get("probe_error"))
        log.error("API probe failed: %s — writing empty complete cursor", notes.get("probe_error"))
        return

    seen = set(cur.get("seen_pmids", []))
    # Cover sentinel + published PMID space used by GenomeCRISPR (db largely 2014–2017 era)
    ranges = [(999999, 999999)] + pmid_ranges(1, 999998, 200000) + pmid_ranges(10_000_000, 40_000_000, 500_000)
    # resume
    start_i = int(cur.get("range_index", 0))
    total_pub = int(cur.get("publications", 0))
    total_scr = int(cur.get("screens", 0))
    total_liv = int(cur.get("liver_screens", 0))

    for i, (lo, hi) in enumerate(ranges):
        if i < start_i:
            continue
        try:
            n_pub, n_scr, n_liv = harvest_range(c, lo, hi, seen)
            total_pub = len(seen)  # unique
            # recompute running screen counts from DB periodically
            if n_pub or n_scr:
                c.commit()
            cur.update({
                "range_index": i + 1,
                "seen_pmids": sorted(seen),
                "publications": len(seen),
                "last_range": [lo, hi],
                "last_range_pubs": n_pub,
                "last_range_screens": n_scr,
            })
            # refresh aggregates
            total_scr = c.execute("SELECT COUNT(*) FROM src_gc_screen").fetchone()[0]
            total_liv = c.execute("SELECT COUNT(*) FROM src_gc_screen WHERE is_liver=1").fetchone()[0]
            cur["screens"] = total_scr
            cur["liver_screens"] = total_liv
            save_cursor(c, SRC, cur)
            if n_pub or (i % 10 == 0):
                log.info(
                    "range[%d/%d] %d-%d pubs_delta=%d screens_db=%d liver=%d seen_pmid=%d",
                    i + 1, len(ranges), lo, hi, n_pub, total_scr, total_liv, len(seen),
                )
        except Exception as e:
            cur["last_error"] = str(e)
            save_cursor(c, SRC, cur, err=str(e))
            log.exception("range %d-%d failed", lo, hi)
            raise

    c.commit()
    n_pub = c.execute("SELECT COUNT(*) FROM src_gc_publication").fetchone()[0]
    n_scr = c.execute("SELECT COUNT(*) FROM src_gc_screen").fetchone()[0]
    n_liv = c.execute("SELECT COUNT(*) FROM src_gc_screen WHERE is_liver=1").fetchone()[0]
    cur.update({
        "complete": True,
        "finished_at": now_iso(),
        "publications": n_pub,
        "screens": n_scr,
        "liver_screens": n_liv,
        "seen_pmids": sorted(seen),
        "note": (
            "POST /api/experiments with pubmed $gte/$lte sweeps (cap~20/call, adaptive split). "
            "GET bibliography 502; sgRNA matrices not fetched. "
            f"Portal {PORTAL}."
        ),
    })
    save_cursor(c, SRC, cur)
    log.info("DONE pubs=%d screens=%d liver=%d db=%s", n_pub, n_scr, n_liv, DB)


if __name__ == "__main__":
    main()
