#!/usr/bin/env python3
"""CNGB CNSA public project catalog -> omics_cnsa_staging.db (never omics_v1.db).

CNSA has no documented public JSON search API (2024 paper: APIs planned).
CNGBdb /search/project/* is an SPA shell without anonymous JSON.
This scanner lists publicly downloadable CNP* project folders from
https://ftp.cngb.org/pub/CNSA/data{1-7}/ HTML directory indexes (no auth),
stores FTP + portal pointers, and optionally enriches titles from
STOmics accessions when omics_stomics_staging.db is present.

Does NOT download raw sequencing data.
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)
import json
from stage_common import *

DB = ROOT / "data" / "omics_cnsa_staging.db"
SRC = "cnsa_scan"
FTP_ROOT = "https://ftp.cngb.org/pub/CNSA"
PORTAL = "https://db.cngb.org/search/project"
STOMICS_DB = ROOT / "data" / "omics_stomics_staging.db"
DATA_DIRS = [f"data{i}" for i in range(1, 8)]
log = get_logger(SRC)
thr = Throttle(0.6)

DDL = """
CREATE TABLE IF NOT EXISTS src_cnsa_project(
  accession TEXT PRIMARY KEY, ftp_dir TEXT, ftp_url TEXT, portal_url TEXT,
  title TEXT, species TEXT, tissues TEXT, disease TEXT, source_note TEXT,
  is_liver INTEGER, liver_tags TEXT, enriched_from TEXT, fetched_at TEXT);
CREATE TABLE IF NOT EXISTS src_cnsa_blocker(
  key TEXT PRIMARY KEY, detail TEXT, fetched_at TEXT);
"""

def list_ftp_dir(subdir):
    url = f"{FTP_ROOT}/{subdir}/"
    r = fetch(url, log, thr)
    return sorted(set(re.findall(r'id="(CNP\d+)"', r.text))), url

def enrich_from_stomics():
    """Map CNP -> title/tissue/disease from STOmics accessions if staging DB exists."""
    out = {}
    if not STOMICS_DB.exists():
        return out
    try:
        sc = sqlite3.connect(STOMICS_DB, timeout=30)
        rows = sc.execute(
            "SELECT dataset_id, title, species, tissues, disease, accessions, is_liver, liver_tags "
            "FROM src_stomics_dataset").fetchall()
        sc.close()
    except Exception as e:
        log.warning("stomics enrich skip: %s", e); return out
    for did, title, species, tissues, disease, accessions, is_liv, ltags in rows:
        try:
            accs = json.loads(accessions or "[]")
        except Exception:
            continue
        for a in accs:
            if not isinstance(a, dict):
                continue
            val = a.get("value") or ""
            if not str(val).startswith("CNP"):
                continue
            out[val] = dict(
                title=title, species=species, tissues=tissues, disease=disease,
                is_liver=is_liv, liver_tags=ltags, enriched_from=f"STOMICS:{did}")
    return out

def tags_for(acc, title, tissues, disease, enrich):
    blob = " ".join([acc or "", title or "", tissues or "", disease or ""])
    t = []
    b = blob.lower()
    if re.search(r"\bhcc\b|hepatocellular", b): t.append("HCC")
    if "cholangi" in b or "bile duct" in b or "biliary" in b: t.append("CHOL")
    if "masld" in b or "nafld" in b or "nash" in b: t.append("MASLD")
    if "hepatocyte" in b: t.append("hepatocyte")
    if re.search(r"\bliver\b", b) or "hepat" in b: t.append("liver")
    if enrich and enrich.get("liver_tags"):
        try:
            t.extend(json.loads(enrich["liver_tags"]))
        except Exception:
            pass
    if is_liver(blob) and not t:
        t.append("liver_related")
    return sorted(set(t))

def main():
    c = connect(DB, DDL)
    cur = load_cursor(c, SRC)
    if cur.get("complete"):
        log.info("already complete; exit"); return

    c.execute("INSERT OR REPLACE INTO src_cnsa_blocker VALUES(?,?,?)", (
        "no_public_json_api",
        "CNSA/CNGBdb project search is SPA-only; no anonymous JSON list API. "
        "Catalog built from public FTP CNP* folders under /pub/CNSA/data1-7. "
        "Project title/sample cards require portal UI or future CNSA API.",
        now_iso()))
    c.commit()

    enrich = enrich_from_stomics()
    log.info("stomics enrich map size=%d", len(enrich))

    done = set(cur.get("done", []))
    dir_stats = cur.get("dir_stats", {})
    for sub in DATA_DIRS:
        if dir_stats.get(sub, {}).get("ok"):
            log.info("skip completed dir %s", sub); continue
        try:
            names, list_url = list_ftp_dir(sub)
        except Exception as e:
            log.warning("dir %s failed: %s", sub, e)
            dir_stats[sub] = {"ok": False, "error": str(e)}
            cur["dir_stats"] = dir_stats
            save_cursor(c, SRC, cur)
            continue
        log.info("%s: %d CNP folders", sub, len(names))
        for acc in names:
            if acc in done:
                continue
            en = enrich.get(acc) or {}
            title = en.get("title")
            species = en.get("species")
            tissues = en.get("tissues")
            disease = en.get("disease")
            ftp_url = f"{FTP_ROOT}/{sub}/{acc}/"
            portal = f"{PORTAL}/{acc}/"
            tags = tags_for(acc, title, tissues, disease, en)
            is_liv = int(bool(tags) or bool(en.get("is_liver")))
            ds = f"CNSA:{acc}"
            upsert_dataset(c, dict(
                canonical_id=ds, title=title or f"CNSA project {acc}",
                organism=None, tissue=None,
                disease_tags=tags, assay="CNSA multi-omics archive",
                platforms=["CNSA"], n_samples=None, access="public",
                license="CC BY 4.0 (CNSA public)", confidence="metadata_only",
                probe_status="ok"))
            add_member(c, ds, "cnsa", acc, "project", portal)
            add_artifact(c, dict(
                artifact_id=f"CNSA:{acc}:ftp", dataset_id=ds, kind="raw",
                role="ftp_project_dir", title=f"{acc}/", filename=acc,
                format="dir", uri=ftp_url, access="public", source="cnsa"))
            add_artifact(c, dict(
                artifact_id=f"CNSA:{acc}:portal", dataset_id=ds, kind="processed",
                role="project_portal", title=acc, format="html",
                uri=portal, access="public", source="cnsa"))
            c.execute(
                "INSERT OR REPLACE INTO src_cnsa_project VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (acc, sub, ftp_url, portal, title, species, tissues, disease,
                 "ftp_listing", is_liv, json.dumps(tags), en.get("enriched_from"), now_iso()))
            done.add(acc)
        dir_stats[sub] = {"ok": True, "n": len(names)}
        cur["done"] = sorted(done)
        cur["dir_stats"] = dir_stats
        save_cursor(c, SRC, cur)
        log.info("after %s stored_unique=%d", sub, len(done))

    cur["complete"] = all(dir_stats.get(d, {}).get("ok") for d in DATA_DIRS)
    cur["blocker"] = "no_public_json_api"
    save_cursor(c, SRC, cur)
    nliv = c.execute("SELECT COUNT(*) FROM src_cnsa_project WHERE is_liver=1").fetchone()[0]
    log.info("DONE projects=%d liver=%d complete=%s", len(done), nliv, cur["complete"])

if __name__ == "__main__":
    main()
