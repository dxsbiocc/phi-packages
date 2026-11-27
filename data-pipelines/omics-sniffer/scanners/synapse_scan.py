#!/usr/bin/env python3
"""Synapse public omics catalogs -> omics_synapse_staging.db (never omics_v1.db).

BLOCKER (2026-10): Synapse REST API (repo-prod.prod.sagebase.org) returns HTTP 403
GEO_RESTRICTION from this region (NIH NOT-OD-25-083). Anonymous search / entity
GET / table query are unavailable here.

Fallback (metadata + URL pointers only):
  1. Sage-Bionetworks/portals GitHub configs — resource table syn IDs + portal names
  2. Public Sage portal sitemaps (ADKP, NF, CancerComplexity, ARK, ELITE, BSMN)
     — study/dataset/tool syn IDs embedded in Explore URLs
  3. Portal + synapse.org entity URL pointers (download still needs Synapse account
     from a non-restricted network)

No bulk file download. Resume via sync_cursors.
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)
import json, re
from urllib.parse import unquote, urlparse, parse_qs
from stage_common import *

DB = ROOT / "data" / "omics_synapse_staging.db"
SRC = "synapse_scan"
GH_RAW = "https://raw.githubusercontent.com/Sage-Bionetworks/portals/master"
GH_API_TREE = "https://api.github.com/repos/Sage-Bionetworks/portals/git/trees/master?recursive=1"
SYN_REPO = "https://repo-prod.prod.sagebase.org/repo/v1"
SYN_WEB = "https://www.synapse.org/Synapse:"
log = get_logger(SRC)
thr = Throttle(0.35)

PORTAL_SITEMAPS = [
    ("adknowledgeportal", "https://adknowledgeportal.synapse.org/sitemap.xml",
     "https://adknowledgeportal.synapse.org", "AD Knowledge Portal (AMP-AD)"),
    ("nf", "https://nf.synapse.org/sitemap.xml",
     "https://nf.synapse.org", "NF Data Portal"),
    ("cancercomplexity", "https://cancercomplexity.synapse.org/sitemap.xml",
     "https://cancercomplexity.synapse.org", "Cancer Complexity Knowledge Portal"),
    ("arkportal", "https://arkportal.synapse.org/sitemap.xml",
     "https://arkportal.synapse.org", "ARK Portal"),
    ("eliteportal", "https://eliteportal.synapse.org/sitemap.xml",
     "https://eliteportal.synapse.org", "ELITE Portal"),
    ("bsmn", "https://bsmn.synapse.org/sitemap.xml",
     "https://bsmn.synapse.org", "BSMN Portal"),
]

SYN_RE = re.compile(r"\bsyn\d+(?:\.\d+)?\b", re.I)

DDL = """
CREATE TABLE IF NOT EXISTS src_synapse_blocker(
  key TEXT PRIMARY KEY, detail TEXT, fetched_at TEXT);
CREATE TABLE IF NOT EXISTS src_synapse_portal_table(
  table_key TEXT PRIMARY KEY, portal TEXT, syn_id TEXT, syn_base TEXT,
  version TEXT, resource_name TEXT, sql_snippet TEXT, config_path TEXT,
  portal_url TEXT, is_liver INTEGER, liver_tags TEXT, fetched_at TEXT);
CREATE TABLE IF NOT EXISTS src_synapse_entity(
  syn_id TEXT PRIMARY KEY, syn_base TEXT, entity_kind TEXT, portal TEXT,
  name_hint TEXT, portal_page_url TEXT, synapse_url TEXT,
  is_liver INTEGER, liver_tags TEXT, raw_json TEXT, fetched_at TEXT);
CREATE INDEX IF NOT EXISTS idx_src_syn_ent_liver ON src_synapse_entity(is_liver);
CREATE INDEX IF NOT EXISTS idx_src_syn_ent_portal ON src_synapse_entity(portal);
"""


def liver_tags(*parts) -> list[str]:
    blob = " ".join(str(x) for x in parts if x)
    if not blob:
        return []
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
    if re.search(r"\bliver\b", b) or "hepatic" in b or "hepat" in b:
        t.append("liver")
    if is_liver(blob) and not t:
        t.append("liver_related")
    return sorted(set(t))


def syn_base(sid: str) -> str:
    return sid.split(".", 1)[0].lower()


def probe_geo_block(c):
    detail = None
    try:
        from http_util import http_get
        thr()
        rr = http_get(f"{SYN_REPO}/entity/syn2580853",
                      headers={"Accept": "application/json"}, timeout=(10, 20))
        detail = f"HTTP {rr.status_code} body={rr.text[:500]}"
    except Exception as e:
        detail = f"entity probe exception: {e}"
    c.execute("INSERT OR REPLACE INTO src_synapse_blocker VALUES(?,?,?)",
              ("geo_restriction",
               "Synapse REST API GEO_RESTRICTION from this region (NIH NOT-OD-25-083). "
               f"Probe: {detail}. Fallback: GitHub portal configs + public portal sitemaps only.",
               now_iso()))
    c.commit()
    log.warning("GEO_RESTRICTION documented: %s", detail[:300])


def ingest_portal_configs(c, cur):
    if cur.get("configs_done"):
        log.info("portal configs already done")
        return
    tree = fetch(GH_API_TREE, log, thr, headers={"Accept": "application/vnd.github+json"}).json()
    paths = [
        t["path"] for t in (tree.get("tree") or [])
        if t["path"].startswith("src/configurations/")
        and t["path"].endswith((".ts", ".tsx", ".js", ".json"))
        and "config-shell" not in t["path"]
    ]
    log.info("portal config files=%d", len(paths))
    n_tables = 0
    for path in paths:
        portal = path.split("/")[2] if len(path.split("/")) > 2 else "unknown"
        text = fetch(f"{GH_RAW}/{path}", log, thr).text
        # resource exports: export const studiesSql = 'SELECT * FROM syn123'
        for m in re.finditer(
            r"(?:export\s+)?(?:const|let|var)\s+(\w+)\s*=\s*[`'\"]([^`'\"]{0,2000})[`'\"]",
            text,
        ):
            name, val = m.group(1), m.group(2)
            for sid in SYN_RE.findall(val):
                base = syn_base(sid)
                ver = sid.split(".", 1)[1] if "." in sid else None
                key = f"{portal}:{name}:{base}"
                tags = liver_tags(portal, name, val, path)
                portal_home = f"https://{portal}.synapse.org" if portal != "adknowledgeportal" else "https://adknowledgeportal.synapse.org"
                # fix known hosts
                host_map = {
                    "adknowledgeportal": "https://adknowledgeportal.synapse.org",
                    "nf": "https://nf.synapse.org",
                    "cancercomplexity": "https://cancercomplexity.synapse.org",
                    "arkportal": "https://arkportal.synapse.org",
                    "bsmn": "https://bsmn.synapse.org",
                    "psychencode": "https://psychencode.synapse.org",
                    "digitalhealth": "https://digitalhealth.synapse.org",
                    "crc-researcher": "https://crc.synapse.org",
                    "stopadportal": "https://stopadportal.synapse.org",
                }
                portal_home = host_map.get(portal, f"https://{portal}.synapse.org")
                c.execute(
                    "INSERT OR REPLACE INTO src_synapse_portal_table VALUES(?,?,?,?,?,?,?,?,?,?,?,?)",
                    (key, portal, sid.lower(), base, ver, name, val[:500], path,
                     portal_home, int(bool(tags)), json.dumps(tags), now_iso()),
                )
                ds = f"SYNAPSE:TABLE:{base}"
                upsert_dataset(c, dict(
                    canonical_id=ds,
                    title=f"Synapse portal table {sid} ({portal}/{name})",
                    organism=None, tissue=None, disease_tags=tags,
                    assay="Synapse portal metadata table",
                    platforms=["Synapse"], n_samples=None,
                    access="public", confidence="metadata_only",
                    probe_status="geo_blocked_entity",
                ))
                add_member(c, ds, "synapse", base, "table", f"{SYN_WEB}{base}")
                add_member(c, ds, "sage_portal", portal, "portal", portal_home)
                add_artifact(c, dict(
                    artifact_id=f"SYNAPSE:{base}:portal_table",
                    dataset_id=ds, kind="processed", role="portal_table_pointer",
                    title=f"{portal} {name}", format="synapse_table",
                    uri=f"{SYN_WEB}{base}", access="public", source="synapse",
                    probe_status="geo_blocked"))
                add_artifact(c, dict(
                    artifact_id=f"SYNAPSE:{base}:github_config:{abs(hash(path+name))%10_000_000}",
                    dataset_id=ds, kind="processed", role="github_portal_config",
                    title=path, format="typescript",
                    uri=f"{GH_RAW}/{path}", access="public", source="synapse"))
                n_tables += 1
        # also harvest any bare syn IDs in config
        for sid in set(SYN_RE.findall(text)):
            base = syn_base(sid)
            tags = liver_tags(portal, path, text[:2000])
            c.execute(
                "INSERT OR IGNORE INTO src_synapse_entity VALUES(?,?,?,?,?,?,?,?,?,?,?)",
                (base, base, "config_reference", portal, None, None,
                 f"{SYN_WEB}{base}", int(bool(tags)), json.dumps(tags),
                 json.dumps({"from": path, "syn": sid.lower()}), now_iso()),
            )
        c.commit()
    cur["configs_done"] = True
    cur["n_config_table_rows"] = n_tables
    save_cursor(c, SRC, cur)
    log.info("portal configs ingested table_rows≈%d", n_tables)


def kind_from_url(url: str) -> str:
    u = url.lower()
    if "/studies" in u or "study=" in u:
        return "study"
    if "/datasets" in u or "dataset=" in u:
        return "dataset"
    if "/projects" in u or "project=" in u or "grant" in u:
        return "project"
    if "/tools" in u:
        return "tool"
    if "/publications" in u:
        return "publication"
    if "/files" in u:
        return "file"
    if "/initiatives" in u:
        return "initiative"
    return "entity"


def name_hint_from_url(url: str) -> str | None:
    # NF: /Explore/Studies/syn123 or path slug
    path = unquote(urlparse(url).path)
    parts = [p for p in path.split("/") if p]
    qs = parse_qs(urlparse(url).query)
    for k in ("Study", "study", "Dataset", "dataset", "Grant Number", "id"):
        if k in qs and qs[k]:
            return qs[k][0]
    if parts:
        last = parts[-1]
        if last.lower().startswith("syn"):
            return last
        if last not in ("DetailsPage", "Explore", "Studies", "Datasets", "Projects"):
            return last[:200]
    return None


def ingest_sitemaps(c, cur):
    done_portals = set(cur.get("sitemaps_done", []))
    for portal, sm_url, home, title in PORTAL_SITEMAPS:
        if portal in done_portals:
            continue
        try:
            r = fetch(sm_url, log, thr)
        except Exception as e:
            log.warning("sitemap fail %s: %s", portal, e)
            c.execute("INSERT OR REPLACE INTO src_synapse_blocker VALUES(?,?,?)",
                      (f"sitemap_{portal}", str(e), now_iso()))
            c.commit()
            done_portals.add(portal)
            cur["sitemaps_done"] = sorted(done_portals)
            save_cursor(c, SRC, cur)
            continue
        locs = re.findall(r"<loc>([^<]+)</loc>", r.text)
        syns_in_sm = set(SYN_RE.findall(r.text))
        log.info("sitemap %s locs=%d syn_mentions=%d", portal, len(locs), len(syns_in_sm))
        # portal catalog dataset
        pds = f"SYNAPSE:PORTAL:{portal}"
        upsert_dataset(c, dict(
            canonical_id=pds, title=title, organism=None, tissue=None,
            disease_tags=liver_tags(portal, title),
            assay="Sage Synapse data portal catalog",
            platforms=["Synapse", portal], n_samples=len(syns_in_sm) or None,
            access="public", confidence="metadata_only", probe_status="ok"))
        add_member(c, pds, "sage_portal", portal, "portal", home)
        add_artifact(c, dict(
            artifact_id=f"SYNAPSE:PORTAL:{portal}:sitemap",
            dataset_id=pds, kind="processed", role="sitemap",
            title=f"{portal} sitemap.xml", format="xml", uri=sm_url,
            access="public", source="synapse"))
        n_new = 0
        for loc in locs:
            found = SYN_RE.findall(loc)
            if not found:
                continue
            kind = kind_from_url(loc)
            hint = name_hint_from_url(loc)
            for sid in found:
                base = syn_base(sid)
                tags = liver_tags(portal, title, loc, hint, kind)
                # prefer richer row via REPLACE if we learn more later — use INSERT OR REPLACE
                c.execute(
                    "INSERT OR REPLACE INTO src_synapse_entity VALUES(?,?,?,?,?,?,?,?,?,?,?)",
                    (base, base, kind, portal, hint, loc, f"{SYN_WEB}{base}",
                     int(bool(tags)), json.dumps(tags),
                     json.dumps({"sitemap_url": loc, "syn": sid.lower()}), now_iso()),
                )
                ds = f"SYNAPSE:{base}"
                upsert_dataset(c, dict(
                    canonical_id=ds,
                    title=hint or f"Synapse {kind} {base}",
                    organism=None, tissue=None, disease_tags=tags,
                    assay=f"Synapse {kind}",
                    platforms=["Synapse", portal], n_samples=None,
                    access="public", confidence="metadata_only",
                    probe_status="geo_blocked_entity",
                ))
                add_member(c, ds, "synapse", base, kind, f"{SYN_WEB}{base}")
                add_member(c, ds, "sage_portal", portal, "portal_page", loc)
                add_artifact(c, dict(
                    artifact_id=f"SYNAPSE:{base}:portal_page",
                    dataset_id=ds, kind="processed", role="portal_detail_page",
                    title=hint or base, format="html", uri=loc,
                    access="public", source="synapse"))
                add_artifact(c, dict(
                    artifact_id=f"SYNAPSE:{base}:entity",
                    dataset_id=ds, kind="processed", role="synapse_entity_pointer",
                    title=base, format="synapse_entity", uri=f"{SYN_WEB}{base}",
                    access="public", source="synapse",
                    probe_status="geo_blocked"))
                n_new += 1
        c.commit()
        done_portals.add(portal)
        cur["sitemaps_done"] = sorted(done_portals)
        cur[f"sitemap_{portal}_locs"] = len(locs)
        cur[f"sitemap_{portal}_syn"] = len(syns_in_sm)
        save_cursor(c, SRC, cur)
        log.info("sitemap %s stored entity_upserts≈%d", portal, n_new)


def main():
    c = connect(DB, DDL)
    cur = load_cursor(c, SRC)
    if cur.get("complete"):
        log.info("already complete; exit")
        return
    cur.setdefault("api", SYN_REPO)
    cur.setdefault("fallback", "github_portals_configs+portal_sitemaps")
    save_cursor(c, SRC, cur)

    probe_geo_block(c)
    ingest_portal_configs(c, cur)
    ingest_sitemaps(c, cur)

    n_ent = c.execute("SELECT COUNT(*) FROM src_synapse_entity").fetchone()[0]
    n_tab = c.execute("SELECT COUNT(*) FROM src_synapse_portal_table").fetchone()[0]
    n_liv = c.execute("SELECT COUNT(*) FROM src_synapse_entity WHERE is_liver=1").fetchone()[0]
    n_ds = c.execute("SELECT COUNT(*) FROM datasets").fetchone()[0]
    cur["complete"] = True
    cur["finished"] = now_iso()
    cur["n_entities"] = n_ent
    cur["n_portal_tables"] = n_tab
    cur["n_liver_entities"] = n_liv
    cur["n_datasets"] = n_ds
    cur["note"] = (
        "Incomplete vs full Synapse public graph: REST API geo-blocked. "
        "Ingested Sage portal resource tables (GitHub) + study/dataset syn IDs "
        "from public portal sitemaps. Entity metadata/files not fetched. "
        "Resume: clear complete flag after geo access restored to enrich via REST."
    )
    save_cursor(c, SRC, cur)
    log.info("DONE entities=%d portal_tables=%d liver=%d datasets=%d",
             n_ent, n_tab, n_liv, n_ds)


if __name__ == "__main__":
    main()
