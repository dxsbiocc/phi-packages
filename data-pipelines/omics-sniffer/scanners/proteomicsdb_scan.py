#!/usr/bin/env python3
"""ProteomicsDB public OData catalog -> omics_proteomicsdb_staging.db (never omics_v1.db).

API: https://www.proteomicsdb.org/proteomicsdb/logic/api_v2/api.xsodata
No auth. Metadata + portal URL pointers only (no expression matrix / spectrum download).

Ingests: Project, Experiment, Tissue (full), liver-ish Sample (filtered), OmicsProject.
File entity has names but no download URLs — store portal experiment pointers instead.
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)
import json, re
from urllib.parse import quote
from stage_common import *

DB = ROOT / "data" / "omics_proteomicsdb_staging.db"
SRC = "proteomicsdb_scan"
ODATA = "https://www.proteomicsdb.org/proteomicsdb/logic/api_v2/api.xsodata"
PORTAL = "https://www.proteomicsdb.org"
log = get_logger(SRC)
thr = Throttle(0.3)
PAGE = 200

DDL = """
CREATE TABLE IF NOT EXISTS src_pdb_project(
  project_id INTEGER PRIMARY KEY, name TEXT, description TEXT, pubmed_id TEXT,
  owner TEXT, g_uid TEXT, is_liver INTEGER, liver_tags TEXT, portal_url TEXT,
  raw_json TEXT, fetched_at TEXT);
CREATE TABLE IF NOT EXISTS src_pdb_experiment(
  experiment_id INTEGER PRIMARY KEY, project_id INTEGER, name TEXT,
  description TEXT, scope TEXT, g_uid TEXT, is_liver INTEGER, liver_tags TEXT,
  portal_url TEXT, raw_json TEXT, fetched_at TEXT);
CREATE INDEX IF NOT EXISTS idx_src_pdb_exp_proj ON src_pdb_experiment(project_id);
CREATE INDEX IF NOT EXISTS idx_src_pdb_exp_liver ON src_pdb_experiment(is_liver);
CREATE TABLE IF NOT EXISTS src_pdb_tissue(
  tissue_id TEXT PRIMARY KEY, name TEXT, definition TEXT, synonym TEXT,
  category TEXT, sap_synonym TEXT, visualize INTEGER, g_uid TEXT,
  is_liver INTEGER, liver_tags TEXT, fetched_at TEXT);
CREATE INDEX IF NOT EXISTS idx_src_pdb_tissue_liver ON src_pdb_tissue(is_liver);
CREATE TABLE IF NOT EXISTS src_pdb_sample(
  sample_id INTEGER PRIMARY KEY, experiment_id INTEGER, name TEXT,
  description TEXT, taxcode INTEGER, tissue_id TEXT, tissue TEXT,
  disease TEXT, treatment TEXT, mass_spectrometer TEXT, protease_name TEXT,
  is_liver INTEGER, liver_tags TEXT, portal_url TEXT, raw_json TEXT, fetched_at TEXT);
CREATE INDEX IF NOT EXISTS idx_src_pdb_sample_liver ON src_pdb_sample(is_liver);
CREATE TABLE IF NOT EXISTS src_pdb_omics_project(
  omics_project_id INTEGER PRIMARY KEY, name TEXT, description TEXT,
  pubmed_id TEXT, owner TEXT, g_uid TEXT, is_liver INTEGER, liver_tags TEXT,
  raw_json TEXT, fetched_at TEXT);
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
    if "cirrhosis" in b:
        t.append("cirrhosis")
    if "hepatitis" in b or re.search(r"\bhbv\b|\bhcv\b", b):
        t.append("hepatitis")
    if re.search(r"hepg2|huh-?7|hep ?3b", b):
        t.append("liver_cell_line")
    if "hepatocyte" in b:
        t.append("hepatocyte")
    if re.search(r"\bliver\b", b) or "hepatic" in b or "hepat" in b:
        t.append("liver")
    if "bile" in b and "liver" not in t:
        t.append("bile")
    if "gall bladder" in b or "gallbladder" in b:
        t.append("gallbladder")
    if is_liver(blob) and not t:
        t.append("liver_related")
    return sorted(set(t))


def odata_page(entity: str, skip: int, top: int = PAGE, extra: str = "") -> dict:
    # Keep $ params in query string; encode filter separately.
    url = f"{ODATA}/{entity}?$format=json&$top={top}&$skip={skip}&$inlinecount=allpages"
    if extra:
        # extra like "$filter=substringof('liver',tolower(TISSUE))"
        if extra.startswith("$filter="):
            filt = extra.split("=", 1)[1]
            url += "&$filter=" + quote(filt, safe="(),'")
        else:
            url += "&" + extra
    return fetch(url, log, thr, headers={"Accept": "application/json"}).json()


def odata_iter(entity: str, cur_key: str, c, cur, extra: str = ""):
    skip = int(cur.get(cur_key, 0))
    while True:
        j = odata_page(entity, skip, PAGE, extra=extra)
        dd = j.get("d") or {}
        rows = dd.get("results") or []
        count = dd.get("__count")
        yield skip, count, rows
        if not rows:
            break
        skip += len(rows)
        cur[cur_key] = skip
        if count is not None and skip >= int(count):
            break
        if len(rows) < PAGE:
            break


def tax_name(code) -> str | None:
    try:
        code = int(code)
    except (TypeError, ValueError):
        return None
    return {9606: "Homo sapiens", 10090: "Mus musculus", 10116: "Rattus norvegicus"}.get(code, f"tax:{code}")


def scan_projects(c, cur):
    if cur.get("projects_done"):
        return
    n = 0
    for skip, count, rows in odata_iter("Project", "project_skip", c, cur):
        cur["project_total"] = count
        for r in rows:
            pid = r.get("PROJECT_ID")
            tags = liver_tags(r.get("NAME"), r.get("DESCRIPTION"))
            portal = f"{PORTAL}/project/{pid}"
            c.execute(
                "INSERT OR REPLACE INTO src_pdb_project VALUES(?,?,?,?,?,?,?,?,?,?,?)",
                (pid, r.get("NAME"), (r.get("DESCRIPTION") or "")[:8000],
                 r.get("PUBMEDID"), r.get("OWNER"), r.get("G_UID"),
                 int(bool(tags)), json.dumps(tags), portal,
                 json.dumps(r)[:20000], now_iso()),
            )
            ds = f"PROTEOMICSDB:PROJECT:{pid}"
            upsert_dataset(c, dict(
                canonical_id=ds, title=r.get("NAME"),
                organism=None, tissue=None, disease_tags=tags,
                assay="ProteomicsDB project (MS proteomics)",
                platforms=["ProteomicsDB"], n_samples=None,
                access="public", confidence="metadata_only", probe_status="ok",
            ))
            add_member(c, ds, "proteomicsdb", str(pid), "project", portal)
            if r.get("PUBMEDID"):
                add_member(c, ds, "pubmed", str(r["PUBMEDID"]), "publication",
                           f"https://pubmed.ncbi.nlm.nih.gov/{r['PUBMEDID']}/")
            add_artifact(c, dict(
                artifact_id=f"PROTEOMICSDB:PROJECT:{pid}:portal",
                dataset_id=ds, kind="processed", role="project_portal",
                title=r.get("NAME"), format="html", uri=portal,
                access="public", source="proteomicsdb"))
            add_artifact(c, dict(
                artifact_id=f"PROTEOMICSDB:PROJECT:{pid}:odata",
                dataset_id=ds, kind="processed", role="odata_project",
                title="OData Project entity", format="json",
                uri=f"{ODATA}/Project({pid})", access="public", source="proteomicsdb"))
            n += 1
        c.commit()
        save_cursor(c, SRC, cur)
        log.info("Project skip=%d count=%s page=%d", skip, count, len(rows))
    cur["projects_done"] = True
    cur["projects_stored"] = n
    save_cursor(c, SRC, cur)


def scan_experiments(c, cur):
    if cur.get("experiments_done"):
        return
    # cache project liver flags for inheritance
    proj_liver = {
        int(r[0]): json.loads(r[1] or "[]")
        for r in c.execute("SELECT project_id, liver_tags FROM src_pdb_project")
    }
    n = 0
    for skip, count, rows in odata_iter("Experiment", "experiment_skip", c, cur):
        cur["experiment_total"] = count
        for r in rows:
            eid = r.get("EXPERIMENT_ID")
            pid = r.get("PROJECT_ID")
            tags = liver_tags(r.get("NAME"), r.get("DESCRIPTION"), r.get("SCOPE"))
            if not tags and pid in proj_liver:
                tags = list(proj_liver.get(pid) or [])
            portal = f"{PORTAL}/project/{pid}/experiment/{eid}" if pid is not None else f"{PORTAL}/project"
            c.execute(
                "INSERT OR REPLACE INTO src_pdb_experiment VALUES(?,?,?,?,?,?,?,?,?,?,?)",
                (eid, pid, r.get("NAME"), (r.get("DESCRIPTION") or "")[:8000],
                 r.get("SCOPE"), r.get("G_UID"), int(bool(tags)), json.dumps(tags),
                 portal, json.dumps(r)[:20000], now_iso()),
            )
            ds = f"PROTEOMICSDB:EXP:{eid}"
            upsert_dataset(c, dict(
                canonical_id=ds, title=r.get("NAME") or f"Experiment {eid}",
                organism=None, tissue=None, disease_tags=tags,
                assay="MS proteomics experiment",
                platforms=["ProteomicsDB"], n_samples=None,
                access="public", confidence="metadata_only", probe_status="ok",
            ))
            add_member(c, ds, "proteomicsdb", str(eid), "experiment", portal)
            if pid is not None:
                add_member(c, ds, "proteomicsdb", str(pid), "project",
                           f"{PORTAL}/project/{pid}")
            add_artifact(c, dict(
                artifact_id=f"PROTEOMICSDB:EXP:{eid}:portal",
                dataset_id=ds, kind="processed", role="experiment_portal",
                title=r.get("NAME"), format="html", uri=portal,
                access="public", source="proteomicsdb"))
            n += 1
        c.commit()
        save_cursor(c, SRC, cur)
        log.info("Experiment skip=%d count=%s page=%d", skip, count, len(rows))
    cur["experiments_done"] = True
    cur["experiments_stored"] = n
    save_cursor(c, SRC, cur)


def scan_tissues(c, cur):
    if cur.get("tissues_done"):
        return
    n = 0
    for skip, count, rows in odata_iter("Tissue", "tissue_skip", c, cur):
        cur["tissue_total"] = count
        for r in rows:
            tid = r.get("TISSUE_ID")
            tags = liver_tags(
                r.get("NAME"), r.get("SYNONYM"), r.get("SAP_SYNONYM"),
                r.get("DEFINITION"), r.get("CATEGORY"),
            )
            c.execute(
                "INSERT OR REPLACE INTO src_pdb_tissue VALUES(?,?,?,?,?,?,?,?,?,?,?)",
                (tid, r.get("NAME"), (r.get("DEFINITION") or "")[:4000],
                 r.get("SYNONYM"), r.get("CATEGORY"), r.get("SAP_SYNONYM"),
                 r.get("VISUALIZE"), r.get("G_UID"),
                 int(bool(tags)), json.dumps(tags), now_iso()),
            )
            n += 1
        c.commit()
        save_cursor(c, SRC, cur)
        log.info("Tissue skip=%d count=%s page=%d", skip, count, len(rows))
    # tissue catalog dataset for liver tissues
    for tid, name, tags_json in c.execute(
        "SELECT tissue_id, name, liver_tags FROM src_pdb_tissue WHERE is_liver=1"
    ):
        tags = json.loads(tags_json or "[]")
        ds = f"PROTEOMICSDB:TISSUE:{tid}"
        upsert_dataset(c, dict(
            canonical_id=ds, title=f"ProteomicsDB tissue: {name}",
            organism="Homo sapiens", tissue=name, disease_tags=tags,
            assay="tissue expression catalog card",
            platforms=["ProteomicsDB"], access="public",
            confidence="metadata_only", probe_status="ok",
        ))
        add_member(c, ds, "proteomicsdb", tid, "tissue",
                   f"{PORTAL}/#/protein")  # expression UI is protein-centric
        add_artifact(c, dict(
            artifact_id=f"PROTEOMICSDB:TISSUE:{tid}:odata",
            dataset_id=ds, kind="processed", role="odata_tissue",
            title=name, format="json",
            uri=f"{ODATA}/Tissue('{tid.replace(':','%3A')}')",
            access="public", source="proteomicsdb"))
    c.commit()
    cur["tissues_done"] = True
    cur["tissues_stored"] = n
    save_cursor(c, SRC, cur)


def scan_liver_samples(c, cur):
    if cur.get("liver_samples_done"):
        return
    # OData substring filter on TISSUE / NAME / DISEASE
    filters = [
        "substringof('liver',tolower(TISSUE))",
        "substringof('liver',tolower(NAME))",
        "substringof('hepat',tolower(TISSUE))",
        "substringof('hepat',tolower(NAME))",
        "substringof('bile',tolower(TISSUE))",
        "substringof('hepg2',tolower(TISSUE))",
        "substringof('hepg2',tolower(NAME))",
    ]
    seen = set(cur.get("liver_sample_ids", []))
    n = 0
    for i, filt in enumerate(filters):
        key = f"liver_sample_skip_{i}"
        # reset skip per filter but dedupe by id
        cur[key] = 0
        for skip, count, rows in odata_iter(
            "Sample", key, c, cur, extra=f"$filter={filt}"
        ):
            cur[f"liver_sample_filter_{i}_total"] = count
            for r in rows:
                sid = r.get("SAMPLE_ID")
                if sid in seen:
                    continue
                eid = r.get("EXPERIMENT_ID")
                tags = liver_tags(
                    r.get("NAME"), r.get("DESCRIPTION"), r.get("TISSUE"),
                    r.get("DISEASE"), r.get("TREATMENT"),
                )
                # resolve project via experiment table
                pid = None
                if eid is not None:
                    row = c.execute(
                        "SELECT project_id FROM src_pdb_experiment WHERE experiment_id=?",
                        (eid,),
                    ).fetchone()
                    if row:
                        pid = row[0]
                portal = (
                    f"{PORTAL}/project/{pid}/experiment/{eid}/sample/{sid}"
                    if pid is not None and eid is not None
                    else f"{PORTAL}/project"
                )
                c.execute(
                    "INSERT OR REPLACE INTO src_pdb_sample VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    (sid, eid, r.get("NAME"), (r.get("DESCRIPTION") or "")[:4000],
                     r.get("TAXCODE"), r.get("TISSUE_ID"), r.get("TISSUE"),
                     r.get("DISEASE"), r.get("TREATMENT"), r.get("MASS_SPECTROMETER"),
                     r.get("PROTEASE_NAME"), int(bool(tags)), json.dumps(tags),
                     portal, json.dumps(r)[:20000], now_iso()),
                )
                # attach sample under experiment dataset if present
                ds = f"PROTEOMICSDB:EXP:{eid}" if eid is not None else f"PROTEOMICSDB:SAMPLE:{sid}"
                if eid is not None:
                    # ensure parent dataset exists
                    if not c.execute("SELECT 1 FROM datasets WHERE canonical_id=?", (ds,)).fetchone():
                        upsert_dataset(c, dict(
                            canonical_id=ds, title=f"Experiment {eid}",
                            assay="MS proteomics experiment",
                            platforms=["ProteomicsDB"], access="public",
                            confidence="metadata_only", probe_status="ok",
                            disease_tags=tags,
                        ))
                else:
                    upsert_dataset(c, dict(
                        canonical_id=ds, title=r.get("NAME") or f"Sample {sid}",
                        organism=tax_name(r.get("TAXCODE")), tissue=r.get("TISSUE"),
                        disease_tags=tags, assay="MS proteomics sample",
                        platforms=["ProteomicsDB"], access="public",
                        confidence="metadata_only", probe_status="ok",
                    ))
                upsert_sample(c, dict(
                    sample_id=f"PROTEOMICSDB:{sid}", dataset_id=ds,
                    title=r.get("NAME"), organism=tax_name(r.get("TAXCODE")),
                    tissue=r.get("TISSUE"), disease=r.get("DISEASE"),
                    condition=r.get("TREATMENT"),
                    attributes={
                        "tissue_id": r.get("TISSUE_ID"),
                        "mass_spectrometer": r.get("MASS_SPECTROMETER"),
                        "protease": r.get("PROTEASE_NAME"),
                        "is_liver": int(bool(tags)),
                        "liver_tags": tags,
                        "portal_url": portal,
                    },
                ))
                add_artifact(c, dict(
                    artifact_id=f"PROTEOMICSDB:SAMPLE:{sid}:portal",
                    dataset_id=ds, kind="processed", role="sample_portal",
                    title=r.get("NAME"), format="html", uri=portal,
                    access="public", source="proteomicsdb"))
                seen.add(sid)
                n += 1
            c.commit()
            cur["liver_sample_ids"] = sorted(seen)[-5000:]  # cap cursor size
            save_cursor(c, SRC, cur)
            log.info("liver Sample filter=%d skip=%d count=%s page=%d", i, skip, count, len(rows))
    cur["liver_samples_done"] = True
    cur["liver_samples_stored"] = len(seen)
    # drop bulky id list from final cursor
    cur.pop("liver_sample_ids", None)
    save_cursor(c, SRC, cur)
    log.info("liver samples DONE unique=%d", len(seen))


def scan_omics_projects(c, cur):
    if cur.get("omics_projects_done"):
        return
    n = 0
    for skip, count, rows in odata_iter("OmicsProject", "omics_project_skip", c, cur):
        cur["omics_project_total"] = count
        for r in rows:
            oid = r.get("OmicsProjectId")
            tags = liver_tags(r.get("Name"), r.get("Description"))
            c.execute(
                "INSERT OR REPLACE INTO src_pdb_omics_project VALUES(?,?,?,?,?,?,?,?,?,?)",
                (oid, r.get("Name"), (r.get("Description") or "")[:8000],
                 r.get("PubMedId"), r.get("Owner"), r.get("G_UID"),
                 int(bool(tags)), json.dumps(tags), json.dumps(r)[:20000], now_iso()),
            )
            ds = f"PROTEOMICSDB:OMICS_PROJECT:{oid}"
            upsert_dataset(c, dict(
                canonical_id=ds, title=r.get("Name"),
                disease_tags=tags, assay="ProteomicsDB omics project",
                platforms=["ProteomicsDB"], access="public",
                confidence="metadata_only", probe_status="ok",
            ))
            add_member(c, ds, "proteomicsdb", str(oid), "omics_project",
                       f"{ODATA}/OmicsProject({oid})")
            if r.get("PubMedId"):
                add_member(c, ds, "pubmed", str(r["PubMedId"]), "publication",
                           f"https://pubmed.ncbi.nlm.nih.gov/{r['PubMedId']}/")
            n += 1
        c.commit()
        save_cursor(c, SRC, cur)
        log.info("OmicsProject skip=%d count=%s page=%d", skip, count, len(rows))
    cur["omics_projects_done"] = True
    cur["omics_projects_stored"] = n
    save_cursor(c, SRC, cur)


def main():
    c = connect(DB, DDL)
    cur = load_cursor(c, SRC)
    if cur.get("complete"):
        log.info("already complete; exit")
        return
    cur["api"] = ODATA
    cur["api_docs"] = f"{PORTAL}/apiv2"
    cur["note_files"] = (
        "File OData entity (~62k) has NAME/G_UID but no download URL; "
        "portal HTML pointers used instead of file bytes."
    )
    save_cursor(c, SRC, cur)
    scan_projects(c, cur)
    scan_experiments(c, cur)
    scan_tissues(c, cur)
    scan_liver_samples(c, cur)
    scan_omics_projects(c, cur)
    n_p = c.execute("SELECT COUNT(*) FROM src_pdb_project").fetchone()[0]
    n_e = c.execute("SELECT COUNT(*) FROM src_pdb_experiment").fetchone()[0]
    n_t = c.execute("SELECT COUNT(*) FROM src_pdb_tissue").fetchone()[0]
    n_tl = c.execute("SELECT COUNT(*) FROM src_pdb_tissue WHERE is_liver=1").fetchone()[0]
    n_s = c.execute("SELECT COUNT(*) FROM src_pdb_sample").fetchone()[0]
    n_el = c.execute("SELECT COUNT(*) FROM src_pdb_experiment WHERE is_liver=1").fetchone()[0]
    cur["complete"] = True
    cur["finished"] = now_iso()
    cur.update(projects=n_p, experiments=n_e, experiment_liver=n_el,
               tissues=n_t, tissue_liver=n_tl, liver_samples=n_s)
    save_cursor(c, SRC, cur)
    log.info("DONE projects=%d experiments=%d (liver=%d) tissues=%d (liver=%d) liver_samples=%d",
             n_p, n_e, n_el, n_t, n_tl, n_s)


if __name__ == "__main__":
    main()
