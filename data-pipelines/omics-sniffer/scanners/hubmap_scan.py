#!/usr/bin/env python3
"""HuBMAP consortium metadata -> omics_hubmap_staging.db (never omics_v1.db).

Public Search API (no Globus token):
  https://search.api.hubmapconsortium.org/v3/search
  https://search.api.hubmapconsortium.org/v3/files/search

Captures published (+ other public) Datasets, Samples, Donors, and non-image
file pointers (assets/portal URLs only). Does NOT download raw assay bytes.
Resumable via sync_cursors (source=hubmap_scan).
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
from http_util import http_get, http_post, TIMEOUT

DB = ROOT / "data" / "omics_hubmap_staging.db"
SRC = "hubmap_scan"
SEARCH = "https://search.api.hubmapconsortium.org/v3/search"
FILES_SEARCH = "https://search.api.hubmapconsortium.org/v3/files/search"
PORTAL = "https://portal.hubmapconsortium.org"
ASSETS = "https://assets.hubmapconsortium.org"
log = get_logger(SRC)
thr = Throttle(0.25)
PAGE = 200

# HuBMAP organ codes (sample_category=organ display_subtype). BD=Blood, LI=Large Intestine — NOT liver/bile.
ORGAN_NAME = {
    "LV": "Liver", "BD": "Blood", "BL": "Bladder", "BM": "Bone Marrow", "BR": "Brain",
    "BV": "Blood Vasculature", "HT": "Heart", "LB": "Main Bronchus (Left)",
    "LE": "Eye (Left)", "LF": "Fallopian Tube (Left)", "LI": "Large Intestine",
    "LK": "Kidney (Left)", "LL": "Lung (Left)", "LN": "Knee (Left)", "LO": "Ovary (Left)",
    "LY": "Lymph Node", "OT": "Other", "PA": "Pancreas", "PL": "Placenta",
    "RB": "Main Bronchus (Right)", "RF": "Fallopian Tube (Right)", "RK": "Kidney (Right)",
    "RL": "Lung (Right)", "RN": "Knee (Right)", "RO": "Ovary (Right)", "SI": "Small Intestine",
    "SK": "Skin", "SP": "Spleen", "TH": "Thymus", "UT": "Uterus",
}

# Useful processed / tabular / sequence-ish extensions (skip ~10M imaging .tif etc.)
USEFUL_EXTS = [
    ".h5ad", ".hdf5", ".h5", ".mtx", ".tsv", ".csv", ".gz", ".bcf", ".bam",
    ".fastq", ".loom", ".zarr", ".rds", ".arrow", ".parquet", ".txt", ".json",
    ".zip", ".xml", ".npz", ".xlsx", ".pkl", ".pickle", ".qs", ".h5mu",
]

DDL = """
CREATE TABLE IF NOT EXISTS src_hubmap_dataset(
  uuid TEXT PRIMARY KEY, hubmap_id TEXT, dataset_type TEXT, status TEXT,
  data_access_level TEXT, title TEXT, description TEXT, group_name TEXT,
  doi_url TEXT, registered_doi TEXT, contains_human_genetic_sequences INTEGER,
  organ_codes TEXT, organ_names TEXT, donor_hubmap_ids TEXT, sample_hubmap_ids TEXT,
  n_files_useful INTEGER, is_liver INTEGER, liver_tags TEXT,
  portal_url TEXT, assets_url TEXT, raw_json TEXT, fetched_at TEXT);
CREATE INDEX IF NOT EXISTS idx_src_hubmap_ds_hid ON src_hubmap_dataset(hubmap_id);
CREATE INDEX IF NOT EXISTS idx_src_hubmap_ds_liver ON src_hubmap_dataset(is_liver);

CREATE TABLE IF NOT EXISTS src_hubmap_sample(
  uuid TEXT PRIMARY KEY, hubmap_id TEXT, organ TEXT, organ_name TEXT,
  sample_category TEXT, display_subtype TEXT, description TEXT,
  donor_hubmap_id TEXT, data_access_level TEXT, group_name TEXT,
  is_liver INTEGER, liver_tags TEXT, portal_url TEXT, raw_json TEXT, fetched_at TEXT);
CREATE INDEX IF NOT EXISTS idx_src_hubmap_sm_hid ON src_hubmap_sample(hubmap_id);
CREATE INDEX IF NOT EXISTS idx_src_hubmap_sm_liver ON src_hubmap_sample(is_liver);

CREATE TABLE IF NOT EXISTS src_hubmap_donor(
  uuid TEXT PRIMARY KEY, hubmap_id TEXT, description TEXT, group_name TEXT,
  data_access_level TEXT, sex TEXT, age TEXT, race TEXT, cause_of_death TEXT,
  is_liver INTEGER, liver_tags TEXT, portal_url TEXT, raw_json TEXT, fetched_at TEXT);
CREATE INDEX IF NOT EXISTS idx_src_hubmap_dn_hid ON src_hubmap_donor(hubmap_id);

CREATE TABLE IF NOT EXISTS src_hubmap_file(
  file_key TEXT PRIMARY KEY, file_uuid TEXT, dataset_uuid TEXT, dataset_hubmap_id TEXT,
  rel_path TEXT, file_extension TEXT, size_bytes INTEGER, data_class TEXT,
  description TEXT, assets_url TEXT, organ_codes TEXT, is_liver INTEGER,
  fetched_at TEXT);
CREATE INDEX IF NOT EXISTS idx_src_hubmap_file_ds ON src_hubmap_file(dataset_hubmap_id);

CREATE TABLE IF NOT EXISTS src_hubmap_blocker(
  key TEXT PRIMARY KEY, detail TEXT, fetched_at TEXT);
"""

LIVER_RE = re.compile(
    r"\bliver\b|hepatic|hepatocyte|\bhcc\b|hepatocellular|cholangi|"
    r"bile[\s-]?duct|biliary|gallbladder|hepatoblastoma",
    re.I,
)


def liver_tags(*parts, organ_codes=None):
    codes = set(organ_codes or [])
    blob = " ".join(str(x) for x in parts if x)
    t = []
    if "LV" in codes or re.search(r"\bliver\b", blob, re.I):
        t.append("liver")
    if re.search(r"fetal\s+liver", blob, re.I):
        t.append("fetal_liver")
    if re.search(r"hepatocyte", blob, re.I):
        t.append("hepatocyte")
    if re.search(r"hepatic", blob, re.I):
        t.append("hepatic")
    if re.search(r"cholangi|bile[\s-]?duct|biliary", blob, re.I):
        t.append("bile_duct")
    if re.search(r"gallbladder", blob, re.I):
        t.append("gallbladder")
    if re.search(r"hepatoblastoma|\bhcc\b|hepatocellular", blob, re.I):
        t.append("HCC_related")
    if LIVER_RE.search(blob) and not t:
        t.append("liver_related")
    # Never treat BD (Blood) or LI (Large Intestine) as liver/bile by code alone.
    return sorted(set(t))


def es_post(url, body, timeout=(10, 120)):
    for a in range(1, MAX_RETRIES + 1):
        thr()
        try:
            r = http_post(
                url,
                data=json.dumps(body),
                headers={"Content-Type": "application/json", "Accept": "application/json"},
                timeout=timeout,
            )
            if r.status_code in (429,) or r.status_code >= 500:
                raise RuntimeError(f"HTTP {r.status_code}")
            if r.status_code == 303:
                loc = r.headers.get("Location")
                if not loc:
                    raise RuntimeError("303 without Location")
                thr()
                g = http_get(loc, timeout=(10, 180))
                g.raise_for_status()
                return g.json()
            r.raise_for_status()
            return r.json()
        except Exception as e:
            log.warning("retry %d/%d %s: %s", a, MAX_RETRIES, url, e)
            time.sleep(min(60, 2 ** a))
    raise RuntimeError(f"giving up {url}")


def search_pages(entity_type, source_fields, cursor_key, process_hit):
    """Paginate entities with search_after on uuid.keyword."""
    cur = load_cursor(c_global, SRC)
    phases = cur.setdefault("phases", {})
    st = phases.setdefault(cursor_key, {})
    if st.get("complete"):
        log.info("%s already complete", cursor_key)
        return
    search_after = st.get("search_after")
    done_n = int(st.get("n", 0))
    total = None
    while True:
        body = {
            "size": PAGE,
            "track_total_hits": True,
            "sort": [{"uuid.keyword": "asc"}],
            "query": {"bool": {"filter": [{"term": {"entity_type.keyword": entity_type}}]}},
            "_source": source_fields,
        }
        if search_after:
            body["search_after"] = search_after
        j = es_post(SEARCH, body)
        hits = (j.get("hits") or {}).get("hits") or []
        if total is None:
            total = ((j.get("hits") or {}).get("total") or {}).get("value")
            log.info("%s total≈%s resume_n=%d", cursor_key, total, done_n)
        if not hits:
            break
        for h in hits:
            process_hit(h.get("_source") or {})
            done_n += 1
        search_after = hits[-1].get("sort")
        st.update(n=done_n, search_after=search_after, total=total)
        phases[cursor_key] = st
        cur["phases"] = phases
        if done_n % 500 < PAGE:
            save_cursor(c_global, SRC, cur)
            c_global.commit()
            log.info("%s progress n=%d/%s", cursor_key, done_n, total)
        if len(hits) < PAGE:
            break
    st["complete"] = True
    st["n"] = done_n
    phases[cursor_key] = st
    cur["phases"] = phases
    save_cursor(c_global, SRC, cur)
    c_global.commit()
    log.info("%s DONE n=%d", cursor_key, done_n)


def donor_meta(donor):
    sex = age = race = cod = None
    md = donor.get("metadata") or {}
    for key in ("organ_donor_data", "living_donor_data"):
        for row in md.get(key) or []:
            term = (row.get("grouping_concept_preferred_term") or row.get("preferred_term") or "").lower()
            val = row.get("data_value")
            units = row.get("units") or ""
            if "sex" in term and not sex:
                sex = val
            elif term == "age" or term.endswith(" age") and not age:
                age = f"{val} {units}".strip() if val is not None else None
            elif "race" in term and not race:
                race = val
            elif "cause of death" in term and not cod:
                cod = val
    return sex, age, race, cod


def organ_label(code, display=None):
    if display:
        return display
    return ORGAN_NAME.get(code or "", code)


def ensure_portal_dataset(c):
    upsert_dataset(c, dict(
        canonical_id="HUBMAP:PORTAL",
        title="HuBMAP Data Portal public catalog",
        organism="Homo sapiens",
        tissue=None,
        disease_tags=[],
        assay="HuBMAP multi-assay atlas",
        platforms=["HuBMAP Search API v3"],
        n_samples=0,
        access="public",
        confidence="metadata_only",
        probe_status="ok",
        license="HuBMAP public / consortium",
    ))
    add_member(c, "HUBMAP:PORTAL", "hubmap_portal", "portal", "catalog",
               f"{PORTAL}/search")


c_global = None  # set in main


def process_donor(src):
    c = c_global
    uuid = src.get("uuid")
    hid = src.get("hubmap_id")
    if not uuid or not hid:
        return
    sex, age, race, cod = donor_meta(src)
    desc = src.get("description") or ""
    tags = liver_tags(desc, src.get("display_subtype"))
    portal = f"{PORTAL}/browse/donor/{uuid}"
    c.execute(
        "INSERT OR REPLACE INTO src_hubmap_donor VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
        (uuid, hid, desc, src.get("group_name"), src.get("data_access_level"),
         sex, age, race, cod, int(bool(tags)), json.dumps(tags), portal,
         json.dumps({k: src.get(k) for k in (
             "uuid", "hubmap_id", "description", "group_name", "data_access_level",
             "metadata", "protocol_url", "created_timestamp", "last_modified_timestamp"
         ) if k in src}, default=str)[:20000], now_iso()),
    )
    # Inventory sample row under portal catalog (insert-once)
    sid = f"HUBMAP:Donor:{hid}"
    if not c.execute("SELECT 1 FROM samples WHERE sample_id=?", (sid,)).fetchone():
        upsert_sample(c, dict(
            sample_id=sid, dataset_id="HUBMAP:PORTAL", title=hid,
            organism="Homo sapiens", tissue=None, disease=(";".join(tags) if tags else None),
            condition="Donor", sex=sex, age=age,
            attributes={"uuid": uuid, "race": race, "cause_of_death": cod,
                        "group_name": src.get("group_name"), "portal_url": portal},
        ))


def process_sample(src):
    c = c_global
    uuid = src.get("uuid")
    hid = src.get("hubmap_id")
    if not uuid or not hid:
        return
    organ = src.get("organ")
    subtype = src.get("display_subtype")
    oname = organ_label(organ, subtype)
    desc = src.get("description") or ""
    donor = src.get("donor") or {}
    donor_hid = donor.get("hubmap_id")
    tags = liver_tags(desc, subtype, oname, organ_codes=[organ] if organ else None)
    portal = f"{PORTAL}/browse/sample/{uuid}"
    c.execute(
        "INSERT OR REPLACE INTO src_hubmap_sample VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
        (uuid, hid, organ, oname, src.get("sample_category"), subtype, desc,
         donor_hid, src.get("data_access_level"), src.get("group_name"),
         int(bool(tags)), json.dumps(tags), portal,
         json.dumps({k: src.get(k) for k in (
             "uuid", "hubmap_id", "organ", "display_subtype", "sample_category",
             "description", "data_access_level", "group_name", "donor",
             "protocol_url", "created_timestamp", "last_modified_timestamp"
         ) if k in src}, default=str)[:20000], now_iso()),
    )
    sid = f"HUBMAP:Sample:{hid}"
    if not c.execute("SELECT 1 FROM samples WHERE sample_id=?", (sid,)).fetchone():
        upsert_sample(c, dict(
            sample_id=sid, dataset_id="HUBMAP:PORTAL", title=hid,
            organism="Homo sapiens", tissue=oname,
            disease=(";".join(tags) if tags else None),
            condition=src.get("sample_category"),
            attributes={"uuid": uuid, "organ_code": organ, "donor_hubmap_id": donor_hid,
                        "group_name": src.get("group_name"), "portal_url": portal},
        ))


def process_dataset(src):
    c = c_global
    uuid = src.get("uuid")
    hid = src.get("hubmap_id")
    if not uuid or not hid:
        return
    origins = src.get("origin_samples") or []
    sources = src.get("source_samples") or []
    organ_codes = []
    organ_names = []
    sample_hids = []
    for s in origins + sources:
        if not isinstance(s, dict):
            continue
        code = s.get("organ")
        if code and code not in organ_codes:
            organ_codes.append(code)
            organ_names.append(organ_label(code, s.get("display_subtype")))
        sh = s.get("hubmap_id")
        if sh and sh not in sample_hids:
            sample_hids.append(sh)
    donor = src.get("donor") or {}
    donor_hids = []
    if isinstance(donor, dict) and donor.get("hubmap_id"):
        donor_hids.append(donor["hubmap_id"])
    title = src.get("title") or f"{src.get('dataset_type')} {hid}"
    desc = src.get("description") or ""
    tags = liver_tags(title, desc, " ".join(organ_names), organ_codes=organ_codes)
    tissue = ", ".join(organ_names) if organ_names else None
    if tags and not tissue:
        tissue = "liver" if "liver" in tags or "fetal_liver" in tags else None
    portal = f"{PORTAL}/browse/dataset/{uuid}"
    assets = f"{ASSETS}/{uuid}/"
    ds = f"HUBMAP:{hid}"
    access = src.get("data_access_level") or "public"
    upsert_dataset(c, dict(
        canonical_id=ds, title=title, organism="Homo sapiens", tissue=tissue,
        disease_tags=tags, assay=src.get("dataset_type"),
        platforms=[src.get("group_name")] if src.get("group_name") else [],
        n_samples=len(sample_hids) or (1 if donor_hids else 0),
        access=access, confidence="metadata_only",
        probe_status=src.get("status"), license="HuBMAP " + access,
    ))
    add_member(c, ds, "hubmap", hid, "dataset", portal)
    if src.get("registered_doi") or src.get("doi_url"):
        doi = src.get("registered_doi") or (src.get("doi_url") or "").replace("https://doi.org/", "")
        try:
            add_member(c, ds, "doi", doi, "doi", src.get("doi_url"))
        except Exception:
            pass  # UNIQUE(source,accession) may collide across datasets
    # nested samples / donor
    for s in origins + sources:
        if not isinstance(s, dict) or not s.get("hubmap_id"):
            continue
        sh = s.get("hubmap_id")
        sid = f"HUBMAP:Sample:{sh}"
        code = s.get("organ")
        oname = organ_label(code, s.get("display_subtype"))
        stags = liver_tags(oname, s.get("description"), organ_codes=[code] if code else None)
        if not c.execute("SELECT 1 FROM samples WHERE sample_id=?", (sid,)).fetchone():
            upsert_sample(c, dict(
                sample_id=sid, dataset_id=ds, title=sh, organism="Homo sapiens",
                tissue=oname, disease=(";".join(stags) if stags else None),
                condition=s.get("sample_category"),
                attributes={"uuid": s.get("uuid"), "organ_code": code,
                            "dataset_hubmap_id": hid},
            ))
        else:
            # re-point inventory row toward a concrete dataset when still on PORTAL
            row = c.execute("SELECT dataset_id FROM samples WHERE sample_id=?", (sid,)).fetchone()
            if row and row[0] == "HUBMAP:PORTAL":
                c.execute("UPDATE samples SET dataset_id=? WHERE sample_id=?", (ds, sid))
    if donor.get("hubmap_id"):
        did = f"HUBMAP:Donor:{donor['hubmap_id']}"
        dtags = liver_tags(donor.get("description"), title, desc)
        if not c.execute("SELECT 1 FROM samples WHERE sample_id=?", (did,)).fetchone():
            upsert_sample(c, dict(
                sample_id=did, dataset_id=ds, title=donor["hubmap_id"],
                organism="Homo sapiens", tissue=tissue,
                disease=(";".join(dtags) if dtags else None), condition="Donor",
                attributes={"uuid": donor.get("uuid"), "dataset_hubmap_id": hid},
            ))
    # dataset-level file pointers (no raw download)
    add_artifact(c, dict(
        artifact_id=f"HUBMAP:{hid}:portal", dataset_id=ds,
        kind="processed", role="portal_page", title=f"HuBMAP portal {hid}",
        format="html", uri=portal, access=access, source="hubmap",
        probe_status=src.get("status"),
    ))
    add_artifact(c, dict(
        artifact_id=f"HUBMAP:{hid}:assets_root", dataset_id=ds,
        kind="processed", role="assets_root", title=f"HuBMAP assets root {hid}",
        format="url", uri=assets, access=access, source="hubmap",
        probe_status=src.get("status"),
    ))
    c.execute(
        "INSERT OR REPLACE INTO src_hubmap_dataset VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
        (uuid, hid, src.get("dataset_type"), src.get("status"), access, title, desc,
         src.get("group_name"), src.get("doi_url"), src.get("registered_doi"),
         int(bool(src.get("contains_human_genetic_sequences"))),
         json.dumps(organ_codes), json.dumps(organ_names), json.dumps(donor_hids),
         json.dumps(sample_hids), 0, int(bool(tags)), json.dumps(tags),
         portal, assets,
         json.dumps({k: src.get(k) for k in (
             "uuid", "hubmap_id", "dataset_type", "status", "data_access_level",
             "title", "description", "group_name", "doi_url", "registered_doi",
             "contains_human_genetic_sequences", "origin_samples", "source_samples",
             "donor", "created_timestamp", "published_timestamp",
             "last_modified_timestamp", "dbgap_study_url", "dbgap_sra_experiment_url",
         ) if k in src}, default=str)[:25000], now_iso()),
    )


def process_files():
    c = c_global
    cur = load_cursor(c, SRC)
    phases = cur.setdefault("phases", {})
    st = phases.setdefault("files_useful", {})
    if st.get("complete"):
        log.info("files_useful already complete")
        return
    search_after = st.get("search_after")
    done_n = int(st.get("n", 0))
    liver_ds = {
        r[0] for r in c.execute(
            "SELECT hubmap_id FROM src_hubmap_dataset WHERE is_liver=1"
        ).fetchall()
    }
    total = None
    while True:
        body = {
            "size": PAGE,
            "track_total_hits": True,
            "sort": [{"file_uuid.keyword": "asc"}],
            "query": {"bool": {"filter": [{"terms": {"file_extension.keyword": USEFUL_EXTS}}]}},
            "_source": [
                "file_uuid", "dataset_uuid", "dataset_hubmap_id", "rel_path",
                "file_extension", "size", "data_class", "description",
                "data_access_level", "organs",
            ],
        }
        if search_after:
            body["search_after"] = search_after
        j = es_post(FILES_SEARCH, body)
        hits = (j.get("hits") or {}).get("hits") or []
        if total is None:
            total = ((j.get("hits") or {}).get("total") or {}).get("value")
            log.info("files_useful total≈%s resume_n=%d", total, done_n)
        if not hits:
            break
        for h in hits:
            src = h.get("_source") or {}
            fu = src.get("file_uuid")
            du = src.get("dataset_uuid")
            dhid = src.get("dataset_hubmap_id")
            rel = src.get("rel_path") or ""
            if not fu or not du or not dhid:
                continue
            ext = src.get("file_extension") or ""
            organs = src.get("organs") or []
            ocodes = []
            for o in organs:
                if isinstance(o, dict) and o.get("code") and o["code"] not in ocodes:
                    ocodes.append(o["code"])
            is_liver = 1 if dhid in liver_ds or "LV" in ocodes else 0
            uri = f"{ASSETS}/{du}/{rel}" if rel else f"{ASSETS}/{du}/"
            fkey = f"{du}:{rel}" if rel else f"{du}:{fu}"
            c.execute(
                "INSERT OR REPLACE INTO src_hubmap_file VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (fkey, fu, du, dhid, rel, ext, src.get("size"), src.get("data_class"),
                 (src.get("description") or "")[:500], uri, json.dumps(ocodes),
                 is_liver, now_iso()),
            )
            ds = f"HUBMAP:{dhid}"
            # ensure dataset exists (may be unpublished edge case)
            if not c.execute("SELECT 1 FROM datasets WHERE canonical_id=?", (ds,)).fetchone():
                upsert_dataset(c, dict(
                    canonical_id=ds, title=dhid, organism="Homo sapiens",
                    tissue=None, disease_tags=[], assay=None, platforms=[],
                    n_samples=0, access=src.get("data_access_level") or "public",
                    confidence="metadata_only", probe_status="from_files_index",
                ))
            add_artifact(c, dict(
                artifact_id=f"HUBMAP:file:{fu}",
                dataset_id=ds,
                kind="processed" if "Processed" in (src.get("data_class") or "") else "raw",
                role=ext.lstrip(".") or "file",
                title=rel or fu,
                filename=rel.split("/")[-1] if rel else fu,
                format=ext.lstrip(".") or None,
                uri=uri,
                size_bytes=src.get("size"),
                access=src.get("data_access_level") or "public",
                source="hubmap",
            ))
            done_n += 1
        search_after = hits[-1].get("sort")
        st.update(n=done_n, search_after=search_after, total=total)
        phases["files_useful"] = st
        cur["phases"] = phases
        save_cursor(c, SRC, cur)
        c.commit()
        if done_n % 2000 < PAGE:
            log.info("files_useful progress n=%d/%s", done_n, total)
        if len(hits) < PAGE:
            break
    # backfill n_files_useful
    c.execute(
        """UPDATE src_hubmap_dataset SET n_files_useful=(
             SELECT COUNT(*) FROM src_hubmap_file f
             WHERE f.dataset_hubmap_id=src_hubmap_dataset.hubmap_id)"""
    )
    st["complete"] = True
    st["n"] = done_n
    phases["files_useful"] = st
    cur["phases"] = phases
    save_cursor(c, SRC, cur)
    c.commit()
    log.info("files_useful DONE n=%d", done_n)


def main():
    global c_global
    assert DB.name != "omics_v1.db"
    c_global = connect(DB, DDL)
    cur = load_cursor(c_global, SRC)
    if cur.get("complete"):
        log.info("already complete; exit")
        return

    c_global.execute(
        "INSERT OR REPLACE INTO src_hubmap_blocker VALUES(?,?,?)",
        (
            "imaging_files_omitted",
            "files index has ~10M public files (mostly .tif imaging). "
            "This scan stores dataset portal/assets roots plus ~115k useful "
            "non-image extensions only (h5ad/csv/tsv/mtx/...). "
            "Consortium-only entities need a Globus token and are not fetched.",
            now_iso(),
        ),
    )
    ensure_portal_dataset(c_global)
    c_global.commit()

    donor_fields = [
        "uuid", "hubmap_id", "description", "group_name", "data_access_level",
        "metadata", "protocol_url", "display_subtype",
        "created_timestamp", "last_modified_timestamp",
    ]
    sample_fields = [
        "uuid", "hubmap_id", "organ", "display_subtype", "sample_category",
        "description", "data_access_level", "group_name", "donor",
        "protocol_url", "created_timestamp", "last_modified_timestamp",
    ]
    dataset_fields = [
        "uuid", "hubmap_id", "dataset_type", "status", "data_access_level",
        "title", "description", "group_name", "doi_url", "registered_doi",
        "contains_human_genetic_sequences", "origin_samples", "source_samples",
        "donor", "created_timestamp", "published_timestamp",
        "last_modified_timestamp", "dbgap_study_url", "dbgap_sra_experiment_url",
    ]

    search_pages("Donor", donor_fields, "donors", process_donor)
    search_pages("Sample", sample_fields, "samples", process_sample)
    search_pages("Dataset", dataset_fields, "datasets", process_dataset)
    process_files()

    cur = load_cursor(c_global, SRC)
    cur["complete"] = True
    cur["blocker"] = (
        "Imaging file inventory (~10M .tif) omitted; consortium-restricted "
        "entities require Globus auth"
    )
    save_cursor(c_global, SRC, cur)
    # summary
    n_ds = c_global.execute("SELECT COUNT(*) FROM src_hubmap_dataset").fetchone()[0]
    n_sm = c_global.execute("SELECT COUNT(*) FROM src_hubmap_sample").fetchone()[0]
    n_dn = c_global.execute("SELECT COUNT(*) FROM src_hubmap_donor").fetchone()[0]
    n_f = c_global.execute("SELECT COUNT(*) FROM src_hubmap_file").fetchone()[0]
    n_liver = c_global.execute(
        "SELECT COUNT(*) FROM src_hubmap_dataset WHERE is_liver=1"
    ).fetchone()[0]
    log.info(
        "DONE datasets=%d samples=%d donors=%d useful_files=%d liver_datasets=%d",
        n_ds, n_sm, n_dn, n_f, n_liver,
    )


if __name__ == "__main__":
    main()
