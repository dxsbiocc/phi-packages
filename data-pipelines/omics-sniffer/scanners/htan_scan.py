#!/usr/bin/env python3
"""HTAN atlas/biospecimen/assay metadata -> omics_htan_staging.db (never omics_v1.db).

Sources (no Synapse token required for this scan):
  - GitHub ncihtan/htan-portal public dumps: syn_data.json.gz, atlases.json,
    syn_metadata.json, research_network.json, entities_v7_0.csv (pointer only)
  - Synapse entity URLs stored as file pointers only (no assay download)

BLOCKER: downloading open-access assay files from Synapse still needs a free
Synapse account; portal ClickHouse / htan-cli portal queries need HTAN team
membership. This scanner uses the published portal metadata dump instead.
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)
import gzip, io, json, re
from collections import defaultdict
from stage_common import *

DB = ROOT / "data" / "omics_htan_staging.db"
SRC = "htan_scan"
GH = "https://raw.githubusercontent.com/ncihtan/htan-portal/master"
log = get_logger(SRC)
thr = Throttle(0.5)

DDL = """
CREATE TABLE IF NOT EXISTS src_htan_atlas(
  htan_id TEXT PRIMARY KEY, name TEXT, title TEXT, lead_institutions TEXT,
  phase TEXT, n_biospecimens INTEGER, n_diagnoses INTEGER, n_assay_files INTEGER,
  components TEXT, is_liver INTEGER, liver_tags TEXT, synapse_meta TEXT,
  raw_json TEXT, fetched_at TEXT);
CREATE TABLE IF NOT EXISTS src_htan_biospecimen(
  biospecimen_id TEXT PRIMARY KEY, atlas_id TEXT, parent_id TEXT,
  biospecimen_type TEXT, site TEXT, tumor_tissue_type TEXT,
  entity_id TEXT, is_liver INTEGER, liver_tags TEXT, raw_json TEXT, fetched_at TEXT);
CREATE TABLE IF NOT EXISTS src_htan_diagnosis(
  diagnosis_key TEXT PRIMARY KEY, atlas_id TEXT, participant_id TEXT,
  primary_diagnosis TEXT, tissue_origin TEXT, site_resection TEXT,
  entity_id TEXT, is_liver INTEGER, liver_tags TEXT, raw_json TEXT, fetched_at TEXT);
CREATE TABLE IF NOT EXISTS src_htan_file(
  file_key TEXT PRIMARY KEY, atlas_id TEXT, component TEXT, filename TEXT,
  file_format TEXT, data_file_id TEXT, parent_file_id TEXT, entity_id TEXT,
  synapse_url TEXT, is_liver INTEGER, liver_tags TEXT, fetched_at TEXT);
CREATE TABLE IF NOT EXISTS src_htan_blocker(
  key TEXT PRIMARY KEY, detail TEXT, fetched_at TEXT);
"""

LIVER_STRICT = re.compile(
    r"\bliver\b|hepatocellular|\bhcc\b|hepatocyte|cholangi|bile duct|biliary|hepatoblastoma|gallbladder",
    re.I,
)

def liver_tags(*parts):
    blob = " ".join(str(x) for x in parts if x)
    if not blob or not LIVER_STRICT.search(blob):
        return []
    t = []
    b = blob.lower()
    if re.search(r"hepatocellular|\bhcc\b", b): t.append("HCC")
    if re.search(r"cholangi|bile duct|biliary", b): t.append("CHOL")
    if "hepatoblastoma" in b: t.append("Hepatoblastoma")
    if "gallbladder" in b: t.append("Gallbladder")
    if re.search(r"\bliver\b|hepatocyte", b): t.append("liver")
    if not t: t.append("liver_related")
    return sorted(set(t))

def rec_dict(cols, rec):
    vals = rec.get("values") or []
    return {cols[i]: (vals[i] if i < len(vals) else None) for i in range(len(cols))}

def syn_url(eid):
    if not eid or not str(eid).startswith("syn"):
        return None
    return f"https://www.synapse.org/#!Synapse:{eid}"

def gh_get(path, binary=False):
    url = f"{GH}/{path}"
    for a in range(1, MAX_RETRIES + 1):
        thr()
        try:
            r = http_get(url, timeout=(10, 300))
            r.raise_for_status()
            return r
        except Exception as e:
            log.warning("retry %d/%d %s: %s", a, MAX_RETRIES, url, e)
            time.sleep(min(60, 2 ** a))
    raise RuntimeError("giving up " + url)

def main():
    c = connect(DB, DDL)
    cur = load_cursor(c, SRC)
    if cur.get("complete"):
        log.info("already complete; exit"); return

    # Document Synapse / portal blockers
    c.execute("INSERT OR REPLACE INTO src_htan_blocker VALUES(?,?,?)", (
        "synapse_assay_download",
        "Open-access assay file download via Synapse needs a free Synapse account; "
        "portal ClickHouse / htan query portal needs HTAN team membership. "
        "This scan uses GitHub-published portal metadata (syn_data.json.gz) — no token.",
        now_iso()))
    c.commit()

    atlases_meta = {a.get("htan_id", "").lower(): a for a in gh_get("data/atlases.json").json()}
    syn_meta = gh_get("data/syn_metadata.json").json()
    try:
        rn = gh_get("data/research_network.json").json()
    except Exception as e:
        log.warning("research_network: %s", e); rn = {}

    log.info("downloading public/syn_data.json.gz ...")
    raw = gh_get("public/syn_data.json.gz").content
    data = json.loads(gzip.decompress(raw))
    atlases = data.get("atlases") or []
    log.info("atlases_in_dump=%d", len(atlases))

    # placeholder portal dataset for catalog-level artifacts (must exist before artifacts FK)
    upsert_dataset(c, dict(
        canonical_id="HTAN:PORTAL", title="HTAN Data Portal metadata dump",
        organism="Homo sapiens", tissue=None, disease_tags=[], assay="HTAN portal catalog",
        platforms=["Synapse", "GitHub htan-portal"], n_samples=0, access="public",
        confidence="metadata_only", probe_status="ok"))
    # entity CSV pointer only (20MB release table)
    add_artifact(c, dict(
        artifact_id="HTAN:entities_v7_0.csv", dataset_id="HTAN:PORTAL",
        kind="processed", role="release_entities", title="HTAN entities_v7_0.csv (pointer)",
        filename="entities_v7_0.csv", format="csv",
        uri=f"{GH}/data/entities_v7_0.csv", access="public", source="htan"))

    clinical_comps = {
        "bts:Biospecimen", "bts:Diagnosis", "bts:Demographics", "bts:Exposure",
        "bts:Therapy", "bts:FollowUp", "bts:FamilyHistory", "bts:MolecularTest",
        "bts:ClinicalDataTier2", "bts:AccessoryManifest",
    }

    done = set(cur.get("done", []))
    for atlas in atlases:
        aid = atlas.get("htan_id") or ""
        if not aid or aid in done:
            continue
        am = atlases_meta.get(aid.lower(), {})
        title = am.get("title")
        if isinstance(title, dict):
            title = title.get("rendered")
        rn_entry = rn.get(aid.lower()) or rn.get(aid) or {}
        name = atlas.get("htan_name") or rn_entry.get("title") or title or aid
        lead = am.get("lead_institutions") or ""
        phase = rn_entry.get("phase")

        bio_block = atlas.get("bts:Biospecimen") or {}
        diag_block = atlas.get("bts:Diagnosis") or {}
        bio_cols = bio_block.get("column_order") or []
        diag_cols = diag_block.get("column_order") or []
        bio_recs = bio_block.get("record_list") or []
        diag_recs = diag_block.get("record_list") or []

        atlas_tags = set()
        liver_participants = set()
        n_files = 0
        components = []

        ds = f"HTAN:{aid}"
        # upsert dataset early so samples FK ok
        upsert_dataset(c, dict(
            canonical_id=ds, title=name, organism="Homo sapiens",
            tissue=None, disease_tags=[], assay="HTAN multi-omics atlas",
            platforms=[], n_samples=len(bio_recs), access="mixed_public_controlled",
            confidence="metadata_only", probe_status="ok",
            license="HTAN / Synapse open + dbGaP controlled"))
        add_member(c, ds, "htan", aid, "atlas", f"https://data.humantumoratlas.org/explore")
        add_member(c, ds, "htan_portal", aid, "atlas", "https://humantumoratlas.org/")

        # diagnoses
        for rec in diag_recs:
            row = rec_dict(diag_cols, rec)
            pid = row.get("HTAN Participant ID") or ""
            tags = liver_tags(
                row.get("Primary Diagnosis"), row.get("Tissue or Organ of Origin"),
                row.get("Site of Resection or Biopsy"), row.get("Metastasis at Diagnosis Site"),
                row.get("Progression or Recurrence Anatomic Site"))
            if tags:
                atlas_tags.update(tags)
                if pid: liver_participants.add(pid)
            key = f"{aid}:{pid}:{row.get('Id') or row.get('entityId') or len(diag_recs)}"
            c.execute(
                "INSERT OR REPLACE INTO src_htan_diagnosis VALUES(?,?,?,?,?,?,?,?,?,?,?)",
                (key, aid, pid, row.get("Primary Diagnosis"), row.get("Tissue or Organ of Origin"),
                 row.get("Site of Resection or Biopsy"), row.get("entityId"),
                 int(bool(tags)), json.dumps(tags), json.dumps(row)[:8000], now_iso()))

        # biospecimens -> samples
        for rec in bio_recs:
            row = rec_dict(bio_cols, rec)
            bid = row.get("HTAN Biospecimen ID") or row.get("Id") or ""
            if not bid:
                continue
            parent = row.get("HTAN Parent ID") or ""
            site = row.get("Site of Resection or Biopsy")
            ttype = row.get("Tumor Tissue Type")
            tags = liver_tags(site, ttype, bid)
            # inherit from participant diagnosis
            root = bid.split("_")
            # participant often HTA*_N
            parts = bid.split("_")
            cand = "_".join(parts[:2]) if len(parts) >= 2 else parent
            if cand in liver_participants or parent in liver_participants:
                tags = sorted(set(tags + ["liver_related"]))
            if tags:
                atlas_tags.update(tags)
            eid = row.get("entityId")
            c.execute(
                "INSERT OR REPLACE INTO src_htan_biospecimen VALUES(?,?,?,?,?,?,?,?,?,?,?)",
                (bid, aid, parent, row.get("Biospecimen Type"), site, ttype, eid,
                 int(bool(tags)), json.dumps(tags), json.dumps(row)[:8000], now_iso()))
            upsert_sample(c, dict(
                sample_id=f"HTAN:{bid}", dataset_id=ds, title=bid,
                organism="Homo sapiens", tissue=site or ttype,
                disease=(";".join(tags) if tags else None),
                condition=row.get("Biospecimen Type"),
                attributes={"parent_id": parent, "entityId": eid, "atlas": aid}))
            if eid:
                add_artifact(c, dict(
                    artifact_id=f"HTAN:{eid}:biospecimen", dataset_id=ds,
                    sample_id=f"HTAN:{bid}", kind="processed", role="biospecimen_metadata",
                    title=f"Biospecimen {bid}", format="csv",
                    uri=syn_url(eid), access="public", source="htan"))

        # assay / other components -> file pointers
        platforms = set()
        for ck, block in atlas.items():
            if not ck.startswith("bts:") or not isinstance(block, dict):
                continue
            components.append(ck)
            cols = block.get("column_order") or []
            recs = block.get("record_list") or []
            if ck in clinical_comps:
                continue
            assay_name = ck.replace("bts:", "")
            platforms.add(assay_name)
            for rec in recs:
                row = rec_dict(cols, rec)
                fname = row.get("Filename") or row.get("File Name")
                fmt = row.get("File Format") or row.get("File format")
                dfid = row.get("HTAN Data File ID") or row.get("Id")
                pfid = row.get("HTAN Parent Data File ID") or row.get("HTAN Parent Biospecimen ID")
                eid = row.get("entityId")
                if not (fname or dfid or eid):
                    continue
                n_files += 1
                tags = liver_tags(fname, dfid, pfid)
                # participant inheritance from file id prefix
                if dfid:
                    pref = "_".join(str(dfid).split("_")[:2])
                    if pref in liver_participants:
                        tags = sorted(set(tags + ["liver_related"]))
                fkey = f"{aid}:{ck}:{dfid or eid or fname}"
                surl = syn_url(eid)
                c.execute(
                    "INSERT OR REPLACE INTO src_htan_file VALUES(?,?,?,?,?,?,?,?,?,?,?,?)",
                    (fkey, aid, ck, fname, fmt, dfid, pfid, eid, surl,
                     int(bool(tags)), json.dumps(tags), now_iso()))
                if eid or fname:
                    add_artifact(c, dict(
                        artifact_id=f"HTAN:{eid or fkey}", dataset_id=ds,
                        kind="processed" if "Level3" in ck or "Level4" in ck else "raw",
                        role=assay_name, title=fname or dfid, filename=fname,
                        format=fmt, uri=surl or f"https://data.humantumoratlas.org/explore",
                        access="public_metadata_synapse_download_may_need_account",
                        source="htan", derived_from=pfid))

        # syn_metadata component -> synapse file pointers
        sm_list = syn_meta.get(aid) or syn_meta.get(aid.upper()) or []
        for item in sm_list:
            sid = item.get("synapseId")
            if not sid:
                continue
            add_artifact(c, dict(
                artifact_id=f"HTAN:meta:{sid}", dataset_id=ds,
                kind="processed", role="metadata_manifest",
                title=f"{item.get('component')} ({item.get('numItems')} items)",
                format="csv", uri=syn_url(sid), access="public", source="htan"))

        # atlas-level liver from name/description too
        atlas_tags.update(liver_tags(name, title, lead, json.dumps(rn_entry)[:2000]))

        upsert_dataset(c, dict(
            canonical_id=ds, title=name, organism="Homo sapiens",
            tissue="liver" if atlas_tags else None,
            disease_tags=sorted(atlas_tags), assay="HTAN multi-omics atlas",
            platforms=sorted(platforms)[:40], n_samples=len(bio_recs),
            access="mixed_public_controlled", confidence="metadata_only",
            probe_status="ok"))
        c.execute(
            "INSERT OR REPLACE INTO src_htan_atlas VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (aid, name, title, lead, str(phase) if phase is not None else None,
             len(bio_recs), len(diag_recs), n_files, json.dumps(components),
             int(bool(atlas_tags)), json.dumps(sorted(atlas_tags)),
             json.dumps(sm_list)[:20000], json.dumps({"htan_id": aid, "htan_name": name})[:5000],
             now_iso()))
        done.add(aid)
        cur["done"] = sorted(done)
        save_cursor(c, SRC, cur)
        log.info("%s biospecimens=%d diagnoses=%d files=%d liver_tags=%s",
                 aid, len(bio_recs), len(diag_recs), n_files, sorted(atlas_tags))

    cur["complete"] = True
    cur["blocker"] = "Synapse account needed for assay file download; metadata scan used GitHub portal dump"
    save_cursor(c, SRC, cur)
    log.info("DONE atlases=%d", len(done))

if __name__ == "__main__":
    main()
