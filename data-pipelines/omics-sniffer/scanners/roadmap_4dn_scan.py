#!/usr/bin/env python3
"""Roadmap Epigenomics + 4D Nucleome -> omics_roadmap_4dn_staging.db (never omics_v1.db).

One staging DB; `source` / src_* tables distinguish roadmap vs 4dn.
Metadata + file URL pointers only (no bulk downloads).

Roadmap: ENCODE portal API award.rfa=Roadmap (also award.project=Roadmap).
4DN: https://data.4dnucleome.org/ search ExperimentSetReplicate (+ file pointers
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)
from embedded experiments_in_set.files / processed_files).
"""
from __future__ import annotations
import json, re
from stage_common import *

DB = ROOT / "data" / "omics_roadmap_4dn_staging.db"
SRC = "roadmap_4dn_scan"
ENCODE = "https://www.encodeproject.org"
ENCODE_SEARCH = f"{ENCODE}/search/"
DN4 = "https://data.4dnucleome.org"
DN4_SEARCH = f"{DN4}/search/"
log = get_logger(SRC)
thr = Throttle(0.6)

DDL = """
CREATE TABLE IF NOT EXISTS src_roadmap_experiment(
  accession TEXT PRIMARY KEY, assay TEXT, organism TEXT, biosample_term TEXT,
  biosample_type TEXT, organ_slims TEXT, status TEXT, lab TEXT, award_rfa TEXT,
  n_files INTEGER, is_liver INTEGER, liver_tags TEXT, raw_json TEXT, fetched_at TEXT);
CREATE INDEX IF NOT EXISTS idx_src_rm_liver ON src_roadmap_experiment(is_liver);

CREATE TABLE IF NOT EXISTS src_4dn_experiment_set(
  accession TEXT PRIMARY KEY, display_title TEXT, description TEXT,
  dataset_label TEXT, experimentset_type TEXT, status TEXT,
  assay_types TEXT, organism TEXT, biosample_types TEXT, tissue_sources TEXT,
  organ_systems TEXT, lab TEXT, award_project TEXT, n_experiments INTEGER,
  n_file_ptrs INTEGER, public_release TEXT, is_liver INTEGER, liver_tags TEXT,
  portal_url TEXT, raw_json TEXT, fetched_at TEXT);
CREATE INDEX IF NOT EXISTS idx_src_4dn_liver ON src_4dn_experiment_set(is_liver);
"""

ENCODE_FIELDS = [
    "accession", "assay_title", "assay_term_name", "status", "lab.title",
    "description", "date_released", "award.rfa", "award.project",
    "biosample_ontology.term_name", "biosample_ontology.classification",
    "biosample_ontology.organ_slims", "biosample_ontology.cell_slims",
    "replicates.library.biosample.organism.scientific_name",
    "replicates.library.biosample.accession",
    "files.accession", "files.href", "files.file_format", "files.output_type",
    "files.file_size", "files.assembly", "files.status", "files.cloud_metadata.url",
]


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
    if re.search(r"hepg2|hep ?3b|huh-?7|plc/prf|sk-hep", b):
        t.append("liver_cell_line")
    if "hepatocyte" in b:
        t.append("hepatocyte")
    if re.search(r"\bliver\b", b) or "hepatic" in b or "hepat" in b:
        t.append("liver")
    if is_liver(blob) and not t:
        t.append("liver_related")
    return sorted(set(t))


def scan_roadmap(c, cur):
    if cur.get("roadmap_done"):
        log.info("roadmap already done")
        return
    # page by assay_title to keep payloads smaller / resumable
    done = set(cur.get("roadmap_assays_done", []))
    fr = fetch(
        ENCODE_SEARCH, log, thr,
        params=[("type", "Experiment"), ("award.rfa", "Roadmap"),
                ("format", "json"), ("limit", 0)],
        headers={"Accept": "application/json"},
    ).json()
    total = fr.get("total") or 0
    assays = [
        t["key"] for fc in fr.get("facets", []) if fc.get("field") == "assay_title"
        for t in fc.get("terms", []) if t.get("doc_count")
    ]
    log.info("Roadmap total=%s assays=%d", total, len(assays))
    cur["roadmap_total"] = total
    cur["roadmap_assays"] = assays
    save_cursor(c, SRC, cur)

    stored = int(cur.get("roadmap_stored", 0))
    for at in assays:
        if at in done:
            continue
        params = [
            ("type", "Experiment"), ("award.rfa", "Roadmap"),
            ("assay_title", at), ("format", "json"), ("limit", "all"),
        ] + [("field", f) for f in ENCODE_FIELDS]
        # large payload — use longer read via raw http_get with retries
        g = None
        for a in range(1, MAX_RETRIES + 1):
            thr()
            try:
                from http_util import http_get
                r = http_get(ENCODE_SEARCH, params=params,
                             headers={"Accept": "application/json"},
                             timeout=(10, 300))
                r.raise_for_status()
                g = r.json().get("@graph", [])
                break
            except Exception as e:
                log.warning("roadmap assay retry %d/%d %s: %s", a, MAX_RETRIES, at, e)
                time.sleep(min(60, 2 ** a))
        if g is None:
            raise RuntimeError(f"roadmap assay give up: {at}")
        for e in g:
            acc = e.get("accession")
            if not acc:
                continue
            bo = e.get("biosample_ontology") or {}
            organs = bo.get("organ_slims") or []
            cells = bo.get("cell_slims") or []
            reps = e.get("replicates") or []
            org = None
            for rp in reps:
                o = ((rp.get("library") or {}).get("biosample") or {}).get("organism") or {}
                if isinstance(o, dict):
                    org = o.get("scientific_name")
                if org:
                    break
            # organism sometimes only in field path as string list — skip if missing
            tags = liver_tags(
                e.get("description"), e.get("assay_title"),
                bo.get("term_name"), organs, cells, acc,
            )
            files = e.get("files") or []
            ds = f"ROADMAP:{acc}"
            upsert_dataset(c, dict(
                canonical_id=ds,
                title=e.get("description") or f"{e.get('assay_title')} {bo.get('term_name')}",
                organism=org, tissue=bo.get("term_name"), disease_tags=tags,
                assay=e.get("assay_title"), platforms=["Roadmap Epigenomics", "ENCODE"],
                n_samples=len(reps) or None, access="public",
                confidence="metadata_only", probe_status=e.get("status"),
            ))
            add_member(c, ds, "roadmap", acc, "experiment",
                       f"{ENCODE}/experiments/{acc}/")
            add_member(c, ds, "encode", acc, "experiment",
                       f"{ENCODE}/experiments/{acc}/")
            for f in files:
                fa = f.get("accession")
                if not fa:
                    continue
                href = f.get("href") or ""
                uri = href if href.startswith("http") else (ENCODE + href if href else None)
                cloud = (f.get("cloud_metadata") or {}).get("url")
                kind = "raw" if f.get("file_format") in ("fastq", "bam", "csfasta") else "processed"
                add_artifact(c, dict(
                    artifact_id=f"ROADMAP:{fa}",
                    dataset_id=ds, kind=kind, role=f.get("output_type"),
                    filename=fa, format=f.get("file_format"), uri=uri,
                    size_bytes=f.get("file_size"), access="public",
                    source="roadmap", probe_status=f.get("status"),
                    derived_from=f.get("assembly"),
                ))
                if cloud:
                    add_artifact(c, dict(
                        artifact_id=f"ROADMAP:{fa}:cloud",
                        dataset_id=ds, kind=kind, role="cloud_url_pointer",
                        filename=fa, format=f.get("file_format"), uri=cloud,
                        size_bytes=f.get("file_size"), access="public",
                        source="roadmap",
                    ))
            for rp in reps:
                b = ((rp.get("library") or {}).get("biosample") or {})
                bacc = b.get("accession") if isinstance(b, dict) else None
                if bacc:
                    upsert_sample(c, dict(
                        sample_id=f"ROADMAP:{bacc}", dataset_id=ds, title=bacc,
                        organism=org, tissue=bo.get("term_name"),
                        attributes={"classification": bo.get("classification"),
                                    "organ_slims": organs, "is_liver": int(bool(tags))},
                    ))
            c.execute(
                "INSERT OR REPLACE INTO src_roadmap_experiment VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (acc, e.get("assay_title"), org, bo.get("term_name"),
                 bo.get("classification"), json.dumps(organs), e.get("status"),
                 (e.get("lab") or {}).get("title"),
                 (e.get("award") or {}).get("rfa") or "Roadmap",
                 len(files), int(bool(tags)), json.dumps(tags),
                 json.dumps(e)[:25000], now_iso()),
            )
            stored += 1
        done.add(at)
        cur["roadmap_assays_done"] = sorted(done)
        cur["roadmap_stored"] = stored
        save_cursor(c, SRC, cur)
        log.info("roadmap assay=%s n=%d stored_total=%d", at, len(g), stored)
    cur["roadmap_done"] = True
    cur["roadmap_finished"] = now_iso()
    save_cursor(c, SRC, cur)
    log.info("roadmap DONE stored=%d", stored)


def _as_list(v):
    if v is None:
        return []
    if isinstance(v, list):
        return v
    return [v]


def _collect_4dn_blob(es: dict) -> str:
    parts = [
        es.get("accession"), es.get("display_title"), es.get("description"),
        es.get("dataset_label"), es.get("condition"),
    ]
    for exp in _as_list(es.get("experiments_in_set")):
        if not isinstance(exp, dict):
            continue
        parts.append(exp.get("display_title"))
        et = exp.get("experiment_type") or {}
        if isinstance(et, dict):
            parts.append(et.get("display_title"))
            parts.append(et.get("experiment_category"))
        else:
            parts.append(et)
        bs = exp.get("biosample") or {}
        if isinstance(bs, dict):
            parts.append(bs.get("display_title"))
            parts.append(bs.get("biosample_type"))
            parts.append(bs.get("biosource_summary"))
            toi = bs.get("tissue_organ_info") or {}
            if isinstance(toi, dict):
                parts.append(toi.get("tissue_source"))
                parts.extend(_as_list(toi.get("organ_system")))
            br = bs.get("biosource") or {}
            if isinstance(br, dict):
                parts.append((br.get("organism") or {}).get("name") if isinstance(br.get("organism"), dict) else br.get("organism"))
                parts.append((br.get("cell_line") or {}).get("display_title") if isinstance(br.get("cell_line"), dict) else None)
    return " ".join(str(x) for x in parts if x)


def _scalar(v):
    if v is None or isinstance(v, (str, int, float)):
        return v
    if isinstance(v, dict):
        for k in ("display_title", "file_format", "file_type", "name", "title", "@id"):
            if v.get(k):
                return str(v.get(k))
        return json.dumps(v, ensure_ascii=False)[:200]
    if isinstance(v, list):
        return ",".join(str(_scalar(x)) for x in v[:5])
    return str(v)


def _file_ptrs(exp: dict) -> list[dict]:
    out = []
    for key, role_prefix in (("files", "raw_or_submitted"), ("processed_files", "processed")):
        for f in _as_list(exp.get(key)):
            if not isinstance(f, dict):
                continue
            href = f.get("href") or f.get("@id") or ""
            if isinstance(href, dict):
                href = href.get("@id") or href.get("href") or ""
            acc = f.get("accession") or (str(href).rstrip("/").split("/")[-1] if href else None)
            if not acc:
                continue
            uri = href if str(href).startswith("http") else (DN4 + str(href) if href else None)
            fmt = _scalar(f.get("file_format") or f.get("file_type"))
            role = _scalar(f.get("file_type")) or role_prefix
            out.append({
                "accession": str(acc),
                "uri": uri,
                "format": fmt,
                "role": role,
                "size": f.get("file_size") or f.get("file_size_bytes"),
                "status": _scalar(f.get("status")),
            })
    return out


def scan_4dn(c, cur):
    if cur.get("4dn_done"):
        log.info("4dn already done")
        return
    page = 100
    offset = int(cur.get("4dn_from", 0))
    # probe total
    probe = fetch(
        DN4_SEARCH, log, thr,
        params={"type": "ExperimentSetReplicate", "limit": 0, "format": "json"},
        headers={"Accept": "application/json"},
    ).json()
    total = int(probe.get("total") or 0)
    cur["4dn_total"] = total
    save_cursor(c, SRC, cur)
    log.info("4DN ExperimentSetReplicate total=%d start_from=%d", total, offset)

    stored = int(cur.get("4dn_stored", 0))
    fields = [
        "accession", "display_title", "description", "dataset_label",
        "experimentset_type", "status", "public_release", "condition",
        "lab.display_title", "award.project", "award.center_title",
        "experiments_in_set.accession",
        "experiments_in_set.display_title",
        "experiments_in_set.experiment_type.display_title",
        "experiments_in_set.experiment_type.experiment_category",
        "experiments_in_set.biosample.accession",
        "experiments_in_set.biosample.display_title",
        "experiments_in_set.biosample.biosample_type",
        "experiments_in_set.biosample.biosource_summary",
        "experiments_in_set.biosample.tissue_organ_info.tissue_source",
        "experiments_in_set.biosample.tissue_organ_info.organ_system",
        "experiments_in_set.biosample.biosource.organism.name",
        "experiments_in_set.files.accession",
        "experiments_in_set.files.href",
        "experiments_in_set.files.file_format",
        "experiments_in_set.files.file_type",
        "experiments_in_set.files.file_size",
        "experiments_in_set.files.status",
        "experiments_in_set.processed_files.accession",
        "experiments_in_set.processed_files.href",
        "experiments_in_set.processed_files.file_format",
        "experiments_in_set.processed_files.file_type",
        "experiments_in_set.processed_files.file_size",
        "experiments_in_set.processed_files.status",
    ]
    while offset < max(total, 1):
        params = [("type", "ExperimentSetReplicate"), ("format", "json"),
                  ("limit", page), ("from", offset)] + [("field", f) for f in fields]
        j = fetch(DN4_SEARCH, log, thr, params=params,
                  headers={"Accept": "application/json"}).json()
        graph = j.get("@graph") or []
        if not graph:
            break
        for es in graph:
            acc = es.get("accession")
            if not acc:
                continue
            blob = _collect_4dn_blob(es)
            tags = liver_tags(blob)
            exps = [e for e in _as_list(es.get("experiments_in_set")) if isinstance(e, dict)]
            assay_types, organisms, btypes, tissues, organs = [], [], [], [], []
            n_files = 0
            for exp in exps:
                et = exp.get("experiment_type") or {}
                if isinstance(et, dict) and et.get("display_title"):
                    assay_types.append(et["display_title"])
                bs = exp.get("biosample") or {}
                if isinstance(bs, dict):
                    if bs.get("biosample_type"):
                        btypes.append(bs["biosample_type"])
                    toi = bs.get("tissue_organ_info") or {}
                    if isinstance(toi, dict):
                        if toi.get("tissue_source"):
                            tissues.append(toi["tissue_source"])
                        organs.extend(_as_list(toi.get("organ_system")))
                    br = bs.get("biosource") or {}
                    if isinstance(br, dict):
                        org = br.get("organism")
                        if isinstance(org, dict) and org.get("name"):
                            organisms.append(org["name"])
                n_files += len(_file_ptrs(exp))
            assay_types = sorted(set(assay_types))
            organisms = sorted(set(organisms))
            btypes = sorted(set(btypes))
            tissues = sorted(set(tissues))
            organs = sorted(set(str(x) for x in organs if x))
            portal = f"{DN4}/experiment-set-replicates/{acc}/"
            ds = f"4DN:{acc}"
            upsert_dataset(c, dict(
                canonical_id=ds,
                title=es.get("display_title") or es.get("dataset_label") or acc,
                organism=", ".join(organisms) if organisms else None,
                tissue=(tissues[0] if len(tissues) == 1 else (", ".join(tissues[:3]) if tissues else None)),
                disease_tags=tags,
                assay="; ".join(assay_types[:8]) if assay_types else "4DN experiment set",
                platforms=["4DN", "4D Nucleome"],
                n_samples=len(exps) or None, access="public",
                confidence="metadata_only", probe_status=es.get("status"),
            ))
            add_member(c, ds, "4dn", acc, "experiment_set", portal)
            for exp in exps:
                ea = exp.get("accession")
                if ea:
                    add_member(c, ds, "4dn", ea, "experiment",
                               f"{DN4}/experiments-seq/{ea}/")
                bs = exp.get("biosample") or {}
                if isinstance(bs, dict) and bs.get("accession"):
                    upsert_sample(c, dict(
                        sample_id=f"4DN:{bs['accession']}", dataset_id=ds,
                        title=bs.get("display_title") or bs["accession"],
                        organism=", ".join(organisms) if organisms else None,
                        tissue=bs.get("biosource_summary") or (tissues[0] if tissues else None),
                        attributes={
                            "biosample_type": bs.get("biosample_type"),
                            "tissue_organ_info": bs.get("tissue_organ_info"),
                            "is_liver": int(bool(tags)),
                            "liver_tags": tags,
                        },
                    ))
                for f in _file_ptrs(exp):
                    add_artifact(c, dict(
                        artifact_id=f"4DN:{f['accession']}",
                        dataset_id=ds,
                        kind="processed" if "processed" in (f.get("role") or "") else "raw",
                        role=f.get("role"), filename=f["accession"],
                        format=f.get("format"), uri=f.get("uri"),
                        size_bytes=f.get("size"), access="public",
                        source="4dn", probe_status=f.get("status"),
                    ))
                    n_files += 0  # already counted
            c.execute(
                "INSERT OR REPLACE INTO src_4dn_experiment_set VALUES("
                "?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (
                    acc, es.get("display_title"), (es.get("description") or "")[:4000],
                    es.get("dataset_label"), es.get("experimentset_type"), es.get("status"),
                    json.dumps(assay_types), json.dumps(organisms), json.dumps(btypes),
                    json.dumps(tissues), json.dumps(organs),
                    (es.get("lab") or {}).get("display_title") if isinstance(es.get("lab"), dict) else None,
                    (es.get("award") or {}).get("project") if isinstance(es.get("award"), dict) else None,
                    len(exps), n_files, es.get("public_release"),
                    int(bool(tags)), json.dumps(tags), portal,
                    json.dumps(es)[:30000], now_iso(),
                ),
            )
            stored += 1
        offset += page
        cur["4dn_from"] = offset
        cur["4dn_stored"] = stored
        save_cursor(c, SRC, cur)
        log.info("4dn from=%d/%d stored=%d page=%d", offset, total, stored, len(graph))
        if len(graph) < page:
            break
    cur["4dn_done"] = True
    cur["4dn_finished"] = now_iso()
    save_cursor(c, SRC, cur)
    log.info("4dn DONE stored=%d", stored)


def main():
    c = connect(DB, DDL)
    cur = load_cursor(c, SRC)
    if cur.get("complete"):
        log.info("already complete; exit")
        return
    cur["encode_api"] = ENCODE_SEARCH
    cur["roadmap_filter"] = "award.rfa=Roadmap"
    cur["4dn_api"] = DN4_SEARCH
    save_cursor(c, SRC, cur)
    scan_roadmap(c, cur)
    scan_4dn(c, cur)
    n_rm = c.execute("SELECT COUNT(*) FROM src_roadmap_experiment").fetchone()[0]
    n_rm_l = c.execute("SELECT COUNT(*) FROM src_roadmap_experiment WHERE is_liver=1").fetchone()[0]
    n_4 = c.execute("SELECT COUNT(*) FROM src_4dn_experiment_set").fetchone()[0]
    n_4_l = c.execute("SELECT COUNT(*) FROM src_4dn_experiment_set WHERE is_liver=1").fetchone()[0]
    cur["complete"] = True
    cur["finished"] = now_iso()
    cur["roadmap_rows"] = n_rm
    cur["roadmap_liver"] = n_rm_l
    cur["4dn_rows"] = n_4
    cur["4dn_liver"] = n_4_l
    save_cursor(c, SRC, cur)
    log.info("DONE roadmap=%d liver=%d | 4dn=%d liver=%d", n_rm, n_rm_l, n_4, n_4_l)


if __name__ == "__main__":
    main()
