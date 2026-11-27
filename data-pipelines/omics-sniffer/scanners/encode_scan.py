#!/usr/bin/env python3
"""ENCODE Experiment metadata (human+mouse) -> omics_encode_staging.db. File pointers only. Resumable by offset."""
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)
import json
from stage_common import *
DB = ROOT / "data" / "omics_encode_staging.db"; SRC = "encode_scan"; BASE = "https://www.encodeproject.org/search/"
from http_util import http_get
def fetch(url, log, thr, params=None, timeout=None, **kw):
    for a in range(1, 9):
        thr()
        try:
            r = http_get(url, params=params, timeout=timeout or (10, 45), **kw); r.raise_for_status(); return r
        except Exception as e:
            log.warning("retry %d/8: %s", a, e); time.sleep(min(60, 2 ** a))
    raise RuntimeError("giving up " + url)
log = get_logger(SRC); thr = Throttle(1.0); BATCH = 500
DDL = """CREATE TABLE IF NOT EXISTS src_encode_experiment(accession TEXT PRIMARY KEY, assay TEXT, organism TEXT, biosample_term TEXT,
 biosample_type TEXT, target TEXT, status TEXT, lab TEXT, n_files INTEGER, is_liver INTEGER, liver_tags TEXT, raw_json TEXT, fetched_at TEXT);"""
FIELDS = ["accession","assay_title","assay_term_name","status","lab.title","target.label","description","date_released",
 "biosample_ontology.term_name","biosample_ontology.classification","biosample_ontology.organ_slims","biosample_ontology.cell_slims",
 "replicates.library.biosample.organism.scientific_name","replicates.library.biosample.accession","files.accession","files.href",
 "files.file_format","files.output_type","files.file_size","files.assembly","files.status"]
LRE = re.compile(r"liver|hepat|hepg2|hep ?3b|huh-?7|plc/prf|sk-hep|kupffer|cholangi|bile", re.I)
def ltags(term, organs, cells):
    blob = " ".join([term or ""] + (organs or []) + (cells or [])); t = []
    if re.search(r"hepg2|hep ?3b|huh-?7|plc/prf|sk-hep", blob, re.I): t.append("liver_cell_line")
    if re.search(r"hepatocyte", blob, re.I): t.append("hepatocyte")
    if re.search(r"\bliver\b", blob, re.I): t.append("liver")
    if LRE.search(blob) and not t: t.append("liver_related")
    return t
def main():
    c = connect(DB, DDL); cur = load_cursor(c, SRC)
    done = set(cur.get("done", []))
    for org in ("Homo sapiens", "Mus musculus"):
        OF = "replicates.library.biosample.donor.organism.scientific_name"
        fr = fetch(BASE, log, thr, params=[("type","Experiment"),("format","json"),("limit",0),(OF,org)], headers={"Accept":"application/json"}).json()
        assays = [t["key"] for fc in fr.get("facets", []) if fc["field"] == "assay_title" for t in fc["terms"] if t.get("doc_count")]
        log.info("%s total=%s assays=%d", org, fr.get("total"), len(assays))
        for at in assays:
            key = f"{org}|{at}"
            if key in done: continue
            p = [("type","Experiment"),("format","json"),("limit","all"),(OF,org),("assay_title",at)] + [("field",f) for f in FIELDS]
            g = fetch(BASE, log, thr, params=p, headers={"Accept":"application/json"}, timeout=(10, 300)).json().get("@graph", [])
            for e in g:
                acc = e["accession"]; bo = e.get("biosample_ontology") or {}
                tag = ltags(bo.get("term_name"), bo.get("organ_slims"), bo.get("cell_slims")); ds = f"ENCODE:{acc}"
                files = e.get("files") or []; reps = e.get("replicates") or []
                upsert_dataset(c, dict(canonical_id=ds, title=e.get("description") or f"{e.get('assay_title')} {bo.get('term_name')}", organism=org,
                    tissue=bo.get("term_name"), disease_tags=tag, assay=e.get("assay_title"), platforms=[], n_samples=len(reps),
                    access="public", confidence="metadata_only", probe_status=e.get("status")))
                add_member(c, ds, "encode", acc, "experiment", f"https://www.encodeproject.org/experiments/{acc}/")
                for rp in reps:
                    b = ((rp.get("library") or {}).get("biosample") or {}).get("accession")
                    if b: upsert_sample(c, dict(sample_id=f"ENCODE:{b}", dataset_id=ds, title=b, organism=org, tissue=bo.get("term_name"),
                        attributes={"classification": bo.get("classification")}))
                for f in files:
                    if not f.get("accession"): continue
                    add_artifact(c, dict(artifact_id=f"ENCODE:{f['accession']}", dataset_id=ds, kind="processed" if f.get("file_format") not in ("fastq","bam") else "raw",
                        role=f.get("output_type"), filename=f["accession"], format=f.get("file_format"), uri="https://www.encodeproject.org"+f.get("href",""),
                        size_bytes=f.get("file_size"), access="public", source="encode", probe_status=f.get("status"), derived_from=f.get("assembly")))
                c.execute("INSERT OR REPLACE INTO src_encode_experiment VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)", (acc, e.get("assay_title"), org, bo.get("term_name"),
                    bo.get("classification"), (e.get("target") or {}).get("label"), e.get("status"), (e.get("lab") or {}).get("title"), len(files),
                    int(bool(tag)), json.dumps(tag), json.dumps(e)[:20000], now_iso()))
            done.add(key); cur["done"] = sorted(done); save_cursor(c, SRC, cur)
            log.info("%s assay=%s n=%d", org, at, len(g))
    cur["complete"] = True; save_cursor(c, SRC, cur); log.info("DONE")
if __name__ == "__main__": main()
