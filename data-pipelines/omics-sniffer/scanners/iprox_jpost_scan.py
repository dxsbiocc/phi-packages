#!/usr/bin/env python3
"""iProX + jPOST public project metadata via ProteomeXchange PROXI (proteomecentral).
Writes ONLY omics_iprox_jpost_staging.db. No raw/mzML downloads (file counts only).
Resumable via sync_cursors (source=iprox_jpost_scan)."""
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)
import json, re
from stage_common import *
DB = ROOT / "data" / "omics_iprox_jpost_staging.db"; SRC = "iprox_jpost_scan"
API = "https://proteomecentral.proteomexchange.org/api/proxi/v0.1/datasets"
log = get_logger(SRC); thr = Throttle(0.5)
DDL = """CREATE TABLE IF NOT EXISTS src_px_dataset(pxd TEXT PRIMARY KEY, repository TEXT, repo_accession TEXT,
 title TEXT, species TEXT, tissue TEXT, disease TEXT, instruments TEXT, pmids TEXT, n_files INTEGER, n_raw INTEGER,
 announce_date TEXT, keywords TEXT, is_liver INTEGER, raw_json TEXT, fetched_at TEXT);"""
REPOS = ("iProX", "jPOST")

def list_repo(repo):
    out, page = [], 1
    while True:
        j = fetch(API, log, thr, params={"pageSize": 500, "pageNumber": page, "resultType": "compact", "repository": repo}).json()
        out += j["datasets"]; rs = j["result_set"]
        if page >= rs["n_available_pages"]: return out
        page += 1

def vals(lst, pat):
    o = []
    for t in lst or []:
        for x in (t.get("terms") if isinstance(t, dict) and "terms" in t else [t]):
            if isinstance(x, dict) and re.search(pat, x.get("name", ""), re.I) and x.get("value"): o.append(x["value"])
            elif isinstance(x, dict) and pat == "ANY" and x.get("name"): o.append(x["name"])
    return o

def main():
    c = connect(DB, DDL); cur = load_cursor(c, SRC); done = set(cur.get("done", []))
    for repo in REPOS:
        rows = list_repo(repo); log.info("%s: %d datasets listed", repo, len(rows))
        for i, r in enumerate(rows):
            pxd = r[0]
            if pxd in done: continue
            try:
                d = fetch(f"{API}/{pxd}", log, thr, ok404=True); d = d.json() if d else {}
            except Exception as e:
                log.error("%s detail failed: %s", pxd, e); save_cursor(c, SRC, {"done": sorted(done)}, str(e)[:300]); continue
            species = "; ".join(sorted(set(vals(d.get("species"), "scientific name|taxonomy")))) or r[3]
            kws = "; ".join(vals(d.get("keywords"), "keyword|submitter") ) or (r[10] or "")
            tissue = "; ".join(vals(d.get("sdrf_metadata") if isinstance(d.get("sdrf_metadata"), list) else [], "tissue|organism part"))
            disease = "; ".join(vals(d.get("sdrf_metadata") if isinstance(d.get("sdrf_metadata"), list) else [], "disease"))
            instr = sorted(set(vals(d.get("instruments"), "ANY"))) or ([r[6]] if r[6] else [])
            pm = sorted(set(vals(d.get("publications"), "PubMed")))
            ids = {x.get("name"): x.get("value") for x in d.get("identifiers", []) if isinstance(x, dict)}
            links = vals(d.get("fullDatasetLinks"), "URI|link")
            racc = next((v for v in ids.values() if v and re.match(r"(IPX|JPST)\d+", str(v))), None)
            if not racc:
                m = re.search(r"(IPX\d+|JPST\d+)", json.dumps(d)); racc = m.group(1) if m else None
            files = d.get("datasetFiles") or []
            n_raw = sum(1 for f in files if "raw" in (f.get("name", "").lower()))
            title = (d.get("title") or r[1] or "").strip()
            desc = (d.get("description") or "")[:2000]
            liver = is_liver(title, desc, kws, tissue, disease)
            slim = {k: d.get(k) for k in ("title", "species", "instruments", "publications", "identifiers", "keywords", "datasetHistory", "fullDatasetLinks")}
            slim["description"] = desc; slim["n_files"] = len(files)
            c.execute("INSERT OR REPLACE INTO src_px_dataset VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (pxd, repo, racc, title, species, tissue, disease, json.dumps(instr), json.dumps(pm), len(files), n_raw, r[9], kws, int(liver), json.dumps(slim, ensure_ascii=False)[:20000], now_iso()))
            tags = ["liver"] if liver else []
            if re.search(r"hepatocellular|\bhcc\b", f"{title} {desc} {kws}", re.I): tags.append("HCC")
            upsert_dataset(c, dict(canonical_id=f"PX:{pxd}", title=title, organism=species, tissue=tissue or None,
                disease_tags=tags, assay="proteomics (MS)", platforms=instr, n_samples=None, access="public",
                confidence="metadata_only", probe_status="ok"))
            add_member(c, f"PX:{pxd}", "proteomexchange", pxd, "dataset", f"https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID={pxd}")
            if racc: add_member(c, f"PX:{pxd}", repo.lower(), racc, "entry",
                f"https://www.iprox.cn/page/project.html?id={racc}" if racc.startswith("IPX") else f"https://repository.jpostdb.org/entry/{racc}")
            add_artifact(c, dict(artifact_id=f"PX:{pxd}:raw", dataset_id=f"PX:{pxd}", kind="raw", role="ms_raw_files",
                title=f"{len(files)} files ({n_raw} raw) - pointer only", uri=links[0] if links else None, access="public", source=repo.lower()))
            done.add(pxd)
            if len(done) % 25 == 0: save_cursor(c, SRC, {"done": sorted(done)}); log.info("%s %d/%d done=%d", repo, i + 1, len(rows), len(done))
        save_cursor(c, SRC, {"done": sorted(done)})
    cur = {"done": sorted(done), "finished": now_iso()}; save_cursor(c, SRC, cur); log.info("FINISHED %d datasets", len(done))
if __name__ == "__main__": main()
