#!/usr/bin/env python3
"""DepMap public release catalog (figshare mirror of latest reachable public release).
samples = cell lines from Model.csv (small annotation file); datasets = one per release data file
with URL/size pointers ONLY (no matrix downloads). Writes ONLY omics_depmap_staging.db."""
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)
import csv, io, json, re
from stage_common import *
DB = ROOT / "data" / "omics_depmap_staging.db"; SRC = "depmap_scan"
log = get_logger(SRC); thr = Throttle(0.5)
RELEASE_ARTICLE = 27993248   # "DepMap 24Q4 Public" (newest DepMap Public release on figshare; portal API is bot-gated)
DDL = """CREATE TABLE IF NOT EXISTS src_depmap_file(file_id TEXT PRIMARY KEY, release TEXT, name TEXT, size_bytes INTEGER,
 url TEXT, md5 TEXT, product TEXT, fetched_at TEXT);
CREATE TABLE IF NOT EXISTS src_depmap_model(model_id TEXT PRIMARY KEY, release TEXT, is_liver INTEGER, raw_json TEXT, fetched_at TEXT);"""
def product(n):
    for pat, p in [(r"^CRISPR|^Achilles|^Avana|^KY|^Humagne|^Screen", "CRISPR gene effect / screens"),
                   (r"^OmicsExpression", "Expression (RNA-seq)"), (r"CN|LoH", "Copy number"),
                   (r"Mutation|\.maf$|Signatures|Microsatellite", "Mutations"), (r"Fusion", "Fusions"),
                   (r"^Model|^Omics.*Profiles|^OmicsProfiles", "Model annotation"), (r"Compound", "Compounds")]:
        if re.search(pat, n): return p
    return "Other"
def main():
    c = connect(DB, DDL)
    art = fetch(f"https://api.figshare.com/v2/articles/{RELEASE_ARTICLE}", log, thr).json()
    rel = re.search(r"\d\dQ\d", art["title"]).group(0); cat = f"DEPMAP:{rel}:models"
    upsert_dataset(c, dict(canonical_id=cat, title=f"DepMap {rel} Public — cell line models (Model.csv)", organism="Homo sapiens",
        assay="cell line annotation", access="public", license="CC BY 4.0", confidence="metadata_only", probe_status="ok"))
    add_member(c, cat, "figshare", str(RELEASE_ARTICLE), "release", art.get("url_public_html"))
    for f in art["files"]:
        p = product(f["name"]); ds = f"DEPMAP:{rel}:{f['name'].rsplit('.',1)[0]}"
        c.execute("INSERT OR REPLACE INTO src_depmap_file VALUES(?,?,?,?,?,?,?,?)", (str(f["id"]), rel, f["name"], f["size"], f["download_url"], f.get("computed_md5"), p, now_iso()))
        upsert_dataset(c, dict(canonical_id=ds, title=f"DepMap {rel} {f['name']}", organism="Homo sapiens", assay=p,
            access="public", license="CC BY 4.0", confidence="metadata_only", probe_status="ok"))
        add_member(c, ds, "depmap_file", f"{rel}/{f['name']}", "file", f["download_url"])
        add_artifact(c, dict(artifact_id=f"DEPMAP:{rel}:file:{f['id']}", dataset_id=ds, kind="processed", role=p, title=f["name"],
            filename=f["name"], format=f["name"].rsplit(".",1)[-1], uri=f["download_url"], size_bytes=f["size"], access="public", source="depmap"))
    c.commit()
    mf = next(f for f in art["files"] if f["name"] == "Model.csv")
    rows = list(csv.DictReader(io.StringIO(fetch(mf["download_url"], log, thr).text)))
    nl = 0
    for r in rows:
        lv = is_liver(r.get("OncotreeLineage"), r.get("OncotreePrimaryDisease"), r.get("OncotreeSubtype")); nl += lv
        c.execute("INSERT OR REPLACE INTO src_depmap_model VALUES(?,?,?,?,?)", (r["ModelID"], rel, int(lv), json.dumps(r, ensure_ascii=False), now_iso()))
        upsert_sample(c, dict(sample_id=f"DEPMAP:{r['ModelID']}", dataset_id=cat, title=r.get("CellLineName"), organism="Homo sapiens",
            tissue=r.get("OncotreeLineage"), disease=r.get("OncotreePrimaryDisease"), condition=r.get("OncotreeSubtype"), sex=r.get("Sex"), age=r.get("Age"),
            attributes={k: r.get(k) for k in ("CCLEName","RRID","OncotreeCode","PrimaryOrMetastasis","SampleCollectionSite","SourceType","ModelType","DepmapModelType","SangerModelID","COSMICID","PatientRace")} | {"is_liver": lv}))
    c.execute("UPDATE datasets SET n_samples=? WHERE canonical_id=?", (len(rows), cat))
    save_cursor(c, SRC, {"release": rel, "article": RELEASE_ARTICLE, "files": len(art["files"]), "models": len(rows), "liver": nl, "finished": now_iso()})
    log.info("FINISHED release=%s files=%d models=%d liver=%d", rel, len(art["files"]), len(rows), nl)
if __name__ == "__main__": main()
