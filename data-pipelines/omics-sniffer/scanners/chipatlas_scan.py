#!/usr/bin/env python3
"""ChIP-Atlas experimentList.tab (streamed metadata) -> omics_chipatlas_staging.db. Pointers to bw/bed only. Idempotent rerun."""
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)
import json
from stage_common import *
from http_util import http_get
DB = ROOT / "data" / "omics_chipatlas_staging.db"; SRC = "chipatlas_scan"; URL = "https://chip-atlas.dbcls.jp/data/metadata/experimentList.tab"
log = get_logger(SRC)
DDL = """CREATE TABLE IF NOT EXISTS src_chipatlas_experiment(srx TEXT, genome TEXT, antigen_class TEXT, antigen TEXT, cell_class TEXT,
 cell_type TEXT, cell_desc TEXT, title TEXT, is_liver INTEGER, raw TEXT, fetched_at TEXT, PRIMARY KEY(srx,genome));"""
ORG = {"hg38":"Homo sapiens","hg19":"Homo sapiens","mm10":"Mus musculus","mm9":"Mus musculus","rn6":"Rattus norvegicus","dm6":"Drosophila melanogaster","dm3":"Drosophila melanogaster","ce11":"Caenorhabditis elegans","ce10":"Caenorhabditis elegans","sacCer3":"Saccharomyces cerevisiae"}
def main():
    c = connect(DB, DDL); cur = load_cursor(c, SRC); n = 0; liv = 0
    r = http_get(URL, stream=True, timeout=(10, 120)); r.raise_for_status()
    for raw in r.iter_lines(decode_unicode=False):
        f = raw.decode("utf-8", "replace").split("\t")
        if len(f) < 9: continue
        srx, g, acl, ant, ccl, ct, cdesc, _logs, title = f[:9]; meta = f[9:]
        lv = int(ccl == "Liver" or bool(re.search(r"hepg2|hep ?3b|huh-?7|hepatocyte|\bliver\b", ct + " " + cdesc, re.I)))
        ds = f"CHIPATLAS:{srx}:{g}"; org = ORG.get(g, g)
        upsert_dataset(c, dict(canonical_id=ds, title=title, organism=org, tissue=f"{ccl} / {ct}", disease_tags=(["liver"] if lv else []),
            assay=f"{acl}" + (f": {ant}" if ant else ""), platforms=[g], n_samples=1, access="public", confidence="metadata_only", probe_status="ok"))
        add_member(c, ds, "sra", f"{srx}@{g}", "experiment", f"https://chip-atlas.org/view?id={srx}")
        base = f"https://chip-atlas.dbcls.jp/data/{g}/eachData"
        add_artifact(c, dict(artifact_id=f"CA:{srx}:{g}:bw", dataset_id=ds, kind="processed", role="coverage", filename=f"{srx}.bw", format="bigWig", uri=f"{base}/bw/{srx}.bw", access="public", source="chipatlas"))
        for q in ("05", "10", "20"):
            add_artifact(c, dict(artifact_id=f"CA:{srx}:{g}:bed{q}", dataset_id=ds, kind="processed", role=f"peaks_q1e-{int(q)}", filename=f"{srx}.{q}.bed", format="bed", uri=f"{base}/bed{q}/{srx}.{q}.bed", access="public", source="chipatlas"))
        c.execute("INSERT OR REPLACE INTO src_chipatlas_experiment VALUES(?,?,?,?,?,?,?,?,?,?,?)", (srx, g, acl, ant, ccl, ct, cdesc, title, lv, json.dumps(meta)[:4000], now_iso()))
        n += 1; liv += lv
        if n % 5000 == 0: cur.update(rows=n, liver=liv); save_cursor(c, SRC, cur); log.info("rows=%d liver=%d", n, liv)
    cur.update(rows=n, liver=liv, complete=True); save_cursor(c, SRC, cur); log.info("DONE rows=%d liver=%d", n, liv)
if __name__ == "__main__": main()
