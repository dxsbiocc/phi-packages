#!/usr/bin/env python3
"""ICGC 25K legacy (release 28) open clinical metadata from the public S3-compatible bucket
s3://icgc25k-open @ object.genomeinformatics.org (official ARGO-documented location).
Fetches ONLY donor / specimen / sample part files (tiny TSV.gz); never exp/ssm/BAM.
datasets = projects, samples = ICGC samples (with donor+specimen fields). access=public_metadata, raw=controlled.
Writes ONLY omics_icgc_staging.db. Resumable by project (source=icgc_scan)."""
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)
import csv, gzip, io, json, re
from stage_common import *
DB = ROOT / "data" / "omics_icgc_staging.db"; SRC = "icgc_scan"
BASE = "https://object.genomeinformatics.org/icgc25k-open"; REL = "release_28"
log = get_logger(SRC); thr = Throttle(0.2)
KINDS = ("donor", "specimen", "sample")
DDL = """CREATE TABLE IF NOT EXISTS src_icgc_project(project_code TEXT PRIMARY KEY, name TEXT, is_liver INTEGER, n_donors INTEGER,
 n_specimens INTEGER, n_samples INTEGER, data_types TEXT, fetched_at TEXT);
CREATE TABLE IF NOT EXISTS src_icgc_donor(icgc_donor_id TEXT PRIMARY KEY, project_code TEXT, raw_json TEXT, fetched_at TEXT);
CREATE TABLE IF NOT EXISTS src_icgc_specimen(icgc_specimen_id TEXT PRIMARY KEY, icgc_donor_id TEXT, project_code TEXT, raw_json TEXT, fetched_at TEXT);
CREATE TABLE IF NOT EXISTS src_icgc_sample(icgc_sample_id TEXT PRIMARY KEY, icgc_specimen_id TEXT, icgc_donor_id TEXT, project_code TEXT, raw_json TEXT, fetched_at TEXT);"""

def s3list(prefix, delim=None):
    keys, pre, tok = [], [], None
    while True:
        p = {"list-type": "2", "prefix": prefix}
        if delim: p["delimiter"] = delim
        if tok: p["continuation-token"] = tok
        t = fetch(BASE, log, thr, params=p).text
        keys += re.findall(r"<Key>([^<]*)</Key>", t); pre += re.findall(r"<Prefix>([^<]*)</Prefix>", t)[1:]
        m = re.search(r"<NextContinuationToken>([^<]*)<", t)
        if not m: return keys, pre
        tok = m.group(1)

def tsv(url, header=None):
    b = gzip.decompress(fetch(url, log, thr).content).decode("utf-8", "replace")
    lines = [l for l in b.splitlines() if l.strip()]
    if header is None: header, lines = lines[0].split("\t"), lines[1:]
    elif lines and lines[0].split("\t")[0] == header[0]: lines = lines[1:]
    return header, [dict(zip(header, l.split("\t"))) for l in lines]

def main():
    c = connect(DB, DDL); cur = load_cursor(c, SRC); done = set(cur.get("done", []))
    hdr = {k: tsv(f"{BASE}/{REL}/headers/{k}.tsv.gz")[0] for k in KINDS}
    names = dict(re.findall(r"- ([A-Z]+-[A-Z]+) (.+)", fetch(f"{BASE}/{REL}/projects_files/README.txt", log, thr).text))
    _, projs = s3list(f"{REL}/data/", "/"); projs = [p.rstrip("/").split("/")[-1] for p in projs]
    projs.sort(key=lambda p: (not is_liver(p, names.get(p)) and not p.startswith("LI"), p))
    log.info("projects=%d done=%d", len(projs), len(done))
    for pc in projs:
        if pc in done: continue
        name = names.get(pc, pc); liver = pc.startswith(("LI", "LIHC")) or is_liver(name)
        keys, _ = s3list(f"{REL}/data/{pc}/")
        dtypes = sorted({k.split("/")[4] for k in keys if k.count("/") >= 5})
        donors, specs, samps = {}, {}, {}
        for k in keys:
            kind = k.split("/")[4] if k.count("/") >= 5 else None
            if kind not in KINDS: continue
            for r in tsv(f"{BASE}/{k}", hdr[kind])[1]:
                if kind == "donor": donors[r["icgc_donor_id"]] = r
                elif kind == "specimen": specs[r["icgc_specimen_id"]] = r
                else: samps[r["icgc_sample_id"]] = r
        ds = f"ICGC:{pc}"; ts = now_iso()
        upsert_dataset(c, dict(canonical_id=ds, title=f"ICGC {pc} — {name}", organism="Homo sapiens",
            disease_tags=(["liver"] if liver else []) + [name], assay="cancer genomics (ICGC 25K release 28)", platforms=dtypes,
            n_samples=len(samps), access="public_metadata", license="ICGC open tier; raw=controlled (DACO/EGA/dbGaP)",
            confidence="metadata_only", probe_status="ok"))
        add_member(c, ds, "icgc", pc, "project", f"{BASE}/{REL}/data/{pc}/")
        add_artifact(c, dict(artifact_id=f"{ds}:raw", dataset_id=ds, kind="raw", role="sequencing", title="Raw/controlled data via DACO (EGA EGAC00001000010 / dbGaP phs000178)",
            uri="https://docs.icgc-argo.org/docs/data-access/icgc-25k-data", access="controlled", source="icgc"))
        for d in donors.values(): c.execute("INSERT OR REPLACE INTO src_icgc_donor VALUES(?,?,?,?)", (d["icgc_donor_id"], pc, json.dumps(d), ts))
        for s in specs.values(): c.execute("INSERT OR REPLACE INTO src_icgc_specimen VALUES(?,?,?,?,?)", (s["icgc_specimen_id"], s.get("icgc_donor_id"), pc, json.dumps(s), ts))
        for sid, s in samps.items():
            sp = specs.get(s.get("icgc_specimen_id"), {}); dn = donors.get(s.get("icgc_donor_id"), {})
            c.execute("INSERT OR REPLACE INTO src_icgc_sample VALUES(?,?,?,?,?,?)", (sid, s.get("icgc_specimen_id"), s.get("icgc_donor_id"), pc, json.dumps(s), ts))
            attrs = {k: v for k, v in {**{f"donor.{k}": v for k, v in dn.items()}, **{f"specimen.{k}": v for k, v in sp.items()}, **s}.items() if v}
            upsert_sample(c, dict(sample_id=f"ICGC:{sid}", dataset_id=ds, title=s.get("submitted_sample_id"), organism="Homo sapiens",
                tissue=sp.get("specimen_type"), disease=name, condition=sp.get("tumour_histological_type") or None,
                sex=dn.get("donor_sex"), age=dn.get("donor_age_at_diagnosis"), attributes=attrs))
        c.execute("INSERT OR REPLACE INTO src_icgc_project VALUES(?,?,?,?,?,?,?,?)", (pc, name, int(liver), len(donors), len(specs), len(samps), json.dumps(dtypes), ts))
        done.add(pc); save_cursor(c, SRC, {"done": sorted(done)})
        log.info("%s liver=%s donors=%d specimens=%d samples=%d (%d/%d)", pc, liver, len(donors), len(specs), len(samps), len(done), len(projs))
    save_cursor(c, SRC, {"done": sorted(done), "finished": now_iso()}); log.info("FINISHED")
if __name__ == "__main__": main()
