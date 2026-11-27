#!/usr/bin/env python3
"""cBioPortal study/sample metadata -> omics_cbioportal_staging.db (never omics_v1.db). Metadata only; resumable."""
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)
import json, sys
from stage_common import *
DB = ROOT / "data" / "omics_cbioportal_staging.db"; SRC = "cbioportal_scan"; API = "https://www.cbioportal.org/api"
log = get_logger(SRC); thr = Throttle(0.3)
DDL = """
CREATE TABLE IF NOT EXISTS src_cbioportal_study(study_id TEXT PRIMARY KEY, name TEXT, cancer_type_id TEXT, cancer_type_name TEXT,
 tissue TEXT, reference_genome TEXT, pmid TEXT, all_sample_count INTEGER, is_liver INTEGER, liver_tags TEXT, overlap_tags TEXT,
 molecular_profiles TEXT, sample_lists TEXT, clinical_attributes TEXT, raw_json TEXT, fetched_at TEXT);
CREATE TABLE IF NOT EXISTS src_cbioportal_cancer_type(cancer_type_id TEXT PRIMARY KEY, name TEXT, tissue TEXT, parent TEXT, raw_json TEXT);
"""
def tags(study, ct):
    blob = " ".join([study["studyId"], study.get("name",""), study.get("description",""), ct.get("name",""), ct.get("tissue","")]).lower()
    t = []
    if any(k in blob for k in ("hepatocellular","lihc","hcc","liver cancer","_liver")): t.append("HCC")
    if any(k in blob for k in ("cholangio","chol_","biliary","bile duct")): t.append("CHOL")
    if "hepatoblastoma" in blob: t.append("Hepatoblastoma")
    if "gallbladder" in blob or "gbc" in study["studyId"]: t.append("Gallbladder")
    if ct.get("tissue") in ("Liver","Biliary Tract") or is_liver(blob): t.append("liver_related")
    o = []
    sid = study["studyId"]
    if "_tcga" in sid and sid.split("_")[0] in ("lihc","chol"): o.append("TCGA-"+sid.split("_")[0].upper())
    if "liri" in blob or "icgc" in blob or "riken" in blob: o.append("ICGC/LIRI-JP?")
    return sorted(set(t)), o
def J(path, **p): return fetch(API+path, log, thr, params=p).json()
def main():
    c = connect(DB, DDL); cur = load_cursor(c, SRC); done = set(cur.get("done", []))
    cts = {x["cancerTypeId"]: x for x in J("/cancer-types")}
    for x in cts.values():
        c.execute("INSERT OR REPLACE INTO src_cbioportal_cancer_type VALUES(?,?,?,?,?)", (x["cancerTypeId"], x.get("name"), x.get("dedicatedColor") and x.get("tissue") or x.get("tissue"), x.get("parent"), json.dumps(x)))
    studies = J("/studies", projection="DETAILED"); log.info("studies=%d done=%d", len(studies), len(done))
    for i, s in enumerate(studies):
        sid = s["studyId"]
        if sid in done: continue
        ct = cts.get(s.get("cancerTypeId"), {}); lt, ov = tags(s, ct)
        mps = J(f"/studies/{sid}/molecular-profiles"); sls = J(f"/studies/{sid}/sample-lists")
        cas = J(f"/studies/{sid}/clinical-attributes"); smp = J(f"/studies/{sid}/samples", pageSize=100000)
        org = "Mus musculus" if "mm10" in (s.get("referenceGenome") or "") or "mouse" in sid else "Homo sapiens"
        ds = f"CBIOPORTAL:{sid}"
        upsert_dataset(c, dict(canonical_id=ds, title=s.get("name"), organism=org, tissue=ct.get("tissue"), disease_tags=lt or [ct.get("name")],
            assay="cancer genomics (cBioPortal)", platforms=sorted({m.get("molecularAlterationType") for m in mps}), n_samples=len(smp),
            access="public", confidence="metadata_only", probe_status="ok"))
        add_member(c, ds, "cbioportal", sid, "study", f"https://www.cbioportal.org/study/summary?id={sid}")
        for m in mps:
            add_artifact(c, dict(artifact_id=f"CBIO:{m['molecularProfileId']}", dataset_id=ds, kind="processed", role=m.get("molecularAlterationType"),
                title=m.get("name"), format=m.get("datatype"), uri=f"{API}/molecular-profiles/{m['molecularProfileId']}", access="public", source="cbioportal"))
        add_artifact(c, dict(artifact_id=f"CBIO:{sid}:datahub", dataset_id=ds, kind="processed", role="study_bundle", title="datahub tarball (pointer only)",
            uri=f"https://cbioportal-datahub.s3.amazonaws.com/{sid}.tar.gz", access="public", source="cbioportal"))
        for x in smp:
            upsert_sample(c, dict(sample_id=f"CBIO:{sid}:{x['sampleId']}", dataset_id=ds, title=x["sampleId"], organism=org, tissue=ct.get("tissue"),
                disease=ct.get("name"), condition=x.get("sampleType"), attributes={"patientId": x.get("patientId")}))
        c.execute("INSERT OR REPLACE INTO src_cbioportal_study VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", (sid, s.get("name"), s.get("cancerTypeId"), ct.get("name"),
            ct.get("tissue"), s.get("referenceGenome"), s.get("pmid"), len(smp), int(bool(lt)), json.dumps(lt), json.dumps(ov),
            json.dumps([m["molecularProfileId"] for m in mps]), json.dumps([{"id":l["sampleListId"],"cat":l.get("category"),"n":l.get("sampleCount")} for l in sls]),
            json.dumps([a["clinicalAttributeId"] for a in cas if not a.get("patientAttribute")]), json.dumps(s), now_iso()))
        done.add(sid); cur["done"] = sorted(done); save_cursor(c, SRC, cur)
        if i % 25 == 0: log.info("%d/%d %s samples=%d liver=%s", i, len(studies), sid, len(smp), lt)
    cur["complete"] = True; save_cursor(c, SRC, cur); log.info("DONE studies=%d", len(done))
if __name__ == "__main__": main()
