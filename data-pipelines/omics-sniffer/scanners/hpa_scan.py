#!/usr/bin/env python3
"""Human Protein Atlas gene/protein metadata -> omics_hpa_staging.db (never omics_v1.db).

Downloads catalog TSV (proteinatlas.tsv.zip) only. Stores gene rows + pointers to
other downloadable TSVs/JSON/XML. Does NOT download image banks or the large
rna_celline matrix (~200MB) — pointer only.
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)
import csv, io, json, zipfile
from stage_common import *

DB = ROOT / "data" / "omics_hpa_staging.db"
SRC = "hpa_scan"
BASE = "https://www.proteinatlas.org"
log = get_logger(SRC)
thr = Throttle(0.5)

DDL = """
CREATE TABLE IF NOT EXISTS src_hpa_gene(
  ensembl TEXT PRIMARY KEY, gene TEXT, gene_synonym TEXT, uniprot TEXT,
  description TEXT, protein_class TEXT, rna_tissue_specificity TEXT,
  rna_tissue_specific_ntpm TEXT, rna_celltype_specificity TEXT,
  rna_celltype_specific_ncpm TEXT, protein_tissue_specificity TEXT,
  hcc_prognostic_tcga TEXT, hcc_prognostic_validation TEXT,
  is_liver INTEGER, liver_tags TEXT, raw_json TEXT, fetched_at TEXT);
CREATE TABLE IF NOT EXISTS src_hpa_download(
  name TEXT PRIMARY KEY, uri TEXT, kind TEXT, note TEXT, fetched_at TEXT);
"""

# Current v25.1 catalog + archived v23 specialized TSVs (pointers; only catalog is fetched)
DOWNLOADS = [
    ("proteinatlas.tsv.zip", f"{BASE}/download/proteinatlas.tsv.zip", "catalog", "fetched into src_hpa_gene"),
    ("proteinatlas.json.gz", f"{BASE}/download/proteinatlas.json.gz", "catalog", "pointer only"),
    ("proteinatlas.xml.gz", f"{BASE}/download/proteinatlas.xml.gz", "catalog", "pointer only; large"),
    ("normal_tissue.tsv.zip", "https://v23.proteinatlas.org/download/normal_tissue.tsv.zip", "tissue_ihc", "v23 archive; not on current /download"),
    ("pathology.tsv.zip", "https://v23.proteinatlas.org/download/pathology.tsv.zip", "pathology", "v23 archive"),
    ("rna_tissue_consensus.tsv.zip", "https://v23.proteinatlas.org/download/rna_tissue_consensus.tsv.zip", "rna_tissue", "v23 archive"),
    ("rna_celline.tsv.zip", "https://v23.proteinatlas.org/download/rna_celline.tsv.zip", "rna_celline", "v23 archive ~200MB; pointer only, NOT downloaded"),
]

def _is_true_prognostic(val):
    """HPA uses 'validated/potential prognostic ...' vs 'unprognostic (...)'."""
    if not val:
        return False
    v = val.strip().lower()
    if v.startswith("unprognostic"):
        return False
    return "prognostic" in v

def gene_tags(row):
    t = []
    rna_ts = row.get("RNA tissue specificity") or ""
    rna_ntpm = row.get("RNA tissue specific nTPM") or ""
    rna_ct = row.get("RNA single cell type specificity") or ""
    rna_ct_ncpm = row.get("RNA single cell type specific nCPM") or ""
    prot_ts = row.get("Protein tissue specificity") or ""
    hcc = row.get("Cancer prognostics - Liver Hepatocellular Carcinoma (TCGA)") or ""
    hcc_v = row.get("Cancer prognostics - Liver Hepatocellular Carcinoma (validation)") or ""
    tissue_blob = " ".join([rna_ts, rna_ntpm, rna_ct, rna_ct_ncpm, prot_ts])
    tl = tissue_blob.lower()
    if re.search(r"\bliver\b", tl):
        t.append("liver")
    if "hepatocyte" in tl:
        t.append("hepatocyte")
    if _is_true_prognostic(hcc):
        t.append("HCC_prognostic")
    if _is_true_prognostic(hcc_v):
        t.append("HCC_prognostic_validation")
    if t and ("HCC_prognostic" in t or "HCC_prognostic_validation" in t):
        t.append("HCC_pathology")
    return sorted(set(t))

def main():
    c = connect(DB, DDL)
    cur = load_cursor(c, SRC)
    if cur.get("complete"):
        log.info("already complete; exit"); return

    upsert_dataset(c, dict(
        canonical_id="HPA:CATALOG", title="Human Protein Atlas downloadable catalogs",
        organism="Homo sapiens", tissue=None, disease_tags=[], assay="HPA resource",
        platforms=["HPA"], n_samples=0, access="public", confidence="metadata_only",
        probe_status="ok"))
    # register download pointers
    for name, uri, kind, note in DOWNLOADS:
        c.execute("INSERT OR REPLACE INTO src_hpa_download VALUES(?,?,?,?,?)",
                  (name, uri, kind, note, now_iso()))
        add_artifact(c, dict(
            artifact_id=f"HPA:download:{name}", dataset_id="HPA:CATALOG",
            kind="processed", role=kind, title=name, filename=name,
            format=name.split(".")[-1], uri=uri, access="public", source="hpa"))
    c.commit()

    log.info("fetching proteinatlas.tsv.zip ...")
    url = f"{BASE}/download/proteinatlas.tsv.zip"
    zbytes = None
    for a in range(1, MAX_RETRIES + 1):
        thr()
        try:
            r = http_get(url, timeout=(10, 300)); r.raise_for_status(); zbytes = r.content; break
        except Exception as e:
            log.warning("retry %d/%d %s: %s", a, MAX_RETRIES, url, e); time.sleep(min(60, 2 ** a))
    if zbytes is None:
        raise RuntimeError("failed to download proteinatlas.tsv.zip")
    zf = zipfile.ZipFile(io.BytesIO(zbytes))
    tsv_name = zf.namelist()[0]
    log.info("zip entry %s", tsv_name)
    n = 0
    n_liver = 0
    with zf.open(tsv_name) as fh:
        text = io.TextIOWrapper(fh, encoding="utf-8")
        reader = csv.DictReader(text, delimiter="\t")
        for row in reader:
            ens = row.get("Ensembl") or ""
            gene = row.get("Gene") or ""
            if not ens:
                continue
            lt = gene_tags(row)
            if lt:
                n_liver += 1
            ds = f"HPA:{ens}"
            upsert_dataset(c, dict(
                canonical_id=ds, title=f"{gene} ({ens})", organism="Homo sapiens",
                tissue="liver" if "liver" in lt or "hepatocyte" in lt else None,
                disease_tags=lt, assay="HPA gene annotation",
                platforms=["IHC", "RNAseq", "MS"], n_samples=1,
                access="public", confidence="metadata_only", probe_status="ok"))
            add_member(c, ds, "hpa", ens, "gene", f"{BASE}/{ens}")
            add_artifact(c, dict(
                artifact_id=f"HPA:{ens}:json", dataset_id=ds, kind="processed",
                role="gene_json", title=f"{ens}.json", format="json",
                uri=f"{BASE}/{ens}.json", access="public", source="hpa"))
            add_artifact(c, dict(
                artifact_id=f"HPA:{ens}:tsv", dataset_id=ds, kind="processed",
                role="gene_tsv", title=f"{ens}.tsv", format="tsv",
                uri=f"{BASE}/{ens}.tsv", access="public", source="hpa"))
            # image bank pointer (do not download)
            add_artifact(c, dict(
                artifact_id=f"HPA:{ens}:entry", dataset_id=ds, kind="processed",
                role="entry_page_images_not_fetched", title=f"{gene} entry",
                format="html", uri=f"{BASE}/{ens}", access="public", source="hpa"))
            upsert_sample(c, dict(
                sample_id=f"HPA:{ens}:gene", dataset_id=ds, title=gene,
                organism="Homo sapiens",
                tissue=row.get("RNA tissue specificity"),
                disease=("HCC" if "HCC_pathology" in lt else None),
                attributes={
                    "uniprot": row.get("Uniprot"),
                    "rna_tissue_specific_nTPM": row.get("RNA tissue specific nTPM"),
                    "rna_celltype": row.get("RNA single cell type specific nCPM"),
                    "hcc_tcga": row.get("Cancer prognostics - Liver Hepatocellular Carcinoma (TCGA)"),
                }))
            slim = {k: row.get(k) for k in (
                "Gene", "Ensembl", "Uniprot", "Gene description", "Protein class",
                "RNA tissue specificity", "RNA tissue specific nTPM",
                "RNA single cell type specificity", "RNA single cell type specific nCPM",
                "Protein tissue specificity",
                "Cancer prognostics - Liver Hepatocellular Carcinoma (TCGA)",
                "Cancer prognostics - Liver Hepatocellular Carcinoma (validation)",
            )}
            c.execute(
                "INSERT OR REPLACE INTO src_hpa_gene VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (ens, gene, row.get("Gene synonym"), row.get("Uniprot"),
                 row.get("Gene description"), row.get("Protein class"),
                 row.get("RNA tissue specificity"), row.get("RNA tissue specific nTPM"),
                 row.get("RNA single cell type specificity"),
                 row.get("RNA single cell type specific nCPM"),
                 row.get("Protein tissue specificity"),
                 row.get("Cancer prognostics - Liver Hepatocellular Carcinoma (TCGA)"),
                 row.get("Cancer prognostics - Liver Hepatocellular Carcinoma (validation)"),
                 int(bool(lt)), json.dumps(lt), json.dumps(slim), now_iso()))
            n += 1
            if n % 2000 == 0:
                save_cursor(c, SRC, {"n": n, "n_liver": n_liver})
                log.info("genes=%d liver_tagged=%d", n, n_liver)
    cur.update({"complete": True, "n": n, "n_liver": n_liver, "version_note": "HPA v25.1 proteinatlas.tsv"})
    save_cursor(c, SRC, cur)
    log.info("DONE genes=%d liver_tagged=%d", n, n_liver)

if __name__ == "__main__":
    main()
