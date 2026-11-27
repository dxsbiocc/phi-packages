#!/usr/bin/env python3
"""PeptideAtlas builds/datasets -> omics_peptideatlas_staging.db (never omics_v1.db).

Sources:
  - buildInfo TSV: https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/buildInfo?output_mode=tsv
  - Per-organism export pages under https://peptideatlas.org/builds/{organism}/
    for public download pointers (FASTA / TSV dump / MySQL dump) — URLs only.

No bulk download of multi-GB dumps. Resume via sync_cursors.
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)
import csv, io, json, re
from urllib.parse import urljoin
from stage_common import *

DB = ROOT / "data" / "omics_peptideatlas_staging.db"
SRC = "peptideatlas_scan"
BUILDINFO_TSV = "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/buildInfo?output_mode=tsv"
BUILD_DETAILS = "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/buildDetails"
BUILDS_ROOT = "https://peptideatlas.org/builds/"
log = get_logger(SRC)
thr = Throttle(0.5)

# organism slug -> builds page path (probed 200s)
ORG_PAGES = [
    # probed HTTP 200 under peptideatlas.org/builds/{slug}/ (candida/toxoplasma/plasma 404)
    "human", "mouse", "yeast", "arabidopsis", "pig", "chicken", "drosophila",
    "zebrafish", "ecoli", "rice", "maize",
]

DDL = """
CREATE TABLE IF NOT EXISTS src_pa_build(
  build_id TEXT PRIMARY KEY, build_name TEXT, organism TEXT,
  peptide_inclusion_cutoff TEXT, n_experiments INTEGER, spectra_searched TEXT,
  peptide_spectrum_matches TEXT, distinct_peptides TEXT, canonical_proteins TEXT,
  description TEXT, is_liver INTEGER, liver_tags TEXT, details_url TEXT,
  raw_json TEXT, fetched_at TEXT);
CREATE INDEX IF NOT EXISTS idx_src_pa_build_liver ON src_pa_build(is_liver);
CREATE INDEX IF NOT EXISTS idx_src_pa_build_org ON src_pa_build(organism);

CREATE TABLE IF NOT EXISTS src_pa_download(
  download_key TEXT PRIMARY KEY, build_id TEXT, organism_slug TEXT,
  rel_path TEXT, filename TEXT, format TEXT, url TEXT, build_name_hint TEXT,
  is_liver INTEGER, fetched_at TEXT);
CREATE INDEX IF NOT EXISTS idx_src_pa_dl_build ON src_pa_download(build_id);
"""


def liver_tags(*parts) -> list[str]:
    blob = " ".join(str(x) for x in parts if x)
    if not blob:
        return []
    b = blob.lower()
    t: list[str] = []
    if re.search(r"\bhcc\b|hepatocellular", b) or "liver cancer" in b:
        t.append("HCC")
    if "cholangi" in b or "bile duct" in b:
        t.append("CHOL")
    if "masld" in b or "nafld" in b or "nash" in b:
        t.append("MASLD")
    if "cirrhosis" in b:
        t.append("cirrhosis")
    if "hepatitis" in b or re.search(r"\bhbv\b|\bhcv\b", b):
        t.append("hepatitis")
    if re.search(r"\bliver\b", b) or "hepatic" in b or "hepat" in b:
        t.append("liver")
    if is_liver(blob) and not t:
        t.append("liver_related")
    return sorted(set(t))


def _to_int(x):
    if x is None or str(x).strip() == "":
        return None
    try:
        return int(str(x).replace(",", "").strip())
    except ValueError:
        return None


def ingest_buildinfo(c, cur):
    if cur.get("buildinfo_done"):
        return
    text = fetch(BUILDINFO_TSV, log, thr).text
    # strip possible HTML wrapper / BOM
    if "<html" in text[:200].lower():
        raise RuntimeError("buildInfo TSV returned HTML")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    rows = list(reader)
    log.info("buildInfo rows=%d cols=%s", len(rows), list(rows[0].keys()) if rows else None)
    n = 0
    for r in rows:
        bid = str(r.get("Build ID") or "").strip()
        if not bid:
            continue
        name = r.get("Build Name")
        org = r.get("Organism")
        desc = r.get("Description")
        tags = liver_tags(name, org, desc)
        details = f"{BUILD_DETAILS}?atlas_build_id={bid}"
        c.execute(
            "INSERT OR REPLACE INTO src_pa_build VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (
                bid, name, org, r.get("Peptide Inclusion Cutoff"),
                _to_int(r.get("Experiments")),
                r.get("spectra searched"), r.get("peptide spectrum matches"),
                r.get("distinct peptides"), r.get("canonical proteins"),
                (desc or "")[:4000], int(bool(tags)), json.dumps(tags),
                details, json.dumps(r, ensure_ascii=False)[:20000], now_iso(),
            ),
        )
        ds = f"PEPTIDEATLAS:BUILD:{bid}"
        upsert_dataset(c, dict(
            canonical_id=ds, title=name or f"PeptideAtlas build {bid}",
            organism=org, tissue=None, disease_tags=tags,
            assay="PeptideAtlas MS build (reprocessed proteomics)",
            platforms=["PeptideAtlas", "TPP"],
            n_samples=_to_int(r.get("Experiments")),
            access="public", confidence="metadata_only", probe_status="ok",
        ))
        add_member(c, ds, "peptideatlas", bid, "build", details)
        add_artifact(c, dict(
            artifact_id=f"PEPTIDEATLAS:BUILD:{bid}:details",
            dataset_id=ds, kind="processed", role="build_details",
            title=name, format="html", uri=details,
            access="public", source="peptideatlas"))
        add_artifact(c, dict(
            artifact_id=f"PEPTIDEATLAS:BUILD:{bid}:buildinfo_tsv",
            dataset_id=ds, kind="processed", role="buildinfo_row",
            title="buildInfo TSV catalog row", format="tsv",
            uri=BUILDINFO_TSV, access="public", source="peptideatlas"))
        n += 1
    c.commit()
    cur["buildinfo_done"] = True
    cur["builds_stored"] = n
    cur["buildinfo_url"] = BUILDINFO_TSV
    save_cursor(c, SRC, cur)
    log.info("buildInfo DONE n=%d", n)


def _guess_format(filename: str) -> str:
    fn = filename.lower()
    for ext in (".fasta.gz", ".fasta", ".tsv.zip", ".tsv", ".mysql.gz", ".xml.gz",
                ".txt.gz", ".gz", ".zip", ".mysql"):
        if fn.endswith(ext):
            return ext.lstrip(".")
    if "." in filename:
        return filename.rsplit(".", 1)[-1]
    return "unknown"


def _build_id_from_path(path: str) -> str | None:
    m = re.search(r"atlas_build_(\d+)", path)
    if m:
        return m.group(1)
    m = re.search(r"atlas_tables_(\d+)", path)
    if m:
        return m.group(1)
    return None


def ingest_org_pages(c, cur):
    done = set(cur.get("org_pages_done", []))
    for slug in ORG_PAGES:
        if slug in done:
            continue
        url = urljoin(BUILDS_ROOT, slug.strip("/") + "/")
        try:
            from http_util import http_get
            thr()
            r = http_get(url, timeout=(10, 45))
            if r.status_code == 404:
                log.info("org page %s HTTP 404 — skip", slug)
                done.add(slug)
                cur["org_pages_done"] = sorted(done)
                cur[f"org_{slug}_error"] = "404"
                save_cursor(c, SRC, cur)
                continue
            r.raise_for_status()
        except Exception as e:
            log.warning("org page %s fail: %s", slug, e)
            done.add(slug)
            cur["org_pages_done"] = sorted(done)
            cur[f"org_{slug}_error"] = str(e)[:300]
            save_cursor(c, SRC, cur)
            continue
        if len(r.content) < 500:
            log.info("org page %s empty/404 size=%d", slug, len(r.content))
            done.add(slug)
            cur["org_pages_done"] = sorted(done)
            save_cursor(c, SRC, cur)
            continue
        html = r.text
        # relative download hrefs like 202601/atlas_build_607.tsv.zip
        hrefs = re.findall(r'href=["\']([^"\']+)["\']', html, re.I)
        n_dl = 0
        for href in hrefs:
            if href.startswith("#") or href.startswith("javascript:"):
                continue
            # keep build product-looking paths
            if not re.search(
                r"(\.fasta|\.tsv|\.mysql|\.xml|\.gz|\.zip|atlas_build|APD_|Homo_sapiens|coordinate)",
                href, re.I,
            ):
                continue
            full = urljoin(url, href)
            # skip navigational peptideatlas links outside builds tree
            if "peptideatlas.org/builds" not in full and not href.startswith(("20", "ncORF", "PAmap", "APD")):
                # allow relative dated dirs
                if not re.match(r"^\d{6}/", href) and not re.match(r"^[A-Za-z0-9_./-]+\.(fasta|tsv|gz|zip|mysql)", href):
                    continue
                if "peptideatlas.org" in full and "/builds/" not in full:
                    continue
            fname = href.rstrip("/").split("/")[-1]
            if not fname or fname.endswith("/"):
                continue
            bid = _build_id_from_path(href) or _build_id_from_path(fname)
            tags = liver_tags(slug, href, fname)
            key = f"{slug}:{href}"
            c.execute(
                "INSERT OR REPLACE INTO src_pa_download VALUES(?,?,?,?,?,?,?,?,?,?)",
                (key, bid, slug, href, fname, _guess_format(fname), full,
                 None, int(bool(tags)), now_iso()),
            )
            # attach to build dataset if known
            if bid and c.execute(
                "SELECT 1 FROM datasets WHERE canonical_id=?",
                (f"PEPTIDEATLAS:BUILD:{bid}",),
            ).fetchone():
                ds = f"PEPTIDEATLAS:BUILD:{bid}"
            else:
                ds = f"PEPTIDEATLAS:EXPORT:{slug}"
                upsert_dataset(c, dict(
                    canonical_id=ds,
                    title=f"PeptideAtlas builds export page: {slug}",
                    organism=slug, disease_tags=liver_tags(slug),
                    assay="PeptideAtlas build exports (pointers)",
                    platforms=["PeptideAtlas"], access="public",
                    confidence="metadata_only", probe_status="ok",
                ))
                add_member(c, ds, "peptideatlas", slug, "export_page", url)
            add_artifact(c, dict(
                artifact_id=f"PEPTIDEATLAS:DL:{abs(hash(key)) % 10_000_000_000}",
                dataset_id=ds, kind="processed", role="build_export_pointer",
                title=fname, filename=fname, format=_guess_format(fname),
                uri=full, access="public", source="peptideatlas",
            ))
            n_dl += 1
        # page pointer artifact on export dataset
        eds = f"PEPTIDEATLAS:EXPORT:{slug}"
        upsert_dataset(c, dict(
            canonical_id=eds,
            title=f"PeptideAtlas builds export page: {slug}",
            organism=slug, disease_tags=liver_tags(slug),
            assay="PeptideAtlas build exports (pointers)",
            platforms=["PeptideAtlas"], access="public",
            confidence="metadata_only", probe_status="ok",
        ))
        add_artifact(c, dict(
            artifact_id=f"PEPTIDEATLAS:EXPORT:{slug}:page",
            dataset_id=eds, kind="processed", role="export_index",
            title=f"{slug} builds page", format="html", uri=url,
            access="public", source="peptideatlas"))
        c.commit()
        done.add(slug)
        cur["org_pages_done"] = sorted(done)
        cur[f"org_{slug}_downloads"] = n_dl
        save_cursor(c, SRC, cur)
        log.info("org page %s downloads=%d", slug, n_dl)


def main():
    c = connect(DB, DDL)
    cur = load_cursor(c, SRC)
    if cur.get("complete"):
        log.info("already complete; exit")
        return
    cur["buildinfo_tsv"] = BUILDINFO_TSV
    cur["builds_root"] = BUILDS_ROOT
    cur["limitation"] = (
        "No dedicated 'liver' PeptideAtlas build in buildInfo; liver peptides live inside "
        "organism-wide Human/Mouse builds. Liver tags only when name/description matches."
    )
    save_cursor(c, SRC, cur)
    ingest_buildinfo(c, cur)
    ingest_org_pages(c, cur)
    n_b = c.execute("SELECT COUNT(*) FROM src_pa_build").fetchone()[0]
    n_bl = c.execute("SELECT COUNT(*) FROM src_pa_build WHERE is_liver=1").fetchone()[0]
    n_d = c.execute("SELECT COUNT(*) FROM src_pa_download").fetchone()[0]
    cur["complete"] = True
    cur["finished"] = now_iso()
    cur["builds"] = n_b
    cur["builds_liver"] = n_bl
    cur["download_ptrs"] = n_d
    save_cursor(c, SRC, cur)
    log.info("DONE builds=%d liver=%d download_ptrs=%d", n_b, n_bl, n_d)


if __name__ == "__main__":
    main()
