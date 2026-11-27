#!/usr/bin/env python3
"""FinnGen public phenotype catalog (latest reachable release, prefer R13)
-> omics_finngen_staging.db (never omics_v1.db).

Metadata + sumstats URL pointers from public manifest only; does NOT download
.gz summary-stat matrices.
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)
import csv
import io
import json
import re
from stage_common import *

DB = ROOT / "data" / "omics_finngen_staging.db"
SRC = "finngen_scan"
log = get_logger(SRC)
thr = Throttle(0.2)

# Prefer DF13 / R13 browser; fall back to R12.
RELEASE_CANDIDATES = [
    {
        "release": "R13",
        "base": "https://r13.finngen.fi",
        "manifest": "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv",
    },
    {
        "release": "R12",
        "base": "https://r12.finngen.fi",
        "manifest": "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv",
    },
]

DDL = """
CREATE TABLE IF NOT EXISTS src_finngen_pheno(
  phenocode TEXT PRIMARY KEY,
  release TEXT,
  phenostring TEXT,
  category TEXT,
  num_cases INTEGER,
  num_controls INTEGER,
  num_gw_significant INTEGER,
  gc_lambda TEXT,
  portal_url TEXT,
  sumstats_https TEXT,
  sumstats_gs TEXT,
  access_note TEXT,
  is_liver INTEGER,
  liver_tags TEXT,
  raw_json TEXT,
  fetched_at TEXT
);
CREATE INDEX IF NOT EXISTS idx_src_finngen_liver ON src_finngen_pheno(is_liver);
CREATE INDEX IF NOT EXISTS idx_src_finngen_release ON src_finngen_pheno(release);
"""


def liver_tags_for(phenostring: str, phenocode: str, category: str) -> list[str]:
    blob = f"{phenostring or ''} {phenocode or ''} {category or ''}"
    b = blob.lower()
    t: list[str] = []
    if re.search(r"\bhcc\b|hepatocellular", b) or "liver cancer" in b or "hepatocellu" in b:
        t.append("HCC")
    if "cholangi" in b or "bile duct" in b or re.search(r"\bbiliary\b", b) or "gallbladder" in b or "bile" in b:
        t.append("CHOL")
    if "masld" in b or "nafld" in b or "nash" in b or "fatty liver" in b:
        t.append("MASLD")
    if "cirrhosis" in b or "cirrhotic" in b or "chirbil" in b or re.search(r"\bcirrh", b):
        t.append("cirrhosis")
    if "hepatitis" in b or re.search(r"\bhbv\b|\bhcv\b", b):
        t.append("hepatitis")
    if "liver transplantation" in b or "liver transplant" in b or (
        "transplant" in b and "liver" in b
    ):
        t.append("transplant-liver")
    if re.search(r"\bliver\b", b) or "hepatic" in b or "hepat" in b:
        t.append("liver")
    if is_liver(blob) and not t:
        t.append("liver_related")
    return sorted(set(t))


def pick_release() -> dict:
    last_err = None
    for cand in RELEASE_CANDIDATES:
        try:
            r = fetch(f"{cand['base']}/api/phenos", log, thr,
                      headers={"Accept": "application/json"})
            phenos = r.json()
            if isinstance(phenos, list) and phenos:
                # try manifest (ok if missing)
                man = None
                try:
                    mr = fetch(cand["manifest"], log, thr, ok404=True)
                    if mr is not None and mr.status_code == 200:
                        man = mr.text
                except Exception as e:
                    log.warning("manifest fetch failed for %s: %s", cand["release"], e)
                return {**cand, "phenos": phenos, "manifest_text": man}
        except Exception as e:
            last_err = e
            log.warning("release %s unavailable: %s", cand["release"], e)
    raise RuntimeError(f"no FinnGen release reachable; last={last_err}")


def parse_manifest(text: str | None) -> dict[str, dict]:
    out: dict[str, dict] = {}
    if not text:
        return out
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    for row in reader:
        code = (row.get("phenocode") or "").strip()
        if code:
            out[code] = row
    return out


def store_pheno(c, release: str, base: str, p: dict, man: dict | None):
    code = p.get("phenocode")
    if not code:
        return False
    phenostring = p.get("phenostring") or (man or {}).get("phenotype") or code
    category = p.get("category") or (man or {}).get("category")
    tags = liver_tags_for(phenostring, code, category or "")
    portal = f"{base}/pheno/{code}"
    sum_https = (man or {}).get("path_https") or None
    sum_gs = (man or {}).get("path_bucket") or None
    # Honest access: browser is public; bulk GCS often needs FinnGen form
    # registration even when object URLs are listed publicly.
    if sum_https:
        access = "registration_for_bulk"
        access_note = (
            "Browser phenotype page is public; GCS sumstats URLs from public "
            "manifest — bulk download may require FinnGen data download registration."
        )
    else:
        access = "public_browser"
        access_note = "Public PheWeb phenotype card; no sumstats URI in manifest."
    n_cases = p.get("num_cases")
    try:
        n_cases_i = int(n_cases) if n_cases is not None else None
    except (TypeError, ValueError):
        n_cases_i = None
    n_ctrl = p.get("num_controls")
    try:
        n_ctrl_i = int(n_ctrl) if n_ctrl is not None else None
    except (TypeError, ValueError):
        n_ctrl_i = None
    n_gw = p.get("num_gw_significant")
    try:
        n_gw_i = int(n_gw) if n_gw is not None else None
    except (TypeError, ValueError):
        n_gw_i = None
    ds = f"FINNGEN:{release}:{code}"
    upsert_dataset(c, dict(
        canonical_id=ds,
        title=f"FinnGen {release} {code}: {phenostring}",
        organism="Homo sapiens",
        tissue=None,
        disease_tags=tags or ([category] if category else []),
        assay=f"GWAS summary statistics (FinnGen {release})",
        platforms=["FinnGen"],
        n_samples=(n_cases_i + n_ctrl_i) if (n_cases_i is not None and n_ctrl_i is not None) else n_cases_i,
        access=access,
        confidence="metadata_only",
        probe_status="ok",
    ))
    add_member(c, ds, "finngen", f"{release}:{code}", "phenotype", portal)
    add_artifact(c, dict(
        artifact_id=f"FINNGEN:{release}:{code}:portal", dataset_id=ds,
        kind="processed", role="phenotype_portal", title=phenostring,
        format="html", uri=portal, access="public", source="finngen"))
    add_artifact(c, dict(
        artifact_id=f"FINNGEN:{release}:{code}:api", dataset_id=ds,
        kind="processed", role="api_pheno_list_record",
        title="FinnGen /api/phenos list record (bulk)",
        format="json", uri=f"{base}/api/phenos", access="public", source="finngen"))
    if sum_https:
        add_artifact(c, dict(
            artifact_id=f"FINNGEN:{release}:{code}:sumstats", dataset_id=ds,
            kind="processed", role="summary_statistics",
            title=f"FinnGen {release} {code} sumstats (pointer; not downloaded)",
            filename=sum_https.rstrip("/").split("/")[-1],
            format="gz", uri=sum_https, access=access, source="finngen"))
    if sum_gs:
        add_artifact(c, dict(
            artifact_id=f"FINNGEN:{release}:{code}:sumstats_gs", dataset_id=ds,
            kind="processed", role="summary_statistics_gs",
            title=f"FinnGen {release} {code} gs:// pointer",
            format="gs", uri=sum_gs, access=access, source="finngen"))
    raw = {k: v for k, v in p.items() if k != "assoc_files"}  # drop local paths
    if man:
        raw["manifest"] = {
            k: man.get(k) for k in ("path_https", "path_bucket", "num_cases", "num_controls", "phenotype", "category")
            if man.get(k) is not None
        }
    c.execute(
        "INSERT OR REPLACE INTO src_finngen_pheno VALUES("
        "?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
        (
            code,
            release,
            phenostring,
            category,
            n_cases_i,
            n_ctrl_i,
            n_gw_i,
            json.dumps(p.get("gc_lambda"), ensure_ascii=False) if p.get("gc_lambda") is not None else None,
            portal,
            sum_https,
            sum_gs,
            access_note,
            int(bool(tags)),
            json.dumps(tags, ensure_ascii=False),
            json.dumps(raw, ensure_ascii=False)[:20000],
            now_iso(),
        ),
    )
    return True


def main():
    c = connect(DB, DDL)
    cur = load_cursor(c, SRC)
    if cur.get("complete"):
        log.info("already complete release=%s; exit", cur.get("release"))
        return

    picked = pick_release()
    release = picked["release"]
    base = picked["base"]
    phenos = picked["phenos"]
    man_map = parse_manifest(picked.get("manifest_text"))
    log.info(
        "release=%s phenos=%d manifest_rows=%d base=%s",
        release, len(phenos), len(man_map), base,
    )
    cur.update({
        "release": release,
        "base": base,
        "manifest": picked["manifest"],
        "total": len(phenos),
        "manifest_rows": len(man_map),
    })
    save_cursor(c, SRC, cur)

    done = set(cur.get("done", []))
    scanned = int(cur.get("scanned", 0))
    liver = int(cur.get("liver", 0))
    errors = int(cur.get("errors", 0))

    for i, p in enumerate(phenos):
        code = p.get("phenocode")
        if not code or code in done:
            continue
        try:
            if store_pheno(c, release, base, p, man_map.get(code)):
                scanned += 1
                if liver_tags_for(
                    p.get("phenostring") or "", code, p.get("category") or ""
                ):
                    liver += 1
                done.add(code)
        except Exception as e:
            errors += 1
            log.warning("pheno %s error: %s", code, e)
        if (i + 1) % 200 == 0:
            cur.update({
                "done": sorted(done),
                "scanned": scanned,
                "liver": liver,
                "errors": errors,
                "index": i + 1,
            })
            save_cursor(c, SRC, cur)
            log.info("progress %d/%d scanned=%d liver=%d errors=%d",
                     i + 1, len(phenos), scanned, liver, errors)

    c.commit()
    cur.update({
        "complete": True,
        "finished_at": now_iso(),
        "done": sorted(done),
        "scanned": scanned,
        "liver": liver,
        "errors": errors,
        "note": (
            f"Used FinnGen {release} public /api/phenos + manifest pointers. "
            "Sumstats .gz files NOT downloaded; access=registration_for_bulk when "
            "GCS HTTPS URI present."
        ),
    })
    save_cursor(c, SRC, cur)
    n_db = c.execute("SELECT COUNT(*) FROM src_finngen_pheno").fetchone()[0]
    n_liv = c.execute("SELECT COUNT(*) FROM src_finngen_pheno WHERE is_liver=1").fetchone()[0]
    n_ss = c.execute(
        "SELECT COUNT(*) FROM src_finngen_pheno WHERE sumstats_https IS NOT NULL"
    ).fetchone()[0]
    log.info(
        "DONE release=%s phenos=%d liver=%d sumstats_ptr=%d errors=%d",
        release, n_db, n_liv, n_ss, errors,
    )


if __name__ == "__main__":
    main()
