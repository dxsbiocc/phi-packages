#!/usr/bin/env python3
"""BioGRID ORCS CRISPR screen catalog -> omics_biogrid_orcs_staging.db (never omics_v1.db).

Pulls SCREEN_INDEX (*.index.tab.txt) from each Latest-Release organism tar.gz by
streaming the archive and stopping after the index member (does NOT retain
genome-wide .screen.tab.txt score matrices). Also records tar download URL
pointers as artifacts.

Release probed: BioGRID ORCS 2.0.18 (Sept 2025) at downloads.thebiogrid.org.
REST webservice (orcsws.thebiogrid.org) requires a free accesskey — not used here.
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_PIPELINE_DIR = str(_Path(__file__).resolve().parent.parent / "pipeline")
if _PIPELINE_DIR not in _sys.path:
    _sys.path.insert(0, _PIPELINE_DIR)
import csv
import gzip
import io
import json
import re
import tarfile
from stage_common import *

DB = ROOT / "data" / "omics_biogrid_orcs_staging.db"
SRC = "biogrid_orcs_scan"
RELEASE = "2.0.18"
LATEST = "https://downloads.thebiogrid.org/Download/BioGRID-ORCS/Latest-Release"
PORTAL = "https://orcs.thebiogrid.org"
DOWNLOADS = "https://downloads.thebiogrid.org/BioGRID-ORCS/Latest-Release/"
log = get_logger(SRC)
thr = Throttle(0.25)

# (label, filename, expected organism hint)
ORGS = [
    ("homo_sapiens", "BIOGRID-ORCS-ALL-homo_sapiens-LATEST.screens.tar.gz"),
    ("mus_musculus", "BIOGRID-ORCS-ALL-mus_musculus-LATEST.screens.tar.gz"),
    ("saccharomyces_cerevisiae_S288C", "BIOGRID-ORCS-ALL-saccharomyces_cerevisiae_S288C-LATEST.screens.tar.gz"),
    ("drosophila_melanogaster", "BIOGRID-ORCS-ALL-drosophila_melanogaster-LATEST.screens.tar.gz"),
    ("chlorocebus_sabaeus", "BIOGRID-ORCS-ALL-chlorocebus_sabaeus-LATEST.screens.tar.gz"),
]

DDL = """
CREATE TABLE IF NOT EXISTS src_orcs_screen(
  screen_id TEXT PRIMARY KEY,
  source_id TEXT,
  source_type TEXT,
  author TEXT,
  screen_name TEXT,
  scores_size INTEGER,
  full_size INTEGER,
  full_size_available TEXT,
  number_of_hits INTEGER,
  analysis TEXT,
  significance_indicator TEXT,
  significance_criteria TEXT,
  throughput TEXT,
  screen_type TEXT,
  screen_format TEXT,
  experimental_setup TEXT,
  duration TEXT,
  condition_name TEXT,
  condition_dosage TEXT,
  moi TEXT,
  library TEXT,
  library_type TEXT,
  library_methodology TEXT,
  enzyme TEXT,
  cell_line TEXT,
  cell_type TEXT,
  phenotype TEXT,
  organism_id TEXT,
  organism_official TEXT,
  notes TEXT,
  screen_rationale TEXT,
  source_db TEXT,
  release TEXT,
  organism_slug TEXT,
  portal_url TEXT,
  score_file_pointer TEXT,
  is_liver INTEGER,
  liver_tags TEXT,
  raw_json TEXT,
  fetched_at TEXT
);
CREATE INDEX IF NOT EXISTS idx_src_orcs_liver ON src_orcs_screen(is_liver);
CREATE INDEX IF NOT EXISTS idx_src_orcs_org ON src_orcs_screen(organism_slug);
CREATE INDEX IF NOT EXISTS idx_src_orcs_lib ON src_orcs_screen(library_type);

CREATE TABLE IF NOT EXISTS src_orcs_release_file(
  filename TEXT PRIMARY KEY,
  organism_slug TEXT,
  url TEXT,
  size_bytes INTEGER,
  index_bytes INTEGER,
  bytes_streamed INTEGER,
  n_screens INTEGER,
  fetched_at TEXT
);
"""

# Cell-line / tissue patterns; avoid bare HCC#### breast/lung cell-line false friends.
LIVER_CELL_RE = re.compile(
    r"hepg2|huh-?7|heparg|hep3b|sk-?hep|plc/?prf|huh-?6|li-?7|egi-?1|"
    r"hepatocyte|hepatoma|hepatoblast|hepatocellular|"
    r"liver[- ]derived|liver cell|\bliver\b|"
    r"cholangiocarcin|\bcholangi|"
    r"bile[- ]duct|biliary",
    re.I,
)


def _dash(v):
    if v is None:
        return None
    s = str(v).strip()
    return None if s in ("", "-") else s


def _int(v):
    s = _dash(v)
    if s is None:
        return None
    try:
        return int(str(s).replace(",", ""))
    except ValueError:
        return None


def liver_tags_for(row: dict) -> list[str]:
    cell = row.get("CELL_LINE") or ""
    ctype = row.get("CELL_TYPE") or ""
    pheno = row.get("PHENOTYPE") or ""
    notes = row.get("NOTES") or ""
    rationale = row.get("SCREEN_RATIONALE") or ""
    cond = row.get("CONDITION_NAME") or ""
    blob = " ".join([cell, ctype, pheno, notes, rationale, cond])
    b = blob.lower()
    t: list[str] = []

    # Reject HCC#### cell lines that are not liver-annotated (e.g. HCC364 lung)
    if re.match(r"^HCC-?\d+", cell.strip(), re.I) and not LIVER_CELL_RE.search(ctype + " " + notes + " " + rationale):
        # still allow explicit liver disease tags elsewhere
        pass
    else:
        if LIVER_CELL_RE.search(cell) or LIVER_CELL_RE.search(ctype):
            t.append("liver")

    if re.search(r"\bhcc\b|hepatocellular", b) or "liver cancer" in b:
        t.append("HCC")
    if "cholangi" in b or "bile duct" in b or re.search(r"\bbiliary\b", b):
        t.append("CHOL")
    if "masld" in b or "nafld" in b or "nash" in b or "fatty liver" in b:
        t.append("MASLD")
    if "cirrhosis" in b or "cirrhotic" in b:
        t.append("cirrhosis")
    if "hepatitis" in b or re.search(r"\bhbv\b|\bhcv\b", b):
        t.append("hepatitis")
    if "hepatoblastoma" in b:
        t.append("hepatoblastoma")
    if re.search(r"\bliver\b", b) or "hepatic" in b or "hepatoma" in b or "hepatocyte" in b:
        t.append("liver")
    # stage_common catch-all, but re-check HCC#### FP
    if is_liver(blob) and not t:
        if re.match(r"^HCC-?\d+", cell.strip(), re.I) and not re.search(
            r"liver|hepat", ctype + " " + notes, re.I
        ):
            pass
        else:
            t.append("liver_related")
    return sorted(set(t))


def stream_index(organism_slug: str, filename: str) -> tuple[str, int, int]:
    """Return (index_text, index_bytes, bytes_streamed). Stops after index member."""
    url = f"{LATEST}/{filename}"
    thr()
    log.info("stream_index org=%s url=%s", organism_slug, url)
    r = fetch(url, log, thr, stream=True)
    # http_util returns Response; use raw stream
    class _Count:
        def __init__(self, f):
            self.f = f
            self.n = 0

        def read(self, size=-1):
            b = self.f.read(size)
            self.n += len(b)
            return b

        def readable(self):
            return True

    # Prefer requests stream with decode_content so Cloudflare gzip of HTTP is handled;
    # the body itself is a .tar.gz so we gunzip again.
    raw = r.raw
    raw.decode_content = True
    cnt = _Count(raw)
    gz = gzip.GzipFile(fileobj=cnt)
    tf = tarfile.open(fileobj=gz, mode="r|")
    index_text = None
    index_bytes = 0
    try:
        for m in tf:
            if m.name.endswith(".index.tab.txt") and m.isreg():
                data = tf.extractfile(m).read()
                index_text = data.decode("utf-8", errors="replace")
                index_bytes = len(data)
                log.info(
                    "index org=%s member=%s index_bytes=%d streamed≈%d",
                    organism_slug, m.name, index_bytes, cnt.n,
                )
                break
            else:
                f = tf.extractfile(m)
                if f:
                    while f.read(1 << 20):
                        pass
        else:
            raise RuntimeError(f"no index.tab.txt in {filename}")
    finally:
        try:
            r.close()
        except Exception:
            pass
    return index_text, index_bytes, cnt.n


def parse_index(text: str) -> list[dict]:
    # Header line starts with #SCREEN_ID
    lines = text.splitlines()
    header = None
    rows = []
    for line in lines:
        if not line.strip():
            continue
        if line.startswith("#SCREEN_ID") or line.startswith("SCREEN_ID"):
            header = line.lstrip("#").split("\t")
            continue
        if line.startswith("#") or header is None:
            continue
        cols = line.split("\t")
        # pad/truncate
        if len(cols) < len(header):
            cols += [""] * (len(header) - len(cols))
        row = dict(zip(header, cols[: len(header)]))
        rows.append(row)
    return rows


def store_screen(c, row: dict, organism_slug: str, filename: str):
    sid = _dash(row.get("SCREEN_ID"))
    if not sid:
        return False
    tags = liver_tags_for(row)
    org_name = _dash(row.get("ORGANISM_OFFICIAL")) or organism_slug.replace("_", " ")
    cell = _dash(row.get("CELL_LINE"))
    ctype = _dash(row.get("CELL_TYPE"))
    pheno = _dash(row.get("PHENOTYPE"))
    lib = _dash(row.get("LIBRARY"))
    lib_type = _dash(row.get("LIBRARY_TYPE"))
    method = _dash(row.get("LIBRARY_METHODOLOGY"))
    author = _dash(row.get("AUTHOR"))
    sname = _dash(row.get("SCREEN_NAME"))
    source_id = _dash(row.get("SOURCE_ID"))
    source_type = _dash(row.get("SOURCE_TYPE")) or "pubmed"
    assay = " / ".join(x for x in [lib_type, method, "CRISPR screen"] if x)
    title = f"ORCS {sid}: {sname or ''} {author or ''} [{cell or ctype or organism_slug}]".strip()
    portal = f"{PORTAL}/Screen/{sid}"
    # Dataset portal often groups by publication; Screen page is per-screen
    score_ptr = (
        f"{LATEST}/{filename}"
        f"#BIOGRID-ORCS-SCREEN_{sid}-{RELEASE}.screen.tab.txt"
    )
    ds = f"ORCS:{sid}"
    upsert_dataset(c, dict(
        canonical_id=ds,
        title=title[:500],
        organism=org_name,
        tissue=ctype,
        disease_tags=tags or [],
        assay=assay,
        platforms=[x for x in [lib, lib_type, _dash(row.get("ENZYME")), _dash(row.get("ANALYSIS"))] if x],
        n_samples=_int(row.get("SCORES_SIZE")),
        access="public",
        license="MIT",
        confidence="metadata_only",
        probe_status="ok",
    ))
    add_member(c, ds, "biogrid_orcs", sid, "screen", portal)
    if source_id and source_type and source_type.lower() == "pubmed":
        add_member(c, ds, "pubmed", str(source_id), "publication",
                   f"https://pubmed.ncbi.nlm.nih.gov/{source_id}/")
    add_artifact(c, dict(
        artifact_id=f"ORCS:{sid}:portal", dataset_id=ds, kind="processed",
        role="screen_portal", title=title[:200], format="html",
        uri=portal, access="public", source="biogrid_orcs"))
    add_artifact(c, dict(
        artifact_id=f"ORCS:{sid}:index", dataset_id=ds, kind="processed",
        role="screen_index_record",
        title=f"ORCS {RELEASE} screen index row (from organism tar)",
        format="tsv",
        uri=f"{LATEST}/{filename}#index",
        access="public", source="biogrid_orcs"))
    add_artifact(c, dict(
        artifact_id=f"ORCS:{sid}:scores_ptr", dataset_id=ds, kind="processed",
        role="score_table_pointer",
        title=f"ORCS screen score table pointer (NOT downloaded)",
        filename=f"BIOGRID-ORCS-SCREEN_{sid}-{RELEASE}.screen.tab.txt",
        format="tsv", uri=score_ptr, access="public", source="biogrid_orcs"))
    # Compact raw (drop huge notes duplication risk — notes already a column)
    raw = {k: _dash(row.get(k)) for k in row}
    c.execute(
        "INSERT OR REPLACE INTO src_orcs_screen VALUES("
        + ",".join(["?"] * 40) + ")",
        (
            sid,
            source_id,
            source_type,
            author,
            sname,
            _int(row.get("SCORES_SIZE")),
            _int(row.get("FULL_SIZE")),
            _dash(row.get("FULL_SIZE_AVAILABLE")),
            _int(row.get("NUMBER_OF_HITS")),
            _dash(row.get("ANALYSIS")),
            _dash(row.get("SIGNIFICANCE_INDICATOR")),
            _dash(row.get("SIGNIFICANCE_CRITERIA")),
            _dash(row.get("THROUGHPUT")),
            _dash(row.get("SCREEN_TYPE")),
            _dash(row.get("SCREEN_FORMAT")),
            _dash(row.get("EXPERIMENTAL_SETUP")),
            _dash(row.get("DURATION")),
            _dash(row.get("CONDITION_NAME")),
            _dash(row.get("CONDITION_DOSAGE")),
            _dash(row.get("MOI")),
            lib,
            lib_type,
            method,
            _dash(row.get("ENZYME")),
            cell,
            ctype,
            pheno,
            _dash(row.get("ORGANISM_ID")),
            org_name,
            (notes := _dash(row.get("NOTES")) or "")[:4000] or None,
            _dash(row.get("SCREEN_RATIONALE")),
            _dash(row.get("SOURCE")),
            RELEASE,
            organism_slug,
            portal,
            score_ptr,
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
        log.info("already complete; exit")
        return

    done_orgs = set(cur.get("done_orgs", []))
    total_screens = int(cur.get("screens", 0))
    liver = int(cur.get("liver", 0))
    errors = int(cur.get("errors", 0))

    for organism_slug, filename in ORGS:
        if organism_slug in done_orgs:
            log.info("skip done org=%s", organism_slug)
            continue
        try:
            text, index_bytes, streamed = stream_index(organism_slug, filename)
            rows = parse_index(text)
            n_new = 0
            n_liv = 0
            for row in rows:
                if store_screen(c, row, organism_slug, filename):
                    n_new += 1
                    if liver_tags_for(row):
                        n_liv += 1
            # HEAD-ish size via Content-Length if present not available after stream;
            # store streamed + pointer
            url = f"{LATEST}/{filename}"
            c.execute(
                "INSERT OR REPLACE INTO src_orcs_release_file VALUES(?,?,?,?,?,?,?,?)",
                (filename, organism_slug, url, None, index_bytes, streamed, n_new, now_iso()),
            )
            # catalog artifact per organism release
            ds_rel = f"ORCS:RELEASE:{RELEASE}:{organism_slug}"
            upsert_dataset(c, dict(
                canonical_id=ds_rel,
                title=f"BioGRID ORCS {RELEASE} screen index — {organism_slug}",
                organism=organism_slug.replace("_", " "),
                assay="CRISPR screen catalog (index only)",
                access="public", license="MIT",
                confidence="metadata_only", probe_status="ok",
                n_samples=n_new,
            ))
            add_member(c, ds_rel, "biogrid_orcs_release", organism_slug, "release", DOWNLOADS)
            add_artifact(c, dict(
                artifact_id=f"ORCS:REL:{RELEASE}:{organism_slug}:tar",
                dataset_id=ds_rel, kind="raw", role="organism_screens_tar_pointer",
                title=filename, filename=filename, format="tar.gz",
                uri=url, access="public", source="biogrid_orcs",
                size_bytes=None,
            ))
            c.commit()
            total_screens += n_new
            liver += n_liv
            done_orgs.add(organism_slug)
            cur.update({
                "done_orgs": sorted(done_orgs),
                "screens": total_screens,
                "liver": liver,
                "errors": errors,
                "release": RELEASE,
                "last_org": organism_slug,
                "last_org_screens": n_new,
                "last_org_liver": n_liv,
                "last_streamed": streamed,
            })
            save_cursor(c, SRC, cur)
            log.info(
                "org=%s screens=%d liver=%d streamed≈%d",
                organism_slug, n_new, n_liv, streamed,
            )
        except Exception as e:
            errors += 1
            cur["errors"] = errors
            cur["last_error_org"] = organism_slug
            save_cursor(c, SRC, cur, err=str(e))
            log.exception("org=%s failed: %s", organism_slug, e)
            raise

    cur.update({
        "complete": True,
        "finished_at": now_iso(),
        "screens": total_screens,
        "liver": liver,
        "errors": errors,
        "note": (
            f"BioGRID ORCS {RELEASE} screen indices streamed from Latest-Release "
            "organism tar.gz (index member only; score matrices NOT downloaded). "
            f"Portal {PORTAL}; REST orcsws.thebiogrid.org needs free accesskey."
        ),
    })
    save_cursor(c, SRC, cur)
    n_db = c.execute("SELECT COUNT(*) FROM src_orcs_screen").fetchone()[0]
    n_liv = c.execute("SELECT COUNT(*) FROM src_orcs_screen WHERE is_liver=1").fetchone()[0]
    log.info("DONE screens=%d liver=%d errors=%d db=%s", n_db, n_liv, errors, DB)


if __name__ == "__main__":
    main()
