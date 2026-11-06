#!/usr/bin/env python3
"""Shared LLM-Wiki helpers for literature-library v0.5.3-draft.

Vault layout (per collection):
  wiki-vault/
    AGENTS.md          # schema / agent conventions (immutable contract)
    raw/               # immutable source snapshots
      publications/<safe_id>/meta.json
      publications/<safe_id>/body.txt   # optional text body / PDF text
    wiki/
      index.md         # theme tree / directory (human + agent editable)
      log.md           # append-only activity log
      sources/<safe_id>.md
      concepts/<slug>.md
      entities/<slug>.md

NO embeddings. Theme tree = index.md headings + optional theme_nodes/theme_pubs in SQLite.
TextIn document-tree is NOT used as KB backbone.

v0.4: default wiki compile is agent-model synthesis (request → agent JSON → apply). No plugin API keys.
Heuristic title+abstract extractors remain as explicit offline fallback only.
v0.4.5: per-source concept contributions + agent/mock concept deep merge (not upsert-link only).
"""
from __future__ import annotations

import json
import re
import shutil
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

VAULT_NAME = "wiki-vault"
SCHEMA_LABEL = "agent-wiki-v0.4.5"

AGENTS_MD = """# Literature Wiki Vault — AGENTS.md

Schema for the LLM Wiki layer beside `library.sqlite` (bibliographic identity / dedup / runs).
Contract marker: `agent-wiki-v0.4.5` (compat: agent-wiki-v0.4).

## Layers (Karpathy/llmwikis style)

1. **raw/** — immutable source snapshots. Never rewrite; only add. Also holds `compile_request.json` / `compile_response.json` for agent compile.
2. **wiki/** — compiled markdown (sources / concepts / entities / index / log). Agents may edit wiki pages; prefer append + link.
3. **This file (AGENTS.md)** — conventions. Do not delete.

## Page types

| Path | Role |
|--|--|
| `wiki/index.md` | Theme tree / table of contents. Headings = branches; links to sources/concepts. |
| `wiki/log.md` | Append-only activity log. |
| `wiki/sources/<id>.md` | One page per publication. Summary / Key claims / Methods from **agent-model compile** (primary). |
| `wiki/concepts/<slug>.md` | Concept / topic / method pages (multi-paper; agent/mock deep merge). |
| `wiki/concepts/_contrib/<slug>/<source_id>.json` | Per-source concept contribution (deterministic; no LLM). |
| `wiki/entities/<slug>.md` | Entity pages (gene/protein symbols). Fragment tokens like `l1`/`ifn` are forbidden. |

## Compile rules (lit-wiki-ingest) — agent-wiki-v0.4

- **Primary:** the **invoking agent's own model** produces structured compile JSON. Scripts only: build grounded prompt + schema → validate → write vault. **No plugin API keys.**
- Flow: `wiki_ingest` / `wiki_compile_request` emits `compile_request.json` → agent fills `compile_response.json` → `wiki_ingest --llm-json` / `wiki_compile_apply`.
- **Offline tests:** `--compile mock` only.
- **Emergency fallback:** `--compile heuristic` or `LIT_WIKI_COMPILE=heuristic` (v0.2 scaffolding; not target quality).
- Schema failure: agent must fix JSON; optional `--allow-heuristic-fallback`.
- Grounding: use only provided title/abstract/metadata; no invented IDs/stats; descriptive biology only.
- Auto-suggest `tags_json` into SQLite when empty/[] so theme-tree suggest works.
- Prefer relative markdown links: `[title](sources/doi_10_1016_….md)`, `[LAYN](../entities/layn.md)`.
- Soft-conflict pubs (`needs_review`) get a banner; do not silently merge wiki pages.
- Suggested theme branches from `lit-wiki-suggest-tree` are JSON proposals only; `lit-wiki-apply-tree --confirm` writes.
- TextIn / commercial PDF tree APIs are optional parsers later — not vault backbone.
- No embeddings.

## Concept deep merge (v0.4.5)

- Each wiki compile writes a **contribution** JSON under `wiki/concepts/_contrib/<slug>/<safe_id>.json` (grounded snippets from that paper only).
- Concept pages keep ## Sources links via upsert; **Definition / Across papers** are updated only via `concept_merge_request` → agent JSON → `concept_merge_apply` (or `--mock`).
- Scripts never invent cross-paper synthesis except `--mock` (deterministic offline).
- Soft-conflict pubs (`needs_review`) must keep a banner on the concept page; do not silent-merge them into synthesis.
- Optional agent fields: `canonical_slug` / `merge_into` → old page becomes a redirect note (no delete).

## Collection binding

This vault lives under `collections/<collection_id>/wiki-vault/`. Collection ≈ one vault root.
"""

INDEX_MD = """# Wiki index

Theme tree for this collection. Edit headings and links freely; papers may appear under multiple branches.

## Unfiled

<!-- Sources land here until hung on a theme branch. -->

## Themes

<!-- Apply confirmed proposals from lit-wiki-suggest-tree here. -->

## Concepts

<!-- High-signal concept pages linked after compile. -->

## Entities

<!-- High-signal entity pages linked after compile. -->
"""

LOG_MD = """# Wiki log

Append-only. Newest entries at the bottom.

"""

# ---------------------------------------------------------------------------
# Lexicons (deterministic, no network)
# ---------------------------------------------------------------------------

METHOD_PATTERNS: list[tuple[str, str]] = [
    # (display name, regex — case insensitive)
    ("scRNA-seq", r"\b(?:scrna[\s-]?seq(?:uencing)?|single[\s-]?cell\s+rna[\s-]?seq(?:uencing)?|single[\s-]?cell\s+sequencing)\b"),
    ("single-cell", r"\bsingle[\s-]?cell\b"),
    ("RNA-seq", r"\brna[\s-]?seq(?:uencing)?\b"),
    ("CITE-seq", r"\bcite[\s-]?seq\b"),
    ("ATAC-seq", r"\batac[\s-]?seq\b"),
    ("ChIP-seq", r"\bchip[\s-]?seq\b"),
    ("flow cytometry", r"\bflow\s+cytometr(?:y|ic)\b"),
    ("CyTOF", r"\bcytof\b"),
    ("mass cytometry", r"\bmass\s+cytometr(?:y|ic)\b"),
    ("CRISPR", r"\bcrispr(?:[\s-]?cas9)?\b"),
    ("limma", r"\blimma\b"),
    ("Seurat", r"\bseurat\b"),
    ("Scanpy", r"\bscanpy\b"),
    ("SC3", r"\bsc3\b"),
    ("TraCeR", r"\btracer\b"),
    ("TCR sequencing", r"\btcr[\s-]?(?:seq(?:uencing)?|repertoire)\b"),
    ("BCR sequencing", r"\bbcr[\s-]?(?:seq(?:uencing)?|repertoire)\b"),
    ("pseudotime", r"\bpseudotime\b"),
    ("trajectory analysis", r"\btrajectory\s+analy\w*\b"),
    ("differential expression", r"\bdifferential(?:ly)?\s+express\w*\b"),
    ("WGCNA", r"\bwgcna\b"),
    ("GSEA", r"\bgsea\b|\bgene\s+set\s+enrichment\b"),
    ("immunohistochemistry", r"\bimmunohistochemistr\w*\b|\bihc\b"),
    ("western blot", r"\bwestern\s+blot\w*\b"),
    ("qPCR", r"\bq(?:rt[\s-]?)?pcr\b"),
    ("spatial transcriptomics", r"\bspatial\s+transcriptom\w*\b"),
    ("multiomics", r"\bmulti[\s-]?omics\b"),
    ("machine learning", r"\bmachine\s+learning\b|\bdeep\s+learning\b"),
    ("LASSO", r"\blasso\b"),
    ("Cox regression", r"\bcox\s+(?:proportional\s+hazards|regression)\b"),
]

DISEASE_TERMS: list[tuple[str, str]] = [
    ("hepatocellular carcinoma", r"\bhepatocellular\s+carcinoma\b|\bHCC\b"),
    ("liver cancer", r"\bliver\s+cancer\b"),
    ("melanoma", r"\bmelanoma\b"),
    ("non-small cell lung cancer", r"\bnon[\s-]?small[\s-]?cell\s+lung\s+cancer\b|\bNSCLC\b"),
    ("lung cancer", r"\blung\s+cancer\b"),
    ("breast cancer", r"\bbreast\s+cancer\b"),
    ("colorectal cancer", r"\bcolorectal\s+cancer\b|\bCRC\b"),
    ("gastric cancer", r"\bgastric\s+cancer\b"),
    ("pancreatic cancer", r"\bpancreatic\s+(?:cancer|ductal\s+adenocarcinoma)\b|\bPDAC\b"),
    ("glioblastoma", r"\bglioblastoma\b|\bGBM\b"),
    ("leukemia", r"\bleukemia\b|\bleukaemia\b"),
    ("lymphoma", r"\blymphoma\b"),
    ("COVID-19", r"\bCOVID[\s-]?19\b|\bSARS[\s-]?CoV[\s-]?2\b"),
]

TOPIC_TERMS: list[tuple[str, str]] = [
    ("T cell exhaustion", r"\bT[\s-]?cell\s+exhaustion\b|\bexhaust(?:ed|ion)\s+(?:CD\d\+?\s+)?T\s+cells?\b|\bexhausted(?:\s+[\w+]+)*\s+T\s+cells?\b|\bexhaustion\b"),
    ("tumor microenvironment", r"\btumor\s+microenvironment\b|\bTME\b|\btumour\s+microenvironment\b"),
    ("immune checkpoint", r"\bimmune\s+checkpoint\b"),
    ("immunotherapy", r"\bimmunotherap\w*\b|\bcheckpoint\s+blockade\b|\bPD[\s-]?1\s+blockade\b"),
    ("infiltrating T cells", r"\binfiltrat\w*\s+T\s+cells?\b|\bTIL\b|\btumor[\s-]?infiltrat\w*\b"),
    ("regulatory T cells", r"\bregulatory\s+T\s+cells?\b|\bTregs?\b"),
    ("CD8 T cells", r"\bCD8\+?\s+T\s+cells?\b"),
    ("CD4 T cells", r"\bCD4\+?\s+T\s+cells?\b"),
    ("antigen presentation", r"\bantigen\s+present\w*\b"),
    ("cytotoxicity", r"\bcytotoxic(?:ity|ic)\b"),
    ("clonal expansion", r"\bclonal\s+expansion\b"),
    ("single-cell atlas", r"\batlas\b"),
]

# High-signal immune / oncology gene-like symbols (boosted even if rare)
KNOWN_GENE_SYMBOLS = {
    "PDCD1", "PD-1", "PD1", "CD274", "PD-L1", "PDL1", "CTLA4", "CTLA-4",
    "HAVCR2", "TIM-3", "TIM3", "LAG3", "LAG-3", "TIGIT", "LAYN", "ENTPD1",
    "CD39", "CD8A", "CD8B", "CD4", "FOXP3", "IL2RA", "CD25", "GZMB", "GZMA",
    "PRF1", "IFNG", "TNF", "IL10", "TGFB1", "TOX", "TCF7", "TCF1", "CXCL13",
    "CXCR5", "CCR7", "SELL", "CD69", "ITGAE", "CD103", "MKI67", "PCNA",
    "TP53", "MYC", "EGFR", "KRAS", "BRAF", "ALK", "MET", "VEGF", "VEGFA",
    "CD3D", "CD3E", "TRAC", "TRBC1", "TRBC2", "NCAM1", "CD56", "NKG7",
    "GNLY", "KLRB1", "FCGR3A", "CD16", "MS4A1", "CD19", "CD79A", "IGKC",
    "TCR", "NKG2A", "NKG2C", "NKG2D", "HLA-E", "HLA-A", "HLA-B", "HLA-C",
    "HLA-DR", "TOP2A", "LDHA", "SLC2A1", "PGK1", "SERPINE1", "CA9",
}

# Short / prefix fragments from hyphen or slash splits — never emit as entities
GENE_FRAGMENT_BLOCK = {
    "L1", "L2", "L3", "D1", "D2", "D3", "A1", "A2", "B1", "B2", "C1", "C2",
    "PD", "IFN", "TGF", "CTL", "CTLA", "TIM", "LAG", "NKG", "HLA", "M1", "M2",
    "TH1", "TH2", "IL", "CD", "DR", "AB", "Ig", "IG",
}

# len==3 allcaps without digit: only keep if allowlisted / known
SHORT_SYMBOL_ALLOW = {
    "TOX", "MYC", "ALK", "MET", "TNF", "TCR",
}

# Entity page slugs that indicate extraction bugs (lint asserts absence)
NOISY_ENTITY_SLUGS = {"l1", "l2", "ifn", "pd", "tgf", "m1", "m2", "cd", "il", "hla"}

# Tokens that look gene-like but are noise
GENE_STOP = {
    "THE", "AND", "FOR", "WITH", "FROM", "THIS", "THAT", "ARE", "WAS", "WERE",
    "HAS", "HAVE", "HAD", "NOT", "BUT", "ALL", "CAN", "MAY", "ITS", "VIA",
    "USING", "USED", "BASED", "ALSO", "THAN", "INTO", "OVER", "AFTER", "BOTH",
    "EACH", "SUCH", "THESE", "THOSE", "THEN", "WHEN", "WHERE", "WHAT", "WHICH",
    "WHILE", "ABOUT", "BETWEEN", "THROUGH", "DURING", "BEFORE", "UNDER", "AGAIN",
    "FURTHER", "ONCE", "HERE", "THERE", "RNA", "DNA", "PCR", "ANOVA", "AUC",
    "ROC", "HR", "OR", "CI", "SD", "SEM", "USA", "UK", "COVID", "HIV", "HCV",
    "HBV", "HCC", "NSCLC", "PDAC", "GBM", "TIL", "TME", "IHC", "FACS", "UMAP",
    "TSNE", "PCA", "DEG", "GO", "KEGG", "GSEA", "WGCNA", "LASSO", "ELISA",
    "IFN", "OTA", "CCK", "LDH", "IAP", "ICC", "OVA", "KPC", "SMAC",
    "TCGA", "ENCODE", "ICGC", "GEO", "SRA", "GTEx", "GTEX", "DBGAP",
    "TIDE", "LIHC", "GSE", "MIES", "MTRG", "ACTR", "CFL",
}




def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def safe_id(publication_id: str) -> str:
    """Filesystem-safe id derived from canonical publication id."""
    s = publication_id.strip()
    s = s.replace(":", "_").replace("/", "_").replace("\\", "_")
    s = re.sub(r"[^\w.\-]+", "_", s, flags=re.UNICODE)
    s = re.sub(r"_+", "_", s).strip("._")
    return (s[:180] or "unknown")


def slugify(text: str) -> str:
    s = text.strip().lower()
    s = s.replace("+", "-plus-")
    s = re.sub(r"[^\w\-]+", "-", s, flags=re.UNICODE)
    s = re.sub(r"-+", "-", s).strip("-")
    return (s[:80] or "untitled")


# Common aliases: short form / synonym → preferred slug (only used if target page exists)
CONCEPT_SLUG_ALIASES: dict[str, str] = {
    "nsclc": "non-small-cell-lung-cancer",
    "non-small-cell-lung-cancer": "non-small-cell-lung-cancer",
    "non small cell lung cancer": "non-small-cell-lung-cancer",
    "scrna-seq": "scrna-seq",
    "sc-rna-seq": "scrna-seq",
    "single-cell-rna-seq": "scrna-seq",
    "single-cell-rna-sequencing": "scrna-seq",
    "single cell rna seq": "scrna-seq",
    "hcc": "hepatocellular-carcinoma",
    "hepatocellular-carcinoma": "hepatocellular-carcinoma",
    "tme": "tumor-microenvironment",
    "tumor-microenvironment": "tumor-microenvironment",
    "tumour-microenvironment": "tumor-microenvironment",
    "treg": "regulatory-t-cells",
    "tregs": "regulatory-t-cells",
    "regulatory-t-cells": "regulatory-t-cells",
    "til": "infiltrating-t-cells",
    "tils": "infiltrating-t-cells",
    "infiltrating-t-cells": "infiltrating-t-cells",
    "exhaustion": "t-cell-exhaustion",
    "t-cell-exhaustion": "t-cell-exhaustion",
    "cd8-t-cells": "cd8-t-cells",
    "cd4-t-cells": "cd4-t-cells",
}


def list_concept_slugs(vault: Path) -> set[str]:
    """Existing concept page stems under wiki/concepts/*.md (excludes _contrib)."""
    d = vault / "wiki" / "concepts"
    if not d.is_dir():
        return set()
    return {p.stem for p in d.glob("*.md") if p.is_file()}


def resolve_related_concept_slug(
    vault: Path,
    name_or_slug: str,
    *,
    existing: Optional[set[str]] = None,
) -> Optional[str]:
    """Map a related_concepts token to an existing concept page slug, or None.

    Resolution order: exact slug → slugify(name) → CONCEPT_SLUG_ALIASES (and
    reverse alias peers) → simple fuzzy (alias target / candidate equals existing).
    Never invents pages; only links to wiki/concepts/<slug>.md that already exist.
    """
    raw = (name_or_slug or "").strip()
    if not raw:
        return None
    pages = existing if existing is not None else list_concept_slugs(vault)
    if not pages:
        return None

    candidates: list[str] = []
    # exact as given (already a slug?)
    candidates.append(raw)
    candidates.append(raw.lower())
    # slugify display name
    candidates.append(slugify(raw))
    # strip backticks / brackets
    cleaned = re.sub(r"[`\[\]]", "", raw).strip()
    if cleaned and cleaned != raw:
        candidates.append(slugify(cleaned))

    # alias table (key may be slug or loose name)
    for key in list(candidates):
        alias = CONCEPT_SLUG_ALIASES.get(key) or CONCEPT_SLUG_ALIASES.get(slugify(key))
        if alias:
            candidates.append(alias)

    # bidirectional peers: if alias maps A→B, also try B's other keys that share value
    for key, val in CONCEPT_SLUG_ALIASES.items():
        if val in candidates or key in candidates:
            candidates.append(val)
            candidates.append(key)

    seen: set[str] = set()
    for c in candidates:
        c = (c or "").strip().lower().strip("-")
        if not c or c in seen:
            continue
        seen.add(c)
        if c in pages:
            return c

    # fuzzy: if slugify token is a compact prefix/suffix of exactly one existing page
    token = slugify(raw)
    if token and len(token) >= 3:
        hits = [p for p in pages if p == token or p.startswith(token + "-") or p.endswith("-" + token)]
        # also: token without hyphens equals page without hyphens
        compact = token.replace("-", "")
        if not hits and compact:
            hits = [p for p in pages if p.replace("-", "") == compact]
        if len(hits) == 1:
            return hits[0]
    return None


def resolve_related_concepts(
    vault: Path,
    related: list[str],
) -> tuple[list[dict[str, str]], list[str]]:
    """Resolve related_concepts → [{display, slug}] for existing pages; list unresolved."""
    pages = list_concept_slugs(vault)
    resolved: list[dict[str, str]] = []
    unresolved: list[str] = []
    seen_slugs: set[str] = set()
    for item in related or []:
        if not isinstance(item, str) or not item.strip():
            continue
        display = item.strip()
        slug = resolve_related_concept_slug(vault, display, existing=pages)
        if slug is None:
            unresolved.append(display)
            continue
        if slug in seen_slugs:
            continue
        seen_slugs.add(slug)
        resolved.append({"display": display, "slug": slug})
    return resolved, unresolved


def collection_dir(project_root: Path, collection_id: str, db_path: Optional[Path] = None) -> Path:
    """Resolve collection path from DB if present, else conventional layout."""
    root = project_root.resolve()
    if db_path is None:
        db_path = root / "library.sqlite"
    if db_path.exists():
        import sqlite3

        conn = sqlite3.connect(str(db_path))
        conn.row_factory = sqlite3.Row
        try:
            row = conn.execute(
                "SELECT path FROM collections WHERE id = ?", (collection_id,)
            ).fetchone()
            if row and row["path"]:
                return Path(row["path"])
        finally:
            conn.close()
    return root / "collections" / collection_id


def vault_root(col_dir: Path) -> Path:
    return col_dir / VAULT_NAME


def init_vault(col_dir: Path, *, collection_id: str, collection_name: str = "") -> dict[str, Any]:
    """Create wiki-vault skeleton. Idempotent for dirs; does not overwrite AGENTS.md/index if present."""
    root = vault_root(col_dir)
    created: list[str] = []
    for sub in (
        root,
        root / "raw" / "publications",
        root / "wiki" / "sources",
        root / "wiki" / "concepts",
        root / "wiki" / "concepts" / "_contrib",
        root / "wiki" / "entities",
    ):
        if not sub.exists():
            sub.mkdir(parents=True, exist_ok=True)
            created.append(str(sub.relative_to(root)) if sub != root else ".")

    agents = root / "AGENTS.md"
    # Refresh conventions when schema bumps (keep user edits only if same schema label present)
    if not agents.exists():
        agents.write_text(AGENTS_MD, encoding="utf-8")
        created.append("AGENTS.md")
    else:
        existing = agents.read_text(encoding="utf-8")
        if "agent-wiki-v0.4.5" not in existing and "Concept deep merge" not in existing:
            agents.write_text(AGENTS_MD, encoding="utf-8")
            created.append("AGENTS.md(updated)")

    index = root / "wiki" / "index.md"
    if not index.exists():
        title = collection_name or collection_id
        body = INDEX_MD.replace("# Wiki index", f"# Wiki index — {title}", 1)
        index.write_text(body, encoding="utf-8")
        created.append("wiki/index.md")
    else:
        # Ensure Concepts / Entities sections exist
        text = index.read_text(encoding="utf-8")
        changed = False
        if "## Concepts" not in text:
            text = text.rstrip() + "\n\n## Concepts\n\n<!-- High-signal concept pages linked after compile. -->\n"
            changed = True
        if "## Entities" not in text:
            text = text.rstrip() + "\n\n## Entities\n\n<!-- High-signal entity pages linked after compile. -->\n"
            changed = True
        if changed:
            index.write_text(text if text.endswith("\n") else text + "\n", encoding="utf-8")

    log = root / "wiki" / "log.md"
    if not log.exists():
        log.write_text(LOG_MD, encoding="utf-8")
        created.append("wiki/log.md")

    meta = {
        "collection_id": collection_id,
        "collection_name": collection_name,
        "created_at": utc_now(),
        "schema": SCHEMA_LABEL,
        "embeddings": False,
    }
    meta_path = root / "vault.json"
    if not meta_path.exists():
        meta_path.write_text(json.dumps(meta, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        created.append("vault.json")
    else:
        try:
            old = json.loads(meta_path.read_text(encoding="utf-8"))
            if old.get("schema") != SCHEMA_LABEL:
                old["schema"] = SCHEMA_LABEL
                old["embeddings"] = False
                old["updated_at"] = utc_now()
                meta_path.write_text(json.dumps(old, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        except json.JSONDecodeError:
            pass

    return {
        "ok": True,
        "vault_root": str(root),
        "created": created,
        "already_present": not bool(created),
        "schema": SCHEMA_LABEL,
    }


def append_log(vault: Path, line: str) -> None:
    log = vault / "wiki" / "log.md"
    log.parent.mkdir(parents=True, exist_ok=True)
    if not log.exists():
        log.write_text(LOG_MD, encoding="utf-8")
    with log.open("a", encoding="utf-8") as f:
        f.write(f"\n- `{utc_now()}` {line}\n")


def _md_escape(text: str) -> str:
    return (text or "").replace("\r\n", "\n").replace("\r", "\n")


def _parse_tags(raw_tags: Any) -> list[str]:
    tags: list[str] = []
    if raw_tags is None:
        return tags
    if isinstance(raw_tags, list):
        return [str(t).strip() for t in raw_tags if str(t).strip()]
    if isinstance(raw_tags, str):
        s = raw_tags.strip()
        if not s or s in ("[]", "null", "None"):
            return []
        try:
            parsed = json.loads(s)
            if isinstance(parsed, list):
                return [str(t).strip() for t in parsed if str(t).strip()]
        except json.JSONDecodeError:
            pass
        return [t.strip() for t in s.split(",") if t.strip()]
    return tags


def tags_are_empty(raw_tags: Any) -> bool:
    return len(_parse_tags(raw_tags)) == 0


# ---------------------------------------------------------------------------
# Deterministic extraction from title + abstract
# ---------------------------------------------------------------------------

_SENT_SPLIT = re.compile(r"(?<=[.!?])\s+(?=[A-Z0-9\"“\(])")


def split_sentences(text: str) -> list[str]:
    text = re.sub(r"\s+", " ", (text or "").strip())
    if not text:
        return []
    parts = _SENT_SPLIT.split(text)
    out: list[str] = []
    for p in parts:
        p = p.strip()
        if len(p) < 20:
            # glue very short fragments to previous if any
            if out:
                out[-1] = (out[-1] + " " + p).strip()
            elif p:
                out.append(p)
            continue
        out.append(p)
    return out


def condense_summary(title: str, abstract: str, *, max_chars: int = 480) -> str:
    """First-paragraph / leading-sentence condensation — not a placeholder stub."""
    abs_text = (abstract or "").strip()
    if not abs_text:
        t = (title or "").strip()
        return t if t else "_No abstract available for condensation._"
    # Prefer first paragraph (double newline) else first 2–3 sentences
    para = abs_text.split("\n\n")[0].strip()
    para = re.sub(r"\s+", " ", para)
    sents = split_sentences(para)
    if not sents:
        chunk = para[:max_chars]
        return chunk + ("…" if len(para) > max_chars else "")
    buf = sents[0]
    for s in sents[1:3]:
        if len(buf) + 1 + len(s) > max_chars:
            break
        buf = buf + " " + s
    if len(buf) > max_chars:
        buf = buf[: max_chars - 1].rstrip() + "…"
    return buf


_CLAIM_CUE = re.compile(
    r"\b(?:we\s+(?:show|found|identify|identified|demonstrate|demonstrated|reveal|revealed|"
    r"observe|observed|report|reported|conclude|propose|provide)|"
    r"(?:results|findings)\s+(?:show|indicate|suggest|demonstrate)|"
    r"(?:is|are|was|were)\s+associated\s+with|"
    r"(?:plays?|played)\s+a\s+(?:critical|key|important)\s+role|"
    r"here\s+we|"
    r"collectively|"
    r"these\s+(?:data|results|findings))\b",
    re.I,
)


def extract_key_claims(abstract: str, *, max_claims: int = 5) -> list[str]:
    sents = split_sentences(abstract or "")
    if not sents:
        return []
    scored: list[tuple[int, str]] = []
    for i, s in enumerate(sents):
        score = 0
        if _CLAIM_CUE.search(s):
            score += 3
        if re.search(r"\b(?:significantly|novel|critical|essential|robustly)\b", s, re.I):
            score += 1
        if 40 <= len(s) <= 320:
            score += 1
        # prefer mid/late abstract slightly (results/conclusion)
        if i >= max(1, len(sents) // 3):
            score += 1
        if score:
            scored.append((score, s))
    scored.sort(key=lambda x: (-x[0],))
    claims = [s for _, s in scored[:max_claims]]
    # Fallback: first 2–3 sentences if no cues
    if not claims:
        claims = sents[: min(3, len(sents))]
    return claims


def extract_methods(text: str) -> list[str]:
    blob = text or ""
    found: list[str] = []
    for name, pat in METHOD_PATTERNS:
        if re.search(pat, blob, flags=re.I):
            if name not in found:
                found.append(name)
    return found


def extract_diseases(text: str) -> list[str]:
    blob = text or ""
    found: list[str] = []
    for name, pat in DISEASE_TERMS:
        if re.search(pat, blob, flags=re.I):
            if name not in found:
                found.append(name)
    return found


def extract_topics(text: str) -> list[str]:
    blob = text or ""
    found: list[str] = []
    for name, pat in TOPIC_TERMS:
        if re.search(pat, blob, flags=re.I):
            if name not in found:
                found.append(name)
    return found


_GENE_TOKEN = re.compile(
    r"\b("
    r"[A-Z][A-Z0-9]{2,7}"  # ALLCAPS len 3–8 (rejects L1)
    r"|[A-Z]{1,3}\d{1,3}[A-Z]?"  # CD8, IL10, CA9
    r"|[A-Z][A-Z0-9]+-[A-Z0-9]+"  # PD-1, CTLA-4, TIM-3, PD-L1, NKG2A-ish
    r"|CD\d{1,3}[A-Z]?"
    r"|IL\d{1,2}[A-Z]?"
    r"|CXC[LRC]\d{1,2}"
    r"|CCR\d{1,2}"
    r"|HLA-[A-Z0-9]+"
    r"|NKG2[A-Z]"
    r")\b"
)

# Hyphenated forms to claim first so PD-L1 does not also yield L1
_HYPHEN_GENE = re.compile(
    r"\b("
    r"PD-L?1|PD-1|CTLA-4|TIM-3|LAG-3|HLA-[A-Z0-9]+|IFN-γ|IFN-gamma|IFN-γ|TGF-β|TGF-beta"
    r"|NKG2[A-Z]|CD\d{1,3}[A-Z]?"
    r")\b",
    re.I,
)


def extract_gene_entities(text: str, *, max_entities: int = 16) -> list[str]:
    """Gene/protein-ish tokens + known symbol boost; drop hyphen/slash fragments."""
    blob = text or ""
    counts: Counter[str] = Counter()
    alias_map = {
        "PD1": "PD-1",
        "PDCD1": "PD-1",
        "PDL1": "PD-L1",
        "CD274": "PD-L1",
        "CTLA4": "CTLA-4",
        "TIM3": "TIM-3",
        "HAVCR2": "TIM-3",
        "LAG3": "LAG-3",
        "TCF1": "TCF7",
        "IFNG": "IFNG",
        "IFN-G": "IFNG",
        "IFN-GAMMA": "IFNG",
        "IFN-Γ": "IFNG",
    }
    known_upper = {k.upper() for k in KNOWN_GENE_SYMBOLS}
    known_compact = {k.upper().replace("-", "") for k in KNOWN_GENE_SYMBOLS}

    # Spans covered by hyphenated / multi-part matches — suppress nested fragments
    covered: list[tuple[int, int]] = []
    for m in _HYPHEN_GENE.finditer(blob):
        covered.append((m.start(), m.end()))
        raw = m.group(1)
        up = raw.upper().replace("Γ", "G").replace("Β", "B")
        tok = raw
        if up in ("IFN-G", "IFN-GAMMA", "IFN-Γ", "IFNG") or raw.lower() in ("ifn-γ", "ifn-gamma"):
            tok = "IFNG"
        elif up in ("TGF-B", "TGF-BETA", "TGF-Β"):
            tok = "TGFB1"
        elif up.replace("-", "") in known_compact or up in known_upper:
            for k in KNOWN_GENE_SYMBOLS:
                if k.upper() == up or k.upper().replace("-", "") == up.replace("-", ""):
                    tok = k
                    break
        # Canonical hyphen display
        for k in ("PD-1", "PD-L1", "CTLA-4", "TIM-3", "LAG-3"):
            if up.replace("-", "") == k.upper().replace("-", ""):
                tok = k
                break
        for k in ("NKG2A", "NKG2C", "NKG2D"):
            if up == k:
                tok = k
                break
        if tok.upper() in GENE_STOP or tok.upper() in GENE_FRAGMENT_BLOCK:
            continue
        counts[tok] += 3

    def _covered(pos: int, end: int) -> bool:
        for a, b in covered:
            # nested or overlapping fragment inside a hyphenated hit
            if pos >= a and end <= b and (pos != a or end != b):
                return True
            # token that starts right after a hyphen inside covered span
            if a <= pos < b:
                return True
        return False

    for m in _GENE_TOKEN.finditer(blob):
        if _covered(m.start(), m.end()):
            continue
        tok = m.group(1)
        up = tok.upper()
        if up in GENE_STOP or up in GENE_FRAGMENT_BLOCK:
            continue
        if len(tok) < 3:
            continue
        if tok in alias_map:
            tok = alias_map[tok]
            up = tok.upper()
        known = up in known_upper or up.replace("-", "") in known_compact or tok in KNOWN_GENE_SYMBOLS
        has_digit = bool(re.search(r"\d", tok))
        hyphenated = "-" in tok
        # len 3 allcaps without digit: allowlist / known only
        if len(tok) == 3 and tok.isupper() and not has_digit and not known:
            if tok not in SHORT_SYMBOL_ALLOW:
                continue
        shaped = has_digit or hyphenated or (tok.isupper() and 4 <= len(tok) <= 6) or known
        if not (known or shaped):
            continue
        # Prefer canonical known display
        for k in KNOWN_GENE_SYMBOLS:
            if k.upper().replace("-", "") == up.replace("-", "") or k.upper() == up:
                if k in ("PD-1", "PD-L1", "CTLA-4", "TIM-3", "LAG-3", "LAYN", "TOX", "TIGIT", "FOXP3", "GZMB", "PRF1", "IFNG", "NKG2A", "NKG2C", "NKG2D"):
                    tok = k
                elif "-" in k:
                    tok = k
                break
        if tok.upper() in GENE_FRAGMENT_BLOCK:
            continue
        # Skip GEO series / consortium acronyms unless known genes
        if re.match(r"^GSE\d+$", up) or re.match(r"^GSM\d+$", up):
            continue
        counts[tok] += 3 if known else 1

    # Synonym / lowercase mentions
    low = blob.lower()
    synonym_hits = {
        "layn": "LAYN",
        "layilin": "LAYN",
        "tox": "TOX",
        "tigit": "TIGIT",
        "foxp3": "FOXP3",
        "havcr2": "TIM-3",
        "tim-3": "TIM-3",
        "tim3": "TIM-3",
        "entpd1": "CD39",
        "cd39": "CD39",
        "gzmb": "GZMB",
        "prf1": "PRF1",
        "tcf7": "TCF7",
        "pd-1": "PD-1",
        "pd1": "PD-1",
        "pd-l1": "PD-L1",
        "ctla-4": "CTLA-4",
        "ctla4": "CTLA-4",
        "lag-3": "LAG-3",
        "lag3": "LAG-3",
        "ifng": "IFNG",
        "ifn-γ": "IFNG",
        "ifn-gamma": "IFNG",
        "interferon-gamma": "IFNG",
        "interferon γ": "IFNG",
        "nkg2a": "NKG2A",
        "nkg2c": "NKG2C",
        "nkg2d": "NKG2D",
    }
    for k, disp in synonym_hits.items():
        if re.search(rf"\b{re.escape(k)}\b", low):
            counts[disp] += 2

    # Final hard filter
    ranked: list[str] = []
    for t, _ in counts.most_common(max_entities * 2):
        up = t.upper()
        if up in GENE_STOP or up in GENE_FRAGMENT_BLOCK:
            continue
        if len(t) < 3:
            continue
        ranked.append(t)
        if len(ranked) >= max_entities:
            break
    return ranked



def compile_extract(title: str, abstract: str) -> dict[str, Any]:
    """Run all deterministic extractors."""
    blob = f"{title or ''}\n{abstract or ''}"
    methods = extract_methods(blob)
    diseases = extract_diseases(blob)
    topics = extract_topics(blob)
    genes = extract_gene_entities(blob)
    claims = extract_key_claims(abstract or "")
    summary = condense_summary(title or "", abstract or "")

    # Concepts = topics + diseases + methods (high-signal); entities = genes + disease short names
    concepts: list[str] = []
    for c in topics + diseases + methods:
        if c not in concepts:
            concepts.append(c)

    entities: list[str] = []
    for e in genes:
        if e not in entities:
            entities.append(e)

    # Tags for SQLite / tree suggest
    tags: list[str] = []
    for t in topics[:4] + diseases[:3] + methods[:4] + genes[:4]:
        if t not in tags:
            tags.append(t)
    tags = tags[:12]

    return {
        "summary": summary,
        "claims": claims,
        "methods": methods,
        "diseases": diseases,
        "topics": topics,
        "concepts": concepts,
        "entities": entities,
        "tags": tags,
    }



def compile_publication(
    pub: dict[str, Any],
    *,
    mode: Optional[str] = None,
    llm_json: Optional[str | Path] = None,
    allow_heuristic_fallback: bool = False,
) -> dict[str, Any]:
    """Compile one publication. Default mode = agent (see wiki_llm.resolve_compile_mode)."""
    from wiki_llm import compile_with_provider, resolve_compile_mode

    resolved = resolve_compile_mode(mode, llm_json=str(llm_json) if llm_json else None)
    ex = compile_with_provider(
        pub,
        mode=resolved,
        llm_json=llm_json,
        allow_heuristic_fallback=allow_heuristic_fallback,
        heuristic_fn=compile_extract,
    )
    backend = ex.get("compile_backend") or resolved
    ex["compile_backend"] = backend
    ex["compile_label"] = f"{backend}-v0.4.0"
    return ex


# ---------------------------------------------------------------------------
# Page writers
# ---------------------------------------------------------------------------

def source_page_markdown(
    pub: dict[str, Any],
    *,
    safe: str,
    relative_raw: str,
    extracted: Optional[dict[str, Any]] = None,
    compile_label: str = "agent-v0.4.0",
) -> str:
    title = pub.get("title") or pub.get("id") or safe
    needs = int(pub.get("needs_review") or 0)
    banner = (
        "\n> **needs_review**: soft ID conflict or conflicting strong IDs — do not auto-merge.\n"
        if needs
        else ""
    )
    tags = pub.get("tags_json") or "[]"
    authors = pub.get("authors_json") or "[]"
    abstract = pub.get("abstract") or ""
    links = []
    if pub.get("doi"):
        links.append(f"- DOI: `{pub['doi']}`")
    if pub.get("pmid"):
        links.append(f"- PMID: `{pub['pmid']}`")
    if pub.get("pmcid"):
        links.append(f"- PMCID: `{pub['pmcid']}`")
    if pub.get("arxiv"):
        links.append(f"- arXiv: `{pub['arxiv']}`")
    link_block = "\n".join(links) if links else "_No strong IDs._"

    ex = extracted or compile_extract(str(title), str(abstract))
    summary = ex.get("summary") or "_No abstract available for condensation._"
    claims = ex.get("claims") or []
    methods = ex.get("methods") or []
    concepts = ex.get("concepts") or []
    entities = ex.get("entities") or []

    claim_block = (
        "\n".join(f"- {_md_escape(c)}" for c in claims)
        if claims
        else "_No claim sentences extracted._"
    )
    method_block = (
        ", ".join(f"`{m}`" for m in methods) if methods else "_None detected._"
    )

    concept_links = []
    for c in concepts[:12]:
        concept_links.append(f"- [`{c}`](../concepts/{slugify(c)}.md)")
    concept_block = "\n".join(concept_links) if concept_links else "_No concepts extracted._"

    entity_links = []
    for e in entities[:16]:
        entity_links.append(f"- [`{e}`](../entities/{slugify(e)}.md)")
    entity_block = "\n".join(entity_links) if entity_links else "_No entities extracted._"

    return f"""# {_md_escape(str(title))}
{banner}
- publication_id: `{pub.get('id')}`
- year: {pub.get('year') or ''}
- venue: {pub.get('venue') or ''}
- type: {pub.get('type') or ''}
- needs_review: {needs}
- tags_json: `{tags}`
- authors_json: `{authors}`
- raw: [`{relative_raw}`](../../raw/{relative_raw})
- compile: `{compile_label}`

## Identifiers

{link_block}

## Summary

{_md_escape(str(summary))}

## Key claims

{claim_block}

## Methods

{method_block}

## Abstract

{_md_escape(str(abstract)) or '_No abstract in SQLite._'}

## Concepts

{concept_block}

## Entities

{entity_block}

## Related sources

<!-- Link other wiki/sources pages here -->
"""


def ensure_section_link(index_path: Path, section: str, link_text: str, rel_target: str) -> bool:
    """Add a markdown bullet under ## <section> if not already linked."""
    text = index_path.read_text(encoding="utf-8") if index_path.exists() else INDEX_MD
    needle = f"]({rel_target})"
    if needle in text or f"](./{rel_target})" in text:
        return False
    bullet = f"- [{link_text}]({rel_target})"
    marker = f"## {section}"
    if marker not in text:
        text = text.rstrip() + f"\n\n{marker}\n\n{bullet}\n"
        index_path.write_text(text, encoding="utf-8")
        return True
    parts = text.split(marker, 1)
    head, rest = parts[0] + marker, parts[1]
    lines = rest.splitlines(keepends=True)
    i = 0
    while i < len(lines) and (lines[i].strip() == "" or lines[i].lstrip().startswith("<!--")):
        i += 1
    new_rest = "".join(lines[:i])
    if not new_rest.endswith("\n"):
        new_rest += "\n"
    new_rest += bullet + "\n" + "".join(lines[i:])
    index_path.write_text(head + new_rest, encoding="utf-8")
    return True


def ensure_index_link(index_path: Path, link_text: str, rel_target: str) -> bool:
    """Add a markdown bullet under ## Unfiled if not already linked. Returns True if wrote."""
    return ensure_section_link(index_path, "Unfiled", link_text, rel_target)


def concept_contrib_dir(vault: Path, slug: str) -> Path:
    return vault / "wiki" / "concepts" / "_contrib" / slug


def list_concept_contributions(vault: Path, slug: str) -> list[dict[str, Any]]:
    """Load per-source contribution JSON files for a concept slug (sorted by source_safe_id)."""
    d = concept_contrib_dir(vault, slug)
    if not d.is_dir():
        return []
    out: list[dict[str, Any]] = []
    for path in sorted(d.glob("*.json")):
        try:
            obj = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if isinstance(obj, dict):
            obj = dict(obj)
            obj["_path"] = str(path.relative_to(vault))
            out.append(obj)
    return out


def write_concept_contribution(
    vault: Path,
    *,
    slug: str,
    concept_name: str,
    source_id: str,
    source_safe_id: str,
    source_title: str,
    source_rel: str,
    summary_snippet: str = "",
    claims: Optional[list[str]] = None,
    evidence: str = "",
    needs_review: bool = False,
) -> Path:
    """Deterministic per-source concept contribution (no LLM)."""
    d = concept_contrib_dir(vault, slug)
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{source_safe_id}.json"
    claim_hits: list[str] = []
    for c in claims or []:
        if not isinstance(c, str) or not c.strip():
            continue
        low = c.lower()
        token = concept_name.lower().split()[0] if concept_name else ""
        if concept_name.lower() in low or (token and token in low):
            claim_hits.append(c.strip())
        elif len(claim_hits) < 2 and evidence and evidence.lower() in low:
            claim_hits.append(c.strip())
    if not claim_hits and claims:
        # keep up to 2 claims as weak grounding when none mention the concept
        claim_hits = [c.strip() for c in claims if isinstance(c, str) and c.strip()][:2]
    payload = {
        "concept_name": concept_name,
        "slug": slug,
        "source_id": source_id,
        "source_safe_id": source_safe_id,
        "source_title": source_title,
        "source_rel": source_rel,
        "summary_snippet": (summary_snippet or "")[:400],
        "claims": claim_hits[:6],
        "evidence": (evidence or "").strip()[:500],
        "needs_review": bool(needs_review),
        "written_at": utc_now(),
    }
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return path


def write_concept_merge_pending_stub(
    vault: Path,
    *,
    slug: str,
    concept_name: str,
    source_ids: list[str],
) -> Path:
    """Soft stub so agent workflows can pick up pending merges (no invented synthesis)."""
    raw = vault / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    path = raw / f"concept_merge_request_{slug}.json"
    stub = {
        "ok": True,
        "action": "concept_merge_pending",
        "schema_version": "wiki-concept-merge-v1",
        "slug": slug,
        "concept_name": concept_name,
        "source_ids": source_ids,
        "n_sources": len(source_ids),
        "instructions": (
            "This concept has contributions from ≥2 sources. "
            "Run concept_merge_request.py for a grounded agent prompt, "
            "fill concept_merge_response JSON with your model, then concept_merge_apply.py. "
            "Offline: concept_merge_apply.py --mock. Scripts do not invent merge text."
        ),
        "next_step": (
            f"python skills/lit-wiki-ingest/scripts/concept_merge_request.py "
            f"--vault {vault} --slug {slug}"
        ),
        "written_at": utc_now(),
    }
    path.write_text(json.dumps(stub, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return path


def concept_page_source_count(path: Path) -> int:
    if not path.exists():
        return 0
    text = path.read_text(encoding="utf-8")
    if "## Sources" not in text:
        return 0
    _, _, rest = text.partition("## Sources")
    # stop at next heading
    body = rest.split("\n## ", 1)[0]
    return sum(1 for line in body.splitlines() if line.strip().startswith("- ["))


def _upsert_linked_page(
    path: Path,
    *,
    title: str,
    slug: str,
    kind: str,
    source_rel: str,
    source_title: str,
    blurb: str,
) -> bool:
    """Create or update a concept/entity page; append source link. Returns True if created."""
    created = False
    if not path.exists():
        if kind == "concept":
            body = (
                f"# {title}\n\n"
                f"Slug: `{slug}`\n\n"
                f"## Definition\n\n"
                f"{blurb or '_Auto-created from wiki compile; refine as needed._'}\n\n"
                f"## Sources\n\n"
                f"- [{source_title}]({source_rel})\n\n"
                f"## Related concepts\n\n"
            )
        else:
            body = (
                f"# {title}\n\n"
                f"Slug: `{slug}`\n\n"
                f"## Notes\n\n"
                f"{blurb or '_Auto-created entity page from wiki compile._'}\n\n"
                f"## Sources\n\n"
                f"- [{source_title}]({source_rel})\n"
            )
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(body, encoding="utf-8")
        return True

    text = path.read_text(encoding="utf-8")
    link_needle = f"]({source_rel})"
    if link_needle not in text:
        marker = "## Sources\n"
        bullet = f"- [{source_title}]({source_rel})"
        if marker in text:
            head, _, rest = text.partition(marker)
            if rest.startswith("\n"):
                rest = rest[1:]
            text = head + marker + "\n" + bullet + "\n" + rest
        else:
            text = text.rstrip() + f"\n\n## Sources\n\n{bullet}\n"
        path.write_text(text, encoding="utf-8")
    return created


def concept_stub(slug: str, name: str) -> str:
    return f"""# {name}

Slug: `{slug}`

## Definition

_Auto-created concept page from wiki compile._

## Sources

## Related concepts

"""


def entity_stub(slug: str, name: str) -> str:
    return f"""# {name}

Slug: `{slug}`

## Notes

_Auto-created entity page from wiki compile._

## Sources

"""


def ingest_publication_to_wiki(
    vault: Path,
    pub: dict[str, Any],
    *,
    raw_file: Optional[Path] = None,
    raw_text: Optional[str] = None,
    extra_tags: Optional[list[str]] = None,
    compile_mode: Optional[str] = None,
    llm_json: Optional[str | Path] = None,
    allow_heuristic_fallback: bool = False,
) -> dict[str, Any]:
    """Copy into raw/, write compiled wiki/sources page, concepts/entities, update index+log.

    Returns suggested_tags for caller to persist into SQLite when tags empty.
    """
    pub_id = str(pub.get("id") or "")
    if not pub_id:
        return {"ok": False, "error": "publication missing id"}
    sid = safe_id(pub_id)
    raw_dir = vault / "raw" / "publications" / sid
    raw_dir.mkdir(parents=True, exist_ok=True)

    title = str(pub.get("title") or pub_id)
    abstract = str(pub.get("abstract") or "")
    extracted = compile_publication(
        pub,
        mode=compile_mode,
        llm_json=llm_json,
        allow_heuristic_fallback=allow_heuristic_fallback,
    )
    compile_label = str(extracted.get("compile_label") or "llm-v0.3.0")

    # Merge existing tags with suggested
    existing_tags = _parse_tags(pub.get("tags_json") or pub.get("tags"))
    suggested_tags = list(extracted.get("tags") or [])
    if extra_tags:
        for t in extra_tags:
            if t and t not in suggested_tags:
                suggested_tags.append(t)
    effective_tags = existing_tags if existing_tags else suggested_tags
    tags_filled = False
    if tags_are_empty(pub.get("tags_json") or pub.get("tags")) and suggested_tags:
        pub = dict(pub)
        pub["tags_json"] = json.dumps(suggested_tags, ensure_ascii=False)
        tags_filled = True
        effective_tags = suggested_tags

    meta = {
        "publication_id": pub_id,
        "safe_id": sid,
        "ingested_at": utc_now(),
        "doi": pub.get("doi"),
        "pmid": pub.get("pmid"),
        "title": pub.get("title"),
        "year": pub.get("year"),
        "needs_review": pub.get("needs_review"),
        "compile": compile_label,
        "compile_backend": extracted.get("compile_backend"),
        "n_concepts": len(extracted.get("concepts") or []),
        "n_entities": len(extracted.get("entities") or []),
        "n_claims": len(extracted.get("claims") or []),
        "tags": effective_tags,
    }
    (raw_dir / "meta.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    body_written = None
    if raw_file and Path(raw_file).is_file():
        dest = raw_dir / Path(raw_file).name
        shutil.copy2(raw_file, dest)
        body_written = dest.name
    elif raw_text is not None:
        (raw_dir / "body.txt").write_text(raw_text, encoding="utf-8")
        body_written = "body.txt"
    else:
        snap = []
        if pub.get("title"):
            snap.append(f"Title: {pub['title']}")
        if pub.get("abstract"):
            snap.append(f"\nAbstract:\n{pub['abstract']}")
        if effective_tags:
            snap.append(f"\nTags: {', '.join(effective_tags)}")
        (raw_dir / "body.txt").write_text("\n".join(snap) + "\n", encoding="utf-8")
        body_written = "body.txt"

    rel_raw = f"publications/{sid}/{body_written}"
    page = vault / "wiki" / "sources" / f"{sid}.md"
    page.write_text(
        source_page_markdown(
            pub, safe=sid, relative_raw=rel_raw, extracted=extracted, compile_label=compile_label
        ),
        encoding="utf-8",
    )

    short_title = title if len(title) <= 80 else title[:77] + "…"
    source_rel_from_concept = f"../sources/{sid}.md"

    concept_pages_created: list[str] = []
    concept_pages_updated: list[str] = []
    concept_merge_pending: list[str] = []
    concept_contrib_paths: list[str] = []
    details_by_name: dict[str, dict[str, Any]] = {}
    for d in extracted.get("concept_details") or []:
        if isinstance(d, dict) and d.get("name"):
            details_by_name[str(d["name"])] = d
    needs_flag = bool(int(pub.get("needs_review") or 0))
    # High-signal: always create for extracted concepts (topics/diseases/methods)
    for cname in (extracted.get("concepts") or [])[:12]:
        cslug = slugify(cname)
        cpath = vault / "wiki" / "concepts" / f"{cslug}.md"
        prior_sources = concept_page_source_count(cpath)
        prior_contribs = list_concept_contributions(vault, cslug)
        blurb = ""
        detail = details_by_name.get(cname) or {}
        evidence = str(detail.get("evidence") or detail.get("note") or "").strip()
        # Prefer a claim that mentions this concept
        for claim in extracted.get("claims") or []:
            if cname.lower().split()[0] in claim.lower() or cname.lower() in claim.lower():
                blurb = claim
                break
        if not blurb and evidence:
            blurb = evidence[:280]
        if not blurb and extracted.get("summary"):
            blurb = str(extracted["summary"])[:280]
        created = _upsert_linked_page(
            cpath,
            title=cname,
            slug=cslug,
            kind="concept",
            source_rel=source_rel_from_concept,
            source_title=short_title,
            blurb=blurb,
        )
        contrib_path = write_concept_contribution(
            vault,
            slug=cslug,
            concept_name=cname,
            source_id=pub_id,
            source_safe_id=sid,
            source_title=short_title,
            source_rel=source_rel_from_concept,
            summary_snippet=str(extracted.get("summary") or "")[:400],
            claims=list(extracted.get("claims") or []),
            evidence=evidence or blurb,
            needs_review=needs_flag,
        )
        concept_contrib_paths.append(str(contrib_path.relative_to(vault)))
        # Pending when this concept already had ≥1 other source/contrib
        other_ids = {
            str(c.get("source_safe_id") or "")
            for c in prior_contribs
            if str(c.get("source_safe_id") or "") and str(c.get("source_safe_id")) != sid
        }
        if prior_sources >= 1 or other_ids:
            if cslug not in concept_merge_pending:
                concept_merge_pending.append(cslug)
            write_concept_merge_pending_stub(
                vault,
                slug=cslug,
                concept_name=cname,
                source_ids=sorted(other_ids | {sid}),
            )
        rel = str(cpath.relative_to(vault))
        if created:
            concept_pages_created.append(rel)
        else:
            concept_pages_updated.append(rel)

    entity_pages_created: list[str] = []
    entity_pages_updated: list[str] = []
    for ename in (extracted.get("entities") or [])[:16]:
        eslug = slugify(ename)
        epath = vault / "wiki" / "entities" / f"{eslug}.md"
        blurb = f"Mentioned in `{short_title}`."
        created = _upsert_linked_page(
            epath,
            title=ename,
            slug=eslug,
            kind="entity",
            source_rel=source_rel_from_concept,
            source_title=short_title,
            blurb=blurb,
        )
        rel = str(epath.relative_to(vault))
        if created:
            entity_pages_created.append(rel)
        else:
            entity_pages_updated.append(rel)

    index_path = vault / "wiki" / "index.md"
    linked = ensure_index_link(index_path, short_title, f"sources/{sid}.md")
    # Also link top concepts/entities into index sections
    for cname in (extracted.get("concepts") or [])[:6]:
        ensure_section_link(index_path, "Concepts", cname, f"concepts/{slugify(cname)}.md")
    for ename in (extracted.get("entities") or [])[:6]:
        ensure_section_link(index_path, "Entities", ename, f"entities/{slugify(ename)}.md")

    append_log(
        vault,
        f"ingest source `{sid}` (`{pub_id}`) raw={body_written} "
        f"claims={len(extracted.get('claims') or [])} "
        f"concepts=+{len(concept_pages_created)} entities=+{len(entity_pages_created)} "
        f"merge_pending={len(concept_merge_pending)} "
        f"tags_filled={tags_filled} index_link={linked} compile={compile_label}",
    )

    return {
        "ok": True,
        "publication_id": pub_id,
        "safe_id": sid,
        "source_page": str(page.relative_to(vault)),
        "raw_dir": str(raw_dir.relative_to(vault)),
        "raw_body": body_written,
        "index_linked": linked,
        "concept_pages_created": concept_pages_created,
        "concept_pages_updated": concept_pages_updated,
        "concept_merge_pending": concept_merge_pending,
        "concept_contrib_paths": concept_contrib_paths,
        "entity_pages_created": entity_pages_created,
        "entity_pages_updated": entity_pages_updated,
        "n_claims": len(extracted.get("claims") or []),
        "n_methods": len(extracted.get("methods") or []),
        "summary_preview": (extracted.get("summary") or "")[:160],
        "suggested_tags": suggested_tags,
        "tags_filled": tags_filled,
        "effective_tags": effective_tags,
        "compile": compile_label,
        "compile_backend": extracted.get("compile_backend"),
    }


def extract_md_links(text: str) -> list[str]:
    # Ignore markdown inside HTML comments (template placeholders)
    cleaned = re.sub(r"<!--.*?-->", "", text, flags=re.S)
    return re.findall(r"\[[^\]]*\]\(([^)]+)\)", cleaned)


def lint_vault(vault: Path) -> dict[str, Any]:
    """Check orphans, missing index entries, broken relative links, needs_review banners."""
    issues: list[dict[str, Any]] = []
    wiki = vault / "wiki"
    sources_dir = wiki / "sources"
    index_path = wiki / "index.md"
    index_text = index_path.read_text(encoding="utf-8") if index_path.exists() else ""

    source_files = sorted(sources_dir.glob("*.md")) if sources_dir.exists() else []
    source_names = {p.name for p in source_files}

    # Missing index entries
    for p in source_files:
        rel = f"sources/{p.name}"
        if rel not in index_text and f"./{rel}" not in index_text:
            issues.append({"code": "missing_index_entry", "path": rel, "severity": "error"})

    # Orphan concepts/entities (no inbound link from any wiki md)
    all_md = list(wiki.rglob("*.md")) if wiki.exists() else []
    inbound: set[str] = set()
    broken: list[dict[str, Any]] = []
    for md in all_md:
        text = md.read_text(encoding="utf-8")
        for link in extract_md_links(text):
            if link.startswith(("http://", "https://", "mailto:")):
                continue
            # strip anchors
            link_path = link.split("#", 1)[0]
            if not link_path:
                continue
            target = (md.parent / link_path).resolve()
            try:
                target.relative_to(vault.resolve())
            except ValueError:
                broken.append({"code": "link_outside_vault", "from": str(md.relative_to(vault)), "link": link})
                continue
            if not target.exists():
                broken.append({"code": "broken_link", "from": str(md.relative_to(vault)), "link": link})
            else:
                try:
                    inbound.add(str(target.relative_to(vault.resolve())))
                except ValueError:
                    pass

    for code_issue in broken:
        issues.append({**code_issue, "severity": "error"})

    for kind in ("concepts", "entities"):
        d = wiki / kind
        if not d.exists():
            continue
        for p in d.glob("*.md"):
            rel = str(p.relative_to(vault))
            if rel not in inbound:
                if not any(p.name in i or rel in i for i in inbound):
                    if rel not in inbound:
                        issues.append({"code": "orphan_page", "path": rel, "severity": "warn"})

    # Noisy entity fragments (PD-L1 → l1, IFN-γ → ifn)
    ent_dir = wiki / "entities"
    if ent_dir.exists():
        for p in ent_dir.glob("*.md"):
            slug = p.stem.lower()
            if slug in NOISY_ENTITY_SLUGS:
                issues.append(
                    {
                        "code": "noisy_entity",
                        "path": str(p.relative_to(vault)),
                        "slug": slug,
                        "severity": "error",
                    }
                )

    # needs_review flags on source pages vs raw meta
    for p in source_files:
        text = p.read_text(encoding="utf-8")
        meta_path = vault / "raw" / "publications" / p.stem / "meta.json"
        nr = 0
        if meta_path.exists():
            try:
                nr = int(json.loads(meta_path.read_text(encoding="utf-8")).get("needs_review") or 0)
            except (json.JSONDecodeError, ValueError):
                issues.append({"code": "bad_raw_meta", "path": str(meta_path.relative_to(vault)), "severity": "warn"})
        if nr and "needs_review" not in text:
            issues.append({"code": "stale_needs_review_banner", "path": str(p.relative_to(vault)), "severity": "warn"})
        # Placeholder summary should not remain after 0.2.1 compile
        if "Summary (placeholder)" in text or "_Stub summary generated by lit-wiki-ingest" in text:
            issues.append({"code": "placeholder_summary", "path": str(p.relative_to(vault)), "severity": "warn"})

    # Raw without wiki source
    raw_pubs = vault / "raw" / "publications"
    if raw_pubs.exists():
        for d in raw_pubs.iterdir():
            if d.is_dir() and f"{d.name}.md" not in source_names:
                issues.append({"code": "raw_without_source_page", "path": str(d.relative_to(vault)), "severity": "error"})

    errors = sum(1 for i in issues if i.get("severity") == "error")
    warns = sum(1 for i in issues if i.get("severity") == "warn")
    return {
        "ok": True,
        "n_issues": len(issues),
        "n_errors": errors,
        "n_warns": warns,
        "n_sources": len(source_files),
        "issues": issues,
    }


def query_wiki_markdown(vault: Path, q: str, *, limit: int = 50) -> list[dict[str, Any]]:
    """Keyword search over wiki markdown (case-insensitive). No vectors."""
    q = (q or "").strip()
    if not q:
        return []
    limit = max(1, min(int(limit), 200))
    wiki = vault / "wiki"
    if not wiki.exists():
        return []
    terms = [t for t in re.split(r"\s+", q.lower()) if t]
    hits: list[dict[str, Any]] = []
    for md in sorted(wiki.rglob("*.md")):
        text = md.read_text(encoding="utf-8")
        low = text.lower()
        if not all(t in low for t in terms):
            continue
        # snippet: first line with a term
        snippet = ""
        for line in text.splitlines():
            if any(t in line.lower() for t in terms):
                snippet = line.strip()[:240]
                break
        if not snippet:
            snippet = text.strip().splitlines()[0][:240] if text.strip() else ""
        hits.append(
            {
                "path": str(md.relative_to(vault)),
                "title": next((ln[2:].strip() for ln in text.splitlines() if ln.startswith("# ")), md.stem),
                "snippet": snippet,
            }
        )
        if len(hits) >= limit:
            break
    return hits


def suggest_theme_tree(
    pubs: list[dict[str, Any]],
    *,
    max_branches: int = 12,
) -> dict[str, Any]:
    """Propose theme branches from tags + title keywords. Proposals only — no writes."""
    tag_counter: Counter[str] = Counter()
    title_words: Counter[str] = Counter()
    stop = {
        "a", "an", "the", "of", "in", "on", "and", "or", "for", "to", "with", "by",
        "from", "as", "at", "is", "are", "be", "this", "that", "into", "via", "using",
        "study", "analysis", "based", "new", "novel",
    }
    pub_ids_by_tag: dict[str, list[str]] = {}

    for pub in pubs:
        pid = str(pub.get("id") or "")
        tags = _parse_tags(pub.get("tags_json") or pub.get("tags") or [])
        for t in tags:
            t_norm = t.strip()
            if not t_norm:
                continue
            tag_counter[t_norm] += 1
            pub_ids_by_tag.setdefault(t_norm, []).append(pid)

        title = str(pub.get("title") or "")
        for w in re.findall(r"[A-Za-z][A-Za-z0-9\-]{2,}", title):
            lw = w.lower()
            if lw in stop:
                continue
            title_words[lw] += 1

    proposals: list[dict[str, Any]] = []
    # Prefer tags
    for tag, n in tag_counter.most_common(max_branches):
        proposals.append(
            {
                "action": "add_branch",
                "name": tag,
                "slug": slugify(tag),
                "parent": "Themes",
                "evidence": {"kind": "tag", "count": n, "publication_ids": pub_ids_by_tag.get(tag, [])[:20]},
                "attach_publication_ids": pub_ids_by_tag.get(tag, [])[:50],
            }
        )

    # Fill with title keywords if needed
    if len(proposals) < max_branches:
        for word, n in title_words.most_common(max_branches * 2):
            if n < 2:
                continue
            name = word
            if any(p["slug"] == slugify(name) for p in proposals):
                continue
            matched = [
                str(p.get("id"))
                for p in pubs
                if word in str(p.get("title") or "").lower()
            ]
            proposals.append(
                {
                    "action": "add_branch",
                    "name": name,
                    "slug": slugify(name),
                    "parent": "Themes",
                    "evidence": {"kind": "title_keyword", "count": n, "publication_ids": matched[:20]},
                    "attach_publication_ids": matched[:50],
                }
            )
            if len(proposals) >= max_branches:
                break

    return {
        "ok": True,
        "n_publications": len(pubs),
        "n_proposals": len(proposals),
        "proposals": proposals,
        "note": "Proposals only. Apply with lit-wiki-apply-tree --confirm after human review.",
    }


def apply_theme_tree(
    vault: Path,
    proposals: list[dict[str, Any]],
    *,
    confirm: bool,
    conn: Any = None,
    collection_id: Optional[str] = None,
) -> dict[str, Any]:
    """Write theme branches into wiki/index.md and optionally theme_nodes/theme_pubs. Requires confirm."""
    if not confirm:
        return {
            "ok": True,
            "wrote": False,
            "reason": "missing --confirm",
            "n_proposals": len(proposals),
            "message": "Dry-run only. Re-run with --confirm to write index + theme_nodes.",
        }

    index_path = vault / "wiki" / "index.md"
    text = index_path.read_text(encoding="utf-8") if index_path.exists() else INDEX_MD
    if "## Themes" not in text:
        text = text.rstrip() + "\n\n## Themes\n\n"

    applied: list[dict[str, Any]] = []
    # Split at ## Themes
    pre, _, post = text.partition("## Themes")
    themes_body = post
    m = re.search(r"\n## [^#]", themes_body)
    if m:
        themes_part = themes_body[: m.start()]
        after = themes_body[m.start() :]
    else:
        themes_part = themes_body
        after = ""

    for prop in proposals:
        if prop.get("action") != "add_branch":
            continue
        name = prop.get("name") or prop.get("slug")
        slug = prop.get("slug") or slugify(str(name))
        heading = f"### {name}"
        if heading in text or f"### {name}\n" in themes_part:
            applied.append({"name": name, "slug": slug, "status": "exists"})
            continue
        lines = [f"\n{heading}\n"]
        for pid in prop.get("attach_publication_ids") or []:
            sid = safe_id(str(pid))
            src = vault / "wiki" / "sources" / f"{sid}.md"
            if src.exists():
                first = src.read_text(encoding="utf-8").splitlines()
                title = next((ln[2:].strip() for ln in first if ln.startswith("# ")), sid)
                short = title if len(title) <= 80 else title[:77] + "…"
                lines.append(f"- [{short}](sources/{sid}.md)\n")
            else:
                lines.append(f"- `{pid}` _(no wiki source page yet)_\n")
        themes_part = themes_part.rstrip() + "\n" + "".join(lines)
        applied.append({"name": name, "slug": slug, "status": "added", "n_links": len(prop.get("attach_publication_ids") or [])})

        # Optional SQLite theme_nodes
        if conn is not None and collection_id:
            node_id = f"theme-{collection_id}-{slug}"[:120]
            parent = f"theme-root-{collection_id}"
            existing = conn.execute("SELECT id FROM theme_nodes WHERE id = ?", (node_id,)).fetchone()
            if not existing:
                conn.execute(
                    "INSERT INTO theme_nodes(id, collection_id, parent_id, name, sort) VALUES (?, ?, ?, ?, ?)",
                    (node_id, collection_id, parent, name, len(applied)),
                )
            for pid in prop.get("attach_publication_ids") or []:
                if conn.execute("SELECT 1 FROM publications WHERE id = ?", (str(pid),)).fetchone():
                    conn.execute(
                        "INSERT OR IGNORE INTO theme_pubs(node_id, publication_id) VALUES (?, ?)",
                        (node_id, str(pid)),
                    )

    new_text = pre + "## Themes" + themes_part
    if after and not after.startswith("\n"):
        new_text += "\n" + after
    else:
        new_text += after
    index_path.write_text(new_text if new_text.endswith("\n") else new_text + "\n", encoding="utf-8")
    append_log(vault, f"apply-tree confirmed branches={len(applied)}")
    if conn is not None:
        conn.commit()

    return {
        "ok": True,
        "wrote": True,
        "n_applied": len(applied),
        "applied": applied,
        "index_path": str(index_path),
    }

# ---------------------------------------------------------------------------
# Identity remap → wiki source rename (0.4.8)
# ---------------------------------------------------------------------------


def _replace_publication_id_in_source_md(text: str, old_id: str, new_id: str, old_safe: str, new_safe: str) -> str:
    """Update publication_id bullet and safe_id path references inside a source page."""
    out = text
    out = out.replace(f"- publication_id: `{old_id}`", f"- publication_id: `{new_id}`")
    out = out.replace(f"publications/{old_safe}/", f"publications/{new_safe}/")
    out = out.replace(f"sources/{old_safe}.md", f"sources/{new_safe}.md")
    out = out.replace(f"`{old_safe}`", f"`{new_safe}`")
    return out


def _replace_links_in_text(text: str, old_safe: str, new_safe: str) -> str:
    if old_safe == new_safe:
        return text
    return (
        text.replace(f"sources/{old_safe}.md", f"sources/{new_safe}.md")
        .replace(f"]({old_safe}.md)", f"]({new_safe}.md)")
        .replace(f"/{old_safe}.md", f"/{new_safe}.md")
        .replace(f"publications/{old_safe}/", f"publications/{new_safe}/")
    )


def wiki_repoint_publication(
    project_root: Path,
    old_id: str,
    new_id: str,
    *,
    collection_id: Optional[str] = None,
    db_path: Optional[Path] = None,
    leave_stub: bool = True,
) -> dict[str, Any]:
    """Rename wiki source/raw paths when a publication canonical id changes.

    For each matching collection vault:
      - move wiki/sources/<old_safe>.md → <new_safe>.md (update publication_id)
      - optionally leave a stub redirect at the old path
      - rename raw/publications/<old_safe>/ when present
      - rename concept contrib JSON keyed by old_safe
      - cheap string replace in index.md / log.md
    """
    root = Path(project_root).resolve()
    if db_path is None:
        db_path = root / "library.sqlite"
    old_safe = safe_id(old_id)
    new_safe = safe_id(new_id)
    if old_id == new_id or old_safe == new_safe:
        return {
            "ok": True,
            "skipped": True,
            "reason": "same_id",
            "old_id": old_id,
            "new_id": new_id,
            "vaults": [],
        }

    col_dirs: list[tuple[str, Path]] = []
    import sqlite3

    if db_path.exists():
        conn = sqlite3.connect(str(db_path))
        conn.row_factory = sqlite3.Row
        try:
            if collection_id:
                rows = conn.execute(
                    "SELECT id, path FROM collections WHERE id=?", (collection_id,)
                ).fetchall()
            else:
                rows = conn.execute("SELECT id, path FROM collections").fetchall()
            for r in rows:
                p = Path(r["path"]) if r["path"] else root / "collections" / r["id"]
                if not p.is_absolute():
                    p = root / p
                col_dirs.append((r["id"], p))
        finally:
            conn.close()
    elif collection_id:
        col_dirs.append((collection_id, root / "collections" / collection_id))

    vault_results: list[dict[str, Any]] = []
    for cid, col_dir in col_dirs:
        vault = vault_root(col_dir)
        if not vault.is_dir():
            continue
        sources_dir = vault / "wiki" / "sources"
        old_page = sources_dir / f"{old_safe}.md"
        new_page = sources_dir / f"{new_safe}.md"
        raw_old = vault / "raw" / "publications" / old_safe
        raw_new = vault / "raw" / "publications" / new_safe
        actions: list[str] = []

        if old_page.is_file():
            body = old_page.read_text(encoding="utf-8")
            body = _replace_publication_id_in_source_md(body, old_id, new_id, old_safe, new_safe)
            sources_dir.mkdir(parents=True, exist_ok=True)
            if new_page.exists() and new_page.resolve() != old_page.resolve():
                actions.append("new_source_exists")
            else:
                new_page.write_text(body, encoding="utf-8")
                actions.append("wrote_new_source")
            if leave_stub:
                stub = (
                    "# Redirect\n\n"
                    "This source moved after identity upgrade/merge.\n\n"
                    f"- old publication_id: `{old_id}`\n"
                    f"- new publication_id: `{new_id}`\n"
                    f"- see: [`{new_safe}.md`]({new_safe}.md)\n"
                )
                old_page.write_text(stub, encoding="utf-8")
                actions.append("left_stub")
            else:
                if old_page.resolve() != new_page.resolve():
                    old_page.unlink()
                    actions.append("removed_old_source")
        else:
            actions.append("no_old_source")

        if raw_old.is_dir():
            raw_new.parent.mkdir(parents=True, exist_ok=True)
            if raw_new.exists():
                for child in raw_old.iterdir():
                    dest = raw_new / child.name
                    if not dest.exists():
                        shutil.move(str(child), str(dest))
                shutil.rmtree(raw_old, ignore_errors=True)
                actions.append("merged_raw")
            else:
                shutil.move(str(raw_old), str(raw_new))
                actions.append("moved_raw")
            meta = raw_new / "meta.json"
            if meta.is_file():
                try:
                    data = json.loads(meta.read_text(encoding="utf-8"))
                    data["publication_id"] = new_id
                    data["safe_id"] = new_safe
                    data["repointed_from"] = old_id
                    meta.write_text(
                        json.dumps(data, ensure_ascii=False, indent=2) + "\n",
                        encoding="utf-8",
                    )
                    actions.append("updated_raw_meta")
                except (json.JSONDecodeError, OSError):
                    actions.append("raw_meta_update_failed")

        contrib_root = vault / "wiki" / "concepts" / "_contrib"
        n_contrib = 0
        if contrib_root.is_dir():
            for slug_dir in contrib_root.iterdir():
                if not slug_dir.is_dir():
                    continue
                old_c = slug_dir / f"{old_safe}.json"
                new_c = slug_dir / f"{new_safe}.json"
                if not old_c.is_file():
                    continue
                try:
                    payload = json.loads(old_c.read_text(encoding="utf-8"))
                except json.JSONDecodeError:
                    payload = {}
                payload["source_safe_id"] = new_safe
                payload["source_id"] = new_id
                payload["publication_id"] = new_id
                if not new_c.exists():
                    new_c.write_text(
                        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
                        encoding="utf-8",
                    )
                old_c.unlink()
                n_contrib += 1
        if n_contrib:
            actions.append(f"repointed_contribs:{n_contrib}")

        for name in ("index.md", "log.md"):
            p = vault / "wiki" / name
            if not p.is_file():
                continue
            before = p.read_text(encoding="utf-8")
            after = _replace_links_in_text(before, old_safe, new_safe)
            if after != before:
                p.write_text(after, encoding="utf-8")
                actions.append(f"updated_{name}")

        if any(a for a in actions if a not in ("no_old_source",)):
            append_log(
                vault,
                f"wiki_repoint {old_id} → {new_id} ({old_safe} → {new_safe})",
            )
            actions.append("logged")

        vault_results.append(
            {
                "collection_id": cid,
                "vault": str(vault),
                "old_safe": old_safe,
                "new_safe": new_safe,
                "actions": actions,
            }
        )

    return {
        "ok": True,
        "old_id": old_id,
        "new_id": new_id,
        "old_safe": old_safe,
        "new_safe": new_safe,
        "n_vaults": len(vault_results),
        "vaults": vault_results,
    }


def wiki_repoint_remaps(
    project_root: Path,
    remaps: list[dict[str, str]],
    *,
    collection_id: Optional[str] = None,
    db_path: Optional[Path] = None,
    leave_stub: bool = True,
) -> dict[str, Any]:
    """Apply wiki_repoint_publication for a list of {old_id,new_id} remaps."""
    results = []
    for r in remaps or []:
        old_id = r.get("old_id") or r.get("from") or ""
        new_id = r.get("new_id") or r.get("to") or ""
        if not old_id or not new_id:
            continue
        results.append(
            wiki_repoint_publication(
                project_root,
                old_id,
                new_id,
                collection_id=collection_id,
                db_path=db_path,
                leave_stub=leave_stub,
            )
        )
    return {
        "ok": True,
        "n": len(results),
        "results": results,
    }
