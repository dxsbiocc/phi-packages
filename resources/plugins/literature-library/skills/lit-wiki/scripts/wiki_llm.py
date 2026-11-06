#!/usr/bin/env python3
"""Wiki compile + concept-merge providers for literature-library (v0.4.5-draft).

Primary path = **agent**: the invoking LiteratureLibrary / host agent uses its
already-running model to fill structured compile JSON. Plugin scripts only
build grounded prompts + schema, validate agent JSON, and write the vault.

Providers:
  - agent (default): emit compile_request.json, or apply --llm-json / response file
  - file: alias of apply-from-path (same JSON shapes)
  - mock: fixture / synthetic offline smoke (tests only)
  - heuristic: emergency offline fallback (explicit flag only)

No plugin-owned LLM API keys. Scripts do not call external LLM HTTP APIs.
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Any, Optional

COMPILE_SCHEMA_VERSION = "wiki-compile-v1"
COMPILE_LABEL_VERSION = "v0.4.0"

SYSTEM_PROMPT = """You compile one scientific publication into a structured wiki page JSON.
Rules:
- Use ONLY the provided title, abstract, and metadata. Do not invent PMIDs, DOIs, sample sizes, p-values, or results absent from the text.
- If a detail is unclear, omit it or note uncertainty in confidence_notes — never fabricate.
- Summary: 2–4 sentences of grounded synthesis (not a verbatim first-paragraph copy).
- key_claims: 3–6 crisp claims supported by the abstract.
- methods / concepts / entities / tags: only terms evidenced in the text.
- entities: gene/protein symbols with kind "gene"|"protein"|"complex"|"other". Map layilin→LAYN. Never emit hyphen fragments (no "l1" from PD-L1; no bare "ifn" from IFN-γ — use IFNG or PD-L1).
- Biology stays descriptive; no clinical advice or treatment recommendations.
- Return a single JSON object matching the schema. No markdown fences.
"""

USER_PROMPT_TEMPLATE = """Compile this publication into JSON with keys:
summary (string), key_claims (string[]), methods (string[]),
concepts ({{name, slug_hint?, evidence?}}[]), entities ({{symbol, name?, kind?}}[]),
tags (string[]), confidence_notes (string, optional).

publication_id: {publication_id}
title: {title}
year: {year}
venue: {venue}
doi: {doi}
pmid: {pmid}

abstract:
{abstract}
"""

RESPONSE_JSON_SCHEMA: dict[str, Any] = {
    "type": "object",
    "required": ["summary", "key_claims"],
    "properties": {
        "summary": {"type": "string", "minLength": 40},
        "key_claims": {"type": "array", "items": {"type": "string"}, "minItems": 1},
        "methods": {"type": "array", "items": {"type": "string"}},
        "concepts": {
            "type": "array",
            "items": {
                "oneOf": [
                    {"type": "string"},
                    {
                        "type": "object",
                        "required": ["name"],
                        "properties": {
                            "name": {"type": "string"},
                            "slug_hint": {"type": "string"},
                            "evidence": {"type": "string"},
                            "note": {"type": "string"},
                        },
                    },
                ]
            },
        },
        "entities": {
            "type": "array",
            "items": {
                "oneOf": [
                    {"type": "string"},
                    {
                        "type": "object",
                        "properties": {
                            "symbol": {"type": "string"},
                            "name": {"type": "string"},
                            "kind": {"type": "string"},
                        },
                    },
                ]
            },
        },
        "tags": {"type": "array", "items": {"type": "string"}},
        "confidence_notes": {"type": "string"},
        "publication_id": {"type": "string"},
    },
}


class CompileRequestNeeded(Exception):
    """Agent mode has no response JSON yet — caller should emit request and stop."""

    def __init__(self, request: dict[str, Any]):
        self.request = request
        super().__init__("agent compile response required")


def resolve_compile_mode(
    cli_mode: Optional[str] = None,
    *,
    llm_json: Optional[str] = None,
) -> str:
    """Resolve compile mode: cli > llm-json implies apply > env > default agent."""
    if cli_mode:
        m = cli_mode.strip().lower()
        # Accept legacy aliases; map discouraged openai/llm HTTP names away from default.
        if m in ("openai", "openai-http", "llm-http"):
            raise ValueError(
                f"compile mode {cli_mode!r} removed: wiki prose is produced by the "
                "invoking agent's model (use --compile agent + --llm-json). "
                "For offline tests use --compile mock; emergency --compile heuristic."
            )
        if m == "llm":
            # Legacy name → agent (no plugin HTTP)
            return "agent"
        if m in ("agent", "heuristic", "mock", "file"):
            return m
        raise ValueError(f"unknown compile mode: {cli_mode}")
    if llm_json:
        return "agent"
    env = (os.environ.get("LIT_WIKI_COMPILE") or "").strip().lower()
    if env in ("openai", "openai-http", "llm-http"):
        raise ValueError(
            "LIT_WIKI_COMPILE=openai/llm-http is not supported; use agent (default)"
        )
    if env == "llm":
        return "agent"
    if env in ("agent", "heuristic", "mock", "file"):
        return env
    return "agent"


def _require_str(v: Any, field: str) -> str:
    if not isinstance(v, str) or not v.strip():
        raise ValueError(f"{field} must be non-empty string")
    return v.strip()


def _require_str_list(v: Any, field: str, *, min_items: int = 0) -> list[str]:
    if not isinstance(v, list):
        raise ValueError(f"{field} must be a list")
    out: list[str] = []
    for i, item in enumerate(v):
        if not isinstance(item, str) or not item.strip():
            raise ValueError(f"{field}[{i}] must be non-empty string")
        out.append(item.strip())
    if len(out) < min_items:
        raise ValueError(f"{field} needs at least {min_items} items")
    return out


def validate_compile_payload(obj: Any) -> dict[str, Any]:
    """Strict schema validation. Returns normalized dict."""
    if not isinstance(obj, dict):
        raise ValueError("compile payload must be an object")
    summary = _require_str(obj.get("summary"), "summary")
    if len(summary) < 40:
        raise ValueError("summary too short")
    key_claims = _require_str_list(obj.get("key_claims"), "key_claims", min_items=1)
    methods = _require_str_list(obj.get("methods") or [], "methods")
    tags = _require_str_list(obj.get("tags") or [], "tags")

    concepts_in = obj.get("concepts") or []
    if not isinstance(concepts_in, list):
        raise ValueError("concepts must be a list")
    concepts: list[dict[str, str]] = []
    for i, c in enumerate(concepts_in):
        if isinstance(c, str):
            name = c.strip()
            if not name:
                raise ValueError(f"concepts[{i}] empty")
            concepts.append({"name": name, "slug_hint": "", "evidence": ""})
        elif isinstance(c, dict):
            name = _require_str(c.get("name"), f"concepts[{i}].name")
            evidence = str(c.get("evidence") or c.get("note") or "").strip()
            concepts.append(
                {
                    "name": name,
                    "slug_hint": str(c.get("slug_hint") or "").strip(),
                    "evidence": evidence,
                }
            )
        else:
            raise ValueError(f"concepts[{i}] must be string or object")

    entities_in = obj.get("entities") or []
    if not isinstance(entities_in, list):
        raise ValueError("entities must be a list")
    entities: list[dict[str, str]] = []
    entity_warnings: list[str] = []
    for i, e in enumerate(entities_in):
        if isinstance(e, str):
            sym = e.strip()
            if not sym:
                raise ValueError(f"entities[{i}] empty")
            name = ""
            kind = "gene"
            note = ""
        elif isinstance(e, dict):
            sym = str(e.get("symbol") or e.get("name") or "").strip()
            if not sym:
                raise ValueError(f"entities[{i}] needs symbol")
            name = str(e.get("name") or "").strip()
            kind = str(e.get("kind") or "gene").strip().lower() or "gene"
            note = str(e.get("note") or e.get("evidence") or "").strip()
        else:
            raise ValueError(f"entities[{i}] must be string or object")
        # Soft-drop speculative / implied entities (do not invent FOXP3-like noise)
        blob = f"{sym} {name} {note}".lower()
        if "implied" in blob or "not claimed" in blob or "not in text" in blob:
            entity_warnings.append(
                f"entities[{i}] stripped (speculative): {sym!r}"
            )
            continue
        entities.append({"symbol": sym, "name": name, "kind": kind})

    notes = obj.get("confidence_notes")
    if notes is not None and not isinstance(notes, str):
        raise ValueError("confidence_notes must be string if present")

    out = {
        "summary": summary,
        "key_claims": key_claims,
        "methods": methods,
        "concepts": concepts,
        "entities": entities,
        "tags": tags,
        "confidence_notes": (notes or "").strip(),
        "schema": COMPILE_SCHEMA_VERSION,
        "publication_id": str(obj.get("publication_id") or "").strip(),
    }
    if entity_warnings:
        out["entity_warnings"] = entity_warnings
    return out


_ENTITY_FRAGMENT_BLOCK = {
    "L1", "L2", "L3", "D1", "D2", "PD", "IFN", "TGF", "CTLA", "TIM", "LAG",
    "NKG", "HLA", "M1", "M2", "IL", "CD", "TH1", "TH2",
}
_ENTITY_ALIAS = {
    "LAYILIN": "LAYN",
    "LAYN": "LAYN",
    "PD1": "PD-1",
    "PDCD1": "PD-1",
    "PD-1": "PD-1",
    "PDL1": "PD-L1",
    "PD-L1": "PD-L1",
    "CD274": "PD-L1",
    "CTLA4": "CTLA-4",
    "CTLA-4": "CTLA-4",
    "TIM3": "TIM-3",
    "TIM-3": "TIM-3",
    "HAVCR2": "TIM-3",
    "LAG3": "LAG-3",
    "LAG-3": "LAG-3",
    "IFNG": "IFNG",
    "IFN-G": "IFNG",
    "IFN-GAMMA": "IFNG",
    "IFN-Γ": "IFNG",
    "NKG2A": "NKG2A",
    "NKG2C": "NKG2C",
    "NKG2D": "NKG2D",
}


def sanitize_entity_symbol(symbol: str) -> Optional[str]:
    s = (symbol or "").strip()
    if not s:
        return None
    up = s.upper().replace("Γ", "G").replace("Β", "B")
    if s.lower() == "layilin":
        return "LAYN"
    if up in _ENTITY_ALIAS:
        return _ENTITY_ALIAS[up]
    compact = up.replace("-", "")
    for k, v in _ENTITY_ALIAS.items():
        if k.replace("-", "") == compact:
            return v
    if up in _ENTITY_FRAGMENT_BLOCK or len(s) < 3:
        return None
    if re.fullmatch(r"GSE\d+|GSM\d+", up):
        return None
    return s


def llm_payload_to_extract(payload: dict[str, Any], *, backend: str = "agent") -> dict[str, Any]:
    """Convert validated payload → wiki_lib extract shape."""
    concepts: list[str] = []
    concept_details: list[dict[str, str]] = []
    for c in payload.get("concepts") or []:
        if isinstance(c, dict):
            name = c.get("name") or ""
            evidence = str(c.get("evidence") or c.get("note") or "").strip()
        else:
            name = str(c)
            evidence = ""
        if name and name not in concepts:
            concepts.append(name)
            concept_details.append({"name": name, "evidence": evidence})

    entities: list[str] = []
    for e in payload.get("entities") or []:
        if isinstance(e, dict):
            sym = sanitize_entity_symbol(str(e.get("symbol") or e.get("name") or ""))
        else:
            sym = sanitize_entity_symbol(str(e))
        if sym and sym not in entities:
            entities.append(sym)

    tags: list[str] = []
    for t in list(payload.get("tags") or []) + concepts[:4] + entities[:4]:
        if t and t not in tags:
            tags.append(t)
    tags = tags[:12]

    return {
        "summary": payload["summary"],
        "claims": list(payload.get("key_claims") or []),
        "methods": list(payload.get("methods") or []),
        "diseases": [],
        "topics": [],
        "concepts": concepts,
        "concept_details": concept_details,
        "entities": entities,
        "tags": tags,
        "confidence_notes": payload.get("confidence_notes") or "",
        "compile_backend": backend,
        "schema": COMPILE_SCHEMA_VERSION,
    }


def build_messages(pub: dict[str, Any]) -> list[dict[str, str]]:
    return [
        {"role": "system", "content": SYSTEM_PROMPT},
        {
            "role": "user",
            "content": USER_PROMPT_TEMPLATE.format(
                publication_id=pub.get("id") or "",
                title=pub.get("title") or "",
                year=pub.get("year") or "",
                venue=pub.get("venue") or "",
                doi=pub.get("doi") or "",
                pmid=pub.get("pmid") or "",
                abstract=pub.get("abstract") or "",
            ),
        },
    ]


def build_compile_request(
    pub: dict[str, Any],
    *,
    response_path: Optional[str] = None,
    request_path: Optional[str] = None,
) -> dict[str, Any]:
    """Build a grounded compile request for the invoking agent to fulfill."""
    messages = build_messages(pub)
    pid = str(pub.get("id") or "")
    abstract = (pub.get("abstract") or "").strip()
    title = (pub.get("title") or "").strip()
    warnings: list[str] = []
    if not abstract:
        warnings.append(
            "abstract_missing: prefer lit-enrich (PMID/DOI fetch) before agent compile; "
            "placeholder titles/empty abstracts yield weak wiki pages."
        )
    # shell-like title
    low = title.lower()
    if not title or any(
        low == str(pub.get(k) or "").lower() or low == f"{k}:{pub.get(k)}".lower()
        for k in ("doi", "pmid", "pmcid", "arxiv", "openalex")
        if pub.get(k)
    ) or low.startswith(("doi:", "pmid:", "pmcid:")):
        warnings.append(
            "shell_or_empty_title: run lit-enrich before compile when only IDs are known."
        )
    out: dict[str, Any] = {
        "ok": True,
        "action": "compile_request",
        "schema_version": COMPILE_SCHEMA_VERSION,
        "publication_id": pid,
        "publication": {
            "id": pid,
            "title": pub.get("title") or "",
            "year": pub.get("year"),
            "venue": pub.get("venue") or "",
            "doi": pub.get("doi") or "",
            "pmid": pub.get("pmid") or "",
            "abstract": pub.get("abstract") or "",
        },
        "system_prompt": SYSTEM_PROMPT,
        "user_prompt": messages[1]["content"],
        "messages": messages,
        "response_schema": RESPONSE_JSON_SCHEMA,
        "instructions": (
            "Use your already-running model to produce one JSON object matching "
            "response_schema from the grounded title/abstract/metadata only. "
            "Write it to response_path (or pass via --llm-json / --apply-compile), "
            "then re-invoke lit-wiki-ingest. Do not set plugin API keys."
            + (
                " If warnings include abstract_missing, run enrich.py first."
                if warnings
                else ""
            )
        ),
        "request_path": request_path,
        "response_path": response_path,
        "next_step": (
            f"Fill {response_path or 'compile_response.json'} then run: "
            "wiki_ingest.py ... --compile agent --llm-json <response.json>"
        ),
    }
    if warnings:
        out["warnings"] = warnings
        out["prefer_enrich_first"] = True
    return out


def write_compile_request(path: Path, request: dict[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(request, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return path


def _index_file_payloads(path: Path) -> dict[str, dict[str, Any]]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    out: dict[str, dict[str, Any]] = {}
    if isinstance(raw, dict):
        if "summary" in raw and "key_claims" in raw:
            pid = str(raw.get("publication_id") or "").strip()
            out[pid or "__default__"] = raw
        else:
            for k, v in raw.items():
                if isinstance(v, dict):
                    out[str(k)] = v
    elif isinstance(raw, list):
        for item in raw:
            if not isinstance(item, dict):
                continue
            pid = str(item.get("publication_id") or item.get("id") or "").strip()
            if not pid:
                raise ValueError("file compile list items need publication_id")
            out[pid] = item
    else:
        raise ValueError("llm-json must be object or list")
    return out


def file_compile(pub: dict[str, Any], llm_json_path: str | Path, *, backend: str = "agent") -> dict[str, Any]:
    path = Path(llm_json_path)
    if not path.is_file():
        raise FileNotFoundError(f"llm-json not found: {path}")
    index = _index_file_payloads(path)
    pid = str(pub.get("id") or "")
    candidates = [pid, pid.replace("doi:", ""), str(pub.get("doi") or ""), "__default__"]
    payload = None
    for c in candidates:
        if c and c in index:
            payload = index[c]
            break
    if payload is None and len(index) == 1:
        payload = next(iter(index.values()))
    if payload is None:
        raise KeyError(f"no compile JSON for publication_id={pid} in {path}")
    validated = validate_compile_payload(payload)
    out = llm_payload_to_extract(validated, backend=backend)
    return out


def mock_compile(pub: dict[str, Any], fixture_path: Optional[str | Path] = None) -> dict[str, Any]:
    """Deterministic mock compile for smoke/unit tests — not a live model call."""
    if fixture_path and Path(fixture_path).is_file():
        return file_compile(pub, fixture_path, backend="mock")
    title = str(pub.get("title") or pub.get("id") or "Untitled")
    abstract = str(pub.get("abstract") or "")
    focus = "tumor-infiltrating T cells"
    if re.search(r"exhaust", abstract + title, re.I):
        focus = "T cell exhaustion in the tumor microenvironment"
    if re.search(r"hepatocellular|HCC|liver cancer", abstract + title, re.I):
        disease = "hepatocellular carcinoma (HCC)"
    else:
        disease = "the studied disease context"
    summary = (
        f"This study examines {focus} in {disease}, based on the provided title and abstract. "
        f"The authors report molecular and cellular observations that clarify immune states "
        f"relevant to immunotherapy research (mock compile; grounded in supplied text only)."
    )
    claims: list[str] = []
    if re.search(r"single[\s-]?cell|scrna", abstract + title, re.I):
        claims.append("Single-cell profiling is used to resolve T cell heterogeneity in the cohort described.")
    if re.search(r"exhaust", abstract, re.I):
        claims.append("Exhausted CD8+ T cell states are highlighted as enriched or functionally relevant.")
    if re.search(r"layilin|LAYN", abstract, re.I):
        claims.append("LAYN (layilin) is linked to immunosuppressive or activated T cell populations.")
    if re.search(r"Treg|regulatory T", abstract, re.I):
        claims.append("Regulatory T cells are reported as preferentially present in tumor tissue.")
    if not claims:
        claims.append(f"Primary findings concern immune or molecular features of {disease} as stated in the abstract.")
    methods: list[str] = []
    for name, pat in [
        ("scRNA-seq", r"single[\s-]?cell|scrna"),
        ("RNA-seq", r"rna[\s-]?seq"),
        ("flow cytometry", r"flow\s+cytometr"),
        ("TCR sequencing", r"\bTCR\b"),
    ]:
        if re.search(pat, abstract + title, re.I):
            methods.append(name)
    concepts = []
    for name, pat in [
        ("T cell exhaustion", r"exhaust"),
        ("hepatocellular carcinoma", r"hepatocellular|\bHCC\b"),
        ("immunotherapy", r"immunotherap|checkpoint"),
        ("infiltrating T cells", r"infiltrat"),
        ("scRNA-seq", r"single[\s-]?cell|scrna"),
    ]:
        if re.search(pat, abstract + title, re.I):
            concepts.append({"name": name, "slug_hint": ""})
    entities = []
    if re.search(r"layilin|\bLAYN\b", abstract, re.I):
        entities.append({"symbol": "LAYN", "name": "layilin", "kind": "gene"})
    if re.search(r"\bCD8\b", abstract + title, re.I):
        entities.append({"symbol": "CD8", "name": "CD8", "kind": "protein"})
    if re.search(r"\bTCR\b", abstract, re.I):
        entities.append({"symbol": "TCR", "name": "T cell receptor", "kind": "complex"})
    if re.search(r"PD-L1|PD-1", abstract, re.I):
        if re.search(r"PD-L1", abstract):
            entities.append({"symbol": "PD-L1", "kind": "protein"})
        if re.search(r"PD-1", abstract):
            entities.append({"symbol": "PD-1", "kind": "protein"})
    if re.search(r"IFN-γ|IFNG|interferon", abstract, re.I):
        entities.append({"symbol": "IFNG", "name": "interferon gamma", "kind": "gene"})
    tags = [c["name"] for c in concepts[:4]] + [e["symbol"] for e in entities[:3]]
    payload = validate_compile_payload(
        {
            "summary": summary,
            "key_claims": claims[:6],
            "methods": methods,
            "concepts": concepts,
            "entities": entities,
            "tags": tags,
            "confidence_notes": "mock provider — synthesised for tests; not a live model call",
            "publication_id": str(pub.get("id") or ""),
        }
    )
    return llm_payload_to_extract(payload, backend="mock")


def compile_with_provider(
    pub: dict[str, Any],
    *,
    mode: str,
    llm_json: Optional[str | Path] = None,
    allow_heuristic_fallback: bool = False,
    heuristic_fn=None,
    emit_request_path: Optional[str | Path] = None,
) -> dict[str, Any]:
    """Dispatch compile. Agent without JSON raises CompileRequestNeeded (or writes request)."""
    mode = (mode or "agent").lower()
    if mode == "file":
        mode = "agent"
    try:
        if mode == "heuristic":
            if heuristic_fn is None:
                raise RuntimeError("heuristic_fn required")
            ex = heuristic_fn(str(pub.get("title") or ""), str(pub.get("abstract") or ""))
            ex = dict(ex)
            ex["compile_backend"] = "heuristic"
            return ex
        if mode == "mock":
            return mock_compile(pub, fixture_path=llm_json)
        if mode == "agent":
            if llm_json:
                return file_compile(pub, llm_json, backend="agent")
            req = build_compile_request(
                pub,
                response_path=str(
                    Path(emit_request_path).with_name("compile_response.json")
                    if emit_request_path
                    else "compile_response.json"
                ),
                request_path=str(emit_request_path) if emit_request_path else None,
            )
            if emit_request_path:
                write_compile_request(Path(emit_request_path), req)
                req["request_path"] = str(Path(emit_request_path))
            raise CompileRequestNeeded(req)
        raise ValueError(f"unknown compile mode: {mode}")
    except CompileRequestNeeded:
        raise
    except Exception as e:
        if allow_heuristic_fallback and heuristic_fn is not None and mode != "heuristic":
            ex = heuristic_fn(str(pub.get("title") or ""), str(pub.get("abstract") or ""))
            ex = dict(ex)
            ex["compile_backend"] = "heuristic_fallback"
            ex["compile_error"] = str(e)
            return ex
        raise


# ---------------------------------------------------------------------------
# Cross-paper concept deep merge (v0.4.5) — agent request → JSON → apply
# ---------------------------------------------------------------------------

CONCEPT_MERGE_SCHEMA_VERSION = "wiki-concept-merge-v1"
CONCEPT_MERGE_LABEL = "concept-merge-v0.4.5"

CONCEPT_MERGE_SYSTEM_PROMPT = """You reconcile one scientific concept page across multiple papers.
Rules:
- Use ONLY the provided Definition, per-source contributions, and Sources. Do not invent papers, IDs, statistics, or mechanisms absent from those texts.
- definition: 2–5 sentences synthesizing the concept as used across the listed sources.
- synthesis: bullet claims; each bullet must cite at least one source_id from the contributions.
- related_concepts / aliases: only if evidenced; otherwise empty arrays.
- If any contribution has needs_review=true, set needs_review_banner to a short warning string (never omit).
- Optional canonical_slug / merge_into: only when this slug should redirect into another existing concept; otherwise omit or null.
- Return a single JSON object matching the schema. No markdown fences.
"""

CONCEPT_MERGE_USER_TEMPLATE = """Merge concept `{slug}` (`{concept_name}`) across papers.

Current Definition:
{current_definition}

Current Sources (markdown):
{current_sources}

Per-source contributions (JSON):
{contributions_json}

Return JSON with keys:
definition (string),
synthesis ({{claim, source_ids[]}}[]),
related_concepts (string[]),
aliases (string[], optional),
needs_review_banner (string|null),
canonical_slug (string|null, optional),
merge_into (string|null, optional),
confidence_notes (string, optional).
"""

CONCEPT_MERGE_RESPONSE_SCHEMA: dict[str, Any] = {
    "type": "object",
    "required": ["definition", "synthesis"],
    "properties": {
        "definition": {"type": "string", "minLength": 40},
        "synthesis": {
            "type": "array",
            "minItems": 1,
            "items": {
                "type": "object",
                "required": ["claim", "source_ids"],
                "properties": {
                    "claim": {"type": "string"},
                    "source_ids": {"type": "array", "items": {"type": "string"}, "minItems": 1},
                },
            },
        },
        "related_concepts": {"type": "array", "items": {"type": "string"}},
        "aliases": {"type": "array", "items": {"type": "string"}},
        "needs_review_banner": {"type": ["string", "null"]},
        "canonical_slug": {"type": ["string", "null"]},
        "merge_into": {"type": ["string", "null"]},
        "confidence_notes": {"type": "string"},
        "slug": {"type": "string"},
    },
}


class ConceptMergeRequestNeeded(Exception):
    def __init__(self, request: dict[str, Any]):
        self.request = request
        super().__init__("agent concept-merge response required")


def _read_section(md: str, heading: str) -> str:
    marker = f"## {heading}"
    if marker not in md:
        return ""
    _, _, rest = md.partition(marker)
    # drop leading newlines after heading line
    rest = rest.lstrip("\n")
    if rest.startswith("\n"):
        rest = rest[1:]
    # heading line may still be empty first line after ## Name\n
    # Actually partition keeps nothing of marker; rest starts after "## Definition"
    # which may be "\n\ntext..."
    body, _, _ = rest.partition("\n## ")
    return body.strip()


def _parse_sources_block(md: str) -> str:
    return _read_section(md, "Sources")


def validate_concept_merge_payload(
    obj: Any,
    *,
    allowed_source_ids: Optional[set[str]] = None,
    require_review_banner: bool = False,
) -> dict[str, Any]:
    if not isinstance(obj, dict):
        raise ValueError("concept merge payload must be an object")
    definition = _require_str(obj.get("definition"), "definition")
    if len(definition) < 40:
        raise ValueError("definition too short")
    synth_in = obj.get("synthesis")
    if not isinstance(synth_in, list) or not synth_in:
        raise ValueError("synthesis must be a non-empty list")
    synthesis: list[dict[str, Any]] = []
    for i, item in enumerate(synth_in):
        if isinstance(item, str):
            claim = item.strip()
            if not claim:
                raise ValueError(f"synthesis[{i}] empty")
            synthesis.append({"claim": claim, "source_ids": []})
            continue
        if not isinstance(item, dict):
            raise ValueError(f"synthesis[{i}] must be object or string")
        claim = _require_str(item.get("claim"), f"synthesis[{i}].claim")
        sids = item.get("source_ids") or item.get("sources") or []
        if not isinstance(sids, list) or not sids:
            raise ValueError(f"synthesis[{i}].source_ids required")
        clean_ids: list[str] = []
        for j, sid in enumerate(sids):
            if not isinstance(sid, str) or not sid.strip():
                raise ValueError(f"synthesis[{i}].source_ids[{j}] invalid")
            sid_s = sid.strip()
            if allowed_source_ids is not None and sid_s not in allowed_source_ids:
                # allow safe_id or full id — check both already in set
                raise ValueError(f"synthesis[{i}] cites unknown source_id {sid_s!r}")
            clean_ids.append(sid_s)
        synthesis.append({"claim": claim, "source_ids": clean_ids})

    related = _require_str_list(obj.get("related_concepts") or [], "related_concepts")
    aliases = _require_str_list(obj.get("aliases") or [], "aliases")
    banner = obj.get("needs_review_banner")
    if banner is not None and not isinstance(banner, str):
        raise ValueError("needs_review_banner must be string or null")
    banner_s = (banner or "").strip()
    if require_review_banner and not banner_s:
        raise ValueError("needs_review_banner required when contributions include needs_review")
    notes = obj.get("confidence_notes")
    if notes is not None and not isinstance(notes, str):
        raise ValueError("confidence_notes must be string if present")
    canonical = obj.get("canonical_slug") or obj.get("merge_into")
    if canonical is not None and not isinstance(canonical, str):
        raise ValueError("canonical_slug/merge_into must be string or null")
    canonical_s = (canonical or "").strip() or None
    return {
        "definition": definition,
        "synthesis": synthesis,
        "related_concepts": related,
        "aliases": aliases,
        "needs_review_banner": banner_s or None,
        "canonical_slug": canonical_s,
        "merge_into": (str(obj.get("merge_into") or "").strip() or None),
        "confidence_notes": (notes or "").strip(),
        "schema": CONCEPT_MERGE_SCHEMA_VERSION,
        "slug": str(obj.get("slug") or "").strip(),
    }


def build_concept_merge_request(
    vault: Path,
    slug: str,
    *,
    request_path: Optional[str] = None,
    response_path: Optional[str] = None,
) -> dict[str, Any]:
    from wiki_lib import list_concept_contributions  # local import to avoid cycles at load

    slug = slug.strip()
    if not slug:
        raise ValueError("slug required")
    page = vault / "wiki" / "concepts" / f"{slug}.md"
    contribs = list_concept_contributions(vault, slug)
    n_page_sources = 0
    if page.exists():
        n_page_sources = _parse_sources_block(page.read_text(encoding="utf-8")).count("- [")
    if len(contribs) < 2 and n_page_sources < 2:
        raise ValueError(
            f"concept {slug!r} needs ≥2 contributions (or ≥2 Sources links); "
            f"have contribs={len(contribs)} page_sources={n_page_sources}"
        )

    md = page.read_text(encoding="utf-8") if page.exists() else ""
    # title from first heading
    concept_name = slug
    for line in md.splitlines():
        if line.startswith("# "):
            concept_name = line[2:].strip() or slug
            break
    if contribs and contribs[0].get("concept_name"):
        concept_name = str(contribs[0]["concept_name"])

    # strip internal _path for prompt
    prompt_contribs = []
    allowed: set[str] = set()
    any_review = False
    for c in contribs:
        row = {k: v for k, v in c.items() if not str(k).startswith("_")}
        prompt_contribs.append(row)
        for key in ("source_id", "source_safe_id"):
            if row.get(key):
                allowed.add(str(row[key]))
        if row.get("needs_review"):
            any_review = True

    user = CONCEPT_MERGE_USER_TEMPLATE.format(
        slug=slug,
        concept_name=concept_name,
        current_definition=_read_section(md, "Definition") or "_No definition yet._",
        current_sources=_parse_sources_block(md) or "_No sources yet._",
        contributions_json=json.dumps(prompt_contribs, ensure_ascii=False, indent=2),
    )
    messages = [
        {"role": "system", "content": CONCEPT_MERGE_SYSTEM_PROMPT},
        {"role": "user", "content": user},
    ]
    return {
        "ok": True,
        "action": "concept_merge_request",
        "schema_version": CONCEPT_MERGE_SCHEMA_VERSION,
        "slug": slug,
        "concept_name": concept_name,
        "n_contributions": len(contribs),
        "allowed_source_ids": sorted(allowed),
        "requires_needs_review_banner": any_review,
        "system_prompt": CONCEPT_MERGE_SYSTEM_PROMPT,
        "user_prompt": user,
        "messages": messages,
        "response_schema": CONCEPT_MERGE_RESPONSE_SCHEMA,
        "contributions": prompt_contribs,
        "current_definition": _read_section(md, "Definition"),
        "current_sources_md": _parse_sources_block(md),
        "instructions": (
            "Use your already-running model to produce one JSON object matching "
            "response_schema from the grounded contributions only. "
            "Write it to response_path then run concept_merge_apply. No plugin API keys."
        ),
        "request_path": request_path,
        "response_path": response_path,
        "next_step": (
            f"Fill {response_path or 'concept_merge_response.json'} then run "
            "concept_merge_apply.py --llm-json <response.json>"
        ),
    }


def mock_concept_merge_payload(
    vault: Path,
    slug: str,
) -> dict[str, Any]:
    """Deterministic offline merge from contribution blurbs — tests only."""
    from wiki_lib import list_concept_contributions

    contribs = list_concept_contributions(vault, slug)
    if not contribs:
        raise ValueError(f"no contributions for {slug}")
    name = str(contribs[0].get("concept_name") or slug)
    titles = [str(c.get("source_title") or c.get("source_id") or "") for c in contribs]
    ids = [str(c.get("source_id") or c.get("source_safe_id") or "") for c in contribs]
    snippets = []
    for c in contribs:
        ev = (c.get("evidence") or "").strip()
        if not ev and c.get("claims"):
            ev = str(c["claims"][0])
        if not ev:
            ev = (c.get("summary_snippet") or "").strip()[:180]
        snippets.append(ev)
    def_parts = [
        f"{name} is discussed across {len(contribs)} compiled sources in this vault "
        f"(mock merge; grounded in contribution snippets only)."
    ]
    for t, snip in zip(titles, snippets):
        if snip:
            def_parts.append(f"From «{t}»: {snip[:160]}")
    definition = " ".join(def_parts)
    if len(definition) < 40:
        definition += " Additional context is limited to the stored contribution records."
    synthesis = []
    for c, sid in zip(contribs, ids):
        claim = ""
        if c.get("claims"):
            claim = str(c["claims"][0])
        elif c.get("evidence"):
            claim = str(c["evidence"])
        else:
            claim = f"Contribution from {c.get('source_title') or sid} mentions {name}."
        synthesis.append({"claim": claim[:400], "source_ids": [sid]})
    any_review = any(bool(c.get("needs_review")) for c in contribs)
    banner = (
        "One or more contributing sources are marked needs_review — treat synthesis cautiously."
        if any_review
        else None
    )
    return validate_concept_merge_payload(
        {
            "slug": slug,
            "definition": definition,
            "synthesis": synthesis,
            "related_concepts": [],
            "aliases": [],
            "needs_review_banner": banner,
            "confidence_notes": "mock concept merge — deterministic; not a live model call",
        },
        allowed_source_ids=set(ids) | {str(c.get("source_safe_id") or "") for c in contribs},
        require_review_banner=any_review,
    )


def render_concept_page(
    *,
    title: str,
    slug: str,
    definition: str,
    synthesis: list[dict[str, Any]],
    sources_md: str,
    related: list[Any],
    aliases: Optional[list[str]] = None,
    needs_review_banner: Optional[str] = None,
    merge_label: str = CONCEPT_MERGE_LABEL,
) -> str:
    banner = ""
    if needs_review_banner:
        banner = f"\n> **needs_review**: {needs_review_banner}\n"
    synth_lines = []
    for item in synthesis:
        claim = item.get("claim") or ""
        sids = item.get("source_ids") or []
        cite = ", ".join(f"`{s}`" for s in sids) if sids else "_uncited_"
        synth_lines.append(f"- {claim} ({cite})")
    synth_block = "\n".join(synth_lines) if synth_lines else "_No synthesis bullets._"
    sources_block = sources_md.strip() if sources_md.strip() else "_No sources._"
    if sources_block and not sources_block.endswith("\n"):
        sources_block += "\n"
    rel_lines = []
    for r in related or []:
        if isinstance(r, dict):
            display = str(r.get("display") or r.get("name") or r.get("slug") or "").strip()
            slug = str(r.get("slug") or "").strip()
            if not slug:
                continue
            if not display:
                display = slug
            rel_lines.append(f"- [`{display}`]({slug}.md)")
        elif isinstance(r, str) and r.strip():
            # legacy: caller already resolved or raw — prefer slugify only as last resort
            s = _slugify_local(r)
            rel_lines.append(f"- [`{r.strip()}`]({s}.md)")
    rel_block = "\n".join(rel_lines) if rel_lines else ""
    alias_block = ""
    if aliases:
        alias_block = "\n## Aliases\n\n" + ", ".join(f"`{a}`" for a in aliases) + "\n"
    return (
        f"# {title}\n"
        f"{banner}\n"
        f"Slug: `{slug}`\n"
        f"- merge: `{merge_label}`\n\n"
        f"## Definition\n\n"
        f"{definition}\n\n"
        f"## Across papers\n\n"
        f"{synth_block}\n\n"
        f"## Sources\n\n"
        f"{sources_block}\n"
        f"## Related concepts\n\n"
        f"{rel_block}\n"
        f"{alias_block}"
    )


def _slugify_local(name: str) -> str:
    s = re.sub(r"[^a-zA-Z0-9]+", "-", (name or "").strip().lower()).strip("-")
    return s or "item"


def _collect_sources_md(vault: Path, slug: str, page_md: str, contribs: list[dict[str, Any]]) -> str:
    """Preserve existing ## Sources bullets; append any contrib links missing."""
    existing = _parse_sources_block(page_md)
    lines = [ln for ln in existing.splitlines() if ln.strip().startswith("- [")]
    seen = set(lines)
    for c in contribs:
        title = str(c.get("source_title") or c.get("source_id") or "source")
        rel = str(c.get("source_rel") or "")
        if not rel:
            continue
        bullet = f"- [{title}]({rel})"
        if not any(rel in ln for ln in lines):
            lines.append(bullet)
            seen.add(bullet)
    return "\n".join(lines) + ("\n" if lines else "")


def apply_concept_merge(
    vault: Path,
    slug: str,
    payload: dict[str, Any],
    *,
    backend: str = "agent",
) -> dict[str, Any]:
    """Validate + rewrite wiki/concepts/<slug>.md; optional redirect on merge_into."""
    from wiki_lib import append_log, list_concept_contributions, slugify

    contribs = list_concept_contributions(vault, slug)
    allowed: set[str] = set()
    any_review = False
    for c in contribs:
        for key in ("source_id", "source_safe_id"):
            if c.get(key):
                allowed.add(str(c[key]))
        if c.get("needs_review"):
            any_review = True
    validated = validate_concept_merge_payload(
        payload,
        allowed_source_ids=allowed or None,
        require_review_banner=any_review,
    )
    page = vault / "wiki" / "concepts" / f"{slug}.md"
    page_md = page.read_text(encoding="utf-8") if page.exists() else ""
    title = slug
    for line in page_md.splitlines():
        if line.startswith("# "):
            title = line[2:].strip() or slug
            break
    if contribs and contribs[0].get("concept_name"):
        title = str(contribs[0]["concept_name"])

    orig_slug = slug
    canonical = validated.get("canonical_slug") or validated.get("merge_into")
    redirect_written = None
    target_slug = slug
    if canonical and slugify(canonical) != slug and canonical != slug:
        target_slug = slugify(canonical) if "/" not in canonical else canonical
        # If merge_into looks like a slug already
        if re.fullmatch(r"[a-z0-9\-]+", canonical):
            target_slug = canonical
        else:
            target_slug = slugify(canonical)
        # Write redirect on old page
        redirect = (
            f"# {title}\n\n"
            f"Slug: `{slug}`\n\n"
            f"> Merged into [`{target_slug}`]({target_slug}.md).\n\n"
            f"See [{target_slug}]({target_slug}.md).\n"
        )
        page.parent.mkdir(parents=True, exist_ok=True)
        page.write_text(redirect, encoding="utf-8")
        redirect_written = str(page.relative_to(vault))
        # Apply body onto canonical page
        page = vault / "wiki" / "concepts" / f"{target_slug}.md"
        page_md = page.read_text(encoding="utf-8") if page.exists() else page_md
        slug = target_slug

    sources_md = _collect_sources_md(vault, slug, page_md, contribs)
    from wiki_lib import resolve_related_concepts
    related_resolved, related_unresolved = resolve_related_concepts(
        vault, list(validated.get("related_concepts") or [])
    )
    body = render_concept_page(
        title=title,
        slug=slug,
        definition=validated["definition"],
        synthesis=validated["synthesis"],
        sources_md=sources_md,
        related=related_resolved,
        aliases=validated.get("aliases") or [],
        needs_review_banner=validated.get("needs_review_banner"),
        merge_label=f"{backend}-{CONCEPT_MERGE_LABEL}",
    )
    page.parent.mkdir(parents=True, exist_ok=True)
    page.write_text(body if body.endswith("\n") else body + "\n", encoding="utf-8")
    append_log(
        vault,
        f"concept_merge `{slug}` backend={backend} synth={len(validated['synthesis'])} "
        f"sources={sources_md.count('- [')} redirect={redirect_written or ''}",
    )
    # Clear soft pending stub if present
    stub = vault / "raw" / f"concept_merge_request_{slug}.json"
    # also clear original slug stub if redirected
    for s in {orig_slug, slug, target_slug}:
        sp = vault / "raw" / f"concept_merge_request_{s}.json"
        if sp.exists():
            try:
                sp.unlink()
            except OSError:
                pass
    return {
        "ok": True,
        "action": "concept_merge_apply",
        "slug": slug,
        "concept_page": str(page.relative_to(vault)),
        "redirect_from": redirect_written,
        "n_synthesis": len(validated["synthesis"]),
        "backend": backend,
        "merge_label": f"{backend}-{CONCEPT_MERGE_LABEL}",
        "needs_review_banner": validated.get("needs_review_banner"),
        "related_resolved": related_resolved,
        "related_unresolved": related_unresolved,
    }
