#!/usr/bin/env python3
"""Ingest publications from hits.jsonl / BibTeX / RIS / ID list / PDF text sidecar / local PDF.

Local PDFs: pdftotext (poppler) for digital text; optional OCR soft-fail
(tesseract or LIT_PDF_OCR_CMD / --ocr_cmd). Never stores PDF bytes in SQLite.
TextIn is not a required parser and is not the KB backbone.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any, Optional

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE / "lib"))
sys.path.insert(0, str(_HERE))

from db import (  # noqa: E402
    apply_id_merges,
    link_to_collection,
    open_db,
    record_pub_file,
    record_run_hit,
    resolve_or_create,
    set_publication_oa,
    utc_now,
)
from formats import latex_to_unicode  # noqa: E402
from oa_verify import verify_record_oa  # noqa: E402


_METADATA_ONLY_SUFFIXES = {".bib", ".ris", ".jsonl", ".ids", ".csv"}


def _parse_bibtex(text: str) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    # Split on @type{key,
    parts = re.split(r"(?=@\w+\s*\{)", text)
    for part in parts:
        part = part.strip()
        if not part.startswith("@"):
            continue
        m = re.match(r"@(\w+)\s*\{\s*([^,]+)\s*,", part)
        if not m:
            continue
        etype, _key = m.group(1), m.group(2)
        fields: dict[str, str] = {}
        for fm in re.finditer(
            r"(\w+)\s*=\s*(\{(?:[^{}]|\{[^{}]*\})*\}|\"[^\"]*\"|[^\s,]+)\s*,?",
            part,
            flags=re.S,
        ):
            k = fm.group(1).lower()
            if k in ("@type",):
                continue
            raw = fm.group(2).strip()
            if raw.startswith("{") and raw.endswith("}"):
                raw = raw[1:-1]
            elif raw.startswith('"') and raw.endswith('"'):
                raw = raw[1:-1]
            fields[k] = raw.strip()
        authors = []
        if "author" in fields:
            authors = [
                latex_to_unicode(a) for a in re.split(r"\s+and\s+", fields["author"]) if a.strip()
            ]
        year = None
        if "year" in fields:
            ym = re.search(r"\d{4}", fields["year"])
            year = ym.group(0) if ym else fields["year"]
        entries.append(
            {
                "source": "bibtex",
                "type": etype.lower(),
                "title": latex_to_unicode(fields.get("title", "").replace("\n", " ")),
                "abstract": fields.get("abstract"),
                "authors": authors,
                "year": year,
                "venue": fields.get("journal") or fields.get("booktitle") or fields.get("publisher"),
                "doi": fields.get("doi"),
                "pmid": fields.get("pmid"),
                "pmcid": fields.get("pmcid"),
                "arxiv": fields.get("arxiv") or fields.get("eprint"),
                "url": fields.get("url"),
            }
        )
    return entries


def _parse_ris(text: str) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    current: dict[str, Any] = {"source": "ris", "authors": []}
    for line in text.splitlines():
        if not line.strip():
            continue
        if line.startswith("ER  -") or line.startswith("ER  - "):
            if current.get("title") or current.get("doi") or current.get("pmid"):
                entries.append(current)
            current = {"source": "ris", "authors": []}
            continue
        m = re.match(r"^([A-Z0-9]{2})\s+-\s+(.*)$", line)
        if not m:
            continue
        tag, val = m.group(1), m.group(2).strip()
        if tag == "TI" or tag == "T1":
            current["title"] = val
        elif tag == "AU" or tag == "A1":
            current.setdefault("authors", []).append(val)
        elif tag == "PY" or tag == "Y1":
            current["year"] = val[:4]
        elif tag == "JO" or tag == "T2" or tag == "JF":
            current["venue"] = val
        elif tag == "DO" or tag == "DOI":
            current["doi"] = val
        elif tag == "AN" and val.isdigit():
            current["pmid"] = val
        elif tag in ("UR", "L1", "L2", "LK"):
            current.setdefault("url", val)
            for k, v in extract_ids_from_text(val).items():
                current.setdefault(k, v)
        elif tag == "AB" or tag == "N2":
            current["abstract"] = val
        elif tag == "TY":
            current["type"] = val.lower()
        elif tag == "ID" and val.upper().startswith("PMC"):
            current["pmcid"] = val
    if current.get("title") or current.get("doi"):
        entries.append(current)
    return entries


_STRONG_ID_KINDS = ("doi", "pmid", "pmcid", "arxiv", "openalex")


def _shell_title_for(draft: dict[str, Any]) -> str:
    return (
        draft.get("doi")
        or draft.get("pmid")
        or draft.get("pmcid")
        or draft.get("arxiv")
        or draft.get("openalex")
        or draft.get("title")
        or "unknown"
    )


def _parse_id_tokens_from_line(line: str) -> dict[str, Any]:
    """Extract strong IDs (and optional title) from one ids-file line.

    Supports:
      doi:10.1016/…
      pmid:28622514
      doi:10.1016/… pmid:28622514
      10.1016/…\\t28622514
      10.1016/…, PMID:28622514
    """
    draft: dict[str, Any] = {"source": "id_list"}
    raw = line.strip()
    if not raw:
        return draft

    labeled = list(
        re.finditer(
            r"(?i)\b(doi|pmid|pmcid|pmc|arxiv|openalex)\s*:\s*(\S+)",
            raw,
        )
    )
    if labeled:
        for m in labeled:
            kind = m.group(1).lower()
            val = m.group(2).rstrip(",;")
            if kind in ("pmc", "pmcid"):
                draft["pmcid"] = val
            elif kind == "doi":
                draft["doi"] = val
            elif kind == "pmid":
                draft["pmid"] = val
            elif kind == "arxiv":
                draft["arxiv"] = val
            elif kind == "openalex":
                draft["openalex"] = val
        if "doi" not in draft:
            m = re.search(r"(?i)doi\.org/(\S+)", raw)
            if m:
                draft["doi"] = m.group(1).rstrip(",;")
        return draft

    parts = [p.strip() for p in re.split(r"[\t,;]+", raw) if p.strip()]
    if len(parts) >= 2:
        for part in parts:
            low = part.lower()
            if re.match(r"10\.\d+/", part) or (low.startswith("10.") and "/" in part):
                draft.setdefault("doi", part)
            elif re.fullmatch(r"\d{5,}", part):
                draft.setdefault("pmid", part)
            elif low.startswith("pmc"):
                draft.setdefault("pmcid", part)
        if any(draft.get(k) for k in _STRONG_ID_KINDS):
            return draft

    low = raw.lower()
    if "doi.org/" in low:
        draft["doi"] = raw.split("doi.org/", 1)[1].strip().rstrip(",;")
    elif re.match(r"10\.\d+/", raw):
        draft["doi"] = raw
    elif re.fullmatch(r"\d{5,}", raw):
        draft["pmid"] = raw
    elif low.startswith("pmc") and re.search(r"\d", raw):
        draft["pmcid"] = raw
    else:
        draft["title"] = raw
    return draft


def _n_strong_ids(d: dict[str, Any]) -> int:
    return sum(1 for k in _STRONG_ID_KINDS if d.get(k))


def _strong_kinds(d: dict[str, Any]) -> set[str]:
    return {k for k in _STRONG_ID_KINDS if d.get(k)}


def _ids_conflict(a: dict[str, Any], b: dict[str, Any]) -> bool:
    for k in _STRONG_ID_KINDS:
        va, vb = a.get(k), b.get(k)
        if va and vb and str(va).lower() != str(vb).lower():
            return True
    return False


def _merge_id_drafts(base: dict[str, Any], extra: dict[str, Any]) -> dict[str, Any]:
    out = dict(base)
    for k in _STRONG_ID_KINDS:
        if extra.get(k) and not out.get(k):
            out[k] = extra[k]
    if extra.get("title") and (
        not out.get("title")
        or out.get("title") in {out.get(x) for x in _STRONG_ID_KINDS}
    ):
        out["title"] = extra["title"]
    if _n_strong_ids(out) >= 2:
        out["co_identity"] = "ids_file"
    return out


def _finalize_id_draft(draft: dict[str, Any]) -> dict[str, Any]:
    out = dict(draft)
    out.setdefault("source", "id_list")
    if not out.get("title"):
        out["title"] = _shell_title_for(out)
    if _n_strong_ids(out) >= 2:
        out["co_identity"] = out.get("co_identity") or "ids_file"
    return out


def _group_ids_complementary(drafts: list[dict[str, Any]]) -> bool:
    if len(drafts) < 2:
        return True
    acc: dict[str, Any] = {}
    for d in drafts:
        if _ids_conflict(acc, d):
            return False
        acc = _merge_id_drafts(acc, d)
    return True


def _try_alternating_pairs(singles: list[dict[str, Any]]) -> list[dict[str, Any]] | None:
    """Partition into (line_i, line_i+1) pairs of different single strong-ID kinds."""
    if len(singles) < 2 or len(singles) % 2 != 0:
        return None
    pairs: list[dict[str, Any]] = []
    for i in range(0, len(singles), 2):
        a, b = singles[i], singles[i + 1]
        if _n_strong_ids(a) != 1 or _n_strong_ids(b) != 1:
            return None
        if _strong_kinds(a) & _strong_kinds(b):
            return None
        if _ids_conflict(a, b):
            return None
        pairs.append(_merge_id_drafts(a, b))
    return pairs


def _emit_id_group(drafts: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Resolve one blank-delimited (or whole-file) group into publication drafts."""
    if not drafts:
        return []
    out: list[dict[str, Any]] = []
    # Lines that already carry 2+ IDs are complete publications
    complete = [d for d in drafts if _n_strong_ids(d) >= 2]
    singles = [d for d in drafts if _n_strong_ids(d) == 1]
    titles = [d for d in drafts if _n_strong_ids(d) == 0 and d.get("title")]

    for d in complete:
        out.append(_finalize_id_draft(d))

    if singles:
        if len(singles) >= 2 and _group_ids_complementary(singles):
            # e.g. doi + pmid (+ pmcid) for the same paper in one group
            merged: dict[str, Any] = {}
            for d in singles:
                merged = _merge_id_drafts(merged, d) if merged else dict(d)
            out.append(_finalize_id_draft(merged))
        else:
            paired = _try_alternating_pairs(singles)
            if paired is not None:
                out.extend(_finalize_id_draft(x) for x in paired)
            else:
                # If any adjacent singles share an ID kind (doi,doi,…), treat as
                # block layout — do NOT cross-pair with later PMIDs.
                adjacent_same_kind = any(
                    _strong_kinds(singles[i]) & _strong_kinds(singles[i + 1])
                    for i in range(len(singles) - 1)
                )
                if adjacent_same_kind:
                    out.extend(_finalize_id_draft(d) for d in singles)
                else:
                    # Greedy adjacent complementary pairs (doi then pmid), else singleton.
                    i = 0
                    while i < len(singles):
                        a = singles[i]
                        if (
                            i + 1 < len(singles)
                            and _n_strong_ids(a) == 1
                            and _n_strong_ids(singles[i + 1]) == 1
                            and not (_strong_kinds(a) & _strong_kinds(singles[i + 1]))
                            and not _ids_conflict(a, singles[i + 1])
                        ):
                            out.append(
                                _finalize_id_draft(_merge_id_drafts(a, singles[i + 1]))
                            )
                            i += 2
                        else:
                            out.append(_finalize_id_draft(a))
                            i += 1

    for d in titles:
        out.append(_finalize_id_draft(d))
    return out


def _parse_id_list(text: str) -> list[dict[str, Any]]:
    """Parse an ID list into publication drafts with parser-level co-identity.

    Co-identity (one publication) when:
      1) Multiple strong IDs on the **same line**, or
      2) Within a blank-line group (or the whole file if no blanks), lines carry
         complementary strong IDs with no same-kind conflicts, or
      3) The group is an even-length alternating dual-ID sequence
         (doi,pmid,doi,pmid,…) / (pmid,doi,…).

    Never merges "any doi + any pmid in the file" across unrelated blocks.
    """
    groups: list[list[dict[str, Any]]] = [[]]
    for line in text.splitlines():
        raw = line.strip()
        if not raw:
            if groups[-1]:
                groups.append([])
            continue
        if raw.startswith("#"):
            continue
        groups[-1].append(_parse_id_tokens_from_line(raw))
    if groups and not groups[-1]:
        groups.pop()

    entries: list[dict[str, Any]] = []
    for g in groups:
        entries.extend(_emit_id_group(g))
    return entries


_DOI_RE = re.compile(r"\b(10\.\d{4,9}/[^\s\"<>\[\]{}]+)", re.I)
_PMID_LABEL_RE = re.compile(r"\b(?:PMID|PubMed\s*ID|PubMed)\s*[:#]?\s*(\d{5,})", re.I)
_PMCID_RE = re.compile(r"\b(PMC\d+)\b", re.I)
_ARXIV_RE = re.compile(
    r"\b(?:arXiv[:\s]*)(\d{4}\.\d{4,5})(?:v\d+)?\b",
    re.I,
)


def _clean_doi_token(raw: str) -> Optional[str]:
    s = (raw or "").strip().rstrip(".,;:)")
    s = re.sub(r"^https?://(dx\.)?doi\.org/", "", s, flags=re.I)
    s = re.sub(r"^doi:\s*", "", s, flags=re.I)
    if re.match(r"10\.\d{4,9}/", s):
        return s
    return None


def extract_ids_from_text(text: str) -> dict[str, Any]:
    """Pull DOI / PMID / PMCID / arXiv from free text (sidecar body or filename)."""
    out: dict[str, Any] = {}
    if not text:
        return out
    # PDF text wraps long DOIs right after the slash: "10.1038/\ns41575-020-0342-4".
    text = re.sub(r"(\b10\.\d{4,9}/)[ \t]*\n[ \t]*(?=\S)", r"\1", text)
    m = _DOI_RE.search(text)
    if m:
        doi = _clean_doi_token(m.group(1))
        if doi:
            out["doi"] = doi
    m = _PMID_LABEL_RE.search(text)
    if m:
        out["pmid"] = m.group(1)
    m = _PMCID_RE.search(text)
    if m:
        out["pmcid"] = m.group(1).upper()
    m = _ARXIV_RE.search(text)
    if m:
        out["arxiv"] = m.group(1)
    return out


def extract_ids_from_filename(name: str) -> dict[str, Any]:
    """Best-effort IDs from a local filename stem.

    Handles ``10.1016_j.cell.2017.05.035.pdf``, ``PMID28622514_…``,
    ``doi-10.1101-2020.01.01.123456.pdf``, ``PMC1234567_…``.
    """
    stem = Path(name).stem
    out: dict[str, Any] = {}
    # Underscore / dash DOI forms → restore slash after 10.xxxx
    m = re.search(r"(10\.\d{4,9})[_-](.+)$", stem, flags=re.I)
    if m:
        candidate = f"{m.group(1)}/{m.group(2).replace('_', '.')}"
        # trim trailing title junk after common DOI end patterns is hard; take until space
        candidate = candidate.split()[0].rstrip(".,;")
        doi = _clean_doi_token(candidate)
        if doi:
            out["doi"] = doi
    m = re.search(r"(?i)(?:pmid|pubmed)[_-]?(\d{5,})", stem)
    if m:
        out["pmid"] = m.group(1)
    m = re.search(r"(?i)\b(PMC\d+)\b", stem)
    if m:
        out["pmcid"] = m.group(1).upper()
    m = re.search(r"(?i)arxiv[_-]?(\d{4}\.\d{4,5})", stem)
    if m:
        out["arxiv"] = m.group(1)
    # Also run free-text extractor on stem with underscores→slash for doi.org style
    soft = extract_ids_from_text(stem.replace("_", "/"))
    for k, v in soft.items():
        out.setdefault(k, v)
    return out


def file_content_hash(path: Path) -> str:
    """SHA-256 of local file bytes (authorized files only; not fetched OA bodies)."""
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _merge_extracted_ids(draft: dict[str, Any], extracted: dict[str, Any]) -> dict[str, Any]:
    out = dict(draft)
    for k in ("doi", "pmid", "pmcid", "arxiv", "openalex"):
        if extracted.get(k) and not out.get(k):
            out[k] = extracted[k]
    return out


def _parse_pdf_text_sidecar(path: Path) -> list[dict[str, Any]]:
    """Expect .txt or .json sidecar next to / instead of PDF binary.

    JSON: {title, abstract|text, doi?, pmid?, authors?, year?, url?}
    TXT: optional YAML-like header --- then body as abstract.
    Also scrapes DOI/PMID from body + filename when header fields are missing.
    """
    fname_ids = extract_ids_from_filename(path.name)
    text = path.read_text(encoding="utf-8", errors="replace")
    if path.suffix.lower() == ".json":
        obj = json.loads(text)
        if isinstance(obj, list):
            rows = []
            for o in obj:
                row = {**o, "source": o.get("source") or "pdf_text"}
                row = _merge_extracted_ids(row, extract_ids_from_text(json.dumps(o, ensure_ascii=False)))
                row = _merge_extracted_ids(row, fname_ids)
                row["_local_file"] = str(path.resolve())
                row["ingest_origin"] = "local_file"
                rows.append(row)
            return rows
        obj.setdefault("source", "pdf_text")
        if "abstract" not in obj and "text" in obj:
            obj["abstract"] = obj["text"]
        blob = " ".join(
            str(obj.get(k) or "")
            for k in ("title", "abstract", "text", "doi", "pmid", "url")
        )
        obj = _merge_extracted_ids(obj, extract_ids_from_text(blob))
        obj = _merge_extracted_ids(obj, fname_ids)
        obj["_local_file"] = str(path.resolve())
        obj["ingest_origin"] = "local_file"
        return [obj]
    # TXT with optional front matter
    title = path.stem
    doi = None
    pmid = None
    pmcid = None
    year = None
    authors: list[str] = []
    url = None
    body = text
    if text.startswith("---"):
        parts = text.split("---", 2)
        if len(parts) >= 3:
            header, body = parts[1], parts[2]
            for line in header.splitlines():
                if ":" not in line:
                    continue
                k, v = line.split(":", 1)
                k, v = k.strip().lower(), v.strip()
                if k == "title":
                    title = v
                elif k == "doi":
                    doi = v
                elif k == "pmid":
                    pmid = v
                elif k == "pmcid" or k == "pmc":
                    pmcid = v
                elif k == "year":
                    year = v
                elif k == "url":
                    url = v
                elif k == "authors":
                    authors = [a.strip() for a in v.split(";") if a.strip()]
    draft: dict[str, Any] = {
        "source": "pdf_text",
        "title": title,
        "abstract": body.strip()[:20000],
        "doi": doi,
        "pmid": pmid,
        "pmcid": pmcid,
        "year": year,
        "authors": authors,
        "url": url,
        "ingest_origin": "local_file",
        "_local_file": str(path.resolve()),
    }
    draft = _merge_extracted_ids(draft, extract_ids_from_text(text))
    draft = _merge_extracted_ids(draft, fname_ids)
    # Prefer a human title when stem was only an ID
    if draft.get("title") == path.stem and (draft.get("doi") or draft.get("pmid")):
        # keep stem as shell title; enrich will replace
        pass
    return [draft]


def _parse_pdf_binary(
    path: Path,
    *,
    try_ocr: bool = True,
    ocr_cmd: Optional[str] = None,
    pdftotext_bin: Optional[str] = None,
    tesseract_bin: Optional[str] = None,
) -> list[dict[str, Any]]:
    """Parse local PDF via sidecar first, else pdftotext (+ optional OCR soft-fail).

    Path-mounts via pub_files + content_hash; does NOT store PDF bytes in DB.
    Missing OCR tools → JSON warnings on the draft, never a hard parse failure.
    """
    from pdf_extract import (  # local sibling
        extract_pdf_text,
        first_nonzero_title_line,
        front_matter_text,
        guess_abstract,
        guess_title,
        pdf_info_title,
    )

    sidecar_txt = path.with_suffix(".txt")
    sidecar_json = path.with_suffix(".json")
    alt_txt = Path(str(path) + ".txt")
    alt_json = Path(str(path) + ".json")
    for sc in (sidecar_txt, sidecar_json, alt_txt, alt_json):
        if sc.is_file():
            rows = _parse_pdf_text_sidecar(sc)
            for r in rows:
                r["_local_file"] = str(path.resolve())
                r["_sidecar_file"] = str(sc.resolve())
                r["ingest_origin"] = "local_file"
                r["extract_method"] = "sidecar"
            return rows

    ids = extract_ids_from_filename(path.name)
    extract = extract_pdf_text(
        path,
        try_ocr=try_ocr,
        ocr_cmd=ocr_cmd,
        pdftotext_bin=pdftotext_bin,
        tesseract_bin=tesseract_bin,
    )
    text = extract.get("text") or ""
    # Only the first two pages: later pages hold the reference list, whose DOIs are other papers.
    front_matter = front_matter_text(path, pdftotext_bin=pdftotext_bin)
    if extract.get("method") != "pdftotext":
        front_matter = "\f".join(text.split("\f")[:2])  # OCR text has no reading-order variant
    text_ids = extract_ids_from_text(front_matter[:12000]) if text else {}
    merged_ids = {**ids, **text_ids}
    has_strong_id = any(merged_ids.get(k) for k in ("doi", "pmid", "pmcid", "arxiv"))
    title = pdf_info_title(path)
    if not title:
        guessed = (
            guess_title(front_matter)
            or (first_nonzero_title_line(text, fallback=path.stem) if text.strip() else path.stem)
        )
        # A header line is a poor title. With a strong id keep the file stem as a placeholder
        # so auto-enrich replaces it with the real title; without one, the guess is all we have.
        title = path.stem if has_strong_id else guessed
    # The abstract is only the text under an "Abstract" heading, never the page body
    # (the body is read from the saved file by lit_get/--with_fulltext).
    abstract = guess_abstract(text)
    draft: dict[str, Any] = {
        "source": "local_pdf",
        "title": title,
        "abstract": abstract,
        "ingest_origin": "local_file",
        "_local_file": str(path.resolve()),
        "extract_method": extract.get("method"),
        "extract_warnings": list(extract.get("warnings") or []),
        "extract_steps": list(extract.get("steps") or []),
        **ids,
    }
    draft = _merge_extracted_ids(draft, text_ids)
    return [draft]




def _tag_local(rows: list[dict[str, Any]], path: Path) -> list[dict[str, Any]]:
    out = []
    for r in rows:
        d = dict(r)
        d.setdefault("ingest_origin", "local_file")
        d.setdefault("_local_file", str(path.resolve()))
        out.append(d)
    return out


def load_records(
    path: Path,
    fmt: Optional[str] = None,
    *,
    try_ocr: bool = True,
    ocr_cmd: Optional[str] = None,
    pdftotext_bin: Optional[str] = None,
    tesseract_bin: Optional[str] = None,
) -> list[dict[str, Any]]:
    fmt = (fmt or "").lower() or None
    suffix = path.suffix.lower()
    if fmt == "pdf" or suffix == ".pdf":
        return _parse_pdf_binary(
            path,
            try_ocr=try_ocr,
            ocr_cmd=ocr_cmd,
            pdftotext_bin=pdftotext_bin,
            tesseract_bin=tesseract_bin,
        )
    text = path.read_text(encoding="utf-8", errors="replace")
    if fmt == "hits" or suffix == ".jsonl":
        rows = []
        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line:
                rows.append(json.loads(line))
        return rows
    if fmt == "bib" or suffix == ".bib":
        return _tag_local(_parse_bibtex(text), path)
    if fmt == "ris" or suffix == ".ris":
        rows = _parse_ris(text)
        # scrape IDs from UR/N1 fields already mapped; also filename
        fname_ids = extract_ids_from_filename(path.name)
        out = []
        for r in rows:
            d = _merge_extracted_ids(r, fname_ids)
            if d.get("url"):
                d = _merge_extracted_ids(d, extract_ids_from_text(str(d.get("url"))))
            out.append(d)
        return _tag_local(out, path)
    if fmt == "pdf_text":
        return _parse_pdf_text_sidecar(path)
    if fmt in ("ids", "id_list") or suffix == ".ids":
        return _parse_id_list(text)
    if suffix in (".txt", ".json"):
        # Heuristic: id list vs pdf text sidecar
        if fmt in ("ids", "id_list") or (
            suffix == ".txt"
            and all(
                (
                    not ln.strip()
                    or ln.strip().startswith("#")
                    or "doi" in ln.lower()
                    or ln.strip()[:2].isdigit()
                    or ln.strip().startswith("10.")
                    or re.match(r"(?i)^(pmid|pmcid|pmc|arxiv|openalex)\s*:", ln.strip())
                )
                for ln in text.splitlines()[:20]
            )
            and not text.startswith("---")
        ):
            # Prefer ids when clearly ID-shaped; else sidecar
            if fmt in ("ids", "id_list") or suffix == ".ids" or (
                not text.startswith("---")
                and "title:" not in text[:500].lower()
                and any(
                    tok in text.lower()
                    for tok in ("doi:", "pmid:", "10.")
                )
            ):
                # If front-matter-ish or long prose → sidecar
                if text.startswith("---") or (len(text) > 400 and " " in text[:200]):
                    return _parse_pdf_text_sidecar(path)
                return _parse_id_list(text)
        return _parse_pdf_text_sidecar(path)
    # Fallback by content
    if "@article" in text.lower() or "@inproceedings" in text.lower():
        return _tag_local(_parse_bibtex(text), path)
    if re.search(r"^TY  -", text, re.M):
        return _tag_local(_parse_ris(text), path)
    return _parse_id_list(text)


def ingest_records(
    *,
    db_path: Path,
    collection_id: str,
    records: list[dict[str, Any]],
    run_id: Optional[str] = None,
    verify_oa: bool = False,
    oa_mock_map: Optional[dict[str, dict[str, Any]]] = None,
) -> dict[str, Any]:
    conn = open_db(db_path)
    counts = {"inserted": 0, "merged": 0, "needs_review": 0, "skipped": 0, "total": len(records)}
    details: list[dict[str, Any]] = []
    files_recorded: list[dict[str, Any]] = []
    oa_checks: list[dict[str, Any]] = []
    merge_out: dict[str, Any] = {"n_applied": 0, "applied": []}
    try:
        col = conn.execute("SELECT id FROM collections WHERE id = ?", (collection_id,)).fetchone()
        if not col:
            return {"ok": False, "error": f"unknown collection_id {collection_id}"}
        for rec in records:
            if not any(rec.get(k) for k in ("doi", "pmid", "pmcid", "arxiv", "openalex", "title")):
                counts["skipped"] += 1
                if run_id:
                    record_run_hit(
                        conn,
                        run_id=run_id,
                        publication_id=None,
                        raw=rec,
                        source=str(rec.get("source") or "unknown"),
                        dedupe_action="skipped",
                    )
                continue
            result = resolve_or_create(conn, rec)
            action = result["action"]
            counts[action] = counts.get(action, 0) + 1
            pub_id = result["publication_id"]
            link_to_collection(conn, collection_id, pub_id, source_run_id=run_id)
            if run_id:
                record_run_hit(
                    conn,
                    run_id=run_id,
                    publication_id=pub_id,
                    raw=rec,
                    source=str(rec.get("source") or "unknown"),
                    dedupe_action=action,
                )
            local_path = rec.get("_local_file") or rec.get("file_path")
            # Reference lists (.bib/.ris/.jsonl/...) describe papers; they are not a paper's file.
            if (
                local_path
                and Path(str(local_path)).suffix.lower() not in _METADATA_ONLY_SUFFIXES
                and str(rec.get("ingest_origin") or "") in ("local_file", "pdf_text", "local_pdf", "")
            ):
                p = Path(str(local_path))
                content_hash = None
                if p.is_file():
                    try:
                        content_hash = file_content_hash(p)
                    except OSError:
                        content_hash = None
                origin = str(rec.get("ingest_origin") or "local_file")
                if origin in ("pdf_text", "local_pdf", "bibtex", "ris", "bib"):
                    origin = "local_file"
                meta = record_pub_file(
                    conn,
                    publication_id=pub_id,
                    file_path=str(p.resolve() if p.exists() else p),
                    content_hash=content_hash,
                    ingest_origin=origin,
                    access_origin=rec.get("access_origin"),
                    access_method=rec.get("access_method"),
                )
                # optional sidecar also recorded
                sc = rec.get("_sidecar_file")
                if sc:
                    sp = Path(str(sc))
                    sc_hash = None
                    if sp.is_file():
                        try:
                            sc_hash = file_content_hash(sp)
                        except OSError:
                            sc_hash = None
                    record_pub_file(
                        conn,
                        publication_id=pub_id,
                        file_path=str(sp.resolve() if sp.exists() else sp),
                        content_hash=sc_hash,
                        ingest_origin="local_file",
                    )
                files_recorded.append(meta)
                result = dict(result)
                result["file"] = meta
            if verify_oa:
                oa_res = verify_record_oa(rec, mock_map=oa_mock_map)
                set_publication_oa(conn, pub_id, oa_res)
                oa_checks.append({"publication_id": pub_id, **oa_res})
                result = dict(result)
                result["oa"] = {
                    "oa_status": oa_res.get("oa_status"),
                    "url": oa_res.get("url"),
                    "stored_pdf_bytes": False,
                }
            details.append(result)
        # Post-pass: merge complementary DOI/PMID shells created across this batch
        # (and any pre-existing high-confidence clusters). Soft conflicts untouched.
        merge_out = apply_id_merges(conn, only_high_confidence=True)
        conn.commit()
    finally:
        conn.close()
    return {
        "ok": True,
        "mode": "live",
        "collection_id": collection_id,
        "run_id": run_id,
        "counts": counts,
        "n_details": len(details),
        "n_files_recorded": len(files_recorded),
        "files_recorded": files_recorded,
        "oa_checks": oa_checks,
        "id_merges": {
            "n_applied": merge_out.get("n_applied", 0),
            "applied": merge_out.get("applied") or [],
        },
    }


def ingest_hits_file(
    *,
    db_path: Path | str,
    collection_id: str,
    run_id: str,
    hits_path: Path | str,
) -> dict[str, Any]:
    path = Path(hits_path)
    records = load_records(path, fmt="hits")
    return ingest_records(
        db_path=Path(db_path),
        collection_id=collection_id,
        records=records,
        run_id=run_id,
    )


def _looks_like_ids_ingest(path: Path, fmt: Optional[str], records: list[dict[str, Any]]) -> bool:
    """True when this ingest came from an ID list (shells expected)."""
    f = (fmt or "").lower()
    if f in ("ids", "id_list"):
        return True
    if path.suffix.lower() == ".ids":
        return True
    if not records:
        return False
    # All rows have ≥1 strong id and no non-shell title/abstract → treat as ids
    strong = ("doi", "pmid", "pmcid", "arxiv", "openalex")
    for r in records:
        if not any(r.get(k) for k in strong):
            return False
        if (r.get("abstract") or "").strip():
            return False
        title = (r.get("title") or "").strip()
        if title and not any(
            title.lower() == str(r.get(k) or "").lower()
            or title.lower() == f"{k}:{r.get(k)}".lower()
            for k in strong
            if r.get(k)
        ):
            return False
    return True


def _should_auto_enrich(path: Path, fmt: Optional[str], records: list[dict[str, Any]]) -> bool:
    """Auto-enrich after ids ingest OR local-file shells missing title/abstract."""
    if _looks_like_ids_ingest(path, fmt, records):
        return True
    f = (fmt or "").lower()
    suffix = path.suffix.lower()
    localish = f in ("pdf_text", "pdf", "bib", "ris") or suffix in (
        ".bib", ".ris", ".pdf", ".txt", ".json"
    )
    if not localish:
        # records tagged local_file also qualify
        localish = any(
            (r.get("ingest_origin") == "local_file") or r.get("_local_file")
            for r in records
        )
    if not localish or not records:
        return False
    from enrich import needs_enrich  # noqa: WPS433 - one rule for "does a registry lookup help?"

    return any(
        needs_enrich({**{k: r.get(k) for k in ("title", "abstract", "doi", "pmid", "arxiv")}, "authors_json": json.dumps(r.get("authors") or [])})
        for r in records
    )


_DIR_SUFFIXES = {".pdf", ".bib", ".ris", ".jsonl"}


def _expand_input(folder: Path) -> list[Path]:
    """Files to import from a folder: PDFs and reference-manager exports, recursively."""
    return sorted(p for p in folder.rglob("*") if p.is_file() and p.suffix.lower() in _DIR_SUFFIXES)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--project_root", required=True)
    p.add_argument("--collection_id", default="")
    p.add_argument("--input_path", required=True)
    p.add_argument(
        "--format",
        default="",
        help="hits|bib|ris|ids|pdf_text|pdf (auto if empty)",
    )
    p.add_argument("--run_id", default="")
    p.add_argument(
        "--access_origin",
        default="",
        help="Label for how the file was obtained (e.g. institutional, purchased, author). "
        "Only for files the user is entitled to hold; recorded, never used to fetch anything.",
    )
    p.add_argument(
        "--no_enrich",
        action="store_true",
        help="Skip auto-enrich of shells after local/ids ingest (default: enrich)",
    )
    p.add_argument(
        "--enrich_mock",
        action="store_true",
        help="When auto-enriching, use offline mock fixture instead of HTTP",
    )
    p.add_argument(
        "--verify_oa",
        action="store_true",
        help="Soft-verify candidate OA/PDF URLs (HEAD/GET); never store PDF bytes",
    )
    p.add_argument(
        "--oa_mock_fixture",
        default="",
        help="JSON map of url→{oa_status|http_status,content_type} for offline OA verify",
    )
    p.add_argument(
        "--no_ocr",
        action="store_true",
        help="Skip OCR fallback when pdftotext yields empty/near-empty text",
    )
    p.add_argument(
        "--ocr_cmd",
        default="",
        help="Optional external OCR command template ({path},{outdir},{out}); "
        "overrides LIT_PDF_OCR_CMD. Soft-fail if it errors.",
    )
    p.add_argument(
        "--pdftotext_bin",
        default="",
        help="Override pdftotext binary path (tests / non-PATH installs)",
    )
    args = p.parse_args()

    root = Path(args.project_root).resolve()
    db_path = root / "library.sqlite"
    from common import open_project_db, resolve_collection  # noqa: E402

    _conn, _err = open_project_db(args.project_root)
    if _err:
        print(json.dumps({"ok": False, "error": _err}))
        return 1
    try:
        _col, _err = resolve_collection(_conn, args.collection_id or None)
    finally:
        _conn.close()
    if _err:
        print(json.dumps({"ok": False, "error": _err}))
        return 1
    args.collection_id = _col["id"]
    path = Path(args.input_path)
    if not path.exists():
        print(json.dumps({"ok": False, "error": f"missing input {path}"}))
        return 1
    oa_mock_map = None
    if args.oa_mock_fixture:
        oa_mock_map = json.loads(Path(args.oa_mock_fixture).read_text(encoding="utf-8"))
    files = _expand_input(path) if path.is_dir() else [path]
    if not files:
        print(json.dumps({"ok": False, "error": f"no .pdf/.bib/.ris/.txt/.jsonl files under {path}"}))
        return 1
    records: list[dict[str, Any]] = []
    try:
        for fp in files:
            records.extend(
                load_records(
                    fp,
                    fmt=(args.format or None) if not path.is_dir() else None,
                    try_ocr=not bool(args.no_ocr),
                    ocr_cmd=args.ocr_cmd or None,
                    pdftotext_bin=args.pdftotext_bin or None,
                )
            )
    except Exception as e:  # noqa: BLE001
        print(json.dumps({"ok": False, "error": f"parse failed: {e}"}))
        return 1
    if args.access_origin:
        for r in records:
            r["access_origin"] = args.access_origin
            r.setdefault("access_method", "user_provided_file")
    out = ingest_records(
        db_path=db_path,
        collection_id=args.collection_id,
        records=records,
        run_id=args.run_id or None,
        verify_oa=bool(args.verify_oa),
        oa_mock_map=oa_mock_map,
    )
    out["input_path"] = str(path)
    out["n_parsed"] = len(records)
    out["verify_oa"] = bool(args.verify_oa)
    extract_warnings = []
    for r in records:
        for w in r.get("extract_warnings") or []:
            extract_warnings.append(w)
        if r.get("extract_method"):
            out.setdefault("extract_methods", [])
            if r["extract_method"] not in out["extract_methods"]:
                out["extract_methods"].append(r["extract_method"])
    if extract_warnings:
        out["extract_warnings"] = extract_warnings

    do_enrich = (
        out.get("ok")
        and not args.no_enrich
        and _should_auto_enrich(path, args.format or None, records)
    )
    if do_enrich:
        from enrich import enrich_publications  # noqa: WPS433

        enrich_out = enrich_publications(
            db_path=db_path,
            collection_id=args.collection_id,
            dry_run=False,
            mock=bool(args.enrich_mock),
            verify_oa=bool(args.verify_oa),
            oa_mock_map=oa_mock_map,
        )
        out["enrich"] = {
            "n_candidates": enrich_out.get("n_candidates"),
            "n_enriched": enrich_out.get("n_enriched"),
            "n_failed": enrich_out.get("n_failed"),
            "n_skipped": enrich_out.get("n_skipped"),
            "mock": enrich_out.get("mock"),
            "oa_checks": enrich_out.get("oa_checks") or [],
            "details": enrich_out.get("details") or [],
        }
        out["enrich_ok"] = bool(enrich_out.get("ok"))
    elif args.no_enrich:
        out["enrich"] = {"skipped": True, "reason": "no_enrich"}

    print(json.dumps(out, ensure_ascii=False))
    return 0 if out.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
