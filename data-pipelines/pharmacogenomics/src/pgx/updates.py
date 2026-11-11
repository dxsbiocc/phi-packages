"""Report which ORCESTRA PharmacoSets are newer than, or missing from, sources.yaml.

The check reads the public ORCESTRA registry and compares it with the ``orcestra_*``
sources. It never downloads a PharmacoSet and never edits sources.yaml.

ORCESTRA is not a clean version feed, so the report is conservative:

* A Zenodo record id is the identity of a deposit. Several configured records are the
  original deposits and ORCESTRA now lists a later re-deposit under a different id. When
  a configured record is not in the registry the status is ``record_not_listed`` and the
  candidates are shown; the check does not claim they are newer.
* ``dateCreated`` is only used to order entries of one dataset family. Some registry
  dates are placeholders, so a ``newer_available`` result is a prompt to look, not proof.
* Variants (unfiltered, u133a, old array) are different cuts of a dataset, not newer
  versions, and are left out of the comparison.
"""

from __future__ import annotations

import json
import re
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from pgx.download import fetch_bytes

REGISTRY_URL = "https://www.orcestra.ca/api/pset/available"
VARIANT_MARKERS = ("unfiltered", "u133a", "oldarray")
# Families that are also served by a non-ORCESTRA source, so an ORCESTRA entry is not a gap.
COVERED_BY = {
    "gdsc": ("gdsc",),
    "gdsc2": ("gdsc_combinations",),
    "prism": ("prism_19q4", "prism_24q2"),
    "nci60": ("nci60_dtp",),
}

_RECORD = re.compile(r"zenodo\.(\d+)|zenodo\.org/records?/(\d+)")

Fetch = Callable[[str], tuple[bytes, str]]


@dataclass(frozen=True)
class Entry:
    name: str
    family: str
    doi: str | None
    record_ids: frozenset[str]
    download_link: str | None
    created: str | None
    variant: bool

    def summary(self) -> dict[str, object]:
        return {"name": self.name, "doi": self.doi, "created": self.created, "download_link": self.download_link}


def _family(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower().replace("²", "2"))


def record_ids(*texts: str | None) -> frozenset[str]:
    found: set[str] = set()
    for text in texts:
        for match in _RECORD.finditer(text or ""):
            found.add(match.group(1) or match.group(2))
    return frozenset(found)


def _created(value: object) -> str | None:
    """A usable ISO date, or None for blanks and the 1970 placeholder."""
    text = str(value or "").strip()
    if not text or text.startswith("1970-"):
        return None
    return text


def parse_registry(payload: object) -> list[Entry]:
    if not isinstance(payload, list):
        raise ValueError("ORCESTRA registry is not a list")
    entries: list[Entry] = []
    for item in payload:
        if not isinstance(item, dict) or not item.get("name"):
            continue
        name = str(item["name"])
        dataset = item.get("dataset") if isinstance(item.get("dataset"), dict) else {}
        family = _family(str(dataset.get("name") or name.split("_")[0]))
        doi = str(item.get("doi") or "").strip() or None
        link = str(item.get("downloadLink") or "").strip() or None
        entries.append(Entry(
            name=name,
            family=family,
            doi=doi,
            record_ids=record_ids(doi, link),
            download_link=link,
            created=_created(item.get("dateCreated")),
            variant=any(marker in _family(name) for marker in VARIANT_MARKERS),
        ))
    return entries


def fetch_registry(fetch: Fetch = fetch_bytes, url: str = REGISTRY_URL) -> list[Entry]:
    payload, _content_type = fetch(url)
    return parse_registry(json.loads(payload))


def load_registry_file(path: Path) -> list[Entry]:
    return parse_registry(json.loads(path.read_text()))


@dataclass(frozen=True)
class Tracked:
    id: str
    family: str
    version: str
    enabled: bool
    planned: bool
    record_ids: frozenset[str]


def tracked_sources(document: dict) -> list[Tracked]:
    found: list[Tracked] = []
    for planned, key in ((False, "sources"), (True, "planned")):
        for source in document.get(key) or []:
            source_id = str(source.get("id") or "")
            if not source_id.startswith("orcestra_"):
                continue
            urls = [str(file.get("url") or "") for file in source.get("files") or []]
            found.append(Tracked(
                id=source_id,
                family=_family(str(source.get("dataset") or source_id.removeprefix("orcestra_"))),
                version=str(source.get("version") or ""),
                enabled=bool(source.get("enabled", True)),
                planned=planned,
                record_ids=record_ids(str(source.get("homepage") or ""), *urls),
            ))
    return found


def _newest_first(entries: list[Entry]) -> list[Entry]:
    return sorted(entries, key=lambda entry: (entry.created or "", entry.name), reverse=True)


def _compare_one(tracked: Tracked, entries: list[Entry]) -> dict[str, object]:
    family = [entry for entry in entries if entry.family == tracked.family and not entry.variant]
    base: dict[str, object] = {"source": tracked.id, "version": tracked.version, "enabled": tracked.enabled}
    matched = [entry for entry in entries if entry.record_ids & tracked.record_ids]
    if tracked.planned:
        return {**base, "status": "planned", "candidates": [entry.summary() for entry in _newest_first(family)]}
    if matched:
        current = max(matched, key=lambda entry: entry.created or "")
        newer = [entry for entry in family if entry.created and current.created and entry.created > current.created]
        return {
            **base,
            "status": "newer_available" if newer else "current",
            "tracked_entry": current.summary(),
            "newer": [entry.summary() for entry in _newest_first(newer)],
        }
    if not family:
        return {**base, "status": "family_not_listed", "candidates": []}
    return {**base, "status": "record_not_listed", "candidates": [entry.summary() for entry in _newest_first(family)]}


def compare(document: dict, entries: list[Entry], since: str | None = None) -> dict[str, object]:
    tracked = tracked_sources(document)
    tracked_families = {item.family for item in tracked}
    enabled_ids = {str(source.get("id")) for source in document.get("sources") or [] if source.get("enabled", True)}
    results = [_compare_one(item, entries) for item in tracked]

    untracked: list[Entry] = []
    covered: list[dict[str, object]] = []
    for entry in _newest_first([entry for entry in entries if not entry.variant]):
        if entry.family in tracked_families:
            continue
        sources = [source for source in COVERED_BY.get(entry.family, ()) if source in enabled_ids]
        if sources:
            covered.append({**entry.summary(), "covered_by": sources})
        elif not since or (entry.created and entry.created[:10] >= since):
            untracked.append(entry)

    counts: dict[str, int] = {}
    for result in results:
        counts[str(result["status"])] = counts.get(str(result["status"]), 0) + 1
    return {
        "registry": {"url": REGISTRY_URL, "entries": len(entries), "variants_ignored": sum(1 for entry in entries if entry.variant)},
        "tracked": results,
        "untracked": [entry.summary() for entry in untracked],
        "covered_elsewhere": covered,
        "summary": {
            **counts,
            "untracked": len(untracked),
            "needs_attention": counts.get("newer_available", 0) + len(untracked),
        },
    }
