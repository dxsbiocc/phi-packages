"""GDSC drug combinations.

A combination name is never a compound. Annotation rows that name two drugs,
and GDSC² anchored-screen COMBI_ID values (`anchor id:library id`), become
combination records. Each component is resolved on its own.
"""

from __future__ import annotations

import re

from pgx.ids import CompoundRecord
from pgx.normalize.compounds import normalize_compound_name

_COMBO_NAME = re.compile(
    r"^(?P<anchor>.+?)\s*[-\s]combo[-\s]*(?P<index>\d+)?\s*\((?P<library>[^)]+)\)\s*$",
    re.IGNORECASE,
)
_DOSE_TAIL = re.compile(
    r"^(?P<name>.+?)\s+(?P<conc>\d+(?:\.\d+)?)\s*(?P<unit>uM|µM|μM|um|nM|mM)\s*$",
    re.IGNORECASE,
)
_ID_PART = re.compile(r"^\d+(?:\|\d+)*$")
_INTERNAL = (
    re.compile(r"^\d+$"),
    re.compile(r"^\d+-D\d+$", re.IGNORECASE),
    re.compile(r"^[A-Za-z][A-Za-z0-9]*_\d+$"),
    re.compile(r"^BAY-[A-Za-z0-9]+_\d+$", re.IGNORECASE),
    re.compile(r"^KIN\d*", re.IGNORECASE),
    re.compile(r"^THZ-", re.IGNORECASE),
    re.compile(r"^QL-", re.IGNORECASE),
    re.compile(r"^TL-\d", re.IGNORECASE),
    re.compile(r"^JW-", re.IGNORECASE),
    re.compile(r"^JQ\d+$", re.IGNORECASE),
    re.compile(r"^VTP-", re.IGNORECASE),
    re.compile(r"^PBD-", re.IGNORECASE),
    re.compile(r"^LMB_", re.IGNORECASE),
    re.compile(r"^THR-", re.IGNORECASE),
    re.compile(r"^CRT\d", re.IGNORECASE),
    re.compile(r"^CCT\d", re.IGNORECASE),
    re.compile(r"^FS\d", re.IGNORECASE),
    re.compile(r"^FY\d", re.IGNORECASE),
    re.compile(r"^ZL\d", re.IGNORECASE),
    re.compile(r"^ICL\d", re.IGNORECASE),
    re.compile(r"^BDP-", re.IGNORECASE),
    re.compile(r"^BPD-", re.IGNORECASE),
    re.compile(r"^BD[OIF]", re.IGNORECASE),
    re.compile(r"^N\d{4}", re.IGNORECASE),
    re.compile(r"^HG-\d", re.IGNORECASE),
    re.compile(r"^XMD", re.IGNORECASE),
    re.compile(r"^Genentech Cpd\b", re.IGNORECASE),
)


def is_annotation_combination(name: str | None) -> bool:
    text = normalize_compound_name(name or "")
    return bool(text and _COMBO_NAME.match(text))


def parse_annotation_combination(drug_id: str, name: str) -> dict[str, str | None] | None:
    """One screened-compound row whose name is a combination, or None.

    The anchor is the named base (`BAY-MPS` in `BAY-MPS-combo-1`). The library
    is the drug in parentheses, with a trailing concentration when present.
    """
    text = normalize_compound_name(name)
    match = _COMBO_NAME.match(text)
    if not match or not drug_id:
        return None
    anchor = match.group("anchor").strip(" -")
    library_name, library_conc, unit = _split_dose(match.group("library").strip())
    return {
        "combination_id": f"gdsc.combination:annotation:{drug_id}",
        "source": "gdsc",
        "screen": "annotation",
        "name": f"{anchor} + {library_name}",
        "anchor_drug_id": None,
        "library_drug_id": None,
        "anchor_name": anchor,
        "library_name": library_name,
        "anchor_concentrations": None,
        "library_concentrations": library_conc,
        "concentration_unit": unit,
        "source_drug_id": drug_id,
    }


def parse_screen_combination(
    anchor_id: str,
    anchor_name: str | None,
    library_id: str,
    library_name: str | None,
    anchor_concentrations: list[float] | None = None,
    library_concentrations: list[float] | None = None,
    screen: str = "gdsc_combination",
) -> dict[str, str | None] | None:
    """One anchored-screen pair. Components stay in anchor, then library, order."""
    if not (_ID_PART.match(str(anchor_id)) and _ID_PART.match(str(library_id))):
        return None
    anchor = normalize_compound_name(anchor_name or "")
    library = normalize_compound_name(library_name or "")
    return {
        "combination_id": f"gdsc.combination:{anchor_id}:{library_id}",
        "source": "gdsc_combinations",
        "screen": screen,
        "name": f"{anchor} + {library}".strip(),
        "anchor_drug_id": anchor_id,
        "library_drug_id": library_id,
        "anchor_name": anchor or None,
        "library_name": library or None,
        "anchor_concentrations": _format_concs(anchor_concentrations or []),
        "library_concentrations": _format_concs(library_concentrations or []),
        "concentration_unit": "uM",
        "source_drug_id": None,
    }


def combination_components(record: dict[str, str | None]) -> list[dict[str, str | None]]:
    """Ordered components. A `|` in a GDSC id is two drugs sharing that role.

    Anchor parts come first, in the source order, then library parts.
    """
    rows: list[dict[str, str | None]] = []
    for role, id_field, name_field, conc_field in (
        ("anchor", "anchor_drug_id", "anchor_name", "anchor_concentrations"),
        ("library", "library_drug_id", "library_name", "library_concentrations"),
    ):
        raw_id = record.get(id_field)
        raw_name = record.get(name_field)
        ids = [part.strip() for part in raw_id.split("|")] if raw_id else [None]
        names = [part.strip() for part in raw_name.split("|")] if raw_name else [None]
        for index in range(max(len(ids), len(names))):
            rows.append({
                "role": role,
                "source_compound_id": ids[index] if index < len(ids) else None,
                "name": names[index] if index < len(names) else None,
                "concentrations": record.get(conc_field),
                "unit": record.get("concentration_unit"),
            })
    return rows


def is_internal_code(name: str | None, synonyms: set[str] | list[str] | None = None) -> bool:
    """Screening codes that are not public compound names.

    Drug-like sponsor codes (GSK, JNJ, AZD, BAY-MPS1, CI-1033, FTY-720) are not
    internal codes. Those still go to PubChem and ChEMBL.
    """
    text = normalize_compound_name(name or "")
    if any(pattern.match(text) for pattern in _INTERNAL):
        return True
    for synonym in synonyms or []:
        folded = synonym.casefold()
        if "wimm synthesis" in folded:
            return True
    return False


def assign_review(record: CompoundRecord) -> str | None:
    """Mark combinations and internal codes so they are not sent to PubChem.

    Returns the review reason when this record should skip name lookup.
    """
    if is_annotation_combination(record.name) or record.status == "combination":
        record.status = "combination"
        record.review_reason = "combination"
        record.match_method = "combination"
        record.confidence = None
        record.pubchem_cid = None
        record.evidence = "combination name is not a single compound"
        return "combination"
    if record.status == "mapped":
        return None
    if record.status == "ambiguous":
        return None
    if record.source == "gdsc" and is_internal_code(record.name, record.synonyms):
        record.status = "unmatched"
        record.review_reason = "internal_code"
        record.match_method = "internal_code"
        record.confidence = None
        record.evidence = "internal screening code; not looked up in PubChem or ChEMBL"
        return "internal_code"
    return None


def chembl_name_candidates(record: CompoundRecord) -> list[str]:
    """Lookup strings for a ChEMBL synonym or preferred-name query."""
    from pgx.normalize.compounds import extra_lookup_names, name_variants, usable_synonym

    found: list[str] = []
    seen: set[str] = set()

    def add(value: str | None) -> None:
        text = normalize_compound_name(value or "")
        key = text.casefold()
        if not text or key in seen:
            return
        seen.add(key)
        found.append(text)

    for original in [record.name, *sorted(record.synonyms)]:
        if not original:
            continue
        if original != record.name and not usable_synonym(original):
            continue
        for _stage, query, _confidence in name_variants(original):
            add(query)
        for extra in extra_lookup_names(original):
            add(extra)
    return found


def _split_dose(text: str) -> tuple[str, str | None, str | None]:
    match = _DOSE_TAIL.match(text.strip())
    if not match:
        return text.strip(), None, None
    unit = match.group("unit")
    if unit.casefold() == "um":
        unit = "uM"
    return match.group("name").strip(), match.group("conc"), unit


def _format_concs(values: list[float]) -> str | None:
    if not values:
        return None
    rendered: list[str] = []
    for value in sorted(set(values)):
        text = f"{value:.6g}"
        if text not in rendered:
            rendered.append(text)
    return ";".join(rendered)
