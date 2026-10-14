"""Compound name normalization and PubChem / ChEMBL resolution.

Salt forms and stereo parentheticals are removed only from the lookup string.
The source name stored on the record is unchanged. A lookup that returns more
than one CID is ambiguous and stops; a looser name is not tried after that.
An NSC number is queried as ``NSC740`` and ``NSC 740`` after a structure and
before the stored name.
"""

from __future__ import annotations

import re
import unicodedata
from collections.abc import Iterable
from typing import Protocol

from pgx.ids import CompoundRecord, apply_pubchem

_NSC = re.compile(r"^NSC\s*(\d+)$", re.IGNORECASE)

SALT_WORDS = (
    "dihydrochloride",
    "dihydrobromide",
    "hydrochloride",
    "hydrobromide",
    "methanesulfonate",
    "methochloride",
    "monohydrate",
    "hemihydrate",
    "dihydrate",
    "disodium",
    "mesylate",
    "mesilate",
    "maleate",
    "sulfate",
    "sulphate",
    "acetate",
    "citrate",
    "phosphate",
    "tartrate",
    "fumarate",
    "besylate",
    "besilate",
    "tosylate",
    "camsylate",
    "succinate",
    "benzoate",
    "oxalate",
    "lactate",
    "nitrate",
    "bromide",
    "chloride",
    "iodide",
    "sodium",
    "potassium",
    "calcium",
    "hydrate",
    "2hcl",
    "2hbr",
    "hcl",
    "hbr",
    "salt",
)
_SALT = "(?:" + "|".join(SALT_WORDS) + ")"
_SALT_TRAIL = re.compile(rf"(?:[\s\-]+|\s*\(){_SALT}\)?\s*$", re.IGNORECASE)
_STEREO = re.compile(
    r"(?:\s*-\s*)?\(\s*(?:[+-]|±|\+/-|R|S|RS|rac|E|Z|cis|trans)\s*\)",
    re.IGNORECASE,
)
_BRACKET = re.compile(r"\s*\[[^\]]*\]\s*$")
_DOSE_PAREN = re.compile(r"\s*\(\s*\d+(?:\.\d+)?\s*(?:um|µm|μm|nm|mm)\s*\)\s*$", re.IGNORECASE)
_PAREN_NAME = re.compile(r"^(.+?)\s*\(([^)]+)\)\s*$")
_STEREO_INNER = re.compile(r"^(?:[+-]|±|\+/-|R|S|RS|rac|E|Z|cis|trans)$", re.IGNORECASE)
_DOSE_INNER = re.compile(r"^\d+(?:\.\d+)?\s*(?:um|µm|μm|nm|mm)$", re.IGNORECASE)


class CompoundClient(Protocol):
    def cids(self, kind: str, query: str) -> list[str]: ...

    def properties_for(self, cids: list[str]) -> dict[str, dict[str, str | None]]: ...

    def chembl_for(self, inchikeys: list[str]) -> dict[str, list[str]]: ...


def normalize_compound_name(name: str) -> str:
    text = unicodedata.normalize("NFKC", name)
    text = text.replace("\u2212", "-").replace("\u2013", "-").replace("\u2014", "-").replace("\u00a0", " ")
    return re.sub(r"\s+", " ", text).strip()


def strip_stereo(name: str) -> str:
    text = name
    previous = None
    while previous != text:
        previous = text
        text = _STEREO.sub(" ", text)
        text = re.sub(r"\s+", " ", text).strip()
    bracket = _BRACKET.search(text)
    if bracket and re.search(r"[RS]", bracket.group(0)):
        text = text[: bracket.start()]
    return text.strip(" -")


def strip_salts(name: str) -> str:
    text = name.strip()
    while True:
        match = _SALT_TRAIL.search(text)
        if not match:
            return text
        nxt = text[: match.start()].strip(" -")
        if len(re.sub(r"[^A-Za-z0-9]", "", nxt)) < 3:
            return text
        text = nxt


def strip_dose(name: str) -> str:
    """Remove a trailing concentration parenthetical from a lookup string."""
    return _DOSE_PAREN.sub("", name).strip()


def name_variants(original: str) -> list[tuple[str, str, str]]:
    """Lookup stages for one source string: `(stage, query, confidence)`.

    Stages are `exact`, `stereo`, `salt`, and `dose`. The original string is not modified.
    """
    base = normalize_compound_name(original)
    if not base:
        return []
    stereo = strip_stereo(base)
    salted = strip_salts(stereo)
    dosed = strip_dose(salted)
    variants = [("exact", base, "medium")]
    if stereo and stereo.casefold() != base.casefold():
        variants.append(("stereo", stereo, "low"))
    if salted and salted.casefold() not in {base.casefold(), stereo.casefold()}:
        variants.append(("salt", salted, "low"))
    seen = {query.casefold() for _stage, query, _confidence in variants}
    if dosed and dosed.casefold() not in seen:
        variants.append(("dose", dosed, "low"))
    return variants


def usable_synonym(text: str) -> bool:
    if not text or text[0] in ")]};,":
        return False
    folded = text.casefold()
    if "back-up solution" in folded or "wimm synthesis" in folded:
        return False
    return text.count("(") == text.count(")") and text.count("[") == text.count("]")


def extra_lookup_names(original: str) -> list[str]:
    """Comma-separated synonym lists and non-dose parenthetical aliases.

    These are extra lookup strings. They are not stored as the compound name.
    A combination name is not split into single-compound queries. A comma
    between digits is a locant (`1,5-`), not a synonym list. A short
    parenthetical counterion such as `TFA` is not a second drug.
    """
    from pgx.normalize.combinations import is_annotation_combination

    base = normalize_compound_name(original)
    if not base or is_annotation_combination(base):
        return []
    extras: list[str] = []
    if "," in base and not re.search(r"\d\s*,\s*\d", base):
        extras.extend(part.strip() for part in base.split(",") if part.strip())
    match = _PAREN_NAME.match(base)
    if match:
        body, inner = match.group(1).strip(), match.group(2).strip()
        if not _STEREO_INNER.match(inner) and not _DOSE_INNER.match(inner):
            extras.append(body)
            if _alias_inner(inner):
                extras.append(inner)
    return [item for item in extras if item and usable_synonym(item)]


def _alias_inner(inner: str) -> bool:
    token = inner.strip()
    if len(token) <= 4 or token.casefold() in SALT_WORDS:
        return False
    return True


def iter_queries(record: CompoundRecord) -> list[tuple[str, str, str, str]]:
    """Ordered `(kind, query, method, confidence)` lookups for one compound."""
    from pgx.normalize.combinations import is_annotation_combination

    if record.status == "combination" or is_annotation_combination(record.name):
        return []
    queries: list[tuple[str, str, str, str]] = []
    seen: set[tuple[str, str]] = set()

    def add(kind: str, query: str, method: str, confidence: str) -> None:
        key = (kind, query.casefold())
        if not query or key in seen:
            return
        seen.add(key)
        queries.append((kind, query, method, confidence))

    if record.inchikey:
        add("inchikey", record.inchikey, "pubchem_inchikey", "high")
    if record.smiles:
        add("smiles", record.smiles, "pubchem_smiles", "high")
    for token in _nsc_tokens(record):
        add("name", token, "pubchem_nsc", "high")
        add("name", token.replace("NSC", "NSC ", 1), "pubchem_nsc", "high")
    if record.name:
        _add_name(add, record.name, "pubchem_name", "pubchem_name_stereo", "pubchem_name_salt", "pubchem_name_dose")
        for extra in extra_lookup_names(record.name):
            _add_name(add, extra, "pubchem_synonym", "pubchem_synonym_stereo", "pubchem_synonym_salt", "pubchem_synonym_dose")
    for synonym in sorted(record.synonyms):
        if not usable_synonym(synonym):
            continue
        _add_name(add, synonym, "pubchem_synonym", "pubchem_synonym_stereo", "pubchem_synonym_salt", "pubchem_synonym_dose")
    if record.brd_id:
        add("name", record.brd_id, "pubchem_brd", "medium")
    return queries


def _nsc_tokens(record: CompoundRecord) -> list[str]:
    """PubChem stores NCI numbers as ``NSC740`` and ``NSC 740`` synonyms."""
    found: list[str] = []
    seen: set[str] = set()
    for raw in (record.source_compound_id, record.name, *sorted(record.synonyms)):
        if not raw:
            continue
        match = _NSC.match(str(raw).strip())
        if not match:
            continue
        token = f"NSC{match.group(1)}"
        if token in seen:
            continue
        seen.add(token)
        found.append(token)
    return found


def _add_name(add, original: str, exact: str, stereo: str, salt: str, dose: str) -> None:
    methods = {"exact": exact, "stereo": stereo, "salt": salt, "dose": dose}
    for stage, query, confidence in name_variants(original):
        add("name", query, methods[stage], confidence)


def resolve_compounds(
    records: list[CompoundRecord],
    client: CompoundClient,
    limit: int | None = None,
) -> int:
    """Resolve unmatched compounds. Returns how many received a single CID."""
    from pgx.normalize.combinations import assign_review

    pending: list[tuple[CompoundRecord, list[tuple[str, str, str, str]]]] = []
    for record in records:
        if assign_review(record):
            continue
        if record.status in {"mapped", "ambiguous", "combination"}:
            continue
        queries = iter_queries(record)
        if queries:
            pending.append((record, queries))
        if limit is not None and len(pending) >= limit:
            break
    active = [(record, queries, 0) for record, queries in pending]
    mapped = 0
    while active:
        if hasattr(client, "prefetch"):
            client.prefetch([(queries[index][0], queries[index][1]) for _record, queries, index in active])
        nxt: list[tuple[CompoundRecord, list[tuple[str, str, str, str]], int]] = []
        for record, queries, index in active:
            kind, query, method, confidence = queries[index]
            cids = client.cids(kind, query)
            if len(cids) > 1:
                apply_pubchem(record, None, None, None, method, f"{query} matched {len(cids)} CIDs")
                record.status = "ambiguous"
                record.confidence = None
            elif len(cids) == 1:
                apply_pubchem(record, cids[0], record.inchikey, record.chembl_id, method, f"{method} {query}", confidence)
                if record.status == "mapped":
                    mapped += 1
            elif index + 1 < len(queries):
                nxt.append((record, queries, index + 1))
            else:
                apply_pubchem(record, None, None, None, "pubchem_miss", "no unique PubChem CID")
        active = nxt
    _resolve_chembl_names(records, client)
    _enrich(records, client)
    _label_structure_gaps(records)
    return mapped


def _resolve_chembl_names(records: list[CompoundRecord], client: CompoundClient) -> None:
    """Map GDSC names that PubChem left unmatched or ambiguous via one ChEMBL synonym.

    A query with more than one ChEMBL molecule stops. A molecule with no standard
    InChIKey is recorded as `no_public_structure` and is not given a CID.
    """
    lookup = getattr(client, "chembl_by_name", None)
    if lookup is None:
        return
    from pgx.normalize.combinations import chembl_name_candidates

    work = [
        record for record in records
        if record.source == "gdsc" and record.status in {"unmatched", "ambiguous"} and not record.review_reason
    ]
    planned: list[tuple[CompoundRecord, list[str]]] = []
    prefetch: list[str] = []
    for record in work:
        names = chembl_name_candidates(record)
        if names:
            planned.append((record, names))
            prefetch.extend(names)
    if hasattr(client, "prefetch_chembl_names"):
        client.prefetch_chembl_names(prefetch)
    for record, names in planned:
        if record.status not in {"unmatched", "ambiguous"} or record.review_reason:
            continue
        gap: tuple[str, str | None, str | None] | None = None
        started = record.status
        for name in names:
            molecules = lookup(name) or []
            if len(molecules) > 1:
                record.status = "ambiguous"
                record.match_method = "chembl_synonym"
                record.confidence = None
                record.pubchem_cid = None
                record.review_reason = None
                record.evidence = f"{name} matched {len(molecules)} ChEMBL molecules"
                gap = None
                break
            if len(molecules) != 1:
                continue
            molecule = molecules[0]
            chembl_id = molecule.get("chembl_id")
            inchikey = molecule.get("inchikey")
            if not inchikey:
                gap = (name, str(chembl_id) if chembl_id else None, None)
                continue
            cids = client.cids("inchikey", str(inchikey))
            if len(cids) > 1:
                record.status = "ambiguous"
                record.inchikey = str(inchikey)
                record.chembl_id = str(chembl_id) if chembl_id else record.chembl_id
                record.match_method = "chembl_synonym"
                record.confidence = None
                record.pubchem_cid = None
                record.review_reason = None
                record.evidence = f"{inchikey} matched {len(cids)} CIDs"
                gap = None
                break
            if len(cids) == 1:
                if record.status == "ambiguous":
                    record.status = "unmatched"
                apply_pubchem(
                    record, cids[0], str(inchikey), str(chembl_id) if chembl_id else None,
                    "chembl_synonym", f"chembl_synonym {name}", "medium",
                )
                record.review_reason = None
                gap = None
                break
            gap = (name, str(chembl_id) if chembl_id else None, str(inchikey))
        else:
            if started == "ambiguous" and record.status != "mapped":
                record.status = "ambiguous"
                record.review_reason = None
                record.pubchem_cid = None
                record.confidence = None
                record.evidence = record.evidence if "matched" in (record.evidence or "") and "ChEMBL" not in (record.evidence or "") else "ChEMBL did not select one PubChem CID"
                continue
            if gap and record.status == "unmatched":
                name, chembl_id, inchikey = gap
                record.status = "unmatched"
                record.review_reason = "no_public_structure"
                record.match_method = "chembl_synonym"
                record.confidence = None
                record.pubchem_cid = None
                record.chembl_id = chembl_id or record.chembl_id
                record.inchikey = inchikey or record.inchikey
                if inchikey:
                    record.evidence = f"{name} had InChIKey {inchikey} and no PubChem CID"
                else:
                    record.evidence = f"{name} matched {chembl_id or 'a ChEMBL molecule'} with no standard InChIKey"


def _label_structure_gaps(records: Iterable[CompoundRecord]) -> None:
    for record in records:
        if record.source != "gdsc" or record.status != "unmatched" or record.review_reason:
            continue
        record.review_reason = "no_public_structure"
        if record.match_method == "pubchem_miss":
            record.evidence = "no unique PubChem CID or ChEMBL structure"


def _enrich(records: Iterable[CompoundRecord], client: CompoundClient) -> None:
    wanted = [record.pubchem_cid for record in records if record.status == "mapped" and record.pubchem_cid]
    props = client.properties_for([cid for cid in wanted if cid]) if wanted else {}
    for record in records:
        if record.status != "mapped" or not record.pubchem_cid:
            continue
        row = props.get(record.pubchem_cid) or {}
        if row.get("inchikey") and not record.inchikey:
            record.inchikey = row["inchikey"]
    keys = [record.inchikey for record in records if record.status == "mapped" and record.inchikey and not record.chembl_id]
    chembl = client.chembl_for([key for key in keys if key]) if keys else {}
    for record in records:
        if record.status != "mapped" or not record.inchikey or record.chembl_id:
            continue
        found = chembl.get(record.inchikey) or []
        if len(found) == 1:
            record.chembl_id = found[0]
