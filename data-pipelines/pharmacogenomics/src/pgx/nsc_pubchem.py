"""Bulk NSC to PubChem CID mapping from the DTP/NCI depositor.

``SID-Map.gz`` rows are ``SID``, source name, source registry id, and the
standardized CID when PubChem has one. Only source ``DTP/NCI`` is used.
Several CIDs for one NSC stay ambiguous. NSCs with no CID fall through to a
synonym file or to concurrent ``NSC####`` name lookups. A multi-name PUG body
is not used: PubChem returns one CID list and drops the name-to-CID pairing.
"""

from __future__ import annotations

import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path

from pgx.ids import CompoundRecord, apply_pubchem
from pgx.nci60 import nsc_id

DTP_SOURCE = "DTP/NCI"
CROSSWALK_NAME = "dtp_nci.tsv"
_NEEDLE = b"DTP/NCI"


def crosswalk_path(cache_dir: Path) -> Path:
    return cache_dir / CROSSWALK_NAME


def nsc_digits(value: object) -> str | None:
    """NSC number as digits, without leading zeros.

    Depositor registry ids may be ``740``, ``NSC740``, ``NSC 740``, or
    ``NSC-740``. The assay files store the unpadded ``NSC{digits}`` id.
    """
    token = nsc_id(value)
    if token:
        return str(int(token[3:]))
    text = "" if value is None else str(value).strip()
    if not text:
        return None
    folded = text.upper().replace(" ", "").replace("-", "").replace("_", "")
    if folded.startswith("NSC") and folded[3:].isdigit():
        return str(int(folded[3:]))
    return None


def _open_binary(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return path.open("rb")


def scan_sid_map(
    path: Path,
    dest: Path | None = None,
    stats: dict | None = None,
) -> dict[str, tuple[set[str], set[str]]]:
    """Return ``NSC digits -> (sids, cids)`` for DTP/NCI substances.

    When ``dest`` is set, the filtered rows are written as
    ``nsc, sid, cid`` so later runs do not rescan the substance map.
    ``stats``, when passed, receives row counts and a sample of registry
    ids that were not NSC numbers.
    """
    sids: dict[str, set[str]] = defaultdict(set)
    cids: dict[str, set[str]] = defaultdict(set)
    writer = None
    handle_out = None
    skipped: list[str] = []
    written = 0
    seen = 0
    if dest is not None:
        dest.parent.mkdir(parents=True, exist_ok=True)
        handle_out = dest.open("w", newline="")
        writer = csv.writer(handle_out)
        writer.writerow(["nsc", "sid", "cid"])
    try:
        with _open_binary(path) as handle:
            for raw in handle:
                seen += 1
                if seen % 5_000_000 == 0:
                    print(f"sid-map lines {seen}", file=sys.stderr, flush=True)
                if _NEEDLE not in raw:
                    continue
                parts = raw.decode("utf-8", "replace").rstrip("\r\n").split("\t")
                if len(parts) < 3 or parts[1] != DTP_SOURCE:
                    continue
                digits = nsc_digits(parts[2])
                if not digits:
                    if len(skipped) < 40:
                        skipped.append(parts[2])
                    continue
                sid = parts[0].strip()
                cid = parts[3].strip() if len(parts) > 3 else ""
                if not cid.isdigit():
                    cid = ""
                if sid:
                    sids[digits].add(sid)
                if cid:
                    cids[digits].add(cid)
                if writer is not None:
                    writer.writerow([digits, sid, cid])
                    written += 1
    finally:
        if handle_out is not None:
            handle_out.close()
    if stats is not None:
        stats["lines"] = seen
        stats["rows"] = written
        stats["nscs"] = len(sids)
        stats["with_cid"] = sum(1 for digits in sids if cids[digits])
        stats["skipped_registry"] = skipped
    return {digits: (sids[digits], cids[digits]) for digits in sids}


def load_crosswalk(path: Path) -> dict[str, tuple[set[str], set[str]]]:
    sids: dict[str, set[str]] = defaultdict(set)
    cids: dict[str, set[str]] = defaultdict(set)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            digits = nsc_digits(row.get("nsc"))
            if not digits:
                continue
            sid = (row.get("sid") or "").strip()
            cid = (row.get("cid") or "").strip()
            if sid:
                sids[digits].add(sid)
            if cid.isdigit():
                cids[digits].add(cid)
    return {digits: (sids[digits], cids[digits]) for digits in sids}


def load_or_scan(sid_map: Path, dest: Path) -> dict[str, tuple[set[str], set[str]]]:
    if dest.exists() and dest.stat().st_size > 0:
        return load_crosswalk(dest)
    return scan_sid_map(sid_map, dest)


def synonym_cids(path: Path, wanted: set[str]) -> dict[str, set[str]]:
    """CID-synonym rows whose synonym is an NSC number in ``wanted`` (digits)."""
    found: dict[str, set[str]] = {digits: set() for digits in wanted}
    if not wanted:
        return found
    with _open_binary(path) as handle:
        for raw in handle:
            if b"NSC" not in raw and b"nsc" not in raw:
                continue
            cid, separator, synonym = raw.decode("utf-8", "replace").rstrip("\r\n").partition("\t")
            if not separator or not cid.strip().isdigit():
                continue
            digits = nsc_digits(synonym.strip().strip('"'))
            if not digits:
                continue
            bucket = found.get(digits)
            if bucket is not None:
                bucket.add(cid.strip())
    return found


def _assign(record: CompoundRecord, cids: set[str], method: str, evidence: str) -> None:
    if len(cids) == 1:
        apply_pubchem(record, next(iter(cids)), record.inchikey, record.chembl_id, method, evidence, "high")
        return
    if len(cids) > 1:
        record.status = "ambiguous"
        record.match_method = method
        record.confidence = None
        record.pubchem_cid = None
        record.evidence = evidence


def apply_deposits(
    records: list[CompoundRecord],
    deposits: dict[str, tuple[set[str], set[str]]],
) -> int:
    """Map NCI-60 records from a DTP/NCI deposit index. Returns how many changed."""
    changed = 0
    for record in records:
        if record.source != "nci60" or record.status in {"mapped", "ambiguous", "combination"}:
            continue
        digits = nsc_digits(record.source_compound_id)
        if not digits:
            continue
        found = deposits.get(digits)
        if not found:
            continue
        _sids, cids = found
        if not cids:
            continue
        before = (record.status, record.pubchem_cid)
        evidence_ids = ",".join(sorted(cids))
        if len(cids) == 1:
            _assign(record, cids, "pubchem_dtp_nci", f"DTP/NCI standardized CID {evidence_ids}")
        else:
            _assign(record, cids, "pubchem_dtp_nci", f"DTP/NCI standardized CIDs {evidence_ids}")
        if (record.status, record.pubchem_cid) != before:
            changed += 1
    return changed


def apply_synonyms(records: list[CompoundRecord], found: dict[str, set[str]]) -> int:
    """Synonym fallback for NSCs the depositor map did not standardize."""
    changed = 0
    for record in records:
        if record.source != "nci60" or record.status in {"mapped", "ambiguous", "combination"}:
            continue
        digits = nsc_digits(record.source_compound_id)
        if not digits:
            continue
        cids = found.get(digits) or set()
        if not cids:
            continue
        before = (record.status, record.pubchem_cid)
        evidence_ids = ",".join(sorted(cids))
        label = f"NSC{digits}"
        if len(cids) == 1:
            _assign(record, cids, "pubchem_nsc", f"PubChem synonym {label} CID {evidence_ids}")
        else:
            _assign(record, cids, "pubchem_nsc", f"PubChem synonym {label} CIDs {evidence_ids}")
        if (record.status, record.pubchem_cid) != before:
            changed += 1
    return changed


def summarize(records: list[CompoundRecord]) -> dict[str, int]:
    summary = {"rows": 0, "mapped": 0, "ambiguous": 0, "unmatched": 0, "pubchem_dtp_nci": 0, "pubchem_nsc": 0}
    for record in records:
        if record.source != "nci60":
            continue
        summary["rows"] += 1
        summary[record.status] = summary.get(record.status, 0) + 1
        if record.match_method in {"pubchem_dtp_nci", "pubchem_nsc"}:
            summary[record.match_method] += 1
    return summary
