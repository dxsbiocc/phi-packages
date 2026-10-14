"""phi-pgx download | process | validate | resolve-compounds"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from pgx.download import download_all, live_enabled, load_sources
from pgx.ids import CompoundRecord, add_compound, compound_rows
from pgx.process import _write_rows, connect, process_all
from pgx.pubchem import resolve_records
from pgx.validate import validate


def _root(args: argparse.Namespace) -> Path:
    return Path(args.root).resolve()


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Pharmacogenomics release pipeline")
    parser.add_argument("--root", default=".", help="Pipeline directory (sources.yaml, raw/, processed/)")
    parser.add_argument("--sources", default=None, help="Path to sources.yaml")
    sub = parser.add_subparsers(dest="command", required=True)

    download = sub.add_parser("download")
    download.add_argument("--live", action="store_true", help="Also download multi-hundred-MB matrices")
    download.add_argument("--resolve-urls", action="store_true", help="Ignore configured URLs and run resolvers")
    download.add_argument("--only", action="append", default=[], help="Source id to download")

    sub.add_parser("process")
    sub.add_parser("validate")

    resolve = sub.add_parser("resolve-compounds")
    resolve.add_argument("--limit", type=int, default=None, help="Maximum unmatched names to send to PubChem")

    args = parser.parse_args(argv)
    root = _root(args)
    sources_path = Path(args.sources) if args.sources else root / "sources.yaml"
    document = load_sources(sources_path)
    raw_root = root / "raw"
    out = root / "processed"

    if args.command == "download":
        results = download_all(document, raw_root, live_enabled(args.live), args.resolve_urls, set(args.only) or None)
        for record in results:
            print(f"{record['status']:10} {record['name']} {record['sha256']}")
        return
    if args.command == "process":
        counts = process_all(document, raw_root, out)
        print(json.dumps(counts, indent=2))
        return
    if args.command == "validate":
        report = validate(out)
        print(json.dumps(report, indent=2))
        return
    if args.command == "resolve-compounds":
        con = connect()
        fetched = con.execute("SELECT * FROM read_parquet(?)", [str(out / "compound_xref" / "part.parquet")]).fetchall()
        columns = [item[0] for item in con.execute("DESCRIBE SELECT * FROM read_parquet(?)", [str(out / "compound_xref" / "part.parquet")]).fetchall()]
        records: dict[tuple[str, str], CompoundRecord] = {}
        for fetched_row in fetched:
            data = dict(zip(columns, fetched_row, strict=True))
            record = CompoundRecord(
                source=data["source"],
                source_compound_id=data["source_compound_id"],
                name=data.get("name"),
                brd_id=data.get("brd_id"),
                gdsc_drug_id=data.get("gdsc_drug_id"),
                pubchem_cid=data.get("pubchem_cid"),
                chembl_id=data.get("chembl_id"),
                inchikey=data.get("inchikey"),
                smiles=data.get("smiles"),
                match_method=data.get("match_method") or "unmatched",
                confidence=data.get("confidence"),
                status=data.get("status") or "unmatched",
                evidence=data.get("evidence") or "",
            )
            if data.get("synonyms"):
                record.synonyms.update(part for part in str(data["synonyms"]).split(";") if part)
            add_compound(records, record)
        updated = resolve_records(list(records.values()), raw_root / "pubchem", limit=args.limit)
        _write_rows(connect(), compound_rows(records), out / "compound_xref", [
            "source", "source_compound_id", "name", "synonyms", "brd_id", "gdsc_drug_id",
            "pubchem_cid", "chembl_id", "inchikey", "smiles", "match_method", "confidence", "status", "evidence",
        ])
        summary: dict[str, dict[str, int]] = {}
        for record in records.values():
            bucket = summary.setdefault(record.source, {"rows": 0, "mapped": 0, "unmatched": 0, "ambiguous": 0})
            bucket["rows"] += 1
            bucket[record.status] = bucket.get(record.status, 0) + 1
        print(json.dumps({"resolved": updated, "compounds": summary}, indent=2))


if __name__ == "__main__":
    main()
