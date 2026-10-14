"""Schema checks, mapping coverage, and the Nutlin-3a / TP53 sanity check."""

from __future__ import annotations

import json
from pathlib import Path

import duckdb

from pgx.ids import SILENT_CONSEQUENCES

REQUIRED = {
    "model_xref": {"source", "source_model_id", "alias", "alias_type", "cvcl", "match_method", "status"},
    "compound_xref": {"source", "source_compound_id", "pubchem_cid", "chembl_id", "inchikey", "status"},
    "releases": {"source", "version", "license", "citation"},
}


class ValidationError(RuntimeError):
    pass


def _connect(out: Path) -> duckdb.DuckDBPyConnection:
    return duckdb.connect()


def _columns(con: duckdb.DuckDBPyConnection, path: Path) -> set[str]:
    described = con.execute("DESCRIBE SELECT * FROM read_parquet(?)", [str(path)]).fetchall()
    return {row[0] for row in described}


def _one(con: duckdb.DuckDBPyConnection, sql: str, params: list[object] | None = None) -> int:
    row = con.execute(sql, params or []).fetchone()
    return int(row[0] if row and row[0] is not None else 0)


def coverage(out: Path) -> dict:
    con = _connect(out)
    report: dict = {"models": {}, "compounds": {}, "tables": {}}
    xref = out / "model_xref" / "part.parquet"
    if not xref.exists():
        raise ValidationError("model_xref parquet is missing")
    missing = REQUIRED["model_xref"] - _columns(con, xref)
    if missing:
        raise ValidationError(f"model_xref is missing columns {sorted(missing)}")
    for source, total, mapped, unmatched, ambiguous in con.execute(
        """
        SELECT source,
               count(*),
               count(*) FILTER (WHERE status = 'mapped'),
               count(*) FILTER (WHERE status = 'unmatched'),
               count(*) FILTER (WHERE status = 'ambiguous')
        FROM read_parquet(?)
        WHERE alias_type = 'primary'
        GROUP BY source
        ORDER BY source
        """,
        [str(xref)],
    ).fetchall():
        report["models"][source] = {
            "primary": int(total),
            "mapped": int(mapped),
            "unmatched": int(unmatched),
            "ambiguous": int(ambiguous),
        }
    if "cellosaurus_flag" in _columns(con, xref):
        report["cellosaurus_flags"] = {
            flag or "": int(count)
            for flag, count in con.execute(
                """
                SELECT coalesce(cellosaurus_flag, ''), count(*)
                FROM read_parquet(?)
                WHERE alias_type = 'primary' AND status = 'mapped'
                GROUP BY 1
                ORDER BY 1
                """,
                [str(xref)],
            ).fetchall()
        }
    gene_xref = out / "gene_xref" / "part.parquet"
    if gene_xref.exists():
        report["genes"] = {
            status: int(count)
            for status, count in con.execute(
                "SELECT status, count(*) FROM read_parquet(?) GROUP BY 1 ORDER BY 1",
                [str(gene_xref)],
            ).fetchall()
        }
    compounds = out / "compound_xref" / "part.parquet"
    if compounds.exists():
        missing = REQUIRED["compound_xref"] - _columns(con, compounds)
        if missing:
            raise ValidationError(f"compound_xref is missing columns {sorted(missing)}")
        for source, total, mapped, unmatched, ambiguous in con.execute(
            """
            SELECT source,
                   count(*),
                   count(*) FILTER (WHERE status = 'mapped'),
                   count(*) FILTER (WHERE status = 'unmatched'),
                   count(*) FILTER (WHERE status = 'ambiguous')
            FROM read_parquet(?)
            GROUP BY source
            ORDER BY source
            """,
            [str(compounds)],
        ).fetchall():
            report["compounds"][source] = {
                "rows": int(total),
                "mapped": int(mapped),
                "unmatched": int(unmatched),
                "ambiguous": int(ambiguous),
            }
    for table in ("gene_effect", "expression", "mutation", "cnv", "drug_response"):
        root = out / table
        parts = list(root.glob("dataset=*/**/*.parquet")) if root.exists() else []
        if not parts:
            report["tables"][table] = {"rows": 0, "with_cvcl": 0}
            continue
        rows = _one(con, "SELECT count(*) FROM read_parquet(?, hive_partitioning=true)", [str(root / "**/*.parquet")])
        with_cvcl = _one(
            con,
            "SELECT count(*) FROM read_parquet(?, hive_partitioning=true) WHERE model_cvcl IS NOT NULL",
            [str(root / "**/*.parquet")],
        )
        report["tables"][table] = {"rows": rows, "with_cvcl": with_cvcl}
    releases = out / "releases" / "part.parquet"
    if not releases.exists():
        raise ValidationError("releases parquet is missing")
    missing = REQUIRED["releases"] - _columns(con, releases)
    if missing:
        raise ValidationError(f"releases is missing columns {sorted(missing)}")
    report["releases"] = con.execute(
        "SELECT source, version, license, citation FROM read_parquet(?) ORDER BY source",
        [str(releases)],
    ).fetchall()
    report["releases"] = [
        {"source": source, "version": version, "license": license, "citation": citation}
        for source, version, license, citation in report["releases"]
    ]
    return report


def nutlin_tp53(out: Path, min_group: int = 3) -> dict:
    """Nutlin-3a should be more effective in TP53 wild-type lines.

    GDSC publishes LN_IC50 and AUC where a lower value means a stronger response.
    A model is TP53-mutant when the mutation table has a non-silent ncbigene:7157
    call. Models absent from the mutation table are excluded.
    """
    response = out / "drug_response"
    mutations = out / "mutation"
    compounds = out / "compound_xref" / "part.parquet"
    missing = [
        name for name, present in (
            ("drug_response", any(response.glob("dataset=*/*.parquet")) if response.exists() else False),
            ("mutation", any(mutations.glob("dataset=*/*.parquet")) if mutations.exists() else False),
            ("compound_xref", compounds.exists()),
        ) if not present
    ]
    if missing:
        return {"status": "skipped", "reason": f"missing {', '.join(missing)}"}
    con = _connect(out)
    result = {"status": "skipped", "reason": "no Nutlin rows", "checks": []}
    frame = con.execute(
        """
        WITH nutlin AS (
          SELECT source_compound_id
          FROM read_parquet(?)
          WHERE source = 'gdsc'
            AND (
              source_compound_id = '1047'
              OR lower(name) LIKE '%nutlin-3a%'
            )
        ),
        response AS (
          SELECT model_cvcl, metric, value
          FROM read_parquet(?, hive_partitioning=true)
          WHERE dataset IN ('gdsc1', 'gdsc2')
            AND source_compound_id IN (SELECT source_compound_id FROM nutlin)
            AND model_cvcl IS NOT NULL
            AND metric IN ('LN_IC50', 'AUC')
            AND value IS NOT NULL
        ),
        mutated AS (
          SELECT DISTINCT model_cvcl
          FROM read_parquet(?, hive_partitioning=true)
          WHERE ncbi_gene_id = '7157'
            AND model_cvcl IS NOT NULL
            AND lower(coalesce(consequence, '')) NOT IN (SELECT unnest(?::VARCHAR[]))
        ),
        sequenced AS (
          SELECT DISTINCT model_cvcl
          FROM read_parquet(?, hive_partitioning=true)
          WHERE model_cvcl IS NOT NULL
        )
        SELECT response.metric,
               response.model_cvcl IN (SELECT model_cvcl FROM mutated) AS mutant,
               median(response.value) AS median_value,
               count(DISTINCT response.model_cvcl) AS models
        FROM response
        WHERE response.model_cvcl IN (SELECT model_cvcl FROM sequenced)
        GROUP BY 1, 2
        ORDER BY 1, 2
        """,
        [
            str(compounds),
            str(response / "**/*.parquet"),
            str(mutations / "**/*.parquet"),
            list(SILENT_CONSEQUENCES),
            str(mutations / "**/*.parquet"),
        ],
    ).fetchall()
    grouped: dict[str, dict[str, tuple[float, int]]] = {}
    for metric, mutant, median_value, models in frame:
        grouped.setdefault(metric, {})["mutant" if mutant else "wild_type"] = (float(median_value), int(models))
    if not grouped:
        return result
    failures = []
    checks = []
    for metric, groups in grouped.items():
        wild = groups.get("wild_type")
        mutant = groups.get("mutant")
        if not wild or not mutant or wild[1] < min_group or mutant[1] < min_group:
            checks.append({"metric": metric, "status": "skipped", "groups": groups})
            continue
        ok = wild[0] < mutant[0]
        checks.append({
            "metric": metric,
            "status": "pass" if ok else "fail",
            "wild_type_median": wild[0],
            "wild_type_models": wild[1],
            "mutant_median": mutant[0],
            "mutant_models": mutant[1],
        })
        if not ok:
            failures.append(metric)
    evaluated = [item for item in checks if item["status"] != "skipped"]
    if not evaluated:
        return {"status": "skipped", "reason": f"fewer than {min_group} models in a TP53 group", "checks": checks}
    return {"status": "fail" if failures else "pass", "checks": checks}


def validate(out: Path, report_dir: Path | None = None) -> dict:
    report = {"coverage": coverage(out), "nutlin_tp53": nutlin_tp53(out)}
    destination = report_dir or (out / "validation")
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    if report["nutlin_tp53"]["status"] == "fail":
        raise ValidationError("Nutlin-3a was not more effective in TP53 wild-type lines")
    return report
