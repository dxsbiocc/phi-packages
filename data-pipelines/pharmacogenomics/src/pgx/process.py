"""Build partitioned Parquet tables from the raw release files."""

from __future__ import annotations

import csv
import shutil
from pathlib import Path

import duckdb

from pgx.ids import (
    add_compound,
    build_model_xref,
    compound_rows,
    compounds_from_gdsc,
    compounds_from_prism_names,
)

GENE_LABEL = r"^(.*) \(([0-9]+)\)$"


def connect() -> duckdb.DuckDBPyConnection:
    return duckdb.connect()


def read_csv_dicts(path: Path) -> list[dict[str, str]]:
    con = connect()
    rows = con.execute(
        "SELECT * FROM read_csv(?, header=true, all_varchar=true, ignore_errors=false)",
        [str(path)],
    ).fetchall()
    columns = [item[0] for item in con.execute(
        "DESCRIBE SELECT * FROM read_csv(?, header=true, all_varchar=true)",
        [str(path)],
    ).fetchall()]
    return [dict(zip(columns, row, strict=True)) for row in rows]


def _write_rows(con: duckdb.DuckDBPyConnection, rows: list[dict[str, str | None]], dest: Path, columns: list[str]) -> None:
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True, exist_ok=True)
    csv_path = dest / "_incoming.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: "" if row.get(column) is None else row.get(column) for column in columns})
    dest_file = dest / "part.parquet"
    con.execute(
        f"COPY (SELECT * FROM read_csv('{csv_path.as_posix()}', header=true, all_varchar=true, nullstr='')) TO '{dest_file.as_posix()}' (FORMAT PARQUET)"
    )
    csv_path.unlink()


def _replace_tree(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def file_for(raw_root: Path, source: str, version: str, name: str) -> Path | None:
    path = raw_root / source / version / name
    return path if path.exists() else None


def specs_by_role(document: dict) -> dict[str, list[tuple[dict, dict]]]:
    found: dict[str, list[tuple[dict, dict]]] = {}
    for source in document.get("sources") or []:
        if not source.get("enabled", True):
            continue
        for spec in source.get("files") or []:
            found.setdefault(spec.get("role") or "", []).append((source, spec))
    return found


def _first_file(raw_root: Path, pairs: list[tuple[dict, dict]]) -> Path | None:
    for source, spec in pairs:
        path = file_for(raw_root, source["id"], str(source["version"]), spec["name"])
        if path:
            return path
    return None


def build_xrefs(document: dict, raw_root: Path, out: Path) -> dict[str, int]:
    roles = specs_by_role(document)
    cmp_path = _first_file(raw_root, roles.get("cmp_models", []))
    depmap_path = _first_file(raw_root, roles.get("depmap_models", []))
    if not cmp_path or not depmap_path:
        raise FileNotFoundError("model harmonization needs CMP model_list and DepMap Model.csv")
    gdsc_rows: list[dict[str, str]] = []
    for source, spec in roles.get("gdsc_dose_response", []):
        path = file_for(raw_root, source["id"], str(source["version"]), spec["name"])
        if not path:
            continue
        con = connect()
        gdsc_rows.extend(dict(zip(
            ["SANGER_MODEL_ID", "COSMIC_ID", "CELL_LINE_NAME"],
            row,
            strict=True,
        )) for row in con.execute(
            """
            SELECT DISTINCT SANGER_MODEL_ID, CAST(COSMIC_ID AS VARCHAR), CELL_LINE_NAME
            FROM read_csv(?, header=true, all_varchar=true)
            """,
            [str(path)],
        ).fetchall())
    prism_rows: list[dict[str, str]] = []
    for source, spec in roles.get("prism_cell_meta", []):
        path = file_for(raw_root, source["id"], str(source["version"]), spec["name"])
        if not path:
            continue
        con = connect()
        columns = [item[0] for item in con.execute("DESCRIBE SELECT * FROM read_csv(?, header=true, all_varchar=true)", [str(path)]).fetchall()]
        tissue = "primary_tissue" if "primary_tissue" in columns else "NULL"
        prism_rows.extend(dict(zip(
            ["depmap_id", "ccle_name", "primary_tissue"],
            row,
            strict=True,
        )) for row in con.execute(
            f"""
            SELECT DISTINCT depmap_id, ccle_name, {tissue}
            FROM read_csv(?, header=true, all_varchar=true)
            WHERE depmap_id IS NOT NULL
            """,
            [str(path)],
        ).fetchall())
    model_rows = build_model_xref(read_csv_dicts(cmp_path), read_csv_dicts(depmap_path), gdsc_rows, prism_rows)
    records = {}
    for source, spec in roles.get("gdsc_compounds", []):
        path = file_for(raw_root, source["id"], str(source["version"]), spec["name"])
        if path:
            for record in compounds_from_gdsc(read_csv_dicts(path)):
                add_compound(records, record)
    for role, id_field in (
        ("prism_primary_compounds", "broad_id"),
        ("prism_secondary_compounds", "broad_id"),
        ("prism_24q2_compounds", "IDs"),
    ):
        for source, spec in roles.get(role, []):
            path = file_for(raw_root, source["id"], str(source["version"]), spec["name"])
            if path:
                for record in compounds_from_prism_names(read_csv_dicts(path), source["id"], id_field):
                    add_compound(records, record)
    con = connect()
    model_columns = list(model_rows[0]) if model_rows else [
        "source", "source_model_id", "alias", "alias_type", "cvcl", "match_method", "status",
        "evidence", "display_name", "tissue", "disease", "species",
    ]
    compound_list = compound_rows(records)
    compound_columns = list(compound_list[0]) if compound_list else [
        "source", "source_compound_id", "name", "synonyms", "brd_id", "gdsc_drug_id",
        "pubchem_cid", "chembl_id", "inchikey", "smiles", "match_method", "status", "evidence",
    ]
    _write_rows(con, model_rows, out / "model_xref", model_columns)
    _write_rows(con, compound_list, out / "compound_xref", compound_columns)
    review_models = [row for row in model_rows if row["status"] != "mapped" and row["alias_type"] == "primary"]
    review_compounds = [row for row in compound_list if row["status"] != "mapped"]
    _write_rows(con, review_models, out / "review" / "models", model_columns)
    _write_rows(con, review_compounds, out / "review" / "compounds", compound_columns)
    releases = []
    for source in document.get("sources") or []:
        if not source.get("enabled", True):
            continue
        releases.append({
            "source": source["id"],
            "version": str(source["version"]),
            "license": source.get("license"),
            "citation": source.get("citation"),
            "homepage": source.get("homepage"),
            "batch": source.get("batch"),
        })
    _write_rows(con, releases, out / "releases", ["source", "version", "license", "citation", "homepage", "batch"])
    mapped = sum(1 for row in model_rows if row["alias_type"] == "primary" and row["status"] == "mapped")
    primary = sum(1 for row in model_rows if row["alias_type"] == "primary")
    return {"model_primary": primary, "model_mapped": mapped, "compounds": len(records)}


def _model_join(source_name: str) -> str:
    return f"""
    LEFT JOIN (
      SELECT source_model_id, any_value(cvcl) AS model_cvcl
      FROM read_parquet('{(Path("__MODEL_XREF__")).as_posix()}')
      WHERE source = '{source_name}' AND alias_type = 'primary' AND status = 'mapped'
      GROUP BY source_model_id
    ) models ON models.source_model_id = src.source_model_id
    """


def melt_gene_matrix(con: duckdb.DuckDBPyConnection, path: Path, dataset: str, metric: str, xref: Path, model_source: str, dest: Path) -> int:
    _replace_tree(dest)
    join = _model_join(model_source).replace(Path("__MODEL_XREF__").as_posix(), (xref / "part.parquet").as_posix())
    con.execute(
        f"""
        COPY (
          WITH wide AS (
            SELECT column0 AS source_model_id, * EXCLUDE (column0)
            FROM read_csv(?, header=true, all_varchar=true)
          ),
          long AS (
            SELECT source_model_id, gene_label, try_cast(nullif(value, 'NA') AS DOUBLE) AS value
            FROM (
              UNPIVOT wide
              ON COLUMNS(* EXCLUDE (source_model_id))
              INTO NAME gene_label VALUE value
            )
          ),
          src AS (
            SELECT
              ? AS dataset,
              source_model_id,
              regexp_extract(gene_label, '{GENE_LABEL}', 2) AS ncbi_gene_id,
              regexp_extract(gene_label, '{GENE_LABEL}', 1) AS symbol,
              value,
              ? AS metric
            FROM long
            WHERE value IS NOT NULL
              AND regexp_matches(gene_label, '{GENE_LABEL}')
          )
          SELECT src.dataset, src.source_model_id, models.model_cvcl, src.ncbi_gene_id, src.symbol, src.value, src.metric
          FROM src
          {join}
        ) TO '{dest.as_posix()}' (FORMAT PARQUET, PARTITION_BY (dataset))
        """,
        [str(path), dataset, metric],
    )
    count = con.execute("SELECT count(*) FROM read_parquet(?, hive_partitioning=true)", [str(dest / "**/*.parquet")]).fetchone()
    return int(count[0] if count else 0)


def load_mutations(con: duckdb.DuckDBPyConnection, path: Path, dataset: str, xref: Path, dest: Path) -> int:
    _replace_tree(dest)
    join = _model_join("depmap").replace(Path("__MODEL_XREF__").as_posix(), (xref / "part.parquet").as_posix())
    con.execute(
        f"""
        COPY (
          WITH src AS (
            SELECT
              ? AS dataset,
              ModelID AS source_model_id,
              nullif(CAST(try_cast(EntrezGeneID AS BIGINT) AS VARCHAR), '0') AS ncbi_gene_id,
              HugoSymbol AS symbol,
              Chrom AS chrom,
              CASE
                WHEN upper(Chrom) IN ('MT', 'M') THEN 'chrM'
                WHEN Chrom ILIKE 'chr%' THEN 'chr' || upper(substr(Chrom, 4))
                WHEN Chrom IS NULL OR Chrom = '' THEN NULL
                ELSE 'chr' || upper(Chrom)
              END AS chrom_norm,
              try_cast(Pos AS BIGINT) AS pos,
              Ref AS ref,
              Alt AS alt,
              ProteinChange AS protein_change,
              MolecularConsequence AS consequence,
              'GRCh38' AS genome_build,
              VepImpact AS vep_impact
            FROM read_csv(?, header=true, all_varchar=true)
            WHERE EntrezGeneID IS NOT NULL AND EntrezGeneID NOT IN ('', 'NA')
          )
          SELECT src.*, models.model_cvcl
          FROM src
          {join}
        ) TO '{dest.as_posix()}' (FORMAT PARQUET, PARTITION_BY (dataset))
        """,
        [dataset, str(path)],
    )
    count = con.execute("SELECT count(*) FROM read_parquet(?, hive_partitioning=true)", [str(dest / "**/*.parquet")]).fetchone()
    return int(count[0] if count else 0)


def load_gdsc_response(con: duckdb.DuckDBPyConnection, path: Path, xref: Path, dest_root: Path) -> int:
    join_cvcl = f"""
      LEFT JOIN (
        SELECT source_model_id, any_value(cvcl) AS model_cvcl
        FROM read_parquet('{(xref / "part.parquet").as_posix()}')
        WHERE source = 'gdsc' AND alias_type = 'primary' AND status = 'mapped'
        GROUP BY source_model_id
      ) models ON models.source_model_id = base.source_model_id
    """
    con.execute(
        f"""
        CREATE OR REPLACE TEMP TABLE gdsc_long AS
        WITH base AS (
          SELECT
            lower(DATASET) AS dataset,
            coalesce(nullif(SANGER_MODEL_ID, ''), 'COSMIC:' || CAST(COSMIC_ID AS VARCHAR), CELL_LINE_NAME) AS source_model_id,
            CAST(DRUG_ID AS VARCHAR) AS source_compound_id,
            DRUG_NAME AS compound_name,
            metric,
            try_cast(metric_value AS DOUBLE) AS value,
            unit,
            NULL::DOUBLE AS dose,
            NULL::VARCHAR AS screen_id,
            NULL::VARCHAR AS qc
          FROM (
            SELECT *, 'AUC' AS metric, CAST(AUC AS VARCHAR) AS metric_value, NULL AS unit FROM read_csv(?, header=true, all_varchar=true)
            UNION ALL
            SELECT *, 'LN_IC50' AS metric, CAST(LN_IC50 AS VARCHAR) AS metric_value, 'ln_uM' AS unit FROM read_csv(?, header=true, all_varchar=true)
          )
        )
        SELECT base.dataset, base.source_model_id, models.model_cvcl, base.source_compound_id,
               NULL::VARCHAR AS pubchem_cid, base.compound_name, base.metric, base.value, base.unit,
               base.dose, base.screen_id, base.qc
        FROM base
        {join_cvcl}
        WHERE base.value IS NOT NULL
        """,
        [str(path), str(path)],
    )
    datasets = [row[0] for row in con.execute("SELECT DISTINCT dataset FROM gdsc_long").fetchall()]
    for dataset in datasets:
        dest = dest_root / f"dataset={dataset}"
        _replace_tree(dest)
        con.execute(
            f"""
            COPY (SELECT * EXCLUDE (dataset) FROM gdsc_long WHERE dataset = ?)
            TO '{dest.as_posix()}/part.parquet' (FORMAT PARQUET)
            """,
            [dataset],
        )
    count = con.execute("SELECT count(*) FROM gdsc_long").fetchone()
    return int(count[0] if count else 0)


def _copy_dataset(con: duckdb.DuckDBPyConnection, sql: str, params: list[object], root: Path, dataset: str) -> int:
    part = root / f"dataset={dataset}"
    _replace_tree(part)
    con.execute(f"COPY ({sql}) TO '{part.as_posix()}/part.parquet' (FORMAT PARQUET)", params)
    count = con.execute("SELECT count(*) FROM read_parquet(?)", [str(part / "part.parquet")]).fetchone()
    return int(count[0] if count else 0)


def load_secondary_curve(con: duckdb.DuckDBPyConnection, path: Path, dataset: str, xref: Path, dest: Path) -> int:
    xref_file = (xref / "part.parquet").as_posix()
    return _copy_dataset(
        con,
        f"""
          WITH base AS (
            SELECT
              depmap_id AS source_model_id,
              broad_id AS source_compound_id,
              name AS compound_name,
              metric,
              try_cast(metric_value AS DOUBLE) AS value,
              screen_id
            FROM (
              SELECT *, 'AUC' AS metric, CAST(auc AS VARCHAR) AS metric_value FROM read_csv(?, header=true, all_varchar=true)
              UNION ALL
              SELECT *, 'IC50' AS metric, CAST(ic50 AS VARCHAR) AS metric_value FROM read_csv(?, header=true, all_varchar=true)
              UNION ALL
              SELECT *, 'EC50' AS metric, CAST(ec50 AS VARCHAR) AS metric_value FROM read_csv(?, header=true, all_varchar=true)
            )
          )
          SELECT
            base.source_model_id,
            models.model_cvcl,
            base.source_compound_id,
            NULL::VARCHAR AS pubchem_cid,
            base.compound_name,
            base.metric,
            base.value,
            NULL::VARCHAR AS unit,
            NULL::DOUBLE AS dose,
            base.screen_id,
            NULL::VARCHAR AS qc
          FROM base
          LEFT JOIN (
            SELECT source_model_id, any_value(cvcl) AS model_cvcl
            FROM read_parquet('{xref_file}')
            WHERE source = 'depmap' AND alias_type = 'primary' AND status = 'mapped'
            GROUP BY source_model_id
          ) models ON models.source_model_id = base.source_model_id
          WHERE base.value IS NOT NULL
        """,
        [str(path), str(path), str(path)],
        dest,
        dataset,
    )


def load_prism24_lfc(con: duckdb.DuckDBPyConnection, path: Path, dataset: str, xref: Path, dest: Path) -> int:
    xref_file = (xref / "part.parquet").as_posix()
    return _copy_dataset(
        con,
        f"""
          WITH base AS (
            SELECT
              regexp_extract(row_id, '^(ACH-[0-9]+)', 1) AS source_model_id,
              regexp_extract(broad_id, '^(BRD-[A-Z0-9-]+)', 1) AS source_compound_id,
              try_cast(dose AS DOUBLE) AS dose,
              screen AS screen_id,
              CASE WHEN broad_id ILIKE '%QC Failure%' THEN 'qc_failure' ELSE NULL END AS qc,
              try_cast(LFC AS DOUBLE) AS value
            FROM read_csv(?, header=true, all_varchar=true)
          )
          SELECT
            base.source_model_id,
            models.model_cvcl,
            base.source_compound_id,
            NULL::VARCHAR AS pubchem_cid,
            NULL::VARCHAR AS compound_name,
            'LFC' AS metric,
            base.value,
            NULL::VARCHAR AS unit,
            base.dose,
            base.screen_id,
            base.qc
          FROM base
          LEFT JOIN (
            SELECT source_model_id, any_value(cvcl) AS model_cvcl
            FROM read_parquet('{xref_file}')
            WHERE source = 'depmap' AND alias_type = 'primary' AND status = 'mapped'
            GROUP BY source_model_id
          ) models ON models.source_model_id = base.source_model_id
          WHERE base.value IS NOT NULL AND base.source_compound_id IS NOT NULL
        """,
        [str(path)],
        dest,
        dataset,
    )


def load_primary_lfc(con: duckdb.DuckDBPyConnection, path: Path, dataset: str, xref: Path, dest: Path) -> int:
    xref_file = (xref / "part.parquet").as_posix()
    return _copy_dataset(
        con,
        f"""
          WITH wide AS (
            SELECT column0 AS source_model_id, * EXCLUDE (column0)
            FROM read_csv(?, header=true, all_varchar=true)
          ),
          long AS (
            SELECT source_model_id, column_name, try_cast(nullif(value, 'NA') AS DOUBLE) AS value
            FROM (
              UNPIVOT wide
              ON COLUMNS(* EXCLUDE (source_model_id))
              INTO NAME column_name VALUE value
            )
          ),
          src AS (
            SELECT
              source_model_id,
              regexp_extract(column_name, '^(BRD-[A-Z0-9-]+)::([0-9.]+)::(.+)$', 1) AS source_compound_id,
              try_cast(regexp_extract(column_name, '^(BRD-[A-Z0-9-]+)::([0-9.]+)::(.+)$', 2) AS DOUBLE) AS dose,
              regexp_extract(column_name, '^(BRD-[A-Z0-9-]+)::([0-9.]+)::(.+)$', 3) AS screen_id,
              value
            FROM long
            WHERE value IS NOT NULL
              AND regexp_matches(column_name, '^(BRD-[A-Z0-9-]+)::')
          )
          SELECT src.source_model_id, models.model_cvcl, src.source_compound_id,
                 NULL::VARCHAR AS pubchem_cid, NULL::VARCHAR AS compound_name, 'LFC' AS metric,
                 src.value, NULL::VARCHAR AS unit, src.dose, src.screen_id, NULL::VARCHAR AS qc
          FROM src
          LEFT JOIN (
            SELECT source_model_id, any_value(cvcl) AS model_cvcl
            FROM read_parquet('{xref_file}')
            WHERE source = 'depmap' AND alias_type = 'primary' AND status = 'mapped'
            GROUP BY source_model_id
          ) models ON models.source_model_id = src.source_model_id
        """,
        [str(path)],
        dest,
        dataset,
    )


def attach_pubchem(con: duckdb.DuckDBPyConnection, response_root: Path, xref: Path) -> None:
    """Fill pubchem_cid on drug_response partitions from a mapped compound_xref."""
    if not response_root.exists():
        return
    xref_file = (xref / "part.parquet").as_posix()
    for part in response_root.glob("dataset=*/part.parquet"):
        temporary = part.with_suffix(".parquet.next")
        con.execute(
            f"""
            COPY (
              SELECT response.* REPLACE (coalesce(response.pubchem_cid, compounds.pubchem_cid) AS pubchem_cid)
              FROM read_parquet('{part.as_posix()}') response
              LEFT JOIN (
                SELECT source_compound_id, any_value(pubchem_cid) AS pubchem_cid
                FROM read_parquet('{xref_file}')
                WHERE status = 'mapped' AND pubchem_cid IS NOT NULL
                GROUP BY source_compound_id
                HAVING count(DISTINCT pubchem_cid) = 1
              ) compounds ON compounds.source_compound_id = response.source_compound_id
            ) TO '{temporary.as_posix()}' (FORMAT PARQUET)
            """
        )
        temporary.replace(part)


def process_all(document: dict, raw_root: Path, out: Path) -> dict[str, int]:
    out.mkdir(parents=True, exist_ok=True)
    counts = build_xrefs(document, raw_root, out)
    roles = specs_by_role(document)
    for name in ("gene_effect", "expression", "cnv", "mutation", "drug_response"):
        existing = out / name
        if existing.exists():
            shutil.rmtree(existing)
    con = connect()
    xref = out / "model_xref"
    matrix_roles = {
        "gene_matrix": melt_gene_matrix,
    }
    for source, spec in roles.get("gene_matrix", []):
        path = file_for(raw_root, source["id"], str(source["version"]), spec["name"])
        if not path:
            continue
        dataset = str(spec["dataset"])
        table = {
            "chronos_gene_effect": "gene_effect",
            "log2_tpm_plus_1": "expression",
            "relative_copy_number": "cnv",
        }.get(str(spec["metric"]), str(spec["metric"]))
        counts[table] = matrix_roles["gene_matrix"](con, path, dataset, str(spec["metric"]), xref, "depmap", out / table)
    for source, spec in roles.get("depmap_mutations", []):
        path = file_for(raw_root, source["id"], str(source["version"]), spec["name"])
        if path:
            counts["mutation"] = load_mutations(con, path, str(spec.get("dataset") or "depmap_mutations_24q2"), xref, out / "mutation")
    response_root = out / "drug_response"
    response_root.mkdir(parents=True, exist_ok=True)
    total_response = 0
    for source, spec in roles.get("gdsc_dose_response", []):
        path = file_for(raw_root, source["id"], str(source["version"]), spec["name"])
        if path:
            total_response += load_gdsc_response(con, path, xref, response_root)
    for source, spec in roles.get("prism_secondary_curve", []):
        path = file_for(raw_root, source["id"], str(source["version"]), spec["name"])
        if path:
            total_response += load_secondary_curve(con, path, str(spec["dataset"]), xref, response_root)
    for source, spec in roles.get("prism_24q2_lfc", []):
        path = file_for(raw_root, source["id"], str(source["version"]), spec["name"])
        if path:
            total_response += load_prism24_lfc(con, path, str(spec["dataset"]), xref, response_root)
    for source, spec in roles.get("prism_primary_lfc", []):
        path = file_for(raw_root, source["id"], str(source["version"]), spec["name"])
        if path:
            total_response += load_primary_lfc(con, path, str(spec["dataset"]), xref, response_root)
    if total_response:
        attach_pubchem(con, response_root, out / "compound_xref")
        counts["drug_response"] = total_response
    return counts
