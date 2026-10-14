"""Build partitioned Parquet tables from the raw release files."""

from __future__ import annotations

import csv
import shutil
from pathlib import Path

import duckdb

from pgx.ids import (
    CompoundRecord,
    add_compound,
    build_model_xref,
    compound_rows,
    compounds_from_gdsc,
    compounds_from_prism_names,
)
from pgx.normalize.cell_lines import CellosaurusIndex
from pgx.normalize.genes import HgncIndex, header_labels, resolve_labels

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
    cellosaurus = None
    cellosaurus_path = _first_file(raw_root, roles.get("cellosaurus", []))
    if cellosaurus_path:
        cellosaurus = CellosaurusIndex.load(cellosaurus_path)
    model_rows = build_model_xref(
        read_csv_dicts(cmp_path), read_csv_dicts(depmap_path), gdsc_rows, prism_rows, cellosaurus=cellosaurus,
    )
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
    _add_batch2_xrefs(raw_root, roles, document, cellosaurus, model_rows, records)
    con = connect()
    model_columns = list(model_rows[0]) if model_rows else [
        "source", "source_model_id", "alias", "alias_type", "cvcl", "match_method", "status",
        "evidence", "display_name", "tissue", "disease", "species", "cellosaurus_flag",
    ]
    compound_list = compound_rows(records)
    compound_columns = list(compound_list[0]) if compound_list else [
        "source", "source_compound_id", "name", "synonyms", "brd_id", "gdsc_drug_id",
        "pubchem_cid", "chembl_id", "inchikey", "smiles", "match_method", "confidence", "status", "evidence",
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


GENE_COLUMNS = [
    "source_label", "source_symbol", "symbol", "ncbigene", "ensembl_gene_id",
    "match_method", "status", "evidence",
]


def _mutation_labels(path: Path) -> list[str]:
    con = connect()
    rows = con.execute(
        """
        SELECT DISTINCT HugoSymbol, nullif(CAST(try_cast(EntrezGeneID AS BIGINT) AS VARCHAR), '0')
        FROM read_csv(?, header=true, all_varchar=true)
        """,
        [str(path)],
    ).fetchall()
    labels: list[str] = []
    for symbol, entrez in rows:
        if symbol and entrez:
            labels.append(f"{symbol} ({entrez})")
        elif symbol:
            labels.append(str(symbol))
    return labels


def build_gene_xref(document: dict, raw_root: Path, out: Path) -> dict[str, int]:
    roles = specs_by_role(document)
    hgnc_path = _first_file(raw_root, roles.get("hgnc_genes", []))
    if not hgnc_path:
        return {}
    labels: list[str] = []
    for source, spec in roles.get("gene_matrix", []):
        path = file_for(raw_root, source["id"], str(source["version"]), spec["name"])
        if path:
            labels.extend(header_labels(path))
    for source, spec in roles.get("depmap_mutations", []):
        path = file_for(raw_root, source["id"], str(source["version"]), spec["name"])
        if path:
            labels.extend(_mutation_labels(path))
    if not labels:
        return {}
    rows = [decision.row() for decision in resolve_labels(labels, HgncIndex.load(hgnc_path))]
    con = connect()
    _write_rows(con, rows, out / "gene_xref", GENE_COLUMNS)
    review = [row for row in rows if row["status"] != "mapped"]
    _write_rows(con, review, out / "review" / "genes", GENE_COLUMNS)
    return {
        "gene_labels": len(rows),
        "gene_mapped": sum(1 for row in rows if row["status"] == "mapped"),
    }


def melt_gene_matrix(
    con: duckdb.DuckDBPyConnection,
    path: Path,
    dataset: str,
    metric: str,
    xref: Path,
    model_source: str,
    dest: Path,
    gene_xref: Path | None = None,
) -> int:
    _replace_tree(dest)
    join = _model_join(model_source).replace(Path("__MODEL_XREF__").as_posix(), (xref / "part.parquet").as_posix())
    gene_file = gene_xref / "part.parquet" if gene_xref and (gene_xref / "part.parquet").exists() else None
    if gene_file:
        gene_sql = f"""
          src AS (
            SELECT
              ? AS dataset,
              long.source_model_id,
              genes.ncbigene AS ncbi_gene_id,
              genes.symbol AS symbol,
              long.value,
              ? AS metric
            FROM long
            JOIN read_parquet('{gene_file.as_posix()}') genes
              ON genes.source_label = long.gene_label AND genes.status = 'mapped'
            WHERE long.value IS NOT NULL
          )
        """
    else:
        gene_sql = f"""
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
        """
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
          {gene_sql}
          SELECT src.dataset, src.source_model_id, models.model_cvcl, src.ncbi_gene_id, src.symbol, src.value, src.metric
          FROM src
          {join}
        ) TO '{dest.as_posix()}' (FORMAT PARQUET, PARTITION_BY (dataset))
        """,
        [str(path), dataset, metric],
    )
    count = con.execute("SELECT count(*) FROM read_parquet(?, hive_partitioning=true)", [str(dest / "**/*.parquet")]).fetchone()
    return int(count[0] if count else 0)


def load_mutations(
    con: duckdb.DuckDBPyConnection,
    path: Path,
    dataset: str,
    xref: Path,
    dest: Path,
    gene_xref: Path | None = None,
) -> int:
    _replace_tree(dest)
    join = _model_join("depmap").replace(Path("__MODEL_XREF__").as_posix(), (xref / "part.parquet").as_posix())
    gene_file = gene_xref / "part.parquet" if gene_xref and (gene_xref / "part.parquet").exists() else None
    if gene_file:
        symbol_sql = f"""
          LEFT JOIN read_parquet('{gene_file.as_posix()}') genes
            ON genes.status = 'mapped'
           AND genes.source_label = src.symbol || ' (' || src.ncbi_gene_id || ')'
        """
        select_sql = "SELECT src.* REPLACE (coalesce(genes.symbol, src.symbol) AS symbol), models.model_cvcl"
    else:
        symbol_sql = ""
        select_sql = "SELECT src.*, models.model_cvcl"
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
          {select_sql}
          FROM src
          {symbol_sql}
          {join}
        ) TO '{dest.as_posix()}' (FORMAT PARQUET, PARTITION_BY (dataset))
        """,
        [dataset, str(path)],
    )
    count = con.execute("SELECT count(*) FROM read_parquet(?, hive_partitioning=true)", [str(dest / "**/*.parquet")]).fetchone()
    return int(count[0] if count else 0)


def load_gdsc_response(
    con: duckdb.DuckDBPyConnection,
    path: Path,
    xref: Path,
    dest_root: Path,
    combination_ids: set[str] | None = None,
) -> int:
    join_cvcl = f"""
      LEFT JOIN (
        SELECT source_model_id, any_value(cvcl) AS model_cvcl
        FROM read_parquet('{(xref / "part.parquet").as_posix()}')
        WHERE source = 'gdsc' AND alias_type = 'primary' AND status = 'mapped'
        GROUP BY source_model_id
      ) models ON models.source_model_id = base.source_model_id
    """
    combo_list = ",".join("'" + drug_id + "'" for drug_id in sorted(combination_ids or []) if drug_id.isdigit())
    if not combo_list:
        combo_list = "''"
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
            NULL::VARCHAR AS qc,
            CASE
              WHEN CAST(DRUG_ID AS VARCHAR) IN ({combo_list})
              THEN 'gdsc.combination:annotation:' || CAST(DRUG_ID AS VARCHAR)
              ELSE NULL
            END AS combination_id
          FROM (
            SELECT *, 'AUC' AS metric, CAST(AUC AS VARCHAR) AS metric_value, NULL AS unit FROM read_csv(?, header=true, all_varchar=true)
            UNION ALL
            SELECT *, 'LN_IC50' AS metric, CAST(LN_IC50 AS VARCHAR) AS metric_value, 'ln_uM' AS unit FROM read_csv(?, header=true, all_varchar=true)
          )
        )
        SELECT base.dataset, base.source_model_id, models.model_cvcl, base.source_compound_id,
               NULL::VARCHAR AS pubchem_cid, base.compound_name, base.metric, base.value, base.unit,
               base.dose, base.screen_id, base.qc,
               base.combination_id,
               CASE WHEN base.combination_id IS NULL THEN 'single' ELSE 'combination' END AS response_class,
               NULL::DOUBLE AS anchor_conc,
               NULL::DOUBLE AS library_conc
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


def write_part(con: duckdb.DuckDBPyConnection, sql: str, params: list[object], root: Path, dataset: str, part_name: str) -> int:
    """Write one parquet file inside a dataset directory without removing sibling parts."""
    directory = root / f"dataset={dataset}"
    directory.mkdir(parents=True, exist_ok=True)
    dest = directory / part_name
    if dest.exists():
        dest.unlink()
    con.execute(f"COPY ({sql}) TO '{dest.as_posix()}' (FORMAT PARQUET)", params)
    count = con.execute("SELECT count(*) FROM read_parquet(?)", [str(dest)]).fetchone()
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
            NULL::VARCHAR AS qc,
            NULL::VARCHAR AS combination_id,
            'single' AS response_class,
            NULL::DOUBLE AS anchor_conc,
            NULL::DOUBLE AS library_conc
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
            base.qc,
            NULL::VARCHAR AS combination_id,
            'single' AS response_class,
            NULL::DOUBLE AS anchor_conc,
            NULL::DOUBLE AS library_conc
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
                 src.value, NULL::VARCHAR AS unit, src.dose, src.screen_id, NULL::VARCHAR AS qc,
                 NULL::VARCHAR AS combination_id, 'single' AS response_class,
                 NULL::DOUBLE AS anchor_conc, NULL::DOUBLE AS library_conc
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


COMBINATION_COLUMNS = [
    "combination_id", "source", "screen", "name",
    "anchor_drug_id", "library_drug_id", "anchor_name", "library_name",
    "anchor_pubchem_cid", "library_pubchem_cid", "anchor_chembl_id", "library_chembl_id",
    "anchor_concentrations", "library_concentrations", "concentration_unit",
    "status", "source_drug_id",
]

_SYNERGY_METRICS = (
    "SYNERGY_DELTA_AUC",
    "SYNERGY_DELTA_EMAX",
    "SYNERGY_DELTA_XMID",
    "SYNERGY_DELTA_XMID_uM",
    "SYNERGY_DELTA_fAUC",
)


def annotation_combination_ids(raw_root: Path, roles: dict[str, list[tuple[dict, dict]]]) -> set[str]:
    from pgx.normalize.combinations import is_annotation_combination

    found: set[str] = set()
    for source, spec in roles.get("gdsc_compounds", []):
        path = file_for(raw_root, source["id"], str(source["version"]), spec["name"])
        if not path:
            continue
        con = connect()
        for drug_id, name in con.execute(
            "SELECT CAST(DRUG_ID AS VARCHAR), DRUG_NAME FROM read_csv(?, header=true, all_varchar=true)",
            [str(path)],
        ).fetchall():
            if drug_id and is_annotation_combination(name):
                found.add(str(drug_id))
    return found


def _blank_combination(**values: str | None) -> dict[str, str | None]:
    row = {column: None for column in COMBINATION_COLUMNS}
    row.update(values)
    row["status"] = row.get("status") or "unmatched"
    return row


def _union_concentrations(left: str | None, right: str | None) -> str | None:
    values: set[float] = set()
    for text in (left, right):
        if not text:
            continue
        for part in str(text).split(";"):
            try:
                values.add(float(part))
            except ValueError:
                continue
    if not values:
        return None
    return ";".join(f"{value:.6g}" for value in sorted(values))


def _merge_combination(rows: dict[str, dict[str, str | None]], incoming: dict[str, str | None]) -> None:
    key = incoming["combination_id"]
    if not key:
        return
    current = rows.get(key)
    if current is None:
        rows[key] = incoming
        return
    current["anchor_concentrations"] = _union_concentrations(current.get("anchor_concentrations"), incoming.get("anchor_concentrations"))
    current["library_concentrations"] = _union_concentrations(current.get("library_concentrations"), incoming.get("library_concentrations"))
    screens = [part for part in f"{current.get('screen') or ''};{incoming.get('screen') or ''}".split(";") if part]
    current["screen"] = ";".join(dict.fromkeys(screens))
    for field in ("anchor_drug_id", "library_drug_id", "anchor_name", "library_name", "concentration_unit", "source_drug_id"):
        if incoming.get(field) and not current.get(field):
            current[field] = incoming[field]


def combinations_from_annotations(raw_root: Path, roles: dict[str, list[tuple[dict, dict]]]) -> dict[str, dict[str, str | None]]:
    from pgx.normalize.combinations import parse_annotation_combination

    rows: dict[str, dict[str, str | None]] = {}
    for source, spec in roles.get("gdsc_compounds", []):
        path = file_for(raw_root, source["id"], str(source["version"]), spec["name"])
        if not path:
            continue
        con = connect()
        for drug_id, name in con.execute(
            "SELECT CAST(DRUG_ID AS VARCHAR), DRUG_NAME FROM read_csv(?, header=true, all_varchar=true)",
            [str(path)],
        ).fetchall():
            parsed = parse_annotation_combination(str(drug_id or ""), str(name or ""))
            if parsed:
                _merge_combination(rows, _blank_combination(**parsed))
    return rows


def _screen_identities(con: duckdb.DuckDBPyConnection, path: Path, screen: str) -> list[dict[str, str | None]]:
    from pgx.normalize.combinations import parse_screen_combination

    fetched = con.execute(
        """
        SELECT
          CAST(ANCHOR_ID AS VARCHAR),
          any_value(ANCHOR_NAME),
          CAST(LIBRARY_ID AS VARCHAR),
          any_value(LIBRARY_NAME),
          list(DISTINCT try_cast(ANCHOR_CONC AS DOUBLE)) FILTER (WHERE try_cast(ANCHOR_CONC AS DOUBLE) IS NOT NULL),
          list(DISTINCT try_cast(LIBRARY_CONC AS DOUBLE)) FILTER (WHERE try_cast(LIBRARY_CONC AS DOUBLE) IS NOT NULL)
        FROM read_csv(?, header=true, all_varchar=true)
        WHERE ANCHOR_ID IS NOT NULL AND LIBRARY_ID IS NOT NULL
        GROUP BY 1, 3
        """,
        [str(path)],
    ).fetchall()
    rows: list[dict[str, str | None]] = []
    for anchor_id, anchor_name, library_id, library_name, anchor_concs, library_concs in fetched:
        parsed = parse_screen_combination(
            str(anchor_id),
            anchor_name,
            str(library_id),
            library_name,
            [float(value) for value in anchor_concs or []],
            [float(value) for value in library_concs or []],
            screen=screen,
        )
        if parsed:
            rows.append(_blank_combination(**parsed))
    return rows


def load_combination_responses(
    con: duckdb.DuckDBPyConnection,
    path: Path,
    dataset: str,
    xref: Path,
    dest: Path,
) -> int:
    xref_file = (xref / "part.parquet").as_posix()
    metrics = ", ".join(_SYNERGY_METRICS)
    return _copy_dataset(
        con,
        f"""
          WITH raw AS (
            SELECT * FROM read_csv(?, header=true, all_varchar=true)
          ),
          long AS (
            SELECT * FROM (
              UNPIVOT raw
              ON {metrics}
              INTO NAME metric VALUE metric_value
            )
          ),
          base AS (
            SELECT
              coalesce(nullif(SIDM, ''), 'COSMIC:' || CAST(COSMIC_ID AS VARCHAR), CELL_LINE_NAME) AS source_model_id,
              CAST(ANCHOR_ID AS VARCHAR) || ':' || CAST(LIBRARY_ID AS VARCHAR) AS source_compound_id,
              ANCHOR_NAME || ' + ' || LIBRARY_NAME AS compound_name,
              'gdsc.combination:' || CAST(ANCHOR_ID AS VARCHAR) || ':' || CAST(LIBRARY_ID AS VARCHAR) AS combination_id,
              try_cast(ANCHOR_CONC AS DOUBLE) AS anchor_conc,
              try_cast(LIBRARY_CONC AS DOUBLE) AS library_conc,
              metric,
              try_cast(metric_value AS DOUBLE) AS value,
              CASE WHEN metric = 'SYNERGY_DELTA_XMID_uM' THEN 'uM' ELSE NULL END AS unit
            FROM long
          )
          SELECT
            base.source_model_id,
            models.model_cvcl,
            base.source_compound_id,
            NULL::VARCHAR AS pubchem_cid,
            base.compound_name,
            base.metric,
            base.value,
            base.unit,
            NULL::DOUBLE AS dose,
            NULL::VARCHAR AS screen_id,
            NULL::VARCHAR AS qc,
            base.combination_id,
            'combination' AS response_class,
            base.anchor_conc,
            base.library_conc
          FROM base
          LEFT JOIN (
            SELECT source_model_id, any_value(cvcl) AS model_cvcl
            FROM read_parquet('{xref_file}')
            WHERE alias_type = 'primary' AND status = 'mapped' AND source IN ('gdsc', 'sanger_cmp')
            GROUP BY source_model_id
            HAVING count(DISTINCT cvcl) = 1
          ) models ON models.source_model_id = base.source_model_id
          WHERE base.value IS NOT NULL AND base.combination_id IS NOT NULL
        """,
        [str(path)],
        dest,
        dataset,
    )


def load_combinations(document: dict, raw_root: Path, out: Path) -> dict[str, int]:
    """Write combination_xref and combination drug_response partitions.

    Single-agent response partitions are left in place. A combination row uses
    response_class `combination` and a combination_id, and its pubchem_cid is null.
    """
    roles = specs_by_role(document)
    rows = combinations_from_annotations(raw_root, roles)
    con = connect()
    response_root = out / "drug_response"
    response_root.mkdir(parents=True, exist_ok=True)
    responses = 0
    for source, spec in roles.get("gdsc_combination_screen", []):
        path = file_for(raw_root, source["id"], str(source["version"]), spec["name"])
        if not path:
            continue
        screen = str(spec.get("dataset") or "gdsc_combination")
        for record in _screen_identities(con, path, screen):
            _merge_combination(rows, record)
        responses += load_combination_responses(con, path, screen, out / "model_xref", response_root)
    if rows:
        _write_rows(con, list(rows.values()), out / "combination_xref", COMBINATION_COLUMNS)
        add_screen_compound_stubs(out, rows)
        refresh_combination_components(out)
    return {"combinations": len(rows), "combination_responses": responses}


def add_screen_compound_stubs(out: Path, rows: dict[str, dict[str, str | None]]) -> None:
    """Add GDSC compound rows for screen components missing from screened_compounds.

    They stay unmatched until resolve-compounds maps each name on its own.
    """
    compounds = out / "compound_xref" / "part.parquet"
    if not compounds.exists():
        return
    con = connect()
    columns = [item[0] for item in con.execute("DESCRIBE SELECT * FROM read_parquet(?)", [str(compounds)]).fetchall()]
    records: dict[tuple[str, str], CompoundRecord] = {}
    for fetched in con.execute(f"SELECT {', '.join(columns)} FROM read_parquet(?)", [str(compounds)]).fetchall():
        data = dict(zip(columns, fetched, strict=True))
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
            review_reason=data.get("review_reason"),
        )
        if data.get("synonyms"):
            record.synonyms.update(part for part in str(data["synonyms"]).split(";") if part)
        add_compound(records, record)
    added = False
    for row in rows.values():
        for role in ("anchor", "library"):
            raw_id = row.get(f"{role}_drug_id")
            raw_name = row.get(f"{role}_name")
            if not raw_id:
                continue
            ids = [part.strip() for part in str(raw_id).split("|")]
            names = [part.strip() for part in str(raw_name).split("|")] if raw_name else []
            for index, drug_id in enumerate(ids):
                if not drug_id or ("gdsc", drug_id) in records:
                    continue
                name = names[index] if index < len(names) else None
                add_compound(records, CompoundRecord(
                    source="gdsc",
                    source_compound_id=drug_id,
                    name=name,
                    gdsc_drug_id=drug_id,
                    match_method="unmatched",
                    status="unmatched",
                    evidence="combination screen component; not in screened_compounds",
                ))
                added = True
    if not added:
        return
    written = compound_rows(records)
    fields = list(written[0]) if written else columns
    _write_rows(con, written, out / "compound_xref", fields)


def _join_parts(values: list[str]) -> str | None:
    if not any(values):
        return None
    return "|".join(values)


def refresh_combination_components(out: Path) -> None:
    """Copy component PubChem and ChEMBL ids from single-agent compound_xref.

    An ambiguous component contributes no CID. A `|` id is two drugs in that
    role, kept in source order. The combination itself is never given a CID.
    Status is mapped when every component has one CID, partial when some do,
    and unmatched otherwise. Several GDSC ids may share one name and one CID.
    """
    combo = out / "combination_xref" / "part.parquet"
    compounds = out / "compound_xref" / "part.parquet"
    if not combo.exists() or not compounds.exists():
        return
    con = connect()
    by_id: dict[str, tuple[str, str | None]] = {}
    by_name: dict[str, tuple[set[str], set[str], set[str | None]]] = {}
    for drug_id, name, cid, chembl in con.execute(
        """
        SELECT source_compound_id, lower(name), pubchem_cid, chembl_id
        FROM read_parquet(?)
        WHERE source = 'gdsc' AND status = 'mapped' AND pubchem_cid IS NOT NULL
        """,
        [str(compounds)],
    ).fetchall():
        current = by_id.get(str(drug_id))
        if current and current[0] != cid:
            by_id[str(drug_id)] = ("", None)
        elif current is None:
            by_id[str(drug_id)] = (str(cid), str(chembl) if chembl else None)
        key = (name or "").casefold()
        if key:
            ids, cids, chembls = by_name.setdefault(key, (set(), set(), set()))
            ids.add(str(drug_id))
            cids.add(str(cid))
            if chembl:
                chembls.add(str(chembl))
    columns = [item[0] for item in con.execute("DESCRIBE SELECT * FROM read_parquet(?)", [str(combo)]).fetchall()]
    rows: list[dict[str, str | None]] = []
    for fetched in con.execute(f"SELECT {', '.join(columns)} FROM read_parquet(?)", [str(combo)]).fetchall():
        row = dict(zip(columns, fetched, strict=True))
        for role in ("anchor", "library"):
            raw_id = row.get(f"{role}_drug_id")
            raw_name = row.get(f"{role}_name")
            if raw_id:
                parts = str(raw_id).split("|")
                cids: list[str] = []
                chembls: list[str] = []
                for part in parts:
                    found = by_id.get(part)
                    if not found or not found[0]:
                        cids.append("")
                        chembls.append("")
                    else:
                        cids.append(found[0])
                        chembls.append(found[1] or "")
                row[f"{role}_pubchem_cid"] = _join_parts(cids)
                row[f"{role}_chembl_id"] = _join_parts(chembls)
            else:
                found_name = by_name.get((raw_name or "").casefold())
                if found_name and len(found_name[1]) == 1:
                    row[f"{role}_pubchem_cid"] = next(iter(found_name[1]))
                    row[f"{role}_chembl_id"] = next(iter(found_name[2])) if len(found_name[2]) == 1 else None
                    if len(found_name[0]) == 1:
                        row[f"{role}_drug_id"] = next(iter(found_name[0]))
                else:
                    row[f"{role}_pubchem_cid"] = None
                    row[f"{role}_chembl_id"] = None
        pieces = []
        for role in ("anchor", "library"):
            joined = row.get(f"{role}_pubchem_cid") or ""
            pieces.extend(joined.split("|") if joined else [""])
        filled = [part for part in pieces if part]
        if filled and len(filled) == len(pieces):
            row["status"] = "mapped"
        elif filled:
            row["status"] = "partial"
        else:
            row["status"] = "unmatched"
        rows.append({column: None if row.get(column) in {None, ""} else str(row.get(column)) for column in COMBINATION_COLUMNS})
    _write_rows(con, rows, out / "combination_xref", COMBINATION_COLUMNS)


def attach_pubchem(con: duckdb.DuckDBPyConnection, response_root: Path, xref: Path) -> None:
    """Fill pubchem_cid on drug_response partitions from a mapped compound_xref."""
    if not response_root.exists():
        return
    xref_file = (xref / "part.parquet").as_posix()
    for part in sorted(response_root.glob("dataset=*/*.parquet")):
        temporary = part.with_suffix(".parquet.next")
        con.execute(
            f"""
            COPY (
              SELECT response.* REPLACE (
                CASE
                  WHEN coalesce(response.response_class, 'single') = 'combination' THEN NULL
                  ELSE coalesce(response.pubchem_cid, compounds.pubchem_cid)
                END AS pubchem_cid
              )
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


def _add_batch2_xrefs(
    raw_root: Path,
    roles: dict[str, list[tuple[dict, dict]]],
    document: dict,
    cellosaurus: object | None,
    model_rows: list[dict[str, str | None]],
    records: dict[tuple[str, str], CompoundRecord],
) -> None:
    from pgx.nci60 import assay_tables, collect_compounds, collect_models
    from pgx.pset import collect_source

    nci_pairs: list[tuple[dict, dict]] = []
    for role in ("nci60_endpoint", "nci60_oneconc", "nci60_doseresp"):
        nci_pairs.extend(roles.get(role, []))
    opened = assay_tables(raw_root, nci_pairs)
    con = connect()
    if opened:
        tables = [path for _source, _spec, path in opened]
        model_rows.extend(collect_models(con, tables, cellosaurus))
        for record in collect_compounds(con, tables):
            add_compound(records, record)
    for source in document.get("sources") or []:
        if not source.get("enabled", True):
            continue
        models, compounds = collect_source(con, raw_root, source, cellosaurus)
        model_rows.extend(models)
        for record in compounds:
            add_compound(records, record)


def _load_batch2_responses(
    con: duckdb.DuckDBPyConnection,
    document: dict,
    raw_root: Path,
    roles: dict[str, list[tuple[dict, dict]]],
    xref: Path,
    response_root: Path,
) -> int:
    from pgx.nci60 import assay_tables, load_responses
    from pgx.pset import load_source

    nci_pairs: list[tuple[dict, dict]] = []
    for role in ("nci60_endpoint", "nci60_oneconc", "nci60_doseresp"):
        nci_pairs.extend(roles.get(role, []))
    opened = assay_tables(raw_root, nci_pairs)
    tables = [(spec["name"], spec.get("role") or "", path) for _source, spec, path in opened]
    total = load_responses(con, tables, xref, response_root) if tables else 0
    for source in document.get("sources") or []:
        if not source.get("enabled", True):
            continue
        total += load_source(con, raw_root, source, xref, response_root)
    return total


def process_all(document: dict, raw_root: Path, out: Path) -> dict[str, int]:
    out.mkdir(parents=True, exist_ok=True)
    counts = build_xrefs(document, raw_root, out)
    counts.update(build_gene_xref(document, raw_root, out))
    roles = specs_by_role(document)
    for name in ("gene_effect", "expression", "cnv", "mutation", "drug_response"):
        existing = out / name
        if existing.exists():
            shutil.rmtree(existing)
    con = connect()
    xref = out / "model_xref"
    gene_xref = out / "gene_xref" if (out / "gene_xref" / "part.parquet").exists() else None
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
        counts[table] = matrix_roles["gene_matrix"](con, path, dataset, str(spec["metric"]), xref, "depmap", out / table, gene_xref)
    for source, spec in roles.get("depmap_mutations", []):
        path = file_for(raw_root, source["id"], str(source["version"]), spec["name"])
        if path:
            counts["mutation"] = load_mutations(con, path, str(spec.get("dataset") or "depmap_mutations_24q2"), xref, out / "mutation", gene_xref)
    response_root = out / "drug_response"
    response_root.mkdir(parents=True, exist_ok=True)
    total_response = 0
    combo_ids = annotation_combination_ids(raw_root, roles)
    for source, spec in roles.get("gdsc_dose_response", []):
        path = file_for(raw_root, source["id"], str(source["version"]), spec["name"])
        if path:
            total_response += load_gdsc_response(con, path, xref, response_root, combo_ids)
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
    total_response += _load_batch2_responses(con, document, raw_root, roles, xref, response_root)
    counts.update(load_combinations(document, raw_root, out))
    total_response += counts.get("combination_responses", 0)
    if total_response:
        attach_pubchem(con, response_root, out / "compound_xref")
        counts["drug_response"] = total_response
    return counts
