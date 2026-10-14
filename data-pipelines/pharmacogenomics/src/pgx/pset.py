"""Frozen ORCESTRA PharmacoSets.

Published TSV extracts are loaded directly. An ``.rds`` file is converted by
``r/extract_pset.R`` when ``Rscript`` is on ``PATH``. Dataset ids are
``ctrpv2``, ``gcsi``, ``fimm``, ``gray``, ``uhnbreast``, and ``ccle``.
Numeric drug ids are prefixed with that dataset id so they cannot collide
with a GDSC drug id.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import duckdb

from pgx.ids import CHEMBL_RE, CompoundRecord, Decision, clean, parse_brd, parse_cvcl, xref_rows

SCRIPT = Path(__file__).resolve().parents[2] / "r" / "extract_pset.R"
DATASETS = {
    "orcestra_ctrpv2": "ctrpv2",
    "orcestra_gcsi": "gcsi",
    "orcestra_fimm": "fimm",
    "orcestra_gray": "gray",
    "orcestra_uhnhbreast": "uhnbreast",
    "orcestra_ccle": "ccle",
    "orcestra_nci60": "nci60_pset",
}
PROFILE_METRICS = {
    "auc": ("AUC", None),
    "area_under_curve": ("AUC", None),
    "ic50": ("IC50", None),
    "ic50_recomputed": ("IC50", None),
    "ln_ic50": ("LN_IC50", "ln_uM"),
    "aac": ("AAC", None),
    "aac_recomputed": ("AAC", None),
    "auc_recomputed": ("AUC", None),
    "ec50": ("EC50", None),
    "ac50_recomputed": ("AC50", None),
}
_CPD = re.compile(r"_CPD_(\d+)", re.IGNORECASE)


def dataset_of(source: dict) -> str:
    declared = clean(source.get("dataset"))
    if declared:
        return declared
    return DATASETS.get(source["id"], source["id"])


def compound_key(dataset: str, value: object) -> str | None:
    text = clean(value)
    if not text:
        return None
    if text.isdigit():
        return f"{dataset}:{text}"
    return text


def extract_rds(rds: Path, dest: Path, script: Path | None = None) -> str:
    """Convert one RDS file. Returns ``extracted`` or ``skipped-no-rscript``."""
    rscript = shutil.which("Rscript")
    if not rscript:
        return "skipped-no-rscript"
    dest.mkdir(parents=True, exist_ok=True)
    subprocess.run([rscript, str(script or SCRIPT), str(rds), str(dest)], check=True)
    if not (dest / "cells.csv").exists() and not (dest / "drugs.csv").exists():
        raise RuntimeError(f"{rds.name} did not produce cells.csv or drugs.csv")
    return "extracted"


def _columns(con: duckdb.DuckDBPyConnection, path: Path) -> list[str]:
    described = con.execute(
        "DESCRIBE SELECT * FROM read_csv(?, header=true, all_varchar=true)",
        [str(path)],
    ).fetchall()
    return [row[0] for row in described]


def _pick(columns: list[str], names: list[str]) -> str | None:
    folded = {column.casefold(): column for column in columns}
    for name in names:
        found = folded.get(name.casefold())
        if found:
            return found
    return None


def _cid_pair(left: object, right: object) -> tuple[str | None, str]:
    """Prefer one source CID. Two different CIDs stay ambiguous."""
    from pgx.ids import _cid

    first = _cid(left)
    second = _cid(right)
    if first and second and first != second:
        return None, "ambiguous"
    chosen = first or second
    if chosen:
        return chosen, "mapped"
    return None, "unmatched"


def _models_from_rows(dataset: str, rows: list[dict[str, str | None]], cellosaurus: object | None) -> list[dict[str, str | None]]:
    out: list[dict[str, str | None]] = []
    seen: set[str] = set()
    for row in rows:
        source_model_id = clean(row.get("source_model_id"))
        if not source_model_id or source_model_id in seen:
            continue
        seen.add(source_model_id)
        name = clean(row.get("cell_name")) or source_model_id
        supplied = parse_cvcl(row.get("cellosaurus"))
        if supplied:
            decision = Decision(supplied, "pset_cellosaurus", "mapped", f"PSet Cellosaurus {supplied}")
            if cellosaurus is not None:
                decision = cellosaurus.harmonize(decision, [name])
        elif cellosaurus is not None:
            decision = cellosaurus.decide_names([name])
        else:
            decision = Decision(None, "unmatched", "unmatched", "no Cellosaurus index")
        out.extend(xref_rows(
            dataset,
            source_model_id,
            decision,
            [("cell_name", name)],
            display_name=name,
            tissue=clean(row.get("tissue")),
            species="Homo sapiens",
        ))
    return out


def _compound_from_row(dataset: str, row: dict[str, str | None]) -> CompoundRecord | None:
    source_compound_id = compound_key(dataset, row.get("source_compound_id"))
    if not source_compound_id:
        return None
    cid, state = _cid_pair(row.get("pubchem_cid"), row.get("pubchem_cid_alt"))
    chembl = clean(row.get("chembl_id"))
    if chembl and not CHEMBL_RE.match(chembl):
        chembl = None
    record = CompoundRecord(
        source=dataset,
        source_compound_id=source_compound_id,
        name=clean(row.get("name")),
        brd_id=parse_brd(row.get("brd_id")),
        pubchem_cid=cid,
        chembl_id=chembl,
        inchikey=clean(row.get("inchikey")),
        smiles=clean(row.get("smiles")),
    )
    if state == "ambiguous":
        record.status = "ambiguous"
        record.match_method = "conflict_pubchem"
        record.confidence = None
        record.evidence = "PSet annotation CIDs disagree"
    elif state == "mapped" and cid:
        record.status = "mapped"
        record.match_method = "pset_pubchem"
        record.confidence = "high"
        record.evidence = "PSet PubChem CID"
    return record


def _read_dicts(con: duckdb.DuckDBPyConnection, path: Path) -> list[dict[str, str | None]]:
    columns = _columns(con, path)
    fetched = con.execute(
        "SELECT * FROM read_csv(?, header=true, all_varchar=true)",
        [str(path)],
    ).fetchall()
    return [dict(zip(columns, row, strict=True)) for row in fetched]


def ctrp_annotations(
    con: duckdb.DuckDBPyConnection,
    sample: Path | None,
    treatment: Path | None,
    screened_cells: set[str],
    screened_drugs: set[str],
) -> tuple[list[dict[str, str | None]], list[dict[str, str | None]]]:
    cells: list[dict[str, str | None]] = []
    drugs: list[dict[str, str | None]] = []
    if sample and sample.exists():
        for row in _read_dicts(con, sample):
            master = clean(row.get("master_ccl_id"))
            if screened_cells and master not in screened_cells:
                continue
            if not screened_cells and clean(row.get("In.CTRPv2.1")) not in {None, "TRUE", "True", "true", "1"}:
                continue
            cells.append({
                "source_model_id": master,
                "cell_name": clean(row.get("CELL_NAME") or row.get("AnnotationDB.CellLine.Name")),
                "cellosaurus": clean(row.get("Cellosaurus.Accession") or row.get("AnnotationDB.Cellosaurus.Accession")),
                "tissue": clean(row.get("ccle_primary_site")),
            })
    if treatment and treatment.exists():
        for row in _read_dicts(con, treatment):
            master = clean(row.get("master_cpd_id"))
            if screened_drugs and master not in screened_drugs:
                continue
            smiles = clean(row.get("AnnotationDB.Canonical.SMILES") or row.get("SMILES") or row.get("AnnotationGx.Connectivity.SMILES"))
            inchikey = clean(row.get("AnnotationGx.InChIKey") or row.get("AnnotationDB.InChIKey"))
            drugs.append({
                "source_compound_id": master,
                "name": clean(row.get("Compound.Name") or row.get("AnnotationGx.Title") or row.get("AnnotationDB.Title")),
                "pubchem_cid": clean(row.get("AnnotationGx.PubChem.CID")),
                "pubchem_cid_alt": clean(row.get("AnnotationDB.PubChem.CID")),
                "smiles": smiles,
                "inchikey": inchikey,
                "brd_id": clean(row.get("Broad.Compound.ID")),
                "chembl_id": clean(row.get("AnnotationDB.ChEMBL.ID")),
            })
    return cells, drugs


def screened_ids(con: duckdb.DuckDBPyConnection, path: Path, cell_names: list[str], drug_names: list[str]) -> tuple[set[str], set[str]]:
    columns = _columns(con, path)
    cell_col = _pick(columns, cell_names)
    drug_col = _pick(columns, drug_names)
    cells: set[str] = set()
    drugs: set[str] = set()
    if cell_col:
        for (value,) in con.execute(
            f'SELECT DISTINCT "{cell_col}" FROM read_csv(?, header=true, all_varchar=true)',
            [str(path)],
        ).fetchall():
            token = clean(value)
            if token:
                cells.add(token)
    if drug_col:
        for (value,) in con.execute(
            f'SELECT DISTINCT "{drug_col}" FROM read_csv(?, header=true, all_varchar=true)',
            [str(path)],
        ).fetchall():
            token = clean(value)
            if token:
                drugs.add(token)
            match = _CPD.search(token or "")
            if match:
                drugs.add(match.group(1))
    return cells, drugs


def _model_join(xref_file: str, dataset: str) -> str:
    return f"""
    LEFT JOIN (
      SELECT source_model_id, any_value(cvcl) AS model_cvcl
      FROM read_parquet('{xref_file}')
      WHERE source = '{dataset}' AND alias_type = 'primary' AND status = 'mapped'
      GROUP BY source_model_id
    ) models ON models.source_model_id = base.source_model_id
    """


def _wrap(body: str, xref_file: str, dataset: str) -> str:
    return f"""
    WITH base AS (
      {body}
    )
    SELECT
      base.source_model_id,
      models.model_cvcl,
      base.source_compound_id,
      NULL::VARCHAR AS pubchem_cid,
      NULL::VARCHAR AS compound_name,
      base.metric,
      base.value,
      base.unit,
      base.dose,
      base.screen_id,
      NULL::VARCHAR AS qc,
      NULL::VARCHAR AS combination_id,
      'single' AS response_class,
      NULL::DOUBLE AS anchor_conc,
      NULL::DOUBLE AS library_conc
    FROM base
    {_model_join(xref_file, dataset)}
    WHERE base.value IS NOT NULL
      AND base.source_model_id IS NOT NULL
      AND base.source_compound_id IS NOT NULL
    """


def auc_sql(dataset: str) -> str:
    return f"""
      SELECT
        CAST(master_ccl_id AS VARCHAR) AS source_model_id,
        '{dataset}:' || regexp_replace(CAST(master_cpd_id AS VARCHAR), '\\.0$', '') AS source_compound_id,
        'AUC' AS metric,
        try_cast(area_under_curve AS DOUBLE) AS value,
        NULL::VARCHAR AS unit,
        NULL::DOUBLE AS dose,
        CAST(experiment_id AS VARCHAR) AS screen_id
      FROM read_csv(?, header=true, all_varchar=true)
    """


def melted_profile_sql(dataset: str) -> str:
    """Profiles already stored one metric per row by extract_pset.R."""
    return f"""
      SELECT
        CAST(source_model_id AS VARCHAR) AS source_model_id,
        CASE
          WHEN regexp_full_match(CAST(source_compound_id AS VARCHAR), '^[0-9]+$')
          THEN '{dataset}:' || regexp_replace(CAST(source_compound_id AS VARCHAR), '\\.0$', '')
          ELSE CAST(source_compound_id AS VARCHAR)
        END AS source_compound_id,
        metric,
        try_cast(value AS DOUBLE) AS value,
        unit,
        NULL::DOUBLE AS dose,
        NULL::VARCHAR AS screen_id
      FROM read_csv(?, header=true, all_varchar=true)
    """


def long_profile_sql(dataset: str, columns: list[str]) -> tuple[str, list[str]]:
    """Melt a profiles table that already has one row per experiment or pair."""
    cell_col = _pick(columns, ["source_model_id", "cellid", "master_ccl_id", "sampleid"])
    drug_col = _pick(columns, ["source_compound_id", "drugid", "master_cpd_id", "treatmentid"])
    if not cell_col or not drug_col:
        raise RuntimeError("profiles table needs a cell id and a drug id")
    pieces = []
    used = []
    for column in columns:
        spec = PROFILE_METRICS.get(column.casefold())
        if not spec:
            continue
        metric, unit = spec
        unit_sql = "NULL" if unit is None else f"'{unit}'"
        drug_expr = (
            f"CASE WHEN regexp_full_match(CAST(\"{drug_col}\" AS VARCHAR), '^[0-9]+$') "
            f"THEN '{dataset}:' || regexp_replace(CAST(\"{drug_col}\" AS VARCHAR), '\\\\.0$', '') "
            f"ELSE CAST(\"{drug_col}\" AS VARCHAR) END"
        )
        pieces.append(f"""
          SELECT
            CAST("{cell_col}" AS VARCHAR) AS source_model_id,
            {drug_expr} AS source_compound_id,
            '{metric}' AS metric,
            try_cast("{column}" AS DOUBLE) AS value,
            {unit_sql}::VARCHAR AS unit,
            NULL::DOUBLE AS dose,
            NULL::VARCHAR AS screen_id
          FROM read_csv(?, header=true, all_varchar=true)
        """)
        used.append(column)
    if not pieces:
        raise RuntimeError("profiles table has no AUC/IC50/AAC/EC50 column")
    return " UNION ALL ".join(pieces), used


def raw_long_sql(dataset: str) -> str:
    return f"""
      SELECT
        CAST(source_model_id AS VARCHAR) AS source_model_id,
        CASE
          WHEN regexp_full_match(CAST(source_compound_id AS VARCHAR), '^[0-9]+$')
          THEN '{dataset}:' || CAST(source_compound_id AS VARCHAR)
          ELSE CAST(source_compound_id AS VARCHAR)
        END AS source_compound_id,
        coalesce(metric, 'VIABILITY') AS metric,
        try_cast(value AS DOUBLE) AS value,
        coalesce(unit, 'percent') AS unit,
        try_cast(dose AS DOUBLE) AS dose,
        NULL::VARCHAR AS screen_id
      FROM read_csv(?, header=true, all_varchar=true)
    """


def ctrp_sensitivity_sql(dataset: str, columns: list[str]) -> str:
    """Unpivot CTRP ``dose_umol_N`` / ``viability_N`` pairs in one scan.

    The drug id is the ``_CPD_`` number inside ``response_id``.
    """
    pairs = []
    for column in columns:
        match = re.fullmatch(r"dose_umol_(\d+)", column)
        if not match:
            continue
        viability = f"viability_{match.group(1)}"
        if viability not in columns:
            continue
        pairs.append((column, viability))
    if not pairs or "master_ccl_id" not in columns or "response_id" not in columns:
        raise RuntimeError("CTRP sensitivity table is missing dose/viability pairs")
    doses = ", ".join(f'try_cast("{dose}" AS DOUBLE)' for dose, _viability in pairs)
    viabilities = ", ".join(f'try_cast("{viability}" AS DOUBLE)' for _dose, viability in pairs)
    return f"""
      SELECT
        CAST(master_ccl_id AS VARCHAR) AS source_model_id,
        '{dataset}:' || regexp_extract(response_id, '_CPD_([0-9]+)', 1) AS source_compound_id,
        'VIABILITY' AS metric,
        unnest([{viabilities}]) AS value,
        'percent' AS unit,
        unnest([{doses}]) AS dose,
        response_id AS screen_id
      FROM read_csv(?, header=true, all_varchar=true)
      WHERE regexp_extract(response_id, '_CPD_([0-9]+)', 1) <> ''
    """


def recomputed_sql(dataset: str, columns: list[str], with_cells: bool) -> str:
    """IC50/AUC from recomputed profiles. Cell id is joined from sensitivity when present."""
    metrics = []
    for column, (metric, unit) in (
        ("ic50_recomputed", ("IC50", "uM")),
        ("auc_recomputed", ("AUC", None)),
        ("aac_recomputed", ("AAC", None)),
        ("ac50_recomputed", ("AC50", "uM")),
        ("EC50", ("EC50", "uM")),
    ):
        if column not in columns:
            continue
        unit_sql = "NULL" if unit is None else f"'{unit}'"
        metrics.append((column, metric, unit_sql))
    if not metrics or "response_id" not in columns:
        raise RuntimeError("recomputed profiles have no response_id metric columns")
    cell_expr = "cells.master_ccl_id" if with_cells else "NULL"
    join = """
      LEFT JOIN (
        SELECT response_id AS response_key, any_value(master_ccl_id) AS master_ccl_id
        FROM read_csv(?, header=true, all_varchar=true)
        GROUP BY response_id
      ) cells ON cells.response_key = profiles.response_id
    """ if with_cells else ""
    pieces = []
    for column, metric, unit_sql in metrics:
        pieces.append(f"""
          SELECT
            CAST({cell_expr} AS VARCHAR) AS source_model_id,
            '{dataset}:' || regexp_extract(profiles.response_id, '_CPD_([0-9]+)', 1) AS source_compound_id,
            '{metric}' AS metric,
            try_cast(profiles."{column}" AS DOUBLE) AS value,
            {unit_sql}::VARCHAR AS unit,
            NULL::DOUBLE AS dose,
            profiles.response_id AS screen_id
          FROM read_csv(?, header=true, all_varchar=true) profiles
          {join}
        """)
    return " UNION ALL ".join(pieces)


def _source_files(raw_root: Path, source: dict) -> dict[str, list[Path]]:
    from pgx.process import file_for

    found: dict[str, list[Path]] = {}
    for spec in source.get("files") or []:
        path = file_for(raw_root, source["id"], str(source["version"]), spec["name"])
        if path:
            found.setdefault(spec.get("role") or "", []).append(path)
    return found


def _extract_dir(raw_root: Path, source: dict, files: dict[str, list[Path]]) -> Path | None:
    dest = raw_root / source["id"] / str(source["version"]) / "extract"
    if (dest / "cells.csv").exists() or (dest / "profiles.csv").exists():
        return dest
    rds_list = files.get("pset_rds") or []
    if not rds_list:
        return None
    if extract_rds(rds_list[0], dest) == "skipped-no-rscript":
        return None
    return dest


def _annotation_rows(
    con: duckdb.DuckDBPyConnection,
    files: dict[str, list[Path]],
    extract: Path | None,
) -> tuple[list[dict[str, str | None]], list[dict[str, str | None]]]:
    profile = (files.get("pset_profiles") or [None])[0]
    cells: set[str] = set()
    drugs: set[str] = set()
    if profile:
        cells, drugs = screened_ids(
            con,
            profile,
            ["master_ccl_id", "source_model_id", "cellid"],
            ["master_cpd_id", "source_compound_id", "drugid"],
        )
    sample = (files.get("pset_cells") or [None])[0]
    treatment = (files.get("pset_drugs") or [None])[0]
    if sample and "master_ccl_id" in _columns(con, sample):
        cell_rows, drug_rows = ctrp_annotations(con, sample, treatment, cells, drugs)
    elif extract and (extract / "cells.csv").exists():
        cell_rows = _read_dicts(con, extract / "cells.csv")
        drug_rows = _read_dicts(con, extract / "drugs.csv") if (extract / "drugs.csv").exists() else []
    elif sample:
        cell_rows = _read_dicts(con, sample)
        drug_rows = _read_dicts(con, treatment) if treatment else []
    else:
        cell_rows, drug_rows = [], []
    known_cells = {clean(row.get("source_model_id")) for row in cell_rows}
    for cell_id in cells:
        if cell_id not in known_cells:
            cell_rows.append({"source_model_id": cell_id, "cell_name": cell_id})
    known_drugs = {clean(row.get("source_compound_id")) for row in drug_rows}
    for drug_id in drugs:
        if drug_id not in known_drugs:
            drug_rows.append({"source_compound_id": drug_id})
    return cell_rows, drug_rows


def collect_source(
    con: duckdb.DuckDBPyConnection,
    raw_root: Path,
    source: dict,
    cellosaurus: object | None,
) -> tuple[list[dict[str, str | None]], list[CompoundRecord]]:
    files = _source_files(raw_root, source)
    extract_ready = (raw_root / source["id"] / str(source["version"]) / "extract" / "cells.csv").exists()
    if not any(role.startswith("pset_") for role in files) and not extract_ready:
        return [], []
    dataset = dataset_of(source)
    extract = _extract_dir(raw_root, source, files)
    cell_rows, drug_rows = _annotation_rows(con, files, extract)
    models = _models_from_rows(dataset, cell_rows, cellosaurus)
    records = [record for row in drug_rows if (record := _compound_from_row(dataset, row))]
    return models, records


def load_source(
    con: duckdb.DuckDBPyConnection,
    raw_root: Path,
    source: dict,
    xref: Path,
    dest_root: Path,
) -> int:
    from pgx.process import write_part

    files = _source_files(raw_root, source)
    extract_ready = (raw_root / source["id"] / str(source["version"]) / "extract" / "cells.csv").exists() or (
        raw_root / source["id"] / str(source["version"]) / "extract" / "profiles.csv"
    ).exists()
    if not any(role.startswith("pset_") for role in files) and not extract_ready:
        return 0
    dataset = dataset_of(source)
    extract = _extract_dir(raw_root, source, files)
    xref_file = (xref / "part.parquet").as_posix()
    directory = dest_root / f"dataset={dataset}"
    if directory.exists():
        shutil.rmtree(directory)
    total = 0
    profile = (files.get("pset_profiles") or [None])[0]
    if profile is None and extract and (extract / "profiles.csv").exists():
        profile = extract / "profiles.csv"
    if profile:
        columns = _columns(con, profile)
        if "area_under_curve" in columns and "master_ccl_id" in columns:
            sql = _wrap(auc_sql(dataset), xref_file, dataset)
            total += write_part(con, sql, [str(profile)], dest_root, dataset, "auc.parquet")
        elif {"metric", "value", "source_model_id", "source_compound_id"} <= set(columns):
            sql = _wrap(melted_profile_sql(dataset), xref_file, dataset)
            total += write_part(con, sql, [str(profile)], dest_root, dataset, "profiles.parquet")
        else:
            body, used = long_profile_sql(dataset, columns)
            sql = _wrap(body, xref_file, dataset)
            total += write_part(con, sql, [str(profile)] * len(used), dest_root, dataset, "profiles.parquet")
    raw = (files.get("pset_raw") or [None])[0]
    if raw is None and extract and (extract / "raw.csv").exists():
        raw = extract / "raw.csv"
    if raw:
        columns = _columns(con, raw)
        if "dose_umol_1" in columns:
            sql = _wrap(ctrp_sensitivity_sql(dataset, columns), xref_file, dataset)
            total += write_part(con, sql, [str(raw)], dest_root, dataset, "raw.parquet")
        elif "dose" in columns and "value" in columns:
            sql = _wrap(raw_long_sql(dataset), xref_file, dataset)
            total += write_part(con, sql, [str(raw)], dest_root, dataset, "raw.parquet")
    recomputed = (files.get("pset_recomputed") or [None])[0]
    if recomputed and raw and "response_id" in _columns(con, raw):
        columns = _columns(con, recomputed)
        body = recomputed_sql(dataset, columns, with_cells=True)
        reads = body.count("read_csv(")
        if reads % 2 != 0:
            raise RuntimeError(f"recomputed SQL has {reads} reads")
        params: list[str] = []
        for _metric in range(reads // 2):
            params.extend([str(recomputed), str(raw)])
        sql = _wrap(body, xref_file, dataset)
        total += write_part(con, sql, params, dest_root, dataset, "recomputed.parquet")
    return total
