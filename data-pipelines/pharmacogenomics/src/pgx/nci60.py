"""NCI-60 DTP growth-inhibition tables.

Compounds are stored as ``NSC{digits}``. A bare NSC number would collide with
a GDSC drug id in the PubChem join, which matches ``source_compound_id`` only.
GI50, TGI, LC50, and IC50 share dataset ``nci60`` and are distinguished by
``metric``. ONECONC and DOSERESP are the same dataset.
"""

from __future__ import annotations

import csv
import shutil
import subprocess
import zipfile
from pathlib import Path

import duckdb

from pgx.ids import CompoundRecord, Decision, clean, xref_rows

NCI_SOURCE = "nci60"
ENDPOINT_METRICS = {"GI50", "TGI", "LC50", "IC50"}
CURATED_CELLS = Path(__file__).resolve().parents[2] / "curated" / "nci60_cell_lines.csv"


def nsc_id(value: object) -> str | None:
    text = clean(value)
    if not text:
        return None
    if text.endswith(".0") and text[:-2].isdigit():
        text = text[:-2]
    if text.isdigit():
        return f"NSC{text}"
    folded = text.upper().replace(" ", "")
    if folded.startswith("NSC") and folded[3:].isdigit():
        return folded
    return None


def open_table(path: Path) -> Path:
    """Return a delimited table, extracting a zip member beside the archive."""
    if path.suffix.lower() != ".zip":
        return path
    dest_dir = path.parent / f"{path.stem}_unzipped"
    dest_dir.mkdir(parents=True, exist_ok=True)
    try:
        with zipfile.ZipFile(path) as archive:
            members = [name for name in archive.namelist() if name.lower().endswith((".csv", ".txt")) and not name.endswith("/")]
            if len(members) != 1:
                raise RuntimeError(f"{path.name} contains {len(members)} tables; expected one")
            target = dest_dir / Path(members[0]).name
            if not target.exists() or target.stat().st_size == 0:
                with archive.open(members[0]) as src, target.open("wb") as dst:
                    shutil.copyfileobj(src, dst)
            return target
    except NotImplementedError as exc:
        extracted = _extract_deflate64(path, dest_dir)
        if extracted:
            return extracted
        raise RuntimeError(
            f"{path.name} uses a compression method the Python zip reader cannot extract. Install 7z and retry."
        ) from exc


def _extract_deflate64(path: Path, dest_dir: Path) -> Path | None:
    seven = shutil.which("7z") or shutil.which("7zz")
    if not seven:
        return None
    subprocess.run([seven, "x", "-y", f"-o{dest_dir}", str(path)], check=True)
    found = [item for item in dest_dir.iterdir() if item.suffix.lower() in {".csv", ".txt"}]
    if len(found) != 1:
        raise RuntimeError(f"{path.name} extracted {len(found)} tables")
    return found[0]


def _columns(con: duckdb.DuckDBPyConnection, path: Path) -> list[str]:
    described = con.execute(
        "DESCRIBE SELECT * FROM read_csv(?, header=true, all_varchar=true)",
        [str(path)],
    ).fetchall()
    return [row[0] for row in described]


def _quote(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def load_curated(path: Path | None = None) -> dict[str, dict[str, str]]:
    """DTP cell-line name to the curated row. Exact `dtp_name` is the lookup key."""
    source = path or CURATED_CELLS
    if not source.exists():
        return {}
    found: dict[str, dict[str, str]] = {}
    with source.open(newline="") as handle:
        for row in csv.DictReader(handle):
            name = clean(row.get("dtp_name"))
            cvcl = clean(row.get("cvcl"))
            if not name or not cvcl:
                continue
            found[name] = {key: (value or "") for key, value in row.items()}
    return found


def curated_decision(row: dict[str, str], cellosaurus: object | None) -> Decision:
    """Curated accession wins. A Cellosaurus problem flag is still attached."""
    cvcl = row["cvcl"]
    notes = row.get("notes") or f"curated {cvcl}"
    flag = None
    if cellosaurus is not None:
        primary = cellosaurus.primary(cvcl)
        if primary:
            if primary != cvcl:
                notes = f"{notes} Secondary {cvcl} -> {primary}."
            cvcl = primary
            flag = cellosaurus.flag(primary)
        else:
            notes = f"{notes} Accession was not in the loaded Cellosaurus file."
    return Decision(cvcl, "curated", "mapped", notes, flag)


def collect_models(
    con: duckdb.DuckDBPyConnection,
    tables: list[Path],
    cellosaurus: object | None,
    curated: dict[str, dict[str, str]] | None = None,
) -> list[dict[str, str | None]]:
    rows: dict[str, str | None] = {}
    for table in tables:
        columns = set(_columns(con, table))
        if "CELL_NAME" not in columns:
            continue
        panel = _quote("PANEL_NAME") if "PANEL_NAME" in columns else "NULL"
        for name, tissue in con.execute(
            f"""
            SELECT DISTINCT CELL_NAME, {panel}
            FROM read_csv(?, header=true, all_varchar=true)
            WHERE CELL_NAME IS NOT NULL AND CELL_NAME <> ''
            """,
            [str(table)],
        ).fetchall():
            token = clean(name)
            if token and token not in rows:
                rows[token] = clean(tissue)
    curated_rows = load_curated() if curated is None else curated
    out: list[dict[str, str | None]] = []
    for name, tissue in rows.items():
        chosen = curated_rows.get(name)
        if chosen:
            decision = curated_decision(chosen, cellosaurus)
        elif cellosaurus is not None:
            decision = cellosaurus.decide_names([name])
        else:
            decision = Decision(None, "unmatched", "unmatched", "no Cellosaurus index")
        out.extend(xref_rows(
            NCI_SOURCE,
            name,
            decision,
            [("cell_name", name)],
            display_name=name,
            tissue=tissue,
            species="Homo sapiens",
        ))
    return out


def collect_compounds(con: duckdb.DuckDBPyConnection, tables: list[Path]) -> list[CompoundRecord]:
    """One record per NSC that appears in an assay table. ``public_nscs.csv`` is only a number list."""
    found: set[str] = set()
    for table in tables:
        columns = set(_columns(con, table))
        if "NSC" not in columns:
            continue
        prefix = ""
        if "PREFIX" in columns:
            prefix = " AND (PREFIX IS NULL OR PREFIX = '' OR PREFIX = 'S')"
        for (value,) in con.execute(
            f"SELECT DISTINCT NSC FROM read_csv(?, header=true, all_varchar=true) WHERE NSC IS NOT NULL{prefix}",
            [str(table)],
        ).fetchall():
            token = nsc_id(value)
            if token:
                found.add(token)
    records: list[CompoundRecord] = []
    for token in sorted(found):
        records.append(CompoundRecord(
            source=NCI_SOURCE,
            source_compound_id=token,
            name=token,
            evidence="NSC number from the DTP assay file; no DTP structure file is pinned",
        ))
    return records


def assay_tables(raw_root: Path, pairs: list[tuple[dict, dict]]) -> list[tuple[dict, dict, Path]]:
    from pgx.process import file_for

    found: list[tuple[dict, dict, Path]] = []
    for source, spec in pairs:
        path = file_for(raw_root, source["id"], str(source["version"]), spec["name"])
        if not path:
            continue
        found.append((source, spec, open_table(path)))
    return found


def _model_join(xref_file: str) -> str:
    return f"""
    LEFT JOIN (
      SELECT source_model_id, any_value(cvcl) AS model_cvcl
      FROM read_parquet('{xref_file}')
      WHERE source = '{NCI_SOURCE}' AND alias_type = 'primary' AND status = 'mapped'
      GROUP BY source_model_id
    ) models ON models.source_model_id = base.source_model_id
    """


def _prefix_filter(columns: set[str]) -> str:
    if "PREFIX" not in columns:
        return ""
    return " AND (PREFIX IS NULL OR PREFIX = '' OR PREFIX = 'S')"


def _response_select(body: str, xref_file: str) -> str:
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
    {_model_join(xref_file)}
    WHERE base.value IS NOT NULL
      AND base.source_model_id IS NOT NULL
      AND base.source_compound_id IS NOT NULL
    """


def _endpoint_sql(metric: str, columns: set[str]) -> str:
    dose = "NULL::DOUBLE"
    unit = "CASE WHEN upper(CONCENTRATION_UNIT) = 'M' THEN 'log10_M' ELSE CONCENTRATION_UNIT END"
    value = "try_cast(AVERAGE AS DOUBLE)"
    return f"""
      SELECT
        CELL_NAME AS source_model_id,
        'NSC' || regexp_replace(CAST(NSC AS VARCHAR), '\\.0$', '') AS source_compound_id,
        '{metric}' AS metric,
        {value} AS value,
        {unit} AS unit,
        {dose} AS dose,
        EXPID AS screen_id
      FROM read_csv(?, header=true, all_varchar=true)
      WHERE CELL_NAME IS NOT NULL {_prefix_filter(columns)}
    """


def _oneconc_sql(columns: set[str]) -> str:
    value_col = "AVERAGE_GIPRCNT" if "AVERAGE_GIPRCNT" in columns else "AVERAGE"
    return f"""
      SELECT
        CELL_NAME AS source_model_id,
        'NSC' || regexp_replace(CAST(NSC AS VARCHAR), '\\.0$', '') AS source_compound_id,
        'ONECONC_GIPRCNT' AS metric,
        try_cast({value_col} AS DOUBLE) AS value,
        'percent' AS unit,
        try_cast(CONCENTRATION AS DOUBLE) AS dose,
        EXPID AS screen_id
      FROM read_csv(?, header=true, all_varchar=true)
      WHERE CELL_NAME IS NOT NULL {_prefix_filter(columns)}
    """


def _doseresp_sql(columns: set[str]) -> str:
    gip = "AVERAGE_GIPRCNT" if "AVERAGE_GIPRCNT" in columns else "GIPRCNT"
    ptc = "AVERAGE_PTC" if "AVERAGE_PTC" in columns else "PTC"
    dose = """
      CASE
        WHEN upper(CONCENTRATION_UNIT) = 'M' THEN pow(10, try_cast(CONCENTRATION AS DOUBLE))
        ELSE try_cast(CONCENTRATION AS DOUBLE)
      END
    """
    pieces = []
    if gip in columns:
        pieces.append(f"""
          SELECT CELL_NAME AS source_model_id,
                 'NSC' || regexp_replace(CAST(NSC AS VARCHAR), '\\.0$', '') AS source_compound_id,
                 'GIPRCNT' AS metric,
                 try_cast({gip} AS DOUBLE) AS value,
                 'percent' AS unit,
                 {dose} AS dose,
                 EXPID AS screen_id
          FROM read_csv(?, header=true, all_varchar=true)
          WHERE CELL_NAME IS NOT NULL {_prefix_filter(columns)}
        """)
    if ptc in columns:
        pieces.append(f"""
          SELECT CELL_NAME AS source_model_id,
                 'NSC' || regexp_replace(CAST(NSC AS VARCHAR), '\\.0$', '') AS source_compound_id,
                 'PTC' AS metric,
                 try_cast({ptc} AS DOUBLE) AS value,
                 'percent' AS unit,
                 {dose} AS dose,
                 EXPID AS screen_id
          FROM read_csv(?, header=true, all_varchar=true)
          WHERE CELL_NAME IS NOT NULL {_prefix_filter(columns)}
        """)
    if not pieces:
        raise RuntimeError("DOSERESP table has neither GIPRCNT nor PTC columns")
    return " UNION ALL ".join(pieces)


def load_responses(
    con: duckdb.DuckDBPyConnection,
    tables: list[tuple[str, str, Path]],
    xref: Path,
    dest_root: Path,
) -> int:
    """Write one parquet part per assay file under ``dataset=nci60``.

    ``tables`` is ``(part_name, role, csv_path)``. ``role`` is
    ``nci60_endpoint``, ``nci60_oneconc``, or ``nci60_doseresp``. Endpoint
    metrics are the part name without a suffix (``GI50``).
    """
    from pgx.process import write_part

    if not tables:
        return 0
    dataset_dir = dest_root / f"dataset={NCI_SOURCE}"
    if dataset_dir.exists():
        shutil.rmtree(dataset_dir)
    xref_file = (xref / "part.parquet").as_posix()
    total = 0
    for part_name, role, path in tables:
        columns = set(_columns(con, path))
        if role == "nci60_endpoint":
            metric = Path(part_name).stem
            if metric not in ENDPOINT_METRICS:
                raise RuntimeError(f"endpoint part {part_name} is not one of {sorted(ENDPOINT_METRICS)}")
            body = _endpoint_sql(metric, columns)
            params = [str(path)]
        elif role == "nci60_oneconc":
            body = _oneconc_sql(columns)
            params = [str(path)]
        elif role == "nci60_doseresp":
            body = _doseresp_sql(columns)
            params = [str(path)] * body.count("read_csv(")
        else:
            raise RuntimeError(f"unknown NCI60 role {role}")
        total += write_part(con, _response_select(body, xref_file), params, dest_root, NCI_SOURCE, f"{Path(part_name).stem}.parquet")
    return total
