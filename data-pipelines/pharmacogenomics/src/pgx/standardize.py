"""Add within-dataset standardized response columns to the drug_response partitions.

Raw metric values are only comparable inside one dataset and metric, and even then "larger"
means stronger for some metrics and weaker for others. After the partitions are written this
step adds these columns so a reader does not have to know which way each metric points:

``response_direction``
    ``lower_is_stronger``, ``higher_is_stronger`` or ``unknown`` (see ``pgx.metrics``).
``sensitivity_percentile``
    Position of the row among all rows of the same dataset, metric, response class and dose,
    from 0 (least sensitive) to 1 (most sensitive). It orders compounds against each other.
``compound_percentile``
    The same, but only among rows for the same compound, so it orders cell lines for one
    drug. It is null when the compound has fewer than ``MIN_GROUP`` rows in the group.

``beyond_tested_range``
    True when an IC50, EC50 or AC50 is larger than the highest dose tested for that model and
    compound, so it was extrapolated from an incomplete curve rather than measured. False when
    it lies inside the tested range; null when no raw doses were loaded or the metric is not one
    of those. The value is never changed or capped, and the row is not hidden.

Percentiles are null where the direction is unknown, so an unverified AUC is never silently
ranked the wrong way round. Ties share the lowest rank, and replicate rows count as rows.
The percentiles are descriptive; they are not a normalisation across datasets.
"""

from __future__ import annotations

from pathlib import Path

import duckdb

from pgx.metrics import HIGHER_IS_STRONGER, LOWER_IS_STRONGER, direction

MIN_GROUP = 10
STANDARDIZED_COLUMNS = ("response_direction", "sensitivity_percentile", "compound_percentile", "beyond_tested_range")
# Concentration-valued metrics reported in the same unit (uM) as the dose of the raw viability rows.
RANGE_METRICS = ("IC50", "EC50", "AC50")
_GROUP_COLUMNS = ("response_class", "dose", "anchor_conc", "library_conc")
_COMPOUND_COLUMNS = ("source_compound_id", "combination_id")


def _columns(con: duckdb.DuckDBPyConnection, table: str) -> set[str]:
    return {str(row[0]) for row in con.execute(f"SELECT column_name FROM (DESCRIBE {table})").fetchall()}


def _partition(columns: set[str], names: tuple[str, ...]) -> list[str]:
    return [f'"{name}"' for name in names if name in columns]


def standardize_dataset(con: duckdb.DuckDBPyConnection, directory: Path) -> int:
    """Standardize every parquet file of one ``dataset=<name>`` directory, in place."""
    files = sorted(directory.glob("*.parquet"))
    if not files:
        return 0
    dataset = directory.name.removeprefix("dataset=")
    names = [file.as_posix() for file in files]
    con.execute(
        "CREATE OR REPLACE TEMP TABLE std_rows AS "
        "SELECT * FROM read_parquet(?, union_by_name=true, filename=true, hive_partitioning=false)",
        [names],
    )
    columns = _columns(con, "std_rows")
    if not {"metric", "value"} <= columns:
        return 0
    metrics = [str(row[0]) for row in con.execute("SELECT DISTINCT metric FROM std_rows WHERE metric IS NOT NULL").fetchall()]
    con.execute("CREATE OR REPLACE TEMP TABLE std_direction(metric VARCHAR, direction VARCHAR)")
    if metrics:
        con.executemany("INSERT INTO std_direction VALUES (?, ?)", [(metric, direction(dataset, metric)) for metric in metrics])

    group = ", ".join(["metric", *_partition(columns, _GROUP_COLUMNS)])
    compound = ", ".join([group, *_partition(columns, _COMPOUND_COLUMNS)])
    stale = [f'"{name}"' for name in STANDARDIZED_COLUMNS if name in columns]
    exclude = f"EXCLUDE ({', '.join(stale)})" if stale else ""
    # strength grows with sensitivity whatever the metric's own direction is.
    strength = (
        f"CASE std_direction.direction WHEN '{LOWER_IS_STRONGER}' THEN -rows.value "
        f"WHEN '{HIGHER_IS_STRONGER}' THEN rows.value END"
    )
    keys = [name for name in ("source_model_id", "source_compound_id", "combination_id") if name in columns]
    has_range = bool(keys) and {"dose", "source_model_id", "source_compound_id"} <= columns
    key_list = ", ".join(f'"{name}"' for name in keys)
    key_join = " AND ".join(f'rows."{name}" IS NOT DISTINCT FROM tested."{name}"' for name in keys)
    metric_list = ", ".join(f"'{name}'" for name in RANGE_METRICS)
    tested_join = (
        f"""LEFT JOIN (
              SELECT {key_list}, max(dose) AS max_dose
              FROM std_rows WHERE metric = 'VIABILITY' AND dose IS NOT NULL GROUP BY {key_list}
            ) tested ON {key_join}"""
        if has_range
        else ""
    )
    computed = (
        f"CASE WHEN upper(rows.metric) IN ({metric_list}) AND isfinite(rows.value) AND tested.max_dose IS NOT NULL "
        "THEN rows.value > tested.max_dose * (1 + 1e-9) END"
        if has_range
        else "NULL::BOOLEAN"
    )
    # Once the raw dose points are dropped a rerun cannot recompute the flag, so it keeps the stored one.
    beyond = f"coalesce({computed}, rows.beyond_tested_range)" if "beyond_tested_range" in columns else computed
    con.execute(
        f"""
        CREATE OR REPLACE TEMP TABLE std_out AS
        WITH scored AS (
          SELECT rows.* {exclude},
                 coalesce(std_direction.direction, 'unknown') AS response_direction,
                 {strength} AS strength,
                 {beyond} AS beyond_tested_range
          FROM std_rows rows
          LEFT JOIN std_direction ON std_direction.metric = rows.metric
          {tested_join}
        ), ranked AS (
          SELECT scored.*,
                 percent_rank() OVER (PARTITION BY {group} ORDER BY strength) AS dataset_rank,
                 count(*) OVER (PARTITION BY {group}) AS dataset_size,
                 percent_rank() OVER (PARTITION BY {compound} ORDER BY strength) AS compound_rank,
                 count(*) OVER (PARTITION BY {compound}) AS compound_size
          FROM scored
        )
        SELECT * EXCLUDE (strength, dataset_rank, dataset_size, compound_rank, compound_size),
               CASE WHEN strength IS NOT NULL AND dataset_size >= {MIN_GROUP} THEN dataset_rank END AS sensitivity_percentile,
               CASE WHEN strength IS NOT NULL AND compound_size >= {MIN_GROUP} THEN compound_rank END AS compound_percentile
        FROM ranked
        """
    )
    for file, name in zip(files, names, strict=True):
        temporary = file.with_suffix(".parquet.next")
        con.execute(
            f"COPY (SELECT * EXCLUDE (filename) FROM std_out WHERE filename = ?) TO '{temporary.as_posix()}' (FORMAT PARQUET)",
            [name],
        )
        temporary.replace(file)
    count = con.execute("SELECT count(*) FROM std_out").fetchone()
    return int(count[0] if count else 0)


RAW_DOSE_FILE = "raw.parquet"


def drop_raw_dose_points(response_root: Path) -> int:
    """Delete the PharmacoSet raw dose-response files, which the served data does not need.

    ``extract_pset`` writes the per-dose viability points to ``raw.parquet`` next to the
    curve-level ``profiles.parquet``. They are about three quarters of the rows. The one thing
    read from them, the highest tested dose, is folded into ``beyond_tested_range`` by
    ``standardize_responses``, so call this after it. Other datasets never write this file name. Returns the number of files removed.
    """
    removed = 0
    for file in sorted(response_root.glob(f"dataset=*/{RAW_DOSE_FILE}")):
        file.unlink()
        removed += 1
    return removed


def standardize_responses(con: duckdb.DuckDBPyConnection, response_root: Path) -> dict[str, int]:
    """Standardize every dataset under ``drug_response/``; returns rows per dataset."""
    done: dict[str, int] = {}
    if not response_root.exists():
        return done
    for directory in sorted(response_root.glob("dataset=*")):
        if directory.is_dir():
            done[directory.name.removeprefix("dataset=")] = standardize_dataset(con, directory)
    return done
