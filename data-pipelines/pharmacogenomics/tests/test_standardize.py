import re
from pathlib import Path

import duckdb
import pytest

from pgx.metrics import POTENCY, VIABILITY_AUC_DATASETS, direction
from pgx.standardize import MIN_GROUP, drop_raw_dose_points, standardize_responses

REPO = Path(__file__).resolve().parents[3]


def write_rows(root, dataset, name, rows):
    """rows: (metric, source_compound_id, model, value, dose)"""
    directory = root / f"dataset={dataset}"
    directory.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute("CREATE TABLE r(metric VARCHAR, source_compound_id VARCHAR, model_cvcl VARCHAR, value DOUBLE, dose DOUBLE, response_class VARCHAR)")
    con.executemany("INSERT INTO r VALUES (?, ?, ?, ?, ?, 'single')", rows)
    con.execute(f"COPY r TO '{(directory / name).as_posix()}' (FORMAT PARQUET)")
    con.close()


def read(root, dataset):
    con = duckdb.connect()
    rows = con.execute(
        "SELECT metric, source_compound_id, model_cvcl, value, response_direction, sensitivity_percentile, compound_percentile "
        "FROM read_parquet(?, hive_partitioning=false) ORDER BY metric, source_compound_id, value",
        [str(root / f"dataset={dataset}" / "*.parquet")],
    ).fetchall()
    con.close()
    return rows


def test_directions_follow_the_documented_semantics():
    assert direction("gdsc1", "LN_IC50") == "lower_is_stronger"
    assert direction("nci60", "gi50") == "lower_is_stronger"
    assert direction("prism_primary_19q4", "LFC") == "lower_is_stronger"
    assert direction("ctrpv2", "AAC") == "higher_is_stronger"
    assert direction("gdsc2", "AUC") == "lower_is_stronger"
    assert direction("ctrpv2", "AUC") == "lower_is_stronger", "CTRPv2.1 AUC is 100 minus AAC"
    assert direction("gcsi", "AUC") == "unknown", "an unchecked PharmacoSet AUC is not assumed to point either way"
    assert direction("gdsc_combination_x", "SYNERGY_DELTA_AUC") == "unknown"
    assert direction(None, None) == "unknown"


def test_python_directions_match_the_mcp_server():
    source = (REPO / "mcp-servers/pharmacogenomics/src/metrics.ts").read_text()

    def names(constant):
        body = re.search(constant + r"\s*=\s*new Set\(\[(.*?)\]\)", source, re.S)
        assert body, constant
        return set(re.findall(r'"([^"]+)"', body.group(1)))

    assert names("POTENCY") == set(POTENCY)
    assert names("VIABILITY_AUC_DATASETS") == set(VIABILITY_AUC_DATASETS)


def test_percentile_one_means_most_sensitive_whichever_way_the_metric_points(tmp_path):
    root = tmp_path / "drug_response"
    count = MIN_GROUP
    # LN_IC50 falls with sensitivity; AAC rises. Both list the same cell line as most sensitive.
    write_rows(root, "gdsc1", "part.parquet", [("LN_IC50", "d1", f"CVCL_{i:04d}", float(i), None) for i in range(count)])
    write_rows(root, "ctrpv2", "part.parquet", [("AAC", "d1", f"CVCL_{i:04d}", float(count - i), None) for i in range(count)])
    con = duckdb.connect()
    assert standardize_responses(con, root) == {"ctrpv2": count, "gdsc1": count}
    for dataset in ("gdsc1", "ctrpv2"):
        rows = {row[2]: row for row in read(root, dataset)}
        assert rows["CVCL_0000"][5] == pytest.approx(1.0), dataset
        assert rows[f"CVCL_{count - 1:04d}"][5] == pytest.approx(0.0), dataset
        assert rows["CVCL_0000"][6] == pytest.approx(1.0), "one compound, so the compound percentile matches"


def test_unknown_direction_is_not_ranked(tmp_path):
    root = tmp_path / "drug_response"
    write_rows(root, "gcsi", "part.parquet", [("AUC", "d1", f"CVCL_{i:04d}", float(i), None) for i in range(MIN_GROUP)])
    standardize_responses(duckdb.connect(), root)
    rows = read(root, "gcsi")
    assert {row[4] for row in rows} == {"unknown"}
    assert all(row[5] is None and row[6] is None for row in rows)


def test_compound_percentile_needs_enough_models_and_dose_groups_stay_apart(tmp_path):
    root = tmp_path / "drug_response"
    big = [("VIABILITY", "big", f"CVCL_{i:04d}", float(i), 1.0) for i in range(MIN_GROUP)]
    small = [("VIABILITY", "small", f"CVCL_{i:04d}", float(i), 1.0) for i in range(3)]
    other_dose = [("VIABILITY", "big", f"CVCL_{i:04d}", float(i) + 100, 10.0) for i in range(MIN_GROUP)]
    write_rows(root, "prism_secondary_19q4", "raw.parquet", big + small + other_dose)
    standardize_responses(duckdb.connect(), root)
    rows = read(root, "prism_secondary_19q4")
    small_rows = [row for row in rows if row[1] == "small"]
    assert all(row[6] is None for row in small_rows), "three models are too few for a compound percentile"
    assert all(row[5] is not None for row in small_rows), "the dataset percentile still applies"
    # The 10 uM rows are all larger, but they are ranked only against each other.
    ten = sorted(row[6] for row in rows if row[3] >= 100)
    assert ten == pytest.approx([i / (MIN_GROUP - 1) for i in range(MIN_GROUP)])


def test_multiple_files_of_one_dataset_are_ranked_together_and_keep_their_layout(tmp_path):
    root = tmp_path / "drug_response"
    write_rows(root, "nci60", "GI50.parquet", [("GI50", "d1", f"CVCL_{i:04d}", float(i), None) for i in range(6)])
    write_rows(root, "nci60", "TGI.parquet", [("TGI", "d1", f"CVCL_{i:04d}", float(i), None) for i in range(MIN_GROUP)])
    standardize_responses(duckdb.connect(), root)
    assert sorted(path.name for path in (root / "dataset=nci60").iterdir()) == ["GI50.parquet", "TGI.parquet"]
    gi50 = [row for row in read(root, "nci60") if row[0] == "GI50"]
    assert len(gi50) == 6
    assert all(row[5] is None for row in gi50), "six rows are below the minimum group size, and TGI rows do not count towards it"


def test_rerunning_does_not_duplicate_columns(tmp_path):
    root = tmp_path / "drug_response"
    write_rows(root, "gdsc1", "part.parquet", [("LN_IC50", "d1", f"CVCL_{i:04d}", float(i), None) for i in range(MIN_GROUP)])
    con = duckdb.connect()
    standardize_responses(con, root)
    before = read(root, "gdsc1")
    standardize_responses(con, root)
    assert read(root, "gdsc1") == before
    columns = [row[0] for row in con.execute("DESCRIBE SELECT * FROM read_parquet(?)", [str(root / "dataset=gdsc1" / "part.parquet")]).fetchall()]
    assert columns.count("sensitivity_percentile") == 1
    assert "filename" not in columns


def test_potency_beyond_the_highest_tested_dose_is_flagged_but_left_alone(tmp_path):
    root = tmp_path / "drug_response"
    directory = root / "dataset=ctrpv2"
    directory.mkdir(parents=True)
    con = duckdb.connect()
    con.execute(
        "CREATE TABLE r(source_model_id VARCHAR, source_compound_id VARCHAR, metric VARCHAR, value DOUBLE, dose DOUBLE, response_class VARCHAR)"
    )
    rows = []
    for model in ("m1", "m2", "m3"):
        for dose in (0.1, 1.0, 10.0):
            rows.append((model, "d1", "VIABILITY", 50.0, dose, "single"))
    rows += [
        ("m1", "d1", "IC50", 2.0, None, "single"),  # inside the tested range
        ("m2", "d1", "IC50", 1e6, None, "single"),  # extrapolated far past 10
        ("m3", "d1", "AUC", 90.0, None, "single"),  # not a concentration
        ("m4", "d1", "IC50", 5.0, None, "single"),  # no raw doses for this model
        ("m2", "d1", "AC50", float("inf"), None, "single"),  # non-finite is not judged
    ]
    con.executemany("INSERT INTO r VALUES (?, ?, ?, ?, ?, ?)", rows)
    con.execute(f"COPY r TO '{(directory / 'part.parquet').as_posix()}' (FORMAT PARQUET)")
    standardize_responses(con, root)
    result = {
        (model, metric): (value, beyond)
        for model, metric, value, beyond in con.execute(
            "SELECT source_model_id, metric, value, beyond_tested_range FROM read_parquet(?) WHERE metric <> 'VIABILITY'",
            [str(directory / "part.parquet")],
        ).fetchall()
    }
    assert result[("m1", "IC50")] == (2.0, False)
    assert result[("m2", "IC50")] == (1e6, True), "the value stays as it was reported"
    assert result[("m3", "AUC")][1] is None
    assert result[("m4", "IC50")][1] is None
    assert result[("m2", "AC50")][1] is None
    assert con.execute(
        "SELECT count(*) FROM read_parquet(?) WHERE metric = 'VIABILITY' AND beyond_tested_range IS NOT NULL", [str(directory / "part.parquet")]
    ).fetchone() == (0,)


def test_raw_dose_points_can_be_dropped_after_the_flag_is_written(tmp_path):
    root = tmp_path / "drug_response"
    directory = root / "dataset=ctrpv2"
    directory.mkdir(parents=True)
    con = duckdb.connect()
    con.execute("CREATE TABLE raw(source_model_id VARCHAR, source_compound_id VARCHAR, metric VARCHAR, value DOUBLE, dose DOUBLE, response_class VARCHAR)")
    con.execute("INSERT INTO raw VALUES ('m1', 'd1', 'VIABILITY', 40.0, 1.0, 'single'), ('m1', 'd1', 'VIABILITY', 60.0, 10.0, 'single')")
    con.execute(f"COPY raw TO '{(directory / 'raw.parquet').as_posix()}' (FORMAT PARQUET)")
    con.execute("CREATE TABLE prof(source_model_id VARCHAR, source_compound_id VARCHAR, metric VARCHAR, value DOUBLE, dose DOUBLE, response_class VARCHAR)")
    con.execute("INSERT INTO prof VALUES ('m1', 'd1', 'IC50', 500.0, NULL, 'single')")
    con.execute(f"COPY prof TO '{(directory / 'profiles.parquet').as_posix()}' (FORMAT PARQUET)")
    standardize_responses(con, root)
    assert drop_raw_dose_points(root) == 1
    assert [path.name for path in directory.iterdir()] == ["profiles.parquet"]
    flag = "SELECT beyond_tested_range FROM read_parquet(?)"
    assert con.execute(flag, [str(directory / "profiles.parquet")]).fetchall() == [(True,)]
    standardize_responses(con, root)  # a rerun has no doses left but keeps what was written
    assert con.execute(flag, [str(directory / "profiles.parquet")]).fetchall() == [(True,)]
