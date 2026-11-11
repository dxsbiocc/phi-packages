import json

import pytest

from pgx.cli import main
from pgx.updates import compare, fetch_registry, parse_registry, record_ids, tracked_sources


def entry(name, dataset, doi="", link="", created="2021-01-01T00:00:00.000Z"):
    return {"name": name, "doi": doi, "downloadLink": link, "dateCreated": created, "dataset": {"name": dataset}}


REGISTRY = [
    # Tracked at the newest record of its family: current.
    entry("CTRPv2.1_2016", "CTRPv2", "10.5281/zenodo.21794049", "https://zenodo.org/records/21794049/files/x.rds", "2026-07-31T14:46:07.030Z"),
    entry("CTRPv2_2015", "CTRPv2", "10.5281/zenodo.7826870", "https://zenodo.org/record/7826870/files/x.rds", "2020-06-24T14:39:26.588Z"),
    # Tracked at an older record than ORCESTRA's newest: newer_available.
    entry("CCLE_2015", "CCLE", "10.5281/zenodo.3905461", "https://zenodo.org/record/3905462/files/CCLE.rds", "2020-06-24T14:39:26.588Z"),
    entry("CCLE_2019", "CCLE", "10.5281/zenodo.18759113", "https://zenodo.org/records/18759113/files/x.RDS", "2025-12-08T14:46:07.030Z"),
    # Variants are not newer versions.
    entry("CCLE_2015_unfiltered", "CCLE", "10.5281/zenodo.4721823", "", "2026-12-01T00:00:00.000Z"),
    # The configured record is not in the registry at all: record_not_listed.
    entry("gCSI_2019", "gCSI", "10.5281/zenodo.7829857", "https://zenodo.org/record/7829857/files/x.rds", "2021-06-11T21:58:16.390Z"),
    entry("gCSI_2018_unfiltered", "gCSI", "", "", "1970-01-01T00:00:00.000Z"),
    # Untracked datasets.
    entry("HMCL", "HMCL", "10.5281/zenodo.23920418", "https://zenodo.org/records/23920419/files/HMCL_PSet.rds", "2026-10-27T14:46:07.030Z"),
    entry("Tavor_2020", "Tavor", "10.5281/zenodo.5979590", "", "2021-03-05T17:05:08.535Z"),
    entry("NCISarcoma_2015", "NCISarcoma", "", "", ""),
    # Covered by non-ORCESTRA sources.
    entry("PRISM_2020", "PRISM", "10.5281/zenodo.7826864", "", "2021-08-18T16:28:45.207Z"),
    entry("GDSC_2020(v2-8.2)", "GDSC", "10.5281/zenodo.5787145", "", "2021-12-16T19:58:28.388Z"),
    entry("GDSC²_Anchor", "GDSC²", "10.5281/zenodo.23921232", "https://zenodo.org/records/23921233/files/a.rds", "2026-10-27T14:46:07.030Z"),
    # A planned source whose family now has a candidate.
    entry("NCI60_2026", "NCI60", "10.5281/zenodo.21043047", "https://zenodo.org/records/21043047/files/NCI60-2026_PSet.rds", "2025-06-26T14:46:07.030Z"),
]

DOCUMENT = {
    "sources": [
        {"id": "depmap", "enabled": True},
        {"id": "gdsc", "enabled": True},
        {"id": "gdsc_combinations", "enabled": True},
        {"id": "prism_19q4", "enabled": True},
        {"id": "orcestra_ctrpv2", "version": "2.1", "dataset": "ctrpv2", "homepage": "https://doi.org/10.5281/zenodo.21794049",
         "files": [{"url": "https://zenodo.org/records/21794049/files/sample.tsv?download=1"}]},
        {"id": "orcestra_ccle", "version": "2015", "dataset": "ccle", "homepage": "https://doi.org/10.5281/zenodo.3905461"},
        {"id": "orcestra_gcsi", "version": "2019", "dataset": "gcsi", "homepage": "https://doi.org/10.5281/zenodo.4737437"},
        {"id": "orcestra_gray", "version": "2017", "dataset": "gray", "homepage": "https://doi.org/10.5281/zenodo.3905458"},
    ],
    "planned": [{"id": "orcestra_nci60", "version": "2026", "enabled": False}],
}


def by_source(report):
    return {item["source"]: item for item in report["tracked"]}


def test_record_ids_read_dois_and_zenodo_urls():
    assert record_ids("https://doi.org/10.5281/zenodo.3905448") == {"3905448"}
    assert record_ids("https://zenodo.org/records/21794049/files/a.rds?download=1") == {"21794049"}
    assert record_ids("https://zenodo.org/record/7826870/files/a.rds", "10.5281/zenodo.7826870") == {"7826870"}
    assert record_ids(None, "https://example.org/a") == frozenset()


def test_registry_parsing_flags_variants_and_placeholder_dates():
    entries = {item.name: item for item in parse_registry(REGISTRY)}
    assert entries["CCLE_2015_unfiltered"].variant is True
    assert entries["CCLE_2019"].variant is False
    assert entries["gCSI_2018_unfiltered"].created is None, "the 1970 placeholder is not a date"
    assert entries["NCISarcoma_2015"].created is None
    assert entries["GDSC²_Anchor"].family == "gdsc2"
    assert entries["CTRPv2.1_2016"].record_ids == {"21794049"}
    assert parse_registry([{"nothing": True}, "x"]) == []
    with pytest.raises(ValueError):
        parse_registry({"not": "a list"})


def test_tracked_sources_include_planned_and_ignore_other_sources():
    tracked = {item.id: item for item in tracked_sources(DOCUMENT)}
    assert set(tracked) == {"orcestra_ctrpv2", "orcestra_ccle", "orcestra_gcsi", "orcestra_gray", "orcestra_nci60"}
    assert tracked["orcestra_ctrpv2"].record_ids == {"21794049"}
    assert tracked["orcestra_nci60"].planned is True
    assert tracked["orcestra_nci60"].family == "nci60"


def test_each_tracked_source_gets_an_honest_status():
    report = compare(DOCUMENT, parse_registry(REGISTRY))
    statuses = by_source(report)
    assert statuses["orcestra_ctrpv2"]["status"] == "current"
    assert statuses["orcestra_ctrpv2"]["newer"] == []

    assert statuses["orcestra_ccle"]["status"] == "newer_available"
    assert [item["name"] for item in statuses["orcestra_ccle"]["newer"]] == ["CCLE_2019"], "the unfiltered variant is not offered as newer"

    gcsi = statuses["orcestra_gcsi"]
    assert gcsi["status"] == "record_not_listed"
    assert [item["name"] for item in gcsi["candidates"]] == ["gCSI_2019"], "variants are left out of the candidates"

    assert statuses["orcestra_gray"]["status"] == "family_not_listed"

    planned = statuses["orcestra_nci60"]
    assert planned["status"] == "planned"
    assert planned["candidates"][0]["name"] == "NCI60_2026"
    assert planned["candidates"][0]["download_link"].endswith("NCI60-2026_PSet.rds")


def test_untracked_and_covered_entries_are_separated():
    report = compare(DOCUMENT, parse_registry(REGISTRY))
    assert [item["name"] for item in report["untracked"]] == ["HMCL", "Tavor_2020", "NCISarcoma_2015"]
    covered = {item["name"]: item["covered_by"] for item in report["covered_elsewhere"]}
    assert covered == {
        "GDSC²_Anchor": ["gdsc_combinations"],
        "PRISM_2020": ["prism_19q4"],
        "GDSC_2020(v2-8.2)": ["gdsc"],
    }
    assert report["registry"] == {"url": "https://www.orcestra.ca/api/pset/available", "entries": 14, "variants_ignored": 2}
    assert report["summary"]["newer_available"] == 1
    assert report["summary"]["untracked"] == 3
    assert report["summary"]["needs_attention"] == 4


def test_since_limits_untracked_entries_by_date():
    report = compare(DOCUMENT, parse_registry(REGISTRY), since="2025-06-01")
    assert [item["name"] for item in report["untracked"]] == ["HMCL"]
    assert report["summary"]["needs_attention"] == 2


def test_fetch_registry_uses_the_injected_fetcher():
    seen = []

    def fetch(url):
        seen.append(url)
        return json.dumps(REGISTRY).encode(), "application/json"

    assert len(fetch_registry(fetch)) == 14
    assert seen == ["https://www.orcestra.ca/api/pset/available"]


@pytest.fixture()
def pipeline(tmp_path):
    import yaml

    (tmp_path / "sources.yaml").write_text(yaml.safe_dump(DOCUMENT))
    (tmp_path / "registry.json").write_text(json.dumps(REGISTRY))
    return tmp_path


def test_cli_prints_the_report(pipeline, capsys):
    main(["--root", str(pipeline), "check-updates", "--registry-file", str(pipeline / "registry.json")])
    report = json.loads(capsys.readouterr().out)
    assert report["summary"]["needs_attention"] == 4


def test_cli_fail_on_updates_sets_the_exit_code(pipeline, capsys):
    registry = str(pipeline / "registry.json")
    with pytest.raises(SystemExit) as raised:
        main(["--root", str(pipeline), "check-updates", "--registry-file", registry, "--fail-on-updates"])
    assert raised.value.code == 1
    capsys.readouterr()
    # Nothing newer and nothing untracked: exit normally.
    quiet = [item for item in REGISTRY if item["dataset"]["name"] in {"CTRPv2", "gCSI", "NCI60"}]
    (pipeline / "quiet.json").write_text(json.dumps(quiet))
    main(["--root", str(pipeline), "check-updates", "--registry-file", str(pipeline / "quiet.json"), "--fail-on-updates"])


def test_cli_reports_an_unreadable_registry_with_exit_code_two(pipeline, capsys):
    with pytest.raises(SystemExit) as raised:
        main(["--root", str(pipeline), "check-updates", "--registry-file", str(pipeline / "missing.json")])
    assert raised.value.code == 2
    assert "could not read the ORCESTRA registry" in capsys.readouterr().err
