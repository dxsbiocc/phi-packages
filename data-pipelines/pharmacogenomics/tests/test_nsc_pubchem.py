from pgx.ids import CompoundRecord
from pgx.nsc_pubchem import apply_deposits, apply_synonyms, nsc_digits, scan_sid_map, summarize, synonym_cids
from pgx.pubchem import resolve_records


SID_MAP = """1\tMolecular Imaging Database (MOLI)\tMOLI000002\t135398523
10\tDTP/NCI\tNSC740\t11433190
11\tDTP/NCI\t741\t1
12\tDTP/NCI\tNSC-741\t2
13\tDTP/NCI\tNSC 742\t
14\tOther Source\tNSC740\t999
15\tDTP/NCI\tnot-an-nsc\t5
16\tDTP/NCI\tNSC000743\t42
17\tNCI-DTP-MPB-TVSL\tNSC740\t777
"""

SYNONYMS = """11433190\tNSC740
8\tNSC741
42\tNSC742
7\tNSC744
9\tNSC-744
"""


def _record(nsc: str) -> CompoundRecord:
    return CompoundRecord(source="nci60", source_compound_id=nsc, name=nsc)


def test_sid_map_keeps_only_dtp_nci_and_ambiguous_cids_stay_ambiguous(tmp_path):
    path = tmp_path / "SID-Map.txt"
    path.write_text(SID_MAP)
    dest = tmp_path / "dtp_nci.tsv"
    stats: dict = {}
    deposits = scan_sid_map(path, dest, stats)
    assert stats["skipped_registry"] == ["not-an-nsc"]
    assert deposits["740"][1] == {"11433190"}
    assert deposits["741"][1] == {"1", "2"}
    assert deposits["742"][0] == {"13"}
    assert deposits["742"][1] == set()
    assert deposits["743"][1] == {"42"}
    assert "777" not in deposits["740"][1]
    assert nsc_digits("NSC-000740") == "740"

    records = [_record("NSC740"), _record("NSC741"), _record("NSC742"), _record("NSC743"), _record("NSC744")]
    records.append(CompoundRecord(source="gdsc", source_compound_id="740", name="other"))
    assert apply_deposits(records, deposits) == 3
    by_id = {record.source_compound_id: record for record in records}
    assert by_id["NSC740"].status == "mapped"
    assert by_id["NSC740"].pubchem_cid == "11433190"
    assert by_id["NSC740"].match_method == "pubchem_dtp_nci"
    assert by_id["NSC740"].confidence == "high"
    assert by_id["NSC741"].status == "ambiguous"
    assert by_id["NSC741"].pubchem_cid is None
    assert by_id["NSC741"].match_method == "pubchem_dtp_nci"
    assert by_id["NSC742"].status == "unmatched"
    assert by_id["NSC743"].pubchem_cid == "42"
    assert by_id["740"].status == "unmatched"

    synonyms = tmp_path / "CID-Synonym-filtered.txt"
    synonyms.write_text(SYNONYMS)
    found = apply_synonyms(records, synonym_cids(synonyms, {"740", "741", "742", "744"}))
    assert found == 2
    assert by_id["NSC741"].status == "ambiguous"
    assert by_id["NSC741"].match_method == "pubchem_dtp_nci"
    assert by_id["NSC742"].status == "mapped"
    assert by_id["NSC742"].pubchem_cid == "42"
    assert by_id["NSC742"].match_method == "pubchem_nsc"
    assert by_id["NSC744"].status == "ambiguous"
    assert by_id["NSC744"].match_method == "pubchem_nsc"
    assert by_id["NSC744"].pubchem_cid is None
    summary = summarize(records)
    assert summary["rows"] == 5
    assert summary["mapped"] == 3
    assert summary["ambiguous"] == 2
    assert summary["pubchem_dtp_nci"] == 3
    assert summary["pubchem_nsc"] == 2


def test_resolve_records_uses_the_crosswalk_before_name_queries(tmp_path, monkeypatch):
    cache = tmp_path / "pubchem"
    cache.mkdir()
    (cache / "dtp_nci.tsv").write_text("nsc,sid,cid\n740,10,11433190\n")
    called = []

    class Stop(Exception):
        pass

    def explode(*_args, **_kwargs):
        called.append(True)
        raise Stop()

    monkeypatch.setattr("pgx.normalize.compounds.resolve_compounds", explode)
    records = [_record("NSC740")]
    try:
        resolve_records(records, cache)
    except Stop:
        pass
    assert records[0].status == "mapped"
    assert records[0].pubchem_cid == "11433190"
    assert called == [True]
