from pathlib import Path
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "skills", "deseq2-qc", "scripts"))
import deqclib as lib
import fixtures


class EvaluateResults(unittest.TestCase):
    def test_eligible_counts(self):
        v = lib.evaluate_results(fixtures.ELIGIBLE_RESULTS, "gene", 0.05, 1.0)
        self.assertEqual(v["disposition"], "eligible")
        self.assertEqual(v["evaluated"], 3)
        self.assertEqual(v["hits"], 2)
        self.assertEqual(v["hit_ids"], {"g1", "g2"})
        self.assertEqual(v["evaluable_ids"], {"g1", "g2", "g3"})

    def test_padj_threshold_is_strict(self):
        v = lib.evaluate_results(fixtures.BOUNDARY_RESULTS, "gene", 0.05, 1.0)
        self.assertEqual(v["evaluated"], 2)
        self.assertEqual(v["hits"], 1)
        self.assertEqual(v["hit_ids"], {"g2"})

    def test_lfc_threshold_is_inclusive(self):
        data = b"gene\tlog2FoldChange\tpadj\ng1\t1\t0.04\ng2\t-1\t0.04\n"
        v = lib.evaluate_results(data, "gene", 0.05, 1.0)
        self.assertEqual(v["hits"], 2)

    def test_all_na_is_invalid(self):
        v = lib.evaluate_results(fixtures.ALLNA_RESULTS, "gene", 0.05, 1.0)
        self.assertEqual(v["disposition"], "invalid")
        self.assertEqual(v["evaluated"], 0)
        self.assertEqual(v["hits"], 0)

    def test_skip_when_no_hits(self):
        v = lib.evaluate_results(fixtures.SKIP_RESULTS, "gene", 0.05, 1.0)
        self.assertEqual(v["disposition"], "skip")
        self.assertEqual(v["evaluated"], 2)
        self.assertEqual(v["hits"], 0)

    def test_duplicate_gene_ids_invalid_without_leaking_ids(self):
        data = b"gene\tlog2FoldChange\tpadj\ng1\t2\t0.001\ng1\t1\t0.01\n"
        v = lib.evaluate_results(data, "gene", 0.05, 1.0)
        self.assertEqual(v["disposition"], "invalid")
        self.assertTrue(v["reasons"])
        for reason in v["reasons"]:
            self.assertNotIn("g1", reason)

    def test_blank_nan_and_out_of_range_cells_invalid(self):
        for cell in (b"", b"NaN", b"nan", b"1.5", b"-0.1", b"abc", b"1e999"):
            data = b"gene\tlog2FoldChange\tpadj\ng1\t2\t" + cell + b"\n"
            v = lib.evaluate_results(data, "gene", 0.05, 1.0)
            self.assertEqual(v["disposition"], "invalid", cell)

    def test_missing_required_column_invalid(self):
        data = b"gene\tlog2FoldChange\ng1\t2\n"
        v = lib.evaluate_results(data, "gene", 0.05, 1.0)
        self.assertEqual(v["disposition"], "invalid")

    def test_missing_gene_id_column_invalid(self):
        v = lib.evaluate_results(fixtures.ELIGIBLE_RESULTS, "nope", 0.05, 1.0)
        self.assertEqual(v["disposition"], "invalid")

    def test_empty_file_invalid(self):
        self.assertEqual(lib.evaluate_results(b"", "gene", 0.05, 1.0)["disposition"], "invalid")

    def test_unnamed_first_column(self):
        v = lib.evaluate_results(fixtures.UNNAMED_RESULTS, None, 0.05, 1.0)
        self.assertEqual(v["disposition"], "eligible")
        self.assertEqual(v["geneIdColumn"], lib.UNNAMED_FIRST)


class ParseCertificate(unittest.TestCase):
    def test_valid(self):
        cert = lib.parse_certificate(fixtures.ELIGIBLE_QC)
        self.assertEqual(cert["disposition"], "eligible")
        self.assertEqual(cert["evaluated"], 3)
        self.assertEqual(cert["hits"], 2)
        self.assertEqual(cert["contrast"], "treated vs control")
        self.assertEqual(cert["resultMd5"], fixtures.md5(fixtures.ELIGIBLE_RESULTS))

    def test_hits_exceeding_evaluated_rejected(self):
        bad = fixtures.certificate(fixtures.ELIGIBLE_RESULTS, "eligible", 2, 3)
        with self.assertRaises(lib.CertificateError):
            lib.parse_certificate(bad)

    def test_disposition_count_mismatch_rejected(self):
        bad = fixtures.certificate(fixtures.ELIGIBLE_RESULTS, "skip", 3, 2)
        with self.assertRaises(lib.CertificateError):
            lib.parse_certificate(bad)

    def test_malformed_rejected(self):
        for bad in (b"", b"garbage", fixtures.ELIGIBLE_QC.replace(b"\t", b","),
                    fixtures.ELIGIBLE_QC + b"extra\trow\n"):
            with self.assertRaises(lib.CertificateError):
                lib.parse_certificate(bad)


class VerifiedEvidence(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.results = os.path.join(self.tmp.name, "results.tsv")
        self.qc = os.path.join(self.tmp.name, "results.deseq2.qc.tsv")

    def write(self, results=None, qc=None):
        if results is not None:
            with open(self.results, "wb") as fh:
                fh.write(results)
        if qc is not None:
            with open(self.qc, "wb") as fh:
                fh.write(qc)

    def test_eligible(self):
        self.write(fixtures.ELIGIBLE_RESULTS, fixtures.ELIGIBLE_QC)
        ev = lib.verified_evidence(self.results, self.qc)
        self.assertEqual(ev["disposition"], "eligible")
        self.assertEqual(ev["evaluated"], 3)
        self.assertEqual(ev["hits"], 2)
        self.assertEqual(ev["hit_ids"], {"g1", "g2"})
        self.assertEqual(ev["contrast"], "treated vs control")
        self.assertEqual(ev["padjThreshold"], 0.05)
        self.assertEqual(ev["log2fcThreshold"], 1.0)
        self.assertEqual(ev["resultMd5"], fixtures.md5(fixtures.ELIGIBLE_RESULTS))

    def test_boundary_strict_padj(self):
        self.write(fixtures.BOUNDARY_RESULTS, fixtures.BOUNDARY_QC)
        ev = lib.verified_evidence(self.results, self.qc)
        self.assertEqual(ev["disposition"], "eligible")
        self.assertEqual(ev["hits"], 1)

    def test_skip_and_invalid_dispositions(self):
        self.write(fixtures.SKIP_RESULTS, fixtures.SKIP_QC)
        self.assertEqual(lib.verified_evidence(self.results, self.qc)["disposition"], "skip")
        self.write(fixtures.ALLNA_RESULTS, fixtures.ALLNA_QC)
        self.assertEqual(lib.verified_evidence(self.results, self.qc)["disposition"], "invalid")

    def test_hash_drift_rejected(self):
        self.write(fixtures.ELIGIBLE_RESULTS + b"g6\t3\t0.001\n", fixtures.ELIGIBLE_QC)
        ev = lib.verified_evidence(self.results, self.qc)
        self.assertEqual(ev["disposition"], "invalid")
        self.assertTrue(any("md5" in r for r in ev["reasons"]))

    def test_count_drift_rejected(self):
        self.write(fixtures.ELIGIBLE_RESULTS, fixtures.certificate(fixtures.ELIGIBLE_RESULTS, "eligible", 99, 2))
        self.assertEqual(lib.verified_evidence(self.results, self.qc)["disposition"], "invalid")

    def test_disposition_drift_rejected(self):
        self.write(fixtures.ELIGIBLE_RESULTS, fixtures.certificate(fixtures.ELIGIBLE_RESULTS, "skip", 3, 0))
        self.assertEqual(lib.verified_evidence(self.results, self.qc)["disposition"], "invalid")

    def test_missing_files_invalid(self):
        ev = lib.verified_evidence(self.results, self.qc)
        self.assertEqual(ev["disposition"], "invalid")
        self.write(fixtures.ELIGIBLE_RESULTS, None)
        self.assertEqual(lib.verified_evidence(self.results, self.qc)["disposition"], "invalid")

    def test_oversized_qc_rejected(self):
        self.write(fixtures.ELIGIBLE_RESULTS, b"x" * (lib.MAX_QC_BYTES + 1))
        self.assertEqual(lib.verified_evidence(self.results, self.qc)["disposition"], "invalid")

    def test_unnamed_gene_column(self):
        self.write(fixtures.UNNAMED_RESULTS, fixtures.UNNAMED_QC)
        ev = lib.verified_evidence(self.results, self.qc)
        self.assertEqual(ev["disposition"], "eligible")
        self.assertEqual(ev["geneIdColumn"], lib.UNNAMED_FIRST)
        self.assertEqual(ev["hit_ids"], {"g1"})

    def test_read_only(self):
        self.write(fixtures.ELIGIBLE_RESULTS, fixtures.ELIGIBLE_QC)
        before = {n: (Path(self.tmp.name) / n).read_bytes() for n in os.listdir(self.tmp.name)}
        lib.verified_evidence(self.results, self.qc)
        after = {n: (Path(self.tmp.name) / n).read_bytes() for n in os.listdir(self.tmp.name)}
        self.assertEqual(before, after)


class StatsHelpers(unittest.TestCase):
    def test_bh_adjust(self):
        self.assertEqual(lib.bh_adjust([]), [])
        q = lib.bh_adjust([0.01, 0.04, 0.03, 0.002])
        for got, want in zip(q, [0.02, 0.04, 0.04, 0.008]):
            self.assertAlmostEqual(got, want)
        self.assertEqual(lib.bh_adjust([0.5, 0.5]), [0.5, 0.5])

    def test_hypergeom_sf_ge(self):
        self.assertAlmostEqual(lib.hypergeom_sf_ge(5, 10, 5, 5), 1.0 / 252.0)
        self.assertEqual(lib.hypergeom_sf_ge(0, 10, 5, 5), 1.0)
        self.assertEqual(lib.hypergeom_sf_ge(6, 10, 5, 5), 0.0)
        self.assertAlmostEqual(lib.hypergeom_sf_ge(1, 10, 5, 5), 1.0 - 1.0 / 252.0)


if __name__ == "__main__":
    unittest.main()