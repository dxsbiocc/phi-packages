from pathlib import Path
import json
import os
import subprocess
import sys
import tempfile
import unittest

import fixtures

SCRIPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "skills", "deseq2-qc", "scripts", "assess.py")


def run(results=None, qc=None, extra=()):
    with tempfile.TemporaryDirectory() as tmp:
        rpath = os.path.join(tmp, "results.tsv")
        qpath = os.path.join(tmp, "results.deseq2.qc.tsv")
        if results is not None:
            with open(rpath, "wb") as fh:
                fh.write(results)
        if qc is not None:
            with open(qpath, "wb") as fh:
                fh.write(qc)
        before = {n: (Path(tmp) / n).read_bytes() for n in os.listdir(tmp)}
        proc = subprocess.run([sys.executable, SCRIPT, "--results", rpath, "--qc", qpath, *extra],
                              capture_output=True, text=True)
        after = {n: (Path(tmp) / n).read_bytes() for n in os.listdir(tmp)}
    return proc, before, after


class Assess(unittest.TestCase):
    def decision(self, results, qc, disposition):
        proc, before, after = run(results, qc)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        obj = json.loads(proc.stdout)
        self.assertEqual(obj["status"], "ok")
        self.assertEqual(obj["disposition"], disposition)
        self.assertEqual(obj["outputs"], [])
        self.assertEqual(before, after)
        return obj, proc.stdout

    def test_eligible(self):
        obj, out = self.decision(fixtures.ELIGIBLE_RESULTS, fixtures.ELIGIBLE_QC, "eligible")
        self.assertEqual(obj["evaluated"], 3)
        self.assertEqual(obj["hits"], 2)
        self.assertEqual(obj["contrast"], "treated vs control")
        self.assertEqual(obj["padjThreshold"], 0.05)
        self.assertEqual(obj["log2fcThreshold"], 1.0)
        self.assertEqual(obj["resultMd5"], fixtures.md5(fixtures.ELIGIBLE_RESULTS))
        self.assertNotIn('"g1"', out)
        self.assertNotIn('"g2"', out)

    def test_strict_padj_boundary(self):
        obj, _ = self.decision(fixtures.BOUNDARY_RESULTS, fixtures.BOUNDARY_QC, "eligible")
        self.assertEqual(obj["evaluated"], 2)
        self.assertEqual(obj["hits"], 1)

    def test_skip(self):
        obj, _ = self.decision(fixtures.SKIP_RESULTS, fixtures.SKIP_QC, "skip")
        self.assertEqual(obj["evaluated"], 2)
        self.assertEqual(obj["hits"], 0)

    def test_invalid_all_na(self):
        obj, _ = self.decision(fixtures.ALLNA_RESULTS, fixtures.ALLNA_QC, "invalid")
        self.assertEqual(obj["evaluated"], 0)
        self.assertTrue(obj["reasons"])

    def test_hash_drift_invalid(self):
        obj, _ = self.decision(fixtures.ELIGIBLE_RESULTS + b"g6\t3\t0.001\n", fixtures.ELIGIBLE_QC, "invalid")
        self.assertTrue(any("md5" in r for r in obj["reasons"]))

    def test_count_drift_invalid(self):
        qc = fixtures.certificate(fixtures.ELIGIBLE_RESULTS, "eligible", 99, 2)
        self.decision(fixtures.ELIGIBLE_RESULTS, qc, "invalid")

    def test_disposition_drift_invalid(self):
        qc = fixtures.certificate(fixtures.ELIGIBLE_RESULTS, "skip", 3, 0)
        self.decision(fixtures.ELIGIBLE_RESULTS, qc, "invalid")

    def test_missing_qc_invalid_exit_zero(self):
        obj, _ = self.decision(fixtures.ELIGIBLE_RESULTS, None, "invalid")
        self.assertIsNone(obj["contrast"])

    def test_missing_results_invalid_exit_zero(self):
        self.decision(None, fixtures.ELIGIBLE_QC, "invalid")

    def test_malformed_qc_invalid(self):
        self.decision(fixtures.ELIGIBLE_RESULTS, b"not a certificate", "invalid")

    def test_read_only_no_sidecar_written(self):
        proc, before, after = run(fixtures.ELIGIBLE_RESULTS, fixtures.ELIGIBLE_QC)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(set(before), {"results.tsv", "results.deseq2.qc.tsv"})
        self.assertEqual(before, after)

    def test_write_options_rejected(self):
        proc, _, _ = run(fixtures.ELIGIBLE_RESULTS, fixtures.ELIGIBLE_QC, extra=("--out", "x.tsv"))
        self.assertEqual(proc.returncode, 2)


if __name__ == "__main__":
    unittest.main()