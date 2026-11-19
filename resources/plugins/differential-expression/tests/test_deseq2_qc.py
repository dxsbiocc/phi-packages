"""Tests for the deseq2-qc skill: gate thresholds, corruption, ORA math.

Run from the plugin root: python -m unittest discover -s tests -v
Fixtures are tiny synthetic tables; no real data or external catalogs.
"""

import json
import math
import os
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction

PLUGIN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(PLUGIN_DIR, "skills", "deseq2-qc", "scripts")
sys.path.insert(0, SCRIPTS)
import deqclib as lib

class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def path(self, name):
        return os.path.join(self.tmp.name, name)

    def write(self, name, text):
        p = self.path(name)
        with open(p, "w", encoding="utf-8", newline="") as fh:
            fh.write(text)
        return p

    def run_script(self, name, **args):
        cmd = [sys.executable, os.path.join(SCRIPTS, name + ".py")]
        for key, value in args.items():
            cmd += ["--" + key, str(value)]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        self.assertIn(proc.returncode, (0, 2), proc.stderr)
        lines = proc.stdout.strip().splitlines()
        self.assertEqual(len(lines), 1, "stdout must be one compact JSON line")
        obj = json.loads(lines[0])
        self.assertNotIn("\n", proc.stdout.strip())
        return proc.returncode, obj


class TestMath(unittest.TestCase):
    def test_bh_known_vector(self):
        q = lib.bh_adjust([0.01, 0.04, 0.03, 0.002])
        for got, want in zip(q, [0.02, 0.04, 0.04, 0.008]):
            self.assertAlmostEqual(got, want, places=12)

    def test_hypergeom_handcomputed(self):
        self.assertAlmostEqual(lib.hypergeom_sf_ge(2, 10, 4, 3), 1 / 3, places=12)
        self.assertAlmostEqual(lib.hypergeom_sf_ge(3, 10, 4, 3), 1 / 30, places=12)
        self.assertEqual(lib.hypergeom_sf_ge(0, 10, 4, 3), 1.0)
        self.assertEqual(lib.hypergeom_sf_ge(4, 10, 4, 3), 0.0)

    def test_hypergeom_exhaustive_small(self):
        for N in range(1, 13):
            for K in range(N + 1):
                for n in range(N + 1):
                    hi = min(K, n)
                    lo = max(0, n - (N - K))
                    for k in range(0, hi + 2):
                        exact = sum(
                            (Fraction(math.comb(K, i) * math.comb(N - K, n - i), math.comb(N, n))
                             for i in range(max(k, lo), hi + 1)),
                            Fraction(0))
                        self.assertAlmostEqual(lib.hypergeom_sf_ge(k, N, K, n),
                                               float(exact), places=9,
                                               msg=(N, K, n, k))

    def test_hypergeom_large_universe_stable(self):
        N, K, n = 200000, 3000, 400
        p8 = lib.hypergeom_sf_ge(8, N, K, n)
        p9 = lib.hypergeom_sf_ge(9, N, K, n)
        p10 = lib.hypergeom_sf_ge(10, N, K, n)
        self.assertTrue(0.0 < p10 <= p9 <= p8 < 1.0)
        self.assertEqual(lib.hypergeom_sf_ge(0, N, K, n), 1.0)
        self.assertEqual(lib.hypergeom_sf_ge(n + 1, N, K, n), 0.0)


class TestORA(Base):
    def test_real_cli_blocks_damage_and_preserves_valid_zero(self):
        import fixtures
        self.gmt = self.write("sets.gmt", "signal\tdesc\tg1\tg2\nnull\tdesc\tg3\n")
        qc = self.path("qc.tsv")
        results = self.path("results.tsv")
        out = self.path("ora.tsv")
        def invoke(data, report):
            with open(results, "wb") as stream:
                stream.write(data)
            with open(qc, "wb") as stream:
                stream.write(report)
            code, decision = self.run_script("ora", results=results, qc=qc, gmt=self.gmt, out=out, min_size=1)
            self.assertEqual(code, 0)
            return decision
        self.assertEqual(invoke(fixtures.SKIP_RESULTS, fixtures.SKIP_QC)["status"], "skip")
        self.assertEqual(invoke(fixtures.SKIP_RESULTS + b"broken", fixtures.SKIP_QC)["status"], "invalid")
        self.assertFalse(os.path.exists(out))
        self.assertEqual(invoke(fixtures.ELIGIBLE_RESULTS, fixtures.ELIGIBLE_QC)["terms_tested"], 2)
        with open(out, encoding="utf-8") as stream:
            rows = stream.read().splitlines()
        self.assertEqual(len(rows), 3)
        self.assertAlmostEqual(float(rows[1].split("\t")[4]), 1 / 3, places=6)
        self.assertAlmostEqual(float(rows[1].split("\t")[5]), 2 / 3, places=6)
        self.write("ora.tsv", "earlier run")
        self.assertEqual(invoke(fixtures.ELIGIBLE_RESULTS, fixtures.ELIGIBLE_QC)["status"], "blocked")
        with open(out) as stream:
            self.assertEqual(stream.read(), "earlier run")
