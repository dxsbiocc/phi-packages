import json
import os
import subprocess
import sys
import tempfile
import unittest

SCRIPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "skills", "deseq2-qc", "scripts", "preflight.py")

INDEP_META = ("sample\tcondition\tdonor\n"
              "s1\ttreated\td1\ns2\ttreated\td2\n"
              "s3\tcontrol\td3\ns4\tcontrol\td4\n")

PAIRED_META = ("sample\tcondition\tdonor\n"
               "s1\ttreated\td1\ns2\tcontrol\td1\n"
               "s3\ttreated\td2\ns4\tcontrol\td2\n"
               "s5\ttreated\td3\ns6\tcontrol\td3\n")

NUM_META = ("sample\tcondition\n"
            "001\t1\n002\t1\n003\t2\n004\t2\n")


def run(metadata, design="independent", reference="control", target="treated",
        donor_col=None, sample_id_col="sample", condition_col="condition"):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "meta.tsv")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(metadata)
        cmd = [sys.executable, SCRIPT, "--metadata", path, "--design", design,
               "--sample_id_col", sample_id_col, "--condition_col", condition_col,
               "--reference_level", reference, "--target_level", target]
        if donor_col is not None:
            cmd += ["--donor_col", donor_col]
        proc = subprocess.run(cmd, capture_output=True, text=True)
    return proc, proc.stdout


class Preflight(unittest.TestCase):
    def ok(self, metadata, **kw):
        proc, out = run(metadata, **kw)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        obj = json.loads(out)
        self.assertEqual(obj["status"], "ok")
        return obj

    def needs(self, metadata, **kw):
        proc, out = run(metadata, **kw)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        obj = json.loads(out)
        self.assertEqual(obj["status"], "needs_input")
        q = obj["question"]
        self.assertTrue(q["header"] and q["question"])
        self.assertGreaterEqual(len(q["options"]), 2)
        for opt in q["options"]:
            self.assertTrue(opt["label"] and opt["description"])
        return obj, out

    def test_independent_ok_wrapper_params(self):
        obj = self.ok(INDEP_META)
        self.assertEqual(obj["wrapper_params"], {
            "study_design": "independent", "sample_id_col": "sample",
            "contrast_variable": "condition", "reference_level": "control",
            "target_level": "treated"})
        self.assertEqual(obj["contrast_counts"], {"treated": 2, "control": 2})

    def test_paired_ok_wrapper_params(self):
        obj = self.ok(PAIRED_META, design="paired", donor_col="donor")
        self.assertEqual(obj["wrapper_params"], {
            "study_design": "paired", "sample_id_col": "sample",
            "contrast_variable": "condition", "reference_level": "control",
            "target_level": "treated", "donor_col": "donor"})
        self.assertEqual(obj["n_pairs"], 3)

    def test_no_invented_wrapper_keys(self):
        for meta, kw in ((INDEP_META, {}), (PAIRED_META, {"design": "paired", "donor_col": "donor"})):
            params = self.ok(meta, **kw)["wrapper_params"]
            self.assertNotIn("design_formula", params)
            self.assertNotIn("donor_subset", params)
            self.assertLessEqual(set(params), {"study_design", "sample_id_col", "contrast_variable",
                                               "reference_level", "target_level", "donor_col"})

    def test_unknown_design_asks_with_mapping_and_descriptive_options(self):
        obj, _ = self.needs(INDEP_META, design="unknown")
        labels = {o["label"] for o in obj["question"]["options"]}
        self.assertIn("Paired samples", labels)
        self.assertIn("Independent samples", labels)
        descs = " ".join(o["description"].lower() for o in obj["question"]["options"])
        self.assertIn("mapping", descs)
        self.assertIn("Descriptive summary only", labels)

    def test_numeric_looking_levels_match_as_exact_strings(self):
        obj = self.ok(NUM_META, reference="2", target="1")
        self.assertEqual(obj["wrapper_params"]["target_level"], "1")
        self.needs(NUM_META, reference="2", target="01")

    def test_missing_group_level_asks(self):
        self.needs(INDEP_META, target="Treated")

    def test_identical_levels_ask(self):
        self.needs(INDEP_META, reference="treated", target="treated")

    def test_missing_column_asks(self):
        self.needs(INDEP_META, condition_col="group")

    def test_missing_values_ask_without_leaking_ids(self):
        meta = INDEP_META + "s5\t\td5\n"
        _, out = self.needs(meta)
        self.assertNotIn("s5", out)

    def test_duplicate_sample_ids_ask_without_leaking_ids(self):
        meta = INDEP_META + "s1\ttreated\td9\n"
        _, out = self.needs(meta)
        self.assertNotIn("s1", out)
        self.assertNotIn("d9", out)

    def test_malformed_row_asks(self):
        self.needs(INDEP_META + "s5\ttreated\n")

    def test_unbalanced_donors_ask_without_exclusion_or_leak(self):
        meta = PAIRED_META + "s7\ttreated\td4\n"
        obj, out = self.needs(meta, design="paired", donor_col="donor")
        self.assertNotIn("d4", out)
        self.assertNotIn("s7", out)
        self.assertNotIn("allow_unpaired", out)
        self.assertNotIn("donor_subset", out)
        self.assertTrue(obj["reasons"])

    def test_paired_without_donor_col_asks(self):
        obj, _ = self.needs(PAIRED_META, design="paired")
        labels = {o["label"] for o in obj["question"]["options"]}
        self.assertEqual(labels, {"Provide donor mapping", "Descriptive summary only"})

    def test_within_group_donor_duplicate_asks(self):
        meta = PAIRED_META + "s7\ttreated\td1\n"
        _, out = self.needs(meta, design="paired", donor_col="donor")
        self.assertNotIn("d1", out)

    def test_insufficient_pairs_asks(self):
        meta = "sample\tcondition\tdonor\ns1\ttreated\td1\ns2\tcontrol\td1\n"
        self.needs(meta, design="paired", donor_col="donor")

    def test_insufficient_independent_replicates_asks(self):
        meta = "sample\tcondition\ns1\ttreated\ns2\tcontrol\ns3\tcontrol\n"
        self.needs(meta)

    def test_independent_donor_conflict_asks(self):
        meta = ("sample\tcondition\tdonor\n"
                "s1\ttreated\td1\ns2\ttreated\td2\n"
                "s3\tcontrol\td1\ns4\tcontrol\td4\n")
        obj, out = self.needs(meta, donor_col="donor")
        self.assertNotIn("d1", out)
        labels = {o["label"] for o in obj["question"]["options"]}
        self.assertIn("Provide donor mapping", labels)

    def test_independent_unique_donors_ok(self):
        obj = self.ok(INDEP_META, donor_col="donor")
        self.assertEqual(obj["wrapper_params"]["study_design"], "independent")

    def test_two_pairs_warns(self):
        meta = "sample\tcondition\tdonor\ns1\ttreated\td1\ns2\tcontrol\td1\ns3\ttreated\td2\ns4\tcontrol\td2\n"
        obj = self.ok(meta, design="paired", donor_col="donor")
        self.assertTrue(obj["warnings"])

    def test_csv_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "meta.csv")
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(INDEP_META.replace("\t", ","))
            proc = subprocess.run(
                [sys.executable, SCRIPT, "--metadata", path, "--design", "independent",
                 "--sample_id_col", "sample", "--condition_col", "condition",
                 "--reference_level", "control", "--target_level", "treated"],
                capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(json.loads(proc.stdout)["status"], "ok")


if __name__ == "__main__":
    unittest.main()