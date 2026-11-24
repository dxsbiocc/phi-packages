"""End-to-end unit and contract tests for the bulk-rna-seq plugin (toolPrefix: rna):
  - rna_inspect (raw_counts vs TPM detection, batch-condition collinearity)
  - rna_qc (library sizes, PCA, correlation heatmap, outlier flagging)
  - rna_compare (multi-contrast DEG integration, shared_up/down/discordant)
  - rna_ora (hypergeometric ORA with explicit background universe)
  - rna_gsea (pre-ranked GSEA ES/NES/FDR + leading-edge + curves)
  - rna_pathway_score (single-sample z-score pathway activity + group tests)
"""

from __future__ import annotations

import json
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
RNA_PLUGIN = REPO_ROOT / "resources" / "plugins" / "bulk-rna-seq"
RNA_SCRIPTS = RNA_PLUGIN / "skills" / "rnaseq-workflow" / "scripts"


def _run_json(cmd: list[str]) -> dict:
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return json.loads(proc.stdout)


class TestBulkRnaSeqPlugin(unittest.TestCase):
    def test_plugin_layout_and_icon(self) -> None:
        self.assertTrue((RNA_PLUGIN / "phi-package.yaml").is_file())
        self.assertTrue((RNA_PLUGIN / "icon.svg").is_file())
        self.assertTrue((RNA_PLUGIN / "agents" / "BulkRnaSeq.md").is_file())
        self.assertTrue((RNA_PLUGIN / "skills" / "rnaseq-workflow" / "SKILL.md").is_file())
        for reserved in ("wrappers", "mcp", "orchestrator"):
            self.assertFalse((RNA_PLUGIN / reserved).exists())

    def test_rna_inspect_and_qc(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            rng = random.Random(101)
            samples = [f"T{i}" for i in range(4)] + [f"C{i}" for i in range(4)]
            lines = ["gene\t" + "\t".join(samples)]
            for g in range(60):
                if g < 15:
                    vals = [rng.randint(400, 600) for _ in range(4)] + [rng.randint(20, 60) for _ in range(4)]
                else:
                    vals = [rng.randint(80, 140) for _ in range(8)]
                lines.append(f"G{g:03d}\t" + "\t".join(str(v) for v in vals))
            mat_path = tmp_dir / "counts.tsv"
            mat_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

            # Confounded batch == condition on purpose to test warning
            meta_lines = ["sample\tcondition\tbatch"] + [
                f"{s}\t{'treated' if i < 4 else 'control'}\t{'B1' if i < 4 else 'B2'}"
                for i, s in enumerate(samples)
            ]
            meta_path = tmp_dir / "meta.tsv"
            meta_path.write_text("\n".join(meta_lines) + "\n", encoding="utf-8")

            insp = _run_json(
                [
                    sys.executable,
                    str(RNA_SCRIPTS / "rna_inspect.py"),
                    "--matrix",
                    str(mat_path),
                    "--metadata",
                    str(meta_path),
                    "--group_col",
                    "condition",
                    "--batch_col",
                    "batch",
                ]
            )
            self.assertEqual(insp["matrix"]["inferred_kind"], "raw_counts")
            self.assertTrue(insp["metadata"]["batch_collinearity"]["confounded"])
            self.assertTrue(any("CONFOUNDED" in w for w in insp["warnings"]))

            qc_dir = tmp_dir / "qc_out"
            qc = _run_json(
                [
                    sys.executable,
                    str(RNA_SCRIPTS / "rna_qc.py"),
                    "--input",
                    str(mat_path),
                    "--metadata",
                    str(meta_path),
                    "--group_col",
                    "condition",
                    "--top_var",
                    "50",
                    "--output_dir",
                    str(qc_dir),
                ]
            )
            self.assertEqual(qc["n_samples"], 8)
            self.assertTrue((qc_dir / "sample-qc.tsv").is_file())
            self.assertTrue((qc_dir / "qc-pca.png").is_file())
            self.assertTrue((qc_dir / "qc-correlation-heatmap.png").is_file())

    def test_rna_compare_ora_gsea_and_pathway_score(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            rng = random.Random(7)

            # Build 80-gene DE tables for two contrasts
            de1_lines = ["gene\tlog2FoldChange\tstat\tpvalue\tpadj"]
            de2_lines = ["gene\tlog2FoldChange\tstat\tpvalue\tpadj"]
            for i in range(80):
                g = f"GENE{i:02d}"
                if i < 12:
                    # Strongly upregulated in both contrasts (Pathway_Up)
                    de1_lines.append(f"{g}\t{2.5 + rng.uniform(0, 1):.3f}\t{6.0 + rng.uniform(0, 2):.3f}\t1e-8\t1e-6")
                    de2_lines.append(f"{g}\t{2.1 + rng.uniform(0, 1):.3f}\t{5.2 + rng.uniform(0, 2):.3f}\t1e-7\t1e-5")
                elif i < 20:
                    # Strongly downregulated in both (Pathway_Down)
                    de1_lines.append(f"{g}\t{-2.4 - rng.uniform(0, 1):.3f}\t{-5.5 - rng.uniform(0, 2):.3f}\t1e-7\t1e-5")
                    de2_lines.append(f"{g}\t{-2.0 - rng.uniform(0, 1):.3f}\t{-4.8 - rng.uniform(0, 2):.3f}\t1e-6\t1e-4")
                elif i == 20:
                    # Discordant gene
                    de1_lines.append(f"{g}\t2.2\t5.0\t1e-6\t1e-4")
                    de2_lines.append(f"{g}\t-2.1\t-4.9\t1e-6\t1e-4")
                else:
                    fc = rng.uniform(-0.4, 0.4)
                    de1_lines.append(f"{g}\t{fc:.3f}\t{fc * 2:.3f}\t0.45\t0.85")
                    de2_lines.append(f"{g}\t{-fc:.3f}\t{-fc * 2:.3f}\t0.50\t0.88")
            de1 = tmp_dir / "drugA.tsv"
            de2 = tmp_dir / "drugB.tsv"
            de1.write_text("\n".join(de1_lines) + "\n", encoding="utf-8")
            de2.write_text("\n".join(de2_lines) + "\n", encoding="utf-8")

            # 1. rna_compare
            cmp_dir = tmp_dir / "cmp_out"
            cmp_res = _run_json(
                [
                    sys.executable,
                    str(RNA_SCRIPTS / "rna_compare.py"),
                    "--inputs",
                    f"DrugA={de1}",
                    f"DrugB={de2}",
                    "--output_dir",
                    str(cmp_dir),
                ]
            )
            self.assertEqual(cmp_res["n_shared_up_all"], 12)
            self.assertEqual(cmp_res["n_shared_down_all"], 8)
            self.assertEqual(cmp_res["n_discordant"], 1)
            self.assertTrue((cmp_dir / "multi-contrast-comparison.png").is_file())

            # GMT with enriched up, enriched down, and background random pathway
            gmt_path = tmp_dir / "hallmark.gmt"
            gmt_path.write_text(
                "HALLMARK_UP\tdesc\t" + "\t".join(f"GENE{i:02d}" for i in range(12)) + "\n"
                "HALLMARK_DOWN\tdesc\t" + "\t".join(f"GENE{i:02d}" for i in range(12, 20)) + "\n"
                "HALLMARK_NOISE\tdesc\t" + "\t".join(f"GENE{i:02d}" for i in range(30, 42)) + "\n",
                encoding="utf-8",
            )

            # 2. rna_ora (upregulated)
            ora_dir = tmp_dir / "ora_out"
            ora_res = _run_json(
                [
                    sys.executable,
                    str(RNA_SCRIPTS / "rna_ora.py"),
                    "--de_table",
                    str(de1),
                    "--gmt",
                    str(gmt_path),
                    "--direction",
                    "up",
                    "--output_dir",
                    str(ora_dir),
                ]
            )
            self.assertEqual(ora_res["n_universe"], 80)
            self.assertGreaterEqual(ora_res["n_significant"], 1)
            self.assertEqual(ora_res["top_pathways"][0]["gene_set"], "HALLMARK_UP")
            self.assertTrue((ora_dir / "ora-dotplot.png").is_file())

            # 3. rna_gsea (pre-ranked)
            gsea_dir = tmp_dir / "gsea_out"
            gsea_res = _run_json(
                [
                    sys.executable,
                    str(RNA_SCRIPTS / "rna_gsea.py"),
                    "--de_table",
                    str(de1),
                    "--gmt",
                    str(gmt_path),
                    "--permutations",
                    "200",
                    "--output_dir",
                    str(gsea_dir),
                ]
            )
            self.assertGreaterEqual(gsea_res["n_significant"], 2)
            by_set = {p["gene_set"]: p for p in gsea_res["top_pathways"]}
            self.assertGreater(by_set["HALLMARK_UP"]["nes"], 1.0)
            self.assertLess(by_set["HALLMARK_DOWN"]["nes"], -1.0)
            self.assertTrue((gsea_dir / "gsea-curves.png").is_file())

            # 4. rna_pathway_score
            samples = [f"T{i}" for i in range(4)] + [f"C{i}" for i in range(4)]
            expr_lines = ["gene\t" + "\t".join(samples)]
            for i in range(50):
                if i < 12:
                    vals = [rng.uniform(10, 12) for _ in range(4)] + [rng.uniform(2, 4) for _ in range(4)]
                elif i < 20:
                    vals = [rng.uniform(2, 4) for _ in range(4)] + [rng.uniform(10, 12) for _ in range(4)]
                else:
                    vals = [rng.uniform(5, 7) for _ in range(8)]
                expr_lines.append(f"GENE{i:02d}\t" + "\t".join(f"{v:.3f}" for v in vals))
            expr_path = tmp_dir / "logexpr.tsv"
            expr_path.write_text("\n".join(expr_lines) + "\n", encoding="utf-8")
            meta_path = tmp_dir / "meta.tsv"
            meta_path.write_text(
                "sample\tcondition\n" + "\n".join(f"{s}\t{'T' if i < 4 else 'C'}" for i, s in enumerate(samples)) + "\n",
                encoding="utf-8",
            )
            ps_dir = tmp_dir / "ps_out"
            ps_res = _run_json(
                [
                    sys.executable,
                    str(RNA_SCRIPTS / "rna_pathway_score.py"),
                    "--input",
                    str(expr_path),
                    "--gmt",
                    str(gmt_path),
                    "--metadata",
                    str(meta_path),
                    "--output_dir",
                    str(ps_dir),
                ]
            )
            self.assertEqual(ps_res["n_pathways_scored"], 3)
            self.assertTrue((ps_dir / "pathway-activity-scores.tsv").is_file())
            self.assertTrue((ps_dir / "pathway-group-tests.tsv").is_file())
            self.assertTrue((ps_dir / "pathway-activity-heatmap.png").is_file())


if __name__ == "__main__":
    unittest.main()
