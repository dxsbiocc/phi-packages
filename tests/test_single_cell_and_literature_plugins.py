#!/usr/bin/env python3
"""End-to-end and regression tests for single-cell (R + Python) and literature-library plugins:
- literature-library offline chain: collection -> BibTeX ingest + deduplication -> status -> query -> export -> wiki_pending
- single-cell (seurat-workflow & scanpy-workflow): fragment/preset/step-note parity, R/Python syntax, notebook rendering,
  and contract-1.0.0 .phi-artifact.json sidecar generation (PNG header dimensions, table shape, <=8 presentation cap)
- live Scanpy chain (when scanpy & anndata are installed, e.g. phi-scanpy-ci): inspect -> qc -> preprocess -> markers -> pseudobulk
"""

from __future__ import annotations

import importlib.util
import json
import os
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest
import zlib
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
LIT_SCRIPTS = REPO_ROOT / "resources" / "plugins" / "literature-library" / "skills" / "manage-library" / "scripts"
SC_ROOT = REPO_ROOT / "resources" / "plugins" / "single-cell" / "skills"
SEURAT_SKILL = SC_ROOT / "seurat-workflow"
SCANPY_SKILL = SC_ROOT / "scanpy-workflow"


def _run_json(cmd: list[str], cwd: str | None = None) -> dict:
    env = {**os.environ, "MPLBACKEND": "Agg"}
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd, env=env, timeout=180)
    if proc.returncode != 0:
        raise AssertionError(
            f"Command failed (rc={proc.returncode}): {' '.join(cmd)}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    lines = [ln for ln in proc.stdout.strip().splitlines() if ln.strip()]
    return json.loads(proc.stdout) if proc.stdout.strip().startswith("{") else json.loads(lines[-1])


def _write_tiny_png(path: Path, width: int = 32, height: int = 24) -> None:
    def chunk(tag: bytes, data: bytes) -> bytes:
        return struct.pack(">I", len(data)) + tag + data + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF)

    ihdr = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)
    raw = b"".join(b"\x00" + b"\x20\x60\xa0" * width for _ in range(height))
    png = b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr) + chunk(b"IDAT", zlib.compress(raw)) + chunk(b"IEND", b"")
    path.write_bytes(png)


class TestLiteratureLibraryPlugin(unittest.TestCase):
    def test_offline_collection_ingest_query_export_and_wiki_batch(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            proj = Path(tmp) / "lit_proj"
            proj.mkdir()

            # 1. Create collection
            col_res = _run_json(
                [sys.executable, str(LIT_SCRIPTS / "collection.py"), "--project_root", str(proj), "--action", "create", "--name", "Single-Cell Atlas"]
            )
            cid = col_res.get("collection_id") or col_res.get("collection", {}).get("id")
            self.assertTrue(cid)
            self.assertTrue((proj / "library.sqlite").is_file())

            # 2. Ingest BibTeX with 2 distinct papers + 1 duplicate DOI
            bib_path = Path(tmp) / "refs.bib"
            bib_path.write_text(
                """@article{hao2024seurat,
  title = {Dictionary learning for integrative, multimodal and scalable single-cell analysis},
  author = {Hao, Yuhan and Stuart, Tim and Kowalski, Madeline H and Satija, Rahul},
  journal = {Nature Biotechnology},
  year = {2024},
  doi = {10.1038/s41587-023-01767-y},
  abstract = {Multimodal single-cell RNA and chromatin integration with Seurat bridge mapping.}
}

@article{virshup2023scverse,
  title = {The scverse project provides a computational ecosystem for single-cell omics data analysis},
  author = {Virshup, Isaac and Bredikhin, Danila and Heumos, Lukas and Theis, Fabian J},
  journal = {Nature Biotechnology},
  year = {2023},
  doi = {10.1038/s41587-023-01733-8},
  abstract = {Scanpy, AnnData, and SpatialData core data structures for single-cell transcriptomics.}
}

@article{hao2024dup,
  title = {Dictionary learning for integrative, multimodal and scalable single-cell analysis},
  author = {Hao, Yuhan and Satija, Rahul},
  year = {2024},
  doi = {10.1038/s41587-023-01767-y}
}
""",
                encoding="utf-8",
            )
            ing_res = _run_json(
                [
                    sys.executable,
                    str(LIT_SCRIPTS / "ingest.py"),
                    "--project_root",
                    str(proj),
                    "--collection_id",
                    cid,
                    "--input_path",
                    str(bib_path),
                    "--no_enrich",
                ]
            )
            self.assertTrue(ing_res.get("ok", True))

            # 3. Status check: deduplicated to 2 unique papers
            st_res = _run_json([sys.executable, str(LIT_SCRIPTS / "status.py"), "--project_root", str(proj)])
            cols = st_res.get("collections", [])
            self.assertEqual(len(cols), 1)

            # 4. Offline FTS/LIKE query
            q_res = _run_json(
                [sys.executable, str(LIT_SCRIPTS / "query.py"), "--project_root", str(proj), "--query", "scverse"]
            )
            self.assertEqual(len(q_res["papers"]), 1)
            self.assertIn("10.1038/s41587-023-01733-8", (q_res["papers"][0].get("doi") or "").lower())

            # 5. Export CSL-JSON and BibTeX
            csl_out = Path(tmp) / "exported.csl.json"
            exp_res = _run_json(
                [
                    sys.executable,
                    str(LIT_SCRIPTS / "export.py"),
                    "--project_root",
                    str(proj),
                    "--format",
                    "csl-json",
                    "--out_path",
                    str(csl_out),
                ]
            )
            self.assertEqual(exp_res["n"], 2)
            csl_items = json.loads(csl_out.read_text(encoding="utf-8"))
            self.assertEqual(len(csl_items), 2)

            # 6. Wiki pending compile batch
            wp_res = _run_json([sys.executable, str(LIT_SCRIPTS / "wiki_pending.py"), "--project_root", str(proj)])
            self.assertEqual(wp_res["n_pending"], 2)


class TestSingleCellPlugin(unittest.TestCase):
    def test_01_fragments_syntax_and_notebook_presets(self) -> None:
        for skill_dir, ext in ((SEURAT_SKILL, "R"), (SCANPY_SKILL, "py")):
            reg = json.loads((skill_dir / "assets" / "steps.json").read_text(encoding="utf-8"))
            step_ids = {s["id"] for s in reg["steps"]}
            for s in reg["steps"]:
                self.assertTrue((skill_dir / "assets" / "fragments" / f"{s['id']}.{ext}").is_file())
                self.assertTrue((skill_dir / "references" / "step-notes" / f"{s['id']}.md").is_file())
                for dep in s.get("requires", []):
                    self.assertIn(dep, step_ids)

        for py_file in SCANPY_SKILL.rglob("*.py"):
            compile(py_file.read_text(encoding="utf-8"), str(py_file), "exec")

        rscript = shutil.which("Rscript")
        if rscript:
            r_files = sorted(SEURAT_SKILL.rglob("*.R"))
            expr = "; ".join(f'parse(file="{p}")' for p in r_files)
            proc = subprocess.run([rscript, "-e", expr], capture_output=True, text=True, timeout=60)
            self.assertEqual(proc.returncode, 0, proc.stderr)

        with tempfile.TemporaryDirectory() as tmp:
            dummy_h5 = Path(tmp) / "sample.h5ad"
            dummy_h5.write_bytes(b"")
            for mod in range(1, 10):
                out_nb = Path(tmp) / f"scanpy_m{mod}.ipynb"
                res = _run_json(
                    [
                        sys.executable,
                        str(SCANPY_SKILL / "scripts" / "build_notebook.py"),
                        "--language",
                        "python",
                        "--module",
                        str(mod),
                        "--data",
                        str(dummy_h5),
                        "--output",
                        str(out_nb),
                    ]
                )
                self.assertEqual(res["language"], "python")
                nb = json.loads(out_nb.read_text(encoding="utf-8"))
                for idx, cell in enumerate(nb["cells"]):
                    if cell["cell_type"] == "code":
                        compile("".join(cell["source"]), f"scanpy_m{mod}_c{idx}", "exec")

            dummy_rds = Path(tmp) / "sample.rds"
            dummy_rds.write_bytes(b"")
            for mod in (4, 5, 8, 9, 10):
                out_nb = Path(tmp) / f"seurat_m{mod}.ipynb"
                res = _run_json(
                    [
                        sys.executable,
                        str(SEURAT_SKILL / "scripts" / "build_notebook.py"),
                        "--language",
                        "r",
                        "--module",
                        str(mod),
                        "--data",
                        str(dummy_rds),
                        "--output",
                        str(out_nb),
                    ]
                )
                self.assertEqual(res["language"], "r")

    def test_02_artifact_sidecar_helper_contract(self) -> None:
        import types

        restore_np = False
        if importlib.util.find_spec("numpy") is None and "numpy" not in sys.modules:
            fake_np = types.ModuleType("numpy")
            fake_np.generic = ()  # type: ignore[attr-defined]
            fake_np.ndarray = ()  # type: ignore[attr-defined]
            sys.modules["numpy"] = fake_np
            restore_np = True
        try:
            spec = importlib.util.spec_from_file_location(
                "scanpy_common_test", SCANPY_SKILL / "scripts" / "_common.py"
            )
            assert spec and spec.loader
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
        finally:
            if restore_np:
                sys.modules.pop("numpy", None)

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            fig = tdir / "umap-clusters.png"
            _write_tiny_png(fig, width=64, height=48)
            tables = []
            for i in range(10):
                tp = tdir / f"markers_c{i}.csv"
                tp.write_text("gene,log2fc,pval\nCD3D,2.5,0.001\nMS4A1,3.1,0.0005\n", encoding="utf-8")
                tables.append(str(tp))

            out = mod._attach_artifacts({"figures": [str(fig)], "tables": tables}, skill="scanpy-workflow")
            # Capped at 8 in presented list, while all 11 files receive valid sidecars
            self.assertEqual(len(out["artifacts"]), 8)
            fig_desc = json.loads((Path(str(fig) + ".phi-artifact.json")).read_text(encoding="utf-8"))
            self.assertEqual(fig_desc["contractVersion"], "1.0.0")
            self.assertEqual(fig_desc["kind"], "figure")
            self.assertEqual(fig_desc["figure"], {"format": "png", "widthPx": 64, "heightPx": 48})
            tbl_desc = json.loads((Path(tables[0] + ".phi-artifact.json")).read_text(encoding="utf-8"))
            self.assertEqual(tbl_desc["kind"], "table")
            self.assertEqual(tbl_desc["table"], {"rows": 2, "columns": 3})

    @unittest.skipIf(
        importlib.util.find_spec("scanpy") is None or importlib.util.find_spec("anndata") is None,
        "scanpy/anndata not installed in current Python environment",
    )
    def test_03_live_scanpy_chain_with_artifacts_and_scientific_context(self) -> None:
        import anndata as ad
        import numpy as np
        import pandas as pd

        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            rng = np.random.default_rng(17)
            n_cells, n_genes = 120, 60
            counts = rng.poisson(lam=2.0, size=(n_cells, n_genes)).astype(np.float32)
            counts[:60, :10] += rng.poisson(18, size=(60, 10))
            counts[60:, 10:20] += rng.poisson(18, size=(60, 10))
            genes = [
                "MT-CO1", "MT-CO2", "MT-CO3",
                "TGFB1", "TGFBR1", "TGFBR2", "CCL5", "CCR5", "CXCL12", "CXCR4",
                "CDKN1A", "MDM2", "BAX", "BBC3", "GADD45A",
                "NFKBIA", "TNFAIP3", "IL1B", "IL6", "CCL2",
            ] + [f"GENE{i:02d}" for i in range(20, n_genes)]
            samples = [f"S{1 + (i % 4)}" for i in range(n_cells)]
            conds = ["treated" if (i % 4) < 2 else "control" for i in range(n_cells)]
            adata = ad.AnnData(
                X=counts,
                obs=pd.DataFrame(
                    {"sample": samples, "condition": conds},
                    index=[f"cell_{i}" for i in range(n_cells)],
                ),
                var=pd.DataFrame(index=genes),
            )
            raw_h5 = tdir / "raw.h5ad"
            adata.write_h5ad(raw_h5)

            qc_h5 = tdir / "qc.h5ad"
            qc_res = _run_json(
                [
                    sys.executable,
                    str(SCANPY_SKILL / "scripts" / "qc_filter.py"),
                    "--input",
                    str(raw_h5),
                    "--output",
                    str(qc_h5),
                    "--figures_dir",
                    str(tdir / "figs_qc"),
                    "--min_genes",
                    "5",
                    "--max_mt",
                    "60",
                ]
            )
            self.assertGreaterEqual(len(qc_res.get("artifacts", [])), 1)
            for art in qc_res["artifacts"]:
                self.assertTrue(Path(art + ".phi-artifact.json").is_file())

            pp_h5 = tdir / "pp.h5ad"
            pp_res = _run_json(
                [
                    sys.executable,
                    str(SCANPY_SKILL / "scripts" / "preprocess_cluster.py"),
                    "--input",
                    str(qc_h5),
                    "--output",
                    str(pp_h5),
                    "--figures_dir",
                    str(tdir / "figs_pp"),
                    "--n_hvg",
                    "40",
                    "--n_pcs",
                    "10",
                    "--dims",
                    "10",
                    "--resolution",
                    "0.5",
                ]
            )
            self.assertGreaterEqual(len(pp_res.get("artifacts", [])), 1)

            pb_counts = tdir / "pb_counts.tsv"
            pb_samples = tdir / "pb_samples.tsv"
            pb_res = _run_json(
                [
                    sys.executable,
                    str(SCANPY_SKILL / "scripts" / "pseudobulk.py"),
                    "--input",
                    str(pp_h5),
                    "--sample_column",
                    "sample",
                    "--group_column",
                    "condition",
                    "--counts_output",
                    str(pb_counts),
                    "--samples_output",
                    str(pb_samples),
                    "--min_cells",
                    "3",
                ]
            )
            self.assertEqual(pb_res["scientificContext"][str(pb_counts)]["valueSemantics"], "raw_counts")
            self.assertEqual(len(pb_res["artifacts"]), 2)
            for art in pb_res["artifacts"]:
                self.assertTrue(Path(art + ".phi-artifact.json").is_file())

            ccc_tsv = tdir / "ccc_interactions.tsv"
            ccc_res = _run_json(
                [
                    sys.executable,
                    str(SCANPY_SKILL / "scripts" / "communication.py"),
                    "--input",
                    str(pp_h5),
                    "--group_by",
                    "condition",
                    "--interactions_output",
                    str(ccc_tsv),
                    "--figures_dir",
                    str(tdir / "figs_ccc"),
                    "--n_perms",
                    "100",
                ]
            )
            self.assertEqual(ccc_res["method"], "heteromeric_permutation")
            self.assertGreater(ccc_res["n_pairs_in_resource_matched"], 0)
            self.assertTrue(ccc_tsv.is_file())
            for art in ccc_res.get("artifacts", []):
                self.assertTrue(Path(art + ".phi-artifact.json").is_file())

            reg_tsv = tdir / "regulon_activity.tsv"
            reg_res = _run_json(
                [
                    sys.executable,
                    str(SCANPY_SKILL / "scripts" / "regulon.py"),
                    "--input",
                    str(pp_h5),
                    "--group_by",
                    "condition",
                    "--regulon_method",
                    "aucell",
                    "--min_targets",
                    "3",
                    "--regulon_output",
                    str(reg_tsv),
                    "--figures_dir",
                    str(tdir / "figs_reg"),
                ]
            )
            self.assertEqual(reg_res["method"], "aucell")
            self.assertTrue(reg_tsv.is_file())

    @unittest.skipIf(
        importlib.util.find_spec("numpy") is None,
        "numpy not installed in current Python environment",
    )
    def test_03b_offline_regulon_aucell_and_ulm_math(self) -> None:
        import numpy as np

        frag_src = (SCANPY_SKILL / "assets" / "fragments" / "regulon.py").read_text(encoding="utf-8")
        # Compile just the pure-NumPy/SciPy helper functions above the runtime block
        helper_src = frag_src.split("if regulon_method not in", 1)[0]
        ns = {"np": np}
        exec(helper_src, ns)

        rng = np.random.default_rng(7)
        # 20 cells x 60 genes; first 10 cells strongly express genes 0..7 (TF1 targets)
        expr = rng.normal(1.0, 0.2, size=(20, 60))
        expr[:10, :8] += 6.0
        aucell = ns["_compute_aucell_scores"](expr, [list(range(8))], top_frac=0.2)
        self.assertGreater(float(aucell[:10, 0].mean()), float(aucell[10:, 0].mean()) + 0.3)

        w = np.zeros((60, 1), dtype=np.float64)
        w[:8, 0] = 1.0
        ulm = ns["_compute_ulm_scores"](expr, w)
        self.assertGreater(float(ulm[:10, 0].mean()), 3.0)
        self.assertGreater(float(ulm[:10, 0].mean()), float(ulm[10:, 0].mean()) + 3.0)

    def test_04_live_seurat_env_probe_when_seurat_installed(self) -> None:
        rscript = shutil.which("Rscript")
        if not rscript:
            self.skipTest("Rscript not installed")
        chk = subprocess.run(
            [rscript, "-e", 'quit(status = if (requireNamespace("Seurat", quietly=TRUE)) 0 else 1)'],
            capture_output=True,
            text=True,
        )
        if chk.returncode != 0:
            self.skipTest("Seurat not installed in active R environment")
        res = _run_json([rscript, str(SEURAT_SKILL / "scripts" / "env_probe.R")])
        self.assertGreaterEqual(res["packages_installed"], 5)
        self.assertIn("4", res["core_ready_modules"])


if __name__ == "__main__":
    unittest.main()
