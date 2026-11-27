#!/usr/bin/env python3
"""Contract, layout, AST, schema, and end-to-end runtime tests for molecular-cloning (0.1.0)."""

from __future__ import annotations

import ast
import importlib.util
import json
import py_compile
import re
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLUGIN_DIR = ROOT / "resources" / "plugins" / "molecular-cloning"
SKILL_DIR = PLUGIN_DIR / "skills" / "cloning-workflow"

HAS_SCIPY_STACK = all(
    importlib.util.find_spec(pkg) is not None
    for pkg in ("numpy", "pandas", "scipy", "sklearn", "matplotlib")
)
HAS_JSONSCHEMA = importlib.util.find_spec("jsonschema") is not None


def _run_json(script_name: str, args: list[str], cwd: Path) -> dict:
    script = SKILL_DIR / "scripts" / script_name
    proc = subprocess.run(
        [sys.executable, "-B", str(script), *args],
        cwd=str(cwd),
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise AssertionError(
            f"{script_name} failed (rc={proc.returncode}):\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return json.loads(proc.stdout)


def _validate_schema_and_sidecars(testcase: unittest.TestCase, payload: dict, schema_name: str) -> None:
    schema_path = SKILL_DIR / "schemas" / schema_name
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    if HAS_JSONSCHEMA:
        import jsonschema

        jsonschema.validate(instance=payload, schema=schema)
    for req in schema.get("required", []):
        testcase.assertIn(req, payload)
    testcase.assertGreater(len(payload["artifacts"]), 0)
    for art in payload["artifacts"]:
        p = Path(art)
        sidecar = p.with_name(p.name + ".phi-artifact.json")
        testcase.assertTrue(p.is_file(), f"missing artifact {art}")
        testcase.assertTrue(sidecar.is_file(), f"missing sidecar {sidecar}")
        meta = json.loads(sidecar.read_text(encoding="utf-8"))
        testcase.assertEqual(meta["contractVersion"], "1.0.0")
        testcase.assertEqual(meta["provenance"]["tool"], payload["tool"])
        testcase.assertEqual(meta["provenance"]["skill"], "cloning-workflow")


def _write_synthetic_ab1(
    path: Path,
    seq: str,
    quals: list[int] | None = None,
    mixed_peak: tuple[int, str, int] | None = None,
) -> None:
    """Synthesize a spec-compliant binary ABIF (.ab1) chromatogram file."""
    seq = seq.upper()
    n_bases = len(seq)
    if quals is None:
        quals = [38] * n_bases
    ploc = [12 + 10 * i for i in range(n_bases)]
    n_pts = 12 + 10 * n_bases + 15

    # Channel order G A T C -> DATA 9, 10, 11, 12
    channels: dict[str, list[int]] = {b: [15] * n_pts for b in "GATC"}
    for idx, b in enumerate(seq):
        c = ploc[idx]
        for dt, amp in ((-2, 260), (-1, 680), (0, 1000), (1, 680), (2, 260)):
            if 0 <= c + dt < n_pts:
                channels[b][c + dt] += amp

    if mixed_peak is not None:
        m_idx, alt_b, alt_amp = mixed_peak
        c = ploc[m_idx]
        for dt, scale in ((-2, 0.26), (-1, 0.68), (0, 1.0), (1, 0.68), (2, 0.26)):
            if 0 <= c + dt < n_pts:
                channels[alt_b][c + dt] += int(round(alt_amp * scale))

    # Build payloads > 4 bytes first starting at offset 128
    cursor = 128
    blobs: list[bytes] = []
    dir_entries: list[bytes] = []

    def add_entry(tag4: bytes, tag_num: int, elem_type: int, elem_size: int, n_elems: int, payload: bytes) -> None:
        nonlocal cursor
        d_size = len(payload)
        if d_size <= 4:
            padded = payload.ljust(4, b"\x00")
            inline_int = struct.unpack(">I", padded)[0]
            dir_entries.append(
                struct.pack(">4sIHHIiIi", tag4, tag_num, elem_type, elem_size, n_elems, d_size, inline_int, 0)
            )
        else:
            d_off = cursor
            blobs.append(payload)
            cursor += d_size
            dir_entries.append(
                struct.pack(">4sIHHIiIi", tag4, tag_num, elem_type, elem_size, n_elems, d_size, d_off, 0)
            )

    add_entry(b"FWO_", 1, 2, 1, 4, b"GATC")
    add_entry(b"PBAS", 2, 2, 1, n_bases, seq.encode("ascii"))
    add_entry(b"PCON", 2, 1, 1, n_bases, bytes(quals))
    add_entry(b"PLOC", 2, 4, 2, n_bases, struct.pack(f">{n_bases}H", *ploc))
    for ch_base, tag_num in (("G", 9), ("A", 10), ("T", 11), ("C", 12)):
        arr = channels[ch_base]
        add_entry(b"DATA", tag_num, 4, 2, n_pts, struct.pack(f">{n_pts}h", *arr))

    dir_offset = cursor
    num_elems = len(dir_entries)
    root_entry = struct.pack(
        ">4sIHHIiIi", b"tdir", 1, 1023, 28, num_elems, num_elems * 28, dir_offset, 0
    )
    header = (b"ABIF" + struct.pack(">H", 101) + root_entry).ljust(128, b"\x00")
    path.write_bytes(header + b"".join(blobs) + b"".join(dir_entries))


class TestMolecularCloningContract(unittest.TestCase):
    def test_plugin_layout_locks_schemas_and_ast(self) -> None:
        manifest = (PLUGIN_DIR / "phi-package.yaml").read_text(encoding="utf-8")
        self.assertIn("id: molecular-cloning", manifest)
        self.assertIn("version: 0.1.0", manifest)
        self.assertIn("toolPrefix: mc", manifest)
        self.assertIn("- agents/MolecularCloning.md", manifest)
        self.assertIn("- skills/cloning-workflow", manifest)

        for req in ("icon.svg", "README.md", "CHANGELOG.md", "agents/MolecularCloning.md"):
            self.assertTrue((PLUGIN_DIR / req).is_file(), f"missing {req}")
        for reserved in ("mcp", "wrappers", "orchestrator"):
            self.assertFalse((PLUGIN_DIR / reserved).exists(), f"reserved dir {reserved} exists")

        for plat in ("darwin-arm64", "darwin-x64", "linux-x64"):
            lock = PLUGIN_DIR / f"environments/mc/locks/{plat}.txt"
            self.assertTrue(lock.is_file(), f"missing lock {plat}")
            ltxt = lock.read_text(encoding="utf-8")
            self.assertIn("@EXPLICIT", ltxt)
            self.assertIn("# spec: resources/plugins/molecular-cloning/environments/mc/environment.yml", ltxt)

        stext = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
        sfm = stext[4 : stext.index("\n---\n", 4)]
        self.assertIn("name: cloning-workflow", sfm)
        self.assertIn("  environment: plugin:mc", sfm)
        self.assertIn("    - MolecularCloning", sfm)
        self.assertIsNone(re.search(r"^\s+toolPrefix:", sfm, re.M))
        for line in sfm.splitlines():
            m = re.match(r"^( +)-\s", line)
            if m:
                self.assertEqual(len(m.group(1)) % 2, 0, f"odd YAML indent: {line}")

        names = re.findall(r"^    - name:\s+(\S+)", sfm, re.M)
        outs = re.findall(r"^      output:\s+(\S+)", sfm, re.M)
        self.assertEqual(names, ["design_primers", "codon_optimize", "assemble_construct", "sanger_verify"])
        self.assertEqual(len(outs), 4)
        for rel_out in outs:
            self.assertTrue((SKILL_DIR / rel_out).is_file())

        exp_scripts = {
            "_common.py",
            "mc_design_primers.py",
            "mc_codon_optimize.py",
            "mc_assemble_construct.py",
            "mc_sanger_verify.py",
        }
        actual_scripts = {p.name for p in (SKILL_DIR / "scripts").glob("*.py")}
        self.assertEqual(actual_scripts, exp_scripts)
        with tempfile.TemporaryDirectory() as pyc_tmp:
            for py in (SKILL_DIR / "scripts").glob("*.py"):
                ast.parse(py.read_text(encoding="utf-8"), filename=str(py))
                py_compile.compile(str(py), cfile=str(Path(pyc_tmp) / f"{py.stem}.pyc"), doraise=True)

        self.assertEqual(list(PLUGIN_DIR.rglob("__pycache__")), [])


@unittest.skipUnless(HAS_SCIPY_STACK, "numpy/pandas/scipy/sklearn/matplotlib required for live E2E")
class TestMolecularCloningRuntimeE2E(unittest.TestCase):
    def test_1_codon_optimize_and_design_primers(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            # Suboptimal DNA CDS with EcoRI (GAATTC), BamHI (GGATCC), BsaI (GGTCTC), TTTT, and rare codons
            raw_cds = (
                "ATG"
                + "AGAAGGCTA"  # R R L (rare E. coli codons)
                + "GAATTCAAG"  # E F K (contains EcoRI GAATTC)
                + "GGATCCCTG"  # G S L (contains BamHI GGATCC)
                + "GGTCTCACC"  # G L T (contains BsaI GGTCTC)
                + "TTTTACGAC"  # F Y D (contains TTTT homopolymer)
                + "GCTGACCTGCTGGAGATGCTGGTGAAGCTGATCGACGCCTCCATCGACCGCATGCTG" * 3
                + "TAA"
            )
            cds_fa = tdir / "raw_cds.fasta"
            cds_fa.write_text(f">TargetCDS\n{raw_cds}\n", encoding="utf-8")

            opt = _run_json(
                "mc_codon_optimize.py",
                [
                    "--sequence_input", str(cds_fa),
                    "--is_dna",
                    "--host", "E_coli_K12",
                    "--add_kozak",
                    "--output_dir", str(tdir / "opt_out"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, opt, "codon-optimize-result.json")
            self.assertTrue(opt["translation_identical"])
            self.assertEqual(opt["remaining_restriction_sites"], 0)
            self.assertEqual(opt["remaining_homopolymers"], 0)
            self.assertGreater(opt["cai_after"], opt["cai_before"])
            self.assertTrue(opt["kozak_prepended"])

            # Use optimized GenBank file directly as input to mc_design_primers!
            gb_tmpl = opt["outputs"]["optimized_construct_gb"]
            p_qpcr = _run_json(
                "mc_design_primers.py",
                [
                    "--template", gb_tmpl,
                    "--mode", "qPCR",
                    "--output_dir", str(tdir / "primers_qpcr"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, p_qpcr, "design-primers-result.json")
            self.assertGreater(p_qpcr["n_candidates"], 0)
            self.assertLessEqual(p_qpcr["top_pair"]["amplicon_len"], 180)
            self.assertGreaterEqual(p_qpcr["top_pair"]["amplicon_len"], 80)

            p_gg = _run_json(
                "mc_design_primers.py",
                [
                    "--template", gb_tmpl,
                    "--mode", "GoldenGate",
                    "--target_start", "7",
                    "--target_end", "180",
                    "--goldengate_enzyme", "BsaI",
                    "--overhang_5p", "AATG",
                    "--overhang_3p", "GCTT",
                    "--output_dir", str(tdir / "primers_gg"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, p_gg, "design-primers-result.json")
            self.assertIn("GGTCTC", p_gg["top_pair"]["fwd_full"])
            self.assertGreater(p_gg["top_pair"]["fwd_full_tm"], p_gg["top_pair"]["fwd_annealing_tm"])

            p_sdm = _run_json(
                "mc_design_primers.py",
                [
                    "--template", gb_tmpl,
                    "--mode", "SiteDirectedMutagenesis",
                    "--mutation_pos", "90",
                    "--mutation_seq", "AAG",
                    "--output_dir", str(tdir / "primers_sdm"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, p_sdm, "design-primers-result.json")
            self.assertGreater(p_sdm["n_candidates"], 0)

    def test_2_assemble_construct_goldengate_gibson_and_restriction(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tdir = Path(tmp)
            b1 = "CGTAGCTAGCTAGCTGATCGATGCTAGCATCGATCGATGCTAGCTGACTGATCGTAGCTG"
            b2 = "TGACTGATGCTAGCTGACTGATGCTAGCTGACTGATGCTAGCTGACTGATCGTAGCTGAC"
            b3 = "ATCGATGCTAGCTGACTGATCGATGCTAGCTGACTGATGCTAGCTGACTGATCGATGCTA"

            # 2a. GoldenGate (BsaI: GGTCTC + 1nt + OH_L + BODY + OH_R + 1nt + GAGACC)
            # Cycle: AATG -> TTCG -> GCTT -> AATG
            vec_gg = tdir / "vec_gg.fa"
            frags_gg = tdir / "frags_gg.fa"
            vec_gg.write_text(f">pVecGG\nAAAGGTCTCAGCTT{b1}AATGTGA" + "GACCAAA\n", encoding="utf-8")
            frags_gg.write_text(
                f">Part1\nAAAGGTCTCAAATG{b2}TTCGTGAGACCAAA\n"
                f">Part2\nAAAGGTCTCATTCG{b3}GCTTTGAGACCAAA\n",
                encoding="utf-8",
            )
            res_gg = _run_json(
                "mc_assemble_construct.py",
                [
                    "--vector", str(vec_gg),
                    "--fragments", str(frags_gg),
                    "--method", "GoldenGate",
                    "--enzyme", "BsaI",
                    "--construct_name", "pGG_Clone",
                    "--output_dir", str(tdir / "asm_gg"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res_gg, "assemble-construct-result.json")
            self.assertEqual(res_gg["topology"], "circular")
            self.assertEqual(res_gg["n_junctions"], 3)
            self.assertEqual(res_gg["retained_type_iis_sites"], 0)
            self.assertEqual(res_gg["assembled_length_bp"], 3 * 4 + len(b1) + len(b2) + len(b3))

            # 2b. Gibson circular assembly (25-bp overlaps)
            ov1 = "GCATGCTAGCTGATCGATGCTAGCA"
            ov2 = "ATCGATCGATGCTAGCTGACTGATC"
            gib_fa = tdir / "gibson_parts.fa"
            gib_fa.write_text(
                f">Backbone\n{ov2}{b1}{ov1}\n>Insert\n{ov1}{b2}{ov2}\n",
                encoding="utf-8",
            )
            res_gib = _run_json(
                "mc_assemble_construct.py",
                [
                    "--fragments", str(gib_fa),
                    "--method", "Gibson",
                    "--construct_name", "pGibson_Clone",
                    "--output_dir", str(tdir / "asm_gib"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res_gib, "assemble-construct-result.json")
            self.assertEqual(res_gib["n_junctions"], 2)
            self.assertEqual(res_gib["assembled_length_bp"], len(b1) + len(b2) + len(ov1) + len(ov2))

            # 2c. RestrictionLigation (EcoRI GAATTC + BamHI GGATCC)
            rl_fa = tdir / "rl_parts.fa"
            rl_fa.write_text(
                f">VecRL\nCCCGGATCC{b1}GAATTCGGG\n>InsRL\nCCCGAATTC{b2}GGATCCGGG\n",
                encoding="utf-8",
            )
            res_rl = _run_json(
                "mc_assemble_construct.py",
                [
                    "--fragments", str(rl_fa),
                    "--method", "RestrictionLigation",
                    "--enzyme", "EcoRI,BamHI",
                    "--construct_name", "pRL_Clone",
                    "--output_dir", str(tdir / "asm_rl"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res_rl, "assemble-construct-result.json")
            self.assertEqual(res_rl["n_junctions"], 2)

            # 2d. Verify Sanger trace (.ab1) against assembled circular GenBank construct (including origin wrap-around!)
            asm_gb = Path(res_gib["outputs"]["assembled_gb"])
            asm_seq = "".join(
                ln.strip()
                for ln in Path(res_gib["outputs"]["assembled_fasta"]).read_text(encoding="utf-8").splitlines()
                if ln.strip() and not ln.startswith(">")
            )
            wrap_read = asm_seq[-45:] + asm_seq[:55]  # 100 bp spanning plasmid origin

            clean_ab1 = tdir / "clone1_clean.ab1"
            _write_synthetic_ab1(clean_ab1, wrap_read)

            res_sv_ok = _run_json(
                "mc_sanger_verify.py",
                [
                    "--reference", str(asm_gb),
                    "--traces", str(clean_ab1),
                    "--output_dir", str(tdir / "sanger_ok"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res_sv_ok, "sanger-verify-result.json")
            self.assertTrue(res_sv_ok["overall_pass"])
            self.assertEqual(res_sv_ok["total_variants"], 0)
            self.assertEqual(res_sv_ok["total_mixed_peaks"], 0)

            # Synthesize a trace carrying both a clean SNV and a 35% secondary peak (mixed colony)
            mut_list = list(asm_seq[10:95])
            mut_list[20] = "A" if mut_list[20] != "A" else "G"
            alt_peak_base = "T" if mut_list[45] != "T" else "C"
            mixed_ab1 = tdir / "clone2_mixed.ab1"
            _write_synthetic_ab1(mixed_ab1, "".join(mut_list), mixed_peak=(45, alt_peak_base, 550))

            res_sv_bad = _run_json(
                "mc_sanger_verify.py",
                [
                    "--reference", str(asm_gb),
                    "--traces", str(mixed_ab1),
                    "--output_dir", str(tdir / "sanger_mixed"),
                ],
                cwd=tdir,
            )
            _validate_schema_and_sidecars(self, res_sv_bad, "sanger-verify-result.json")
            self.assertFalse(res_sv_bad["overall_pass"])
            self.assertGreaterEqual(res_sv_bad["total_variants"], 1)
            self.assertGreaterEqual(res_sv_bad["total_mixed_peaks"], 1)


if __name__ == "__main__":
    unittest.main()
