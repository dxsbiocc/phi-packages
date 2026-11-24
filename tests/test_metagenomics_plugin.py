#!/usr/bin/env python3
"""Repo-level tests for the metagenomics plugin (outside the plugin directory, so
nothing here ships in the plugin archive — package contract §2 takes git-tracked
files at the package root only).

Run:  python3 -m pytest tests/test_metagenomics_plugin.py
  or: python3 -m unittest tests.test_metagenomics_plugin

Structural tests need pyyaml + jsonschema. Functional tests additionally need the
plugin environment's packages (pandas/numpy/scipy/matplotlib) and are skipped with
a reason when those are absent. Fixtures are synthetic and generated deterministically
into a temp dir; no real data or network is used.

Checks mirror docs/contracts/{plugin,skill,environment,agent,artifact}.md and the
literature-library tests/test_plugin_contract.py precedent.
"""
from __future__ import annotations

import gzip
import json
import random
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
PLUGIN = ROOT / "resources/plugins/metagenomics"
ARTIFACT_SCHEMA = json.loads((ROOT / "docs/contracts/artifact.schema.json").read_text())

try:
    import jsonschema
except ImportError:  # pragma: no cover
    jsonschema = None

NAME_RE = re.compile(r"^[a-z0-9][a-z0-9-]{0,63}$")
TOOL_RE = re.compile(r"^[a-z][a-z0-9_]{0,31}$")
ARG_RE = re.compile(r"^[a-z][a-z0-9_]*$")
PATH_FORMATS = {"input-path", "project-path"}
TYPES = {"string", "number", "integer", "boolean", "array"}
PLATFORMS = {"darwin-arm64": "osx-arm64", "darwin-x64": "osx-64", "linux-x64": "linux-64"}


def frontmatter(path: Path):
    text = path.read_text(encoding="utf-8")
    assert text.startswith("---\n"), f"{path} needs frontmatter"
    head, body = text.split("\n---\n", 1)
    return yaml.safe_load(head[4:]), body


def markdown_links(text: str):
    for link in re.findall(r"\[[^\]]+\]\(([^)\s]+)\)", text):
        if not link.startswith(("http://", "https://", "#")):
            yield link.split("#")[0]


def script_tools():
    """Yield (skill_dir, tool_dict) for every declared script tool."""
    for sd in sorted((PLUGIN / "skills").iterdir()):
        fm, _ = frontmatter(sd / "SKILL.md")
        for tool in fm.get("phi", {}).get("scripts", []):
            yield sd, tool


@unittest.skipIf(jsonschema is None, "jsonschema not installed")
class PluginStructureTests(unittest.TestCase):
    def test_manifest_fields(self):
        m = yaml.safe_load((PLUGIN / "phi-package.yaml").read_text())
        self.assertEqual(m["schemaVersion"], 1)
        self.assertRegex(m["id"], r"^[a-z][a-z0-9-]{1,63}$")
        self.assertEqual(m["id"], PLUGIN.name)
        self.assertEqual(m["type"], "plugin")
        self.assertRegex(m["version"], r"^\d+\.\d+\.\d+(-[0-9A-Za-z.-]+)?$")
        self.assertTrue(1 <= len(m["title"]) <= 80)
        self.assertTrue(1 <= len(m["summary"].strip()) <= 300)
        self.assertRegex(m["toolPrefix"], r"^[a-z][a-z0-9]{1,11}$")
        self.assertNotIn(m["toolPrefix"], {"skill", "env", "http", "wrapper", "agent", "db", "mcp"})

    def test_tool_prefix_unique(self):
        m = yaml.safe_load((PLUGIN / "phi-package.yaml").read_text())
        prefix = m["toolPrefix"]
        for other in (ROOT / "resources/plugins").iterdir():
            if other.name == PLUGIN.name or not (other / "phi-package.yaml").exists():
                continue
            om = yaml.safe_load((other / "phi-package.yaml").read_text())
            self.assertNotEqual(om.get("toolPrefix"), prefix, other.name)
        for sk in (ROOT / "resources/skills").iterdir():
            if (sk / "SKILL.md").exists():
                sfm, _ = frontmatter(sk / "SKILL.md")
                self.assertNotEqual(sfm.get("phi", {}).get("toolPrefix"), prefix, sk.name)

    def test_no_reserved_dirs(self):
        for reserved in ("mcp", "wrappers", "orchestrator"):
            self.assertFalse((PLUGIN / reserved).exists(), reserved)

    def test_components_complete(self):
        m = yaml.safe_load((PLUGIN / "phi-package.yaml").read_text())
        agents = set(m["components"]["agents"])
        skills = set(m["components"]["skills"])
        for rel in agents | skills:
            self.assertTrue((PLUGIN / rel).exists(), rel)
        self.assertEqual({str(p.relative_to(PLUGIN)) for p in (PLUGIN / "agents").glob("*.md")}, agents)
        self.assertEqual({str(p.relative_to(PLUGIN)) for p in (PLUGIN / "skills").iterdir() if p.is_dir()}, skills)

    def test_environments_and_locks(self):
        m = yaml.safe_load((PLUGIN / "phi-package.yaml").read_text())
        for name, e in m.get("environments", {}).items():
            spec_path = PLUGIN / e["spec"]
            spec = yaml.safe_load(spec_path.read_text())
            self.assertEqual(spec["name"], name)
            self.assertRegex(spec["name"], r"^[a-z][a-z0-9-]{0,62}$")
            self.assertTrue(set(spec) <= {"name", "channels", "dependencies", "description", "host", "sourcePackages"})
            self.assertTrue(spec["channels"])
            self.assertGreaterEqual(len(spec["dependencies"]), 1)
            for plat, sub in PLATFORMS.items():
                lock = spec_path.parent / "locks" / f"{plat}.txt"
                self.assertTrue(lock.exists(), f"{name}/{plat}")
                lines = [l.strip() for l in lock.read_text().splitlines()
                         if l.strip() and not l.strip().startswith("#")]
                self.assertEqual(lines[0], "@EXPLICIT", plat)
                for l in lines[1:]:
                    self.assertRegex(l, r"^https://\S+#[0-9a-f]{32}$", plat)
                    self.assertTrue(f"/{sub}/" in l or "/noarch/" in l, f"{plat}: {l[:80]}")

    def test_agent_contract(self):
        names = set()
        for other in (ROOT / "resources/plugins").iterdir():
            for af in (other / "agents").glob("*.md"):
                names.add(frontmatter(af)[0]["name"])
        for af in (PLUGIN / "agents").glob("*.md"):
            a, _ = frontmatter(af)
            self.assertRegex(a["name"], r"^[A-Z][A-Za-z0-9]*$")
            self.assertFalse(a["name"].endswith("Agent"))
            self.assertEqual(a["name"], af.stem)
            self.assertTrue(1 <= len(a["description"].strip()) <= 1024)
            self.assertTrue(a["tools"])
            for reserved in ("spawns", "outputSchema"):
                self.assertNotIn(reserved, a)
            self.assertNotEqual(a.get("visibility", "entry"), "internal")
            self.assertIn(a.get("delegationMode", "optional"), ("required-first", "preferred", "optional"))
            if "fallback" in a:
                fb = a["fallback"]
                self.assertIsInstance(fb["afterFailures"], int)
                self.assertGreater(fb["afterFailures"], 0)
            for s in a.get("skills", []):
                self.assertTrue((PLUGIN / "skills" / s / "SKILL.md").exists(), s)
            names.discard(a["name"])
        # plugin agent names must not clash with other plugins' agents
        for af in (PLUGIN / "agents").glob("*.md"):
            self.assertNotIn(frontmatter(af)[0]["name"], names)

    def test_skill_frontmatter_and_tools(self):
        final_names = set()
        for sd, tool in script_tools():
            fm, _ = frontmatter(sd / "SKILL.md")
            self.assertEqual(fm["name"], sd.name)
            self.assertRegex(fm["name"], NAME_RE)
            self.assertTrue(1 <= len(fm["description"].strip()) <= 1024)
            phi = fm["phi"]
            self.assertNotIn("toolPrefix", phi)
            self.assertEqual(phi.get("environment"), "plugin:metagenomics")
            self.assertTrue(set(phi.get("attachTo", ["main"])) <= {"main", "Metagenomics"})
            tn = tool["name"]
            self.assertRegex(tn, TOOL_RE)
            final = f"mg_{tn}"
            self.assertNotIn(final, final_names)
            final_names.add(final)
            self.assertTrue(1 <= len(tool["description"].strip()) <= 1024)
            self.assertEqual(tool["run"][0], "python")
            script = (sd / tool["run"][1]).resolve()
            self.assertTrue(script.is_file())
            self.assertTrue(str(script).startswith(str(sd.resolve())))
            args = tool["args"]
            self.assertEqual(args["type"], "object")
            self.assertIs(args["additionalProperties"], False)
            for pn, ps in args["properties"].items():
                self.assertRegex(pn, ARG_RE)
                self.assertIn(ps["type"], TYPES)
                if ps["type"] == "array":
                    self.assertIn(ps["items"]["type"], TYPES - {"array"})
                if "format" in ps:
                    self.assertIn(ps["format"], PATH_FORMATS)
                    self.assertEqual(ps["type"], "string")
            for r in args.get("required", []):
                self.assertIn(r, args["properties"])
            self.assertIn(tool["approval"], ("read", "write"))
            if tool["approval"] == "read":
                self.assertTrue(all(ps.get("format") != "project-path"
                                    for ps in args["properties"].values()), final)
            self.assertTrue(1 <= tool.get("timeoutSeconds", 600) <= 86400)
            if "output" in tool:
                sp = (sd / tool["output"]).resolve()
                self.assertTrue(sp.is_file())
                jsonschema.Draft7Validator.check_schema(json.loads(sp.read_text()))
        self.assertEqual(len(final_names), 19)

    def test_markdown_links_resolve(self):
        for sd in sorted((PLUGIN / "skills").iterdir()):
            for md in [sd / "SKILL.md"] + list((sd / "references").rglob("*.md")):
                for link in markdown_links(md.read_text()):
                    self.assertTrue((md.parent / link).resolve().exists(),
                                    f"{md.relative_to(PLUGIN)}: {link}")


def _deps_available():
    try:
        import pandas, numpy, scipy, matplotlib  # noqa: F401
        return True
    except ImportError:
        return False


def make_fixtures(root: Path):
    """Synthetic per-sample reports + MAG tables, deterministic (seed 7)."""
    inp = root / "in"
    inp.mkdir(parents=True, exist_ok=True)
    rng = random.Random(7)
    samples = ["s1", "s2", "s3", "s4"]
    species = ["Escherichia coli", "Bacteroides fragilis", "Lactobacillus acidophilus",
               "Staphylococcus aureus", "Pseudomonas aeruginosa", "Bifidobacterium longum"]
    taxids = [562, 817, 1579, 1280, 287, 216816]
    for i, s in enumerate(samples):
        lines = ["100.00\t50000\t0\tD\t2\tBacteria"]
        for j, (sp, tid) in enumerate(zip(species, taxids)):
            n = 8000 - j * 1000 + (3000 if (i < 2 and j == 0) else 0) + rng.randint(-200, 200)
            lines.append(f"{n/500:.2f}\t{n}\t{n}\tS\t{tid}\t  {sp}")
        (inp / f"{s}.k2report").write_text("\n".join(lines) + "\n")
    for i, s in enumerate(samples):
        lines = ["name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads"]
        for j, (sp, tid) in enumerate(zip(species, taxids)):
            n = 7000 - j * 900 + (2500 if (i < 2 and j == 1) else 0)
            lines.append(f"{sp}\t{tid}\tS\t{n-100}\t100\t{n}\t{n/50000:.4f}")
        (inp / f"{s}.bracken").write_text("\n".join(lines) + "\n")
    for i, s in enumerate(samples):
        lines = ["#mpa_vOct22_CHOCOPhlAnSGB_202403", "#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species"]
        vals = [30.0 - j * 4 + (10 if (i >= 2 and j == 2) else 0) + rng.uniform(-1, 1) for j in range(len(species))]
        total = sum(vals)
        for j, (sp, tid) in enumerate(zip(species, taxids)):
            clade = f"k__Bacteria|p__Prot|c__Cls|o__Ord|f__Fam|g__Gen{j}|s__{sp.replace(' ', '_')}"
            lines.append(f"{clade}\t{tid}\t{vals[j]/total*100:.4f}\t")
            lines.append(f"{clade}|t__SGB{j}\t\t{vals[j]/total*100:.4f}\t")
        (inp / f"{s}.profile").write_text("\n".join(lines) + "\n")
    for s in samples:
        lines = ["#SampleID\tMetaphlan_Analysis"]
        for j, sp in enumerate(species):
            lines.append(f"k__Bacteria|p__Prot|g__Gen{j}|s__{sp.replace(' ', '_')}\t{5000 - j*500}")
        (inp / f"{s}.mpa").write_text("\n".join(lines) + "\n")
    (inp / "metadata.tsv").write_text("sample\tgroup\ns1\tcase\ns2\tcase\ns3\tcontrol\ns4\tcontrol\n")
    gf = ["# Gene Family\ts1\ts2\ts3\ts4"]
    for k in ["UniRef90_A", "UniRef90_B", "UniRef90_C"]:
        gf.append(f"{k}\t{rng.randint(100,900)}\t{rng.randint(100,900)}\t{rng.randint(100,900)}\t{rng.randint(100,900)}")
        gf.append(f"{k}|g__Gen0.s__Escherichia_coli\t{rng.randint(10,90)}\t{rng.randint(10,90)}\t{rng.randint(10,90)}\t{rng.randint(10,90)}")
    gf.append("UNMAPPED\t500\t600\t550\t520")
    (inp / "genefamilies.tsv").write_text("\n".join(gf) + "\n")
    (inp / "mapping.tsv").write_text("UniRef90_A\tK00001\nUniRef90_B\tK00002\nUniRef90_C\tK00001\n")
    (inp / "report.tsv").write_text(
        "Assembly\tmegahit_s1\tspades_s1\n# contigs (>= 0 bp)\t12000\t9000\n"
        "# contigs (>= 1000 bp)\t4000\t3500\nTotal length (>= 0 bp)\t45000000\t48000000\n"
        "Largest contig\t150000\t280000\nN50\t5200\t8100\nN75\t2100\t3400\n"
        "L50\t2600\t1700\nL75\t6100\t4000\nGC (%)\t48.5\t49.1\n")
    bins = [("bin.1", 97.2, 1.1), ("bin.2", 93.5, 4.2), ("bin.3", 78.0, 6.5),
            ("bin.4", 55.1, 9.8), ("bin.5", 40.0, 3.0), ("bin.6", 88.0, 15.0)]
    q = ["Name\tCompleteness\tContamination\tGenome_Size\tContig_N50\tGC_Content\tTotal_Contigs"]
    for b, c, x in bins:
        q.append(f"{b}\t{c}\t{x}\t{rng.randint(2,5)*1000000}\t{rng.randint(5,50)*1000}\t{rng.uniform(35,60):.2f}\t{rng.randint(50,400)}")
    (inp / "quality_report.tsv").write_text("\n".join(q) + "\n")
    qa = ["Bin Id\tMarker lineage\t# genomes\t# markers\t# marker sets\t0\t1\t2\t3\t4\t5+\tCompleteness\tContamination\tStrain heterogeneity"]
    for b, c, x in bins[:3]:
        qa.append(f"{b}\tk__Bacteria (UID203)\t5449\t104\t58\t0\t100\t4\t0\t0\t0\t{c}\t{x}\t0.00")
    (inp / "checkm_qa.tsv").write_text("\n".join(qa) + "\n")
    bd = inp / "bins"
    bd.mkdir(exist_ok=True)
    for b, _, _ in bins:
        (bd / f"{b}.fa").write_text(f">{b}_ctg1\n" + "ACGT" * 250 + "\n")
    g = ["user_genome\tclassification\tfastani_reference"]
    cls = {
        "bin.1": "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli",
        "bin.2": "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides fragilis",
        "bin.3": "d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus acidophilus",
        "bin.4": "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__",
    }
    for b, c in cls.items():
        g.append(f"{b}\t{c}\tGCF_000005845.2")
    (inp / "gtdbtk.bac120.summary.tsv").write_text("\n".join(g) + "\n")
    cv = ["Genome\tMean\tCovered_fraction"]
    for b, _, _ in bins:
        cv.append(f"{b}\t{rng.uniform(1,80):.2f}\t{rng.uniform(0.2,0.99):.3f}")
    (inp / "coverm.tsv").write_text("\n".join(cv) + "\n")
    (inp / "Cdb.csv").write_text(
        "genome,secondary_cluster,primary_cluster\nbin.1.fa,1_1,1\nbin.2.fa,2_1,2\nbin.3.fa,3_1,3\nbin.4.fa,1_1,1\n")
    pk = inp / "prokka"
    pk.mkdir()
    for b, cds, trna in [("bin.1", 1200, 20), ("bin.2", 900, 15), ("bin.3", 700, 10)]:
        d = pk / b
        d.mkdir()
        (d / f"{b}.txt").write_text(
            f"organism: Genus species\ncontigs: 50\nbases: 3000000\nCDS: {cds}\n"
            f"gene: {cds + 23}\nrRNA: 3\ntRNA: {trna}\ntmRNA: 1\n")
    tsv = ("locus_tag\tftype\tlength_bp\tgene\tEC_number\tproduct\n"
           "B1_00001\tCDS\t900\tgA\t-\thypothetical protein\n")
    for label in ("5S", "16S", "23S"):
        tsv += f"B1_r{label}\trRNA\t100\t-\t-\t{label} ribosomal RNA\n"
    tsv += "B1_t1\ttRNA\t76\t-\t-\ttRNA-Ala\n"
    (pk / "bin.1" / "bin.1.tsv").write_text(tsv)
    (pk / "bin.2" / "bin.2.tsv").write_text(tsv.replace("23S ribosomal RNA", "16S ribosomal RNA"))
    pg = inp / "prodigal"
    pg.mkdir()
    faa = (">bin4_ctg1_1 # 3 # 300 # 1 # ID=1_1;partial=00\nMKVL\n"
           ">bin4_ctg1_2 # 400 # 700 # 1 # ID=1_2\nMABC\n"
           ">bin4_ctg2_1 # 1 # 200 # -1 # ID=2_1\nMDEF\n")
    with gzip.open(pg / "MEGAHIT-bin.4_prodigal.faa.gz", "wt") as fh:
        fh.write(faa)
    (pg / "SPAdes-bin.5_prodigal.gff").write_text(
        "##gff-version 3\nctg1\tProdigal\tCDS\t3\t300\t.\t+\t0\tID=1_1\n"
        "ctg1\tProdigal\tCDS\t400\t700\t.\t+\t0\tID=1_2\n")
    eg = inp / "eggnog"
    eg.mkdir()
    hdr = ("#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\t"
           "Description\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\t"
           "KEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n")
    for b, rows in [("bin.1", [("q1", "K", "gyrA", "ko:K00001,ko:K00002"),
                               ("q2", "KT", "-", "ko:K00001"), ("q3", "-", "-", "-")]),
                    ("bin.2", [("q1", "C", "dnaA", "ko:K00003"), ("q2", "-", "-", "-")])]:
        lines = [hdr]
        for q, cog, name, ko in rows:
            lines.append(f"{q}\t562.123\t1e-50\t500\tCOG123@1|root\tBacteria\t{cog}\tdesc\t"
                         f"{name}\t-\t-\t{ko}\t-\t-\t-\t-\t-\t-\t-\t-\t-\n")
        (eg / f"{b}.emapper.annotations").write_text("".join(lines))
    (inp / "coverm_multi.tsv").write_text(
        "Genome\ts1 Mean\ts1 Covered_fraction\ts2 Mean\ts2 Covered_fraction\n"
        "bin.1\t10\t0.9\t20\t0.8\nbin.2\t30\t0.2\t40\t0.9\n")
    (inp / "q_small.tsv").write_text(
        "bin\tcompleteness\tcontamination\tgenome_size\nbin.1\t95\t2\t1000000\nbin.2\t90\t3\t3000000\n")
    (inp / "classify.tree").write_text(
        "((bin.1:0.1,bin.2:0.2):0.3,(bin.3:0.4,RS_GCF_000005845.2:0.5):0.6,GB_GCA_000001.1:0.7);\n")
    (inp / "noref.tree").write_text("(bin.1:0.1,bin.2:0.2);\n")
    return inp


@unittest.skipUnless(_deps_available() and jsonschema is not None,
                     "plugin environment packages (pandas/numpy/scipy/matplotlib) or jsonschema missing")
class ScriptFunctionalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory(prefix="mg-test-")
        cls.root = Path(cls._tmp.name)
        cls.inp = make_fixtures(cls.root)
        cls.out = cls.root / "out"
        cls.out.mkdir()
        cls.tools = {}
        for sd, tool in script_tools():
            cls.tools[tool["name"]] = (sd, tool)
        cls.artifact_validator = jsonschema.Draft202012Validator(ARTIFACT_SCHEMA)

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def script(self, name):
        sd, tool = self.tools[name]
        return (sd / tool["run"][1]).resolve()

    def run_tool(self, name, argv, expect_ok=True):
        sd, tool = self.tools[name]
        r = subprocess.run([sys.executable, str(self.script(name))] + argv,
                           capture_output=True, text=True, cwd=str(self.root))
        if not expect_ok:
            self.assertNotEqual(r.returncode, 0, f"mg_{name} should fail: {r.stdout[:200]}")
            return r
        self.assertEqual(r.returncode, 0, f"mg_{name} failed: {r.stderr[-400:]}")
        obj = json.loads(r.stdout)  # exactly one JSON object, nothing else (§5.5)
        if "output" in tool:
            schema = json.loads((sd / tool["output"]).read_text())
            jsonschema.Draft7Validator(schema).validate(obj)
        self.check_artifacts(obj, name, sd.name)
        return obj

    def check_artifacts(self, obj, name, skill):
        for path in obj.get("artifacts", []):
            p = Path(path)
            self.assertTrue(p.is_file(), f"artifact missing: {path}")
            desc = p.with_name(p.name + ".phi-artifact.json")
            self.assertTrue(desc.is_file(), f"descriptor missing: {desc}")
            d = json.loads(desc.read_text())
            errs = list(self.artifact_validator.iter_errors(d))
            self.assertFalse(errs, f"{desc}: {errs[0].message if errs else ''}")
            self.assertEqual(d["file"], p.name)
            self.assertEqual(d["provenance"]["tool"], f"mg_{name}")
            self.assertEqual(d["provenance"]["skill"], skill)

    def o(self, *p):
        return str(self.out.joinpath(*p))

    def k2(self):
        return [str(self.inp / f"s{i}.k2report") for i in (1, 2, 3, 4)]

    def mp(self):
        return [str(self.inp / f"s{i}.profile") for i in (1, 2, 3, 4)]

    def meta(self):
        return str(self.inp / "metadata.tsv")

    def merged_k2(self):
        path = self.o("kraken2.tsv")
        self.run_tool("merge", ["--inputs", *self.k2(), "--kind", "kraken2", "--rank", "S",
                                "--output", path, "--metadata", self.meta()])
        return path

    def merged_mp(self):
        path = self.o("metaphlan.tsv")
        self.run_tool("merge", ["--inputs", *self.mp(), "--kind", "metaphlan", "--rank", "S",
                                "--output", path])
        return path

    def test_flag_parity(self):
        for name, (sd, tool) in self.tools.items():
            r = subprocess.run([sys.executable, str(self.script(name)), "--help"],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, f"mg_{name} --help: {r.stderr[-200:]}")
            flags = {f.lstrip("-") for f in re.findall(r"(?<![\w-])(--[a-z][a-z0-9_-]*)", r.stdout)} - {"help"}
            missing = set(tool["args"]["properties"]) - flags
            self.assertFalse(missing, f"mg_{name}: declared args not flags: {missing}")

    def test_inspect(self):
        obj = self.run_tool("inspect", ["--input", str(self.inp)])
        self.assertIn("mg_merge", obj["ready_steps"])
        kinds = {f["kind"] for f in obj["files"]}
        self.assertTrue({"eggnog_annotations", "prokka_summary", "prokka_features",
                         "prodigal_faa", "phylo_tree"} <= kinds, kinds)
        self.assertIn("mg_mag_function (mag-workflow)", obj["ready_steps"])
        self.assertIn("mg_mag_tree (mag-workflow)", obj["ready_steps"])
        self.run_tool("inspect", ["--input", str(self.inp / "quality_report.tsv")])

    def test_merge_value_kinds(self):
        k2 = self.run_tool("merge", ["--inputs", *self.k2(), "--kind", "kraken2", "--rank", "S",
                                     "--output", self.o("k2.tsv")])
        self.assertEqual(k2["value_kind"], "counts")
        mp = self.run_tool("merge", ["--inputs", *self.mp(), "--kind", "metaphlan", "--rank", "S",
                                     "--output", self.o("mp.tsv")])
        self.assertEqual(mp["value_kind"], "relative")

    def test_rarefaction_and_determinism(self):
        mat = self.merged_k2()
        a = self.run_tool("rarefaction", ["--input", mat, "--output_dir", self.o("rare"),
                                          "--metadata", self.meta(), "--group_col", "group"])
        self.run_tool("rarefaction", ["--input", mat, "--output_dir", self.o("rare2"),
                                      "--metadata", self.meta(), "--group_col", "group"])
        t1 = (self.out / "rare" / "rarefaction.tsv").read_bytes()
        t2 = (self.out / "rare2" / "rarefaction.tsv").read_bytes()
        self.assertEqual(t1, t2, "same-seed rarefaction not byte-identical")
        self.assertEqual(a["metric"], "observed")

    def test_rarefaction_rejects_relative(self):
        self.run_tool("rarefaction", ["--input", self.merged_mp(), "--output_dir", self.o("rare-bad")],
                      expect_ok=False)

    def test_diversity_counts_and_relative(self):
        counts = self.run_tool("diversity", ["--input", self.merged_k2(), "--values", "counts",
                                             "--output_dir", self.o("div"), "--metadata", self.meta()])
        self.assertIsNotNone(counts["alpha"]["s1"]["observed"])
        rel = self.run_tool("diversity", ["--input", self.merged_mp(), "--values", "relative",
                                          "--output_dir", self.o("divmp"), "--metadata", self.meta()])
        for row in rel["alpha"].values():
            self.assertIsNone(row["observed"])
            self.assertIsNone(row["chao1"])
            self.assertIsNotNone(row["shannon"])

    def test_differential(self):
        obj = self.run_tool("differential", ["--input", self.merged_k2(), "--metadata", self.meta(),
                                             "--output_dir", self.o("diff"), "--values", "counts",
                                             "--case", "case", "--control", "control"])
        self.assertEqual(obj["test"], "mannwhitneyu")
        self.assertIn("n_significant", obj)

    def test_function_actions(self):
        gf = str(self.inp / "genefamilies.tsv")
        self.run_tool("function", ["--input", gf, "--action", "normalize", "--output", self.o("gf-rel.tsv"), "--unit", "relative"])
        self.run_tool("function", ["--input", gf, "--action", "regroup", "--output", self.o("gf-ko.tsv"), "--mapping", str(self.inp / "mapping.tsv")])
        unstrat = self.run_tool("function", ["--input", gf, "--action", "unstratified", "--output", self.o("gf-unstrat.tsv")])
        strat = self.run_tool("function", ["--input", gf, "--action", "stratified", "--output", self.o("gf-strat.tsv")])
        self.assertGreater(unstrat["rows_out"], 0)
        self.assertGreater(strat["rows_out"], 0)

    def test_plot_barplot_and_clustered_heatmap(self):
        mat = self.merged_k2()
        self.run_tool("plot", ["--input", mat, "--output_dir", self.o("plot-bar"), "--kind", "barplot", "--metadata", self.meta()])
        hm = self.run_tool("plot", ["--input", mat, "--output_dir", self.o("plot-hm"),
                                    "--kind", "heatmap", "--trans", "log10", "--cluster", "both"])
        self.assertTrue(hm["cluster"]["applied"])

    def test_plot_rejects_cluster_on_barplot(self):
        self.run_tool("plot", ["--input", self.merged_k2(), "--output_dir", self.o("plot-bad"),
                               "--kind", "barplot", "--cluster", "rows"], expect_ok=False)

    def test_assembly_stats(self):
        obj = self.run_tool("assembly_stats", ["--input", str(self.inp / "report.tsv")])
        self.assertEqual(obj["tool"], "mg_assembly_stats")

    def test_mag_chain(self):
        qc = self.run_tool("bin_qc", ["--input", str(self.inp / "quality_report.tsv"), "--output_dir", self.o("binqc")])
        self.assertEqual(qc["tier_counts"], {"HQ": 2, "MQ": 2, "LQ": 2})
        self.run_tool("bin_qc", ["--input", str(self.inp / "checkm_qa.tsv"), "--output_dir", self.o("binqc-checkm")])
        qctsv = qc["outputs"]["table"]
        filt = self.run_tool("bin_filter", ["--bins_dir", str(self.inp / "bins"), "--quality", qctsv, "--output_dir", self.o("clean")])
        self.assertEqual(filt["n_passed"], 4)
        table = self.run_tool("mag_table", ["--quality", qctsv, "--gtdb_summary", str(self.inp / "gtdbtk.bac120.summary.tsv"),
                                            "--coverm", str(self.inp / "coverm.tsv"), "--drep_clusters", str(self.inp / "Cdb.csv"),
                                            "--output_dir", self.o("magtable")])
        self.assertEqual(table["joined"], ["quality", "gtdbtk", "coverm", "drep"])

    def test_mag_table_normalize_multisample(self):
        obj = self.run_tool("mag_table", ["--quality", str(self.inp / "q_small.tsv"),
                                          "--coverm", str(self.inp / "coverm_multi.tsv"),
                                          "--normalize", "--output_dir", self.o("mt-multi")])
        self.assertEqual(obj["coverm_samples"], ["s1", "s2"])
        self.assertEqual(obj["n_masked_mean"], 1)  # bin.2 s1 fraction 0.2 < 0.3
        self.assertTrue(obj["normalized"])
        self.assertEqual(obj["normalization_basis"], "genome_size_weighted")
        import csv
        with open(self.o("mt-multi", "mag-table.tsv")) as fh:
            rows = {r["bin"]: r for r in csv.DictReader(fh, delimiter="\t")}
        self.assertEqual(rows["bin.1"]["rel_abund_s1"], "1.0")
        self.assertAlmostEqual(float(rows["bin.1"]["rel_abund_s2"]), 2 / 14, places=6)
        self.assertAlmostEqual(float(rows["bin.2"]["rel_abund_s2"]), 12 / 14, places=6)
        self.assertEqual(rows["bin.2"]["s1 Mean"], "")  # masked

    def test_mag_table_normalize_requires_coverm(self):
        self.run_tool("mag_table", ["--quality", str(self.inp / "q_small.tsv"),
                                    "--normalize", "--output_dir", self.o("mt-bad")],
                      expect_ok=False)

    def test_mag_function(self):
        obj = self.run_tool("mag_function", [
            "--prokka_dir", str(self.inp / "prokka"),
            "--prodigal", str(self.inp / "prodigal/MEGAHIT-bin.4_prodigal.faa.gz"),
            str(self.inp / "prodigal/SPAdes-bin.5_prodigal.gff"),
            "--eggnog", str(self.inp / "eggnog/bin.1.emapper.annotations"),
            str(self.inp / "eggnog/bin.2.emapper.annotations"),
            "--output_dir", self.o("magfunc")])
        self.assertEqual(obj["n_bins"], 5)
        self.assertEqual(obj["mimag_structural"], {"pass": 1, "fail": 1, "unknown": 3})
        self.assertEqual(obj["n_kos"], 3)
        self.assertEqual(obj["n_annotated_queries"], 3)
        self.assertEqual(set(obj["sources"]), {"prokka", "prodigal", "eggnog"})
        import csv
        with open(obj["outputs"]["gene_counts"]) as fh:
            counts = {r["bin"]: r for r in csv.DictReader(fh, delimiter="\t")}
        self.assertEqual(counts["bin.1"]["n_cds"], "1200")  # prokka CDS count
        self.assertEqual(counts["bin.4"]["n_cds"], "3")
        self.assertEqual(counts["bin.5"]["n_cds"], "2")
        self.assertEqual(counts["bin.1"]["mimag_structural_pass"], "True")
        self.assertEqual(counts["bin.2"]["mimag_structural_pass"], "False")
        with open(obj["outputs"]["ko_matrix"]) as fh:
            ko = {r["bin"]: r for r in csv.DictReader(fh, delimiter="\t")}
        self.assertEqual(ko["bin.1"]["K00001"], "2")
        self.assertEqual(ko["bin.2"]["K00003"], "1")

    def test_mag_function_bin_map_and_errors(self):
        bm = self.o("binmap.tsv")
        Path(bm).write_text("file\tbin\nMEGAHIT-bin.4_prodigal.faa.gz\tbin.9\n")
        obj = self.run_tool("mag_function", [
            "--prodigal", str(self.inp / "prodigal/MEGAHIT-bin.4_prodigal.faa.gz"),
            "--bin_map", bm, "--output_dir", self.o("magfunc-map")])
        self.assertEqual(obj["n_bins"], 1)
        self.run_tool("mag_function", ["--output_dir", self.o("magfunc-bad")], expect_ok=False)

    def test_mag_tree(self):
        obj = self.run_tool("mag_tree", ["--tree", str(self.inp / "classify.tree"),
                                         "--summary", str(self.inp / "gtdbtk.bac120.summary.tsv"),
                                         "--output_dir", self.o("magtree")])
        self.assertEqual(obj["n_user_tips"], 3)
        self.assertEqual(obj["n_reference_tips"], 2)
        self.assertEqual(obj["user_identification"], "summary")
        import csv
        with open(obj["outputs"]["distance_matrix"]) as fh:
            dm = {r["bin"]: r for r in csv.DictReader(fh, delimiter="\t")}
        self.assertAlmostEqual(float(dm["bin.1"]["bin.2"]), 0.3, places=9)
        self.assertAlmostEqual(float(dm["bin.1"]["bin.3"]), 1.4, places=9)
        self.assertAlmostEqual(float(dm["bin.2"]["bin.3"]), 1.5, places=9)
        with open(obj["outputs"]["nearest_reference"]) as fh:
            nr = {r["bin"]: r for r in csv.DictReader(fh, delimiter="\t")}
        self.assertEqual(nr["bin.1"]["nearest_reference"], "GB_GCA_000001.1")
        self.assertAlmostEqual(float(nr["bin.1"]["patristic_distance"]), 1.1, places=9)
        self.assertEqual(nr["bin.3"]["nearest_reference"], "RS_GCF_000005845.2")
        self.assertAlmostEqual(float(nr["bin.3"]["patristic_distance"]), 0.9, places=9)
        heur = self.run_tool("mag_tree", ["--tree", str(self.inp / "classify.tree"),
                                          "--output_dir", self.o("magtree-heur")])
        self.assertEqual(heur["user_identification"], "heuristic")
        self.assertEqual(Path(heur["outputs"]["nearest_reference"]).read_text(),
                         Path(obj["outputs"]["nearest_reference"]).read_text())
        noref = self.run_tool("mag_tree", ["--tree", str(self.inp / "noref.tree"),
                                           "--output_dir", self.o("magtree-noref")])
        self.assertEqual(noref["n_reference_tips"], 0)
        self.assertNotIn("nearest_reference", noref["outputs"])

    def test_merge_sgb(self):
        obj = self.run_tool("merge", ["--inputs", *self.mp(), "--kind", "metaphlan",
                                      "--rank", "T", "--output", self.o("sgb.tsv")])
        self.assertEqual(obj["value_kind"], "relative")
        self.assertEqual(obj["taxa_after"], 6)
        self.assertEqual(obj["sgb_species"]["SGB0"], "Escherichia_coli")
        self.run_tool("merge", ["--inputs", *self.k2(), "--kind", "kraken2", "--rank", "T",
                                "--output", self.o("sgb-bad.tsv")], expect_ok=False)

    def test_function_differential_chain(self):
        norm = self.o("gf-rel.tsv")
        self.run_tool("function", ["--input", str(self.inp / "genefamilies.tsv"),
                                   "--action", "normalize", "--output", norm, "--unit", "relative"])
        obj = self.run_tool("differential", ["--input", norm, "--metadata", self.meta(),
                                             "--values", "relative", "--case", "case",
                                             "--control", "control", "--output_dir", self.o("fdiff")])
        self.assertEqual(obj["test"], "mannwhitneyu")
        ko = self.o("gf-ko.tsv")
        self.run_tool("function", ["--input", str(self.inp / "genefamilies.tsv"),
                                   "--action", "regroup", "--output", ko,
                                   "--mapping", str(self.inp / "mapping.tsv")])
        self.run_tool("differential", ["--input", ko, "--metadata", self.meta(),
                                       "--values", "relative", "--output_dir", self.o("fdiff-ko")])


if __name__ == "__main__":
    unittest.main()