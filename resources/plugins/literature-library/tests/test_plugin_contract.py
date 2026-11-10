#!/usr/bin/env python3
"""Structure checks: manifest, skill frontmatter (skill contract 1.1.0), and that every
declared tool argument is a real flag of the script it runs."""
from __future__ import annotations

import re
import subprocess
import sys
import unittest
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
SKILLS = ROOT / "skills"
NAME_RE = re.compile(r"^[a-z0-9][a-z0-9-]{0,63}$")
TOOL_RE = re.compile(r"^[a-z][a-z0-9_]{0,31}$")
ARG_RE = re.compile(r"^[a-z][a-z0-9_]*$")
PATH_FORMATS = {"input-path", "project-path"}
TYPES = {"string", "number", "integer", "boolean", "array"}


def frontmatter(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    assert text.startswith("---\n"), f"{path} needs frontmatter"
    return yaml.safe_load(text.split("\n---\n", 1)[0][4:])


def script_flags(script: Path) -> set[str]:
    out = subprocess.run([sys.executable, str(script), "--help"], capture_output=True, text=True, check=False)
    assert out.returncode == 0, f"{script.name} --help failed: {out.stderr[-300:]}"
    return {f.lstrip("-") for f in re.findall(r"(?<![\w-])(--[a-z][a-z0-9_-]*)", out.stdout)}


class ManifestTests(unittest.TestCase):
    def setUp(self):
        self.manifest = yaml.safe_load((ROOT / "phi-package.yaml").read_text(encoding="utf-8"))

    def test_components_exist(self):
        self.assertEqual(self.manifest["type"], "plugin")
        for rel in self.manifest["components"]["skills"]:
            self.assertTrue((ROOT / rel / "SKILL.md").is_file(), rel)
        for rel in self.manifest["components"]["agents"]:
            self.assertTrue((ROOT / rel).is_file(), rel)

    def test_agent_lists_the_manifest_skills(self):
        agent = frontmatter(ROOT / self.manifest["components"]["agents"][0])
        in_manifest = {Path(s).name for s in self.manifest["components"]["skills"]}
        self.assertEqual(set(agent["skills"]), in_manifest)

    def test_no_orphan_skill_directories(self):
        in_manifest = {Path(s).name for s in self.manifest["components"]["skills"]}
        on_disk = {p.name for p in SKILLS.iterdir() if p.is_dir()}
        self.assertEqual(on_disk, in_manifest)


class SkillTests(unittest.TestCase):
    def skills(self):
        return sorted(p for p in SKILLS.iterdir() if p.is_dir())

    def test_frontmatter(self):
        for d in self.skills():
            fm = frontmatter(d / "SKILL.md")
            self.assertEqual(fm["name"], d.name)
            self.assertRegex(fm["name"], NAME_RE)
            self.assertTrue(1 <= len(fm["description"].strip()) <= 1024, d.name)

    def test_scripts_need_an_environment(self):
        for d in self.skills():
            fm = frontmatter(d / "SKILL.md")
            has_scripts = any(d.glob("scripts/**/*.py"))
            self.assertEqual(has_scripts, bool(fm.get("phi", {}).get("environment")), d.name)

    def test_no_script_imports_another_skill(self):
        for d in self.skills():
            for py in d.glob("scripts/**/*.py"):
                text = py.read_text(encoding="utf-8")
                self.assertNotRegex(text, r"parents\[\d\][^\n]*skills|/ \"lit-", str(py))

    def test_tools_match_their_scripts(self):
        seen = set()
        for d in self.skills():
            fm = frontmatter(d / "SKILL.md")
            for tool in fm.get("phi", {}).get("scripts", []):
                name = tool["name"]
                self.assertRegex(name, TOOL_RE)
                self.assertNotIn(name, seen, f"duplicate tool {name}")
                seen.add(name)
                self.assertTrue(1 <= len(tool["description"].strip()) <= 1024, name)
                self.assertIn(tool["approval"], ("read", "write"))
                self.assertEqual(tool["run"][0], "python")
                script = d / tool["run"][1]
                self.assertTrue(script.is_file(), f"{name}: {script}")
                args = tool["args"]
                self.assertEqual(args["type"], "object")
                self.assertIs(args["additionalProperties"], False)
                props = args["properties"]
                flags = script_flags(script)
                for prop, spec in props.items():
                    self.assertRegex(prop, ARG_RE)
                    self.assertIn(spec["type"], TYPES, f"{name}.{prop}")
                    if "format" in spec:
                        self.assertIn(spec["format"], PATH_FORMATS)
                        self.assertEqual(spec["type"], "string")
                    self.assertIn(prop, flags, f"{name}: argument {prop} is not a flag of {script.name}")
                for req in args.get("required", []):
                    self.assertIn(req, props)
                if tool["approval"] == "read":
                    for prop, spec in props.items():
                        self.assertNotEqual(spec.get("format"), "project-path", f"{name}.{prop}")
        self.assertEqual(len(seen), 12, sorted(seen))

    def test_skill_text_mentions_every_tool(self):
        fm = frontmatter(SKILLS / "manage-library" / "SKILL.md")
        text = (SKILLS / "manage-library" / "SKILL.md").read_text(encoding="utf-8").split("\n---\n", 2)[-1]
        for tool in fm["phi"]["scripts"]:
            self.assertIn(f"lit_{tool['name']}", text + (SKILLS / "compile-wiki" / "SKILL.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
