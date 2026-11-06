"""Template catalog lookup for the omics visualization CLI.

Selection text (title, aliases, use_when, avoid_when, input_shape) comes from
references/plots.yaml and references/catalog/*.yaml, the same files the app
gallery reads. Input roles come from references/template_contracts.json.
Lookups are exact: nothing here ranks or guesses which template fits.
"""

from __future__ import annotations

import json
import os
import re
import sys
from typing import Any

sys.dont_write_bytecode = True

from viz_common import SKILL_ROOT
from viz_templates import installed_preview

REFERENCES = os.path.join(SKILL_ROOT, "references")
CONTRACTS = os.path.join(REFERENCES, "template_contracts.json")
ROLE_TYPES = ("number", "numbers", "any", "table")

_LIST_ITEM = re.compile(r"^-(?:\s+|$)")
_MAP_KEY = re.compile(r"^([A-Za-z_][A-Za-z0-9_-]*):(?:\s+(.*))?$")


class YamlSubsetError(ValueError):
    pass


def _strip_comment(text: str) -> str:
    index = text.find(" #")
    return text if index < 0 else text[:index].rstrip()


def _scalar(text: str, where: str) -> str:
    text = text.strip()
    if text.startswith('"'):
        end = 1
        while end < len(text):
            if text[end] == "\\":
                end += 2
                continue
            if text[end] == '"':
                break
            end += 1
        if end >= len(text):
            raise YamlSubsetError(f"{where}: unterminated double-quoted string")
        rest = text[end + 1 :].strip()
        if rest and not rest.startswith("#"):
            raise YamlSubsetError(f"{where}: unexpected text after a quoted string")
        return json.loads(text[: end + 1])
    if text.startswith(("'", "[", "{", "&", "*", "!", "|", ">")):
        raise YamlSubsetError(f"{where}: unsupported YAML syntax {text[:1]!r}")
    return _strip_comment(text)


class _Parser:
    """Block mappings, block lists, double-quoted or plain scalars, and `>`/`|` blocks."""

    def __init__(self, text: str, name: str) -> None:
        self.raw = text.splitlines()
        self.name = name
        self.index = 0

    def where(self, index: int | None = None) -> str:
        return f"{self.name}:{(self.index if index is None else index) + 1}"

    def peek(self) -> tuple[int, str] | None:
        while self.index < len(self.raw):
            line = self.raw[self.index]
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                if "\t" in line[: len(line) - len(line.lstrip())]:
                    raise YamlSubsetError(f"{self.where()}: tab indentation")
                return len(line) - len(line.lstrip(" ")), stripped
            self.index += 1
        return None

    def parse(self) -> Any:
        first = self.peek()
        if first is None:
            return None
        value = self.node(first[0])
        if self.peek() is not None:
            raise YamlSubsetError(f"{self.where()}: unexpected indentation")
        return value

    def node(self, indent: int) -> Any:
        line = self.peek()
        if line is None or line[0] < indent:
            return None
        if _LIST_ITEM.match(line[1]):
            return self.sequence(line[0])
        return self.mapping(line[0])

    def mapping(self, indent: int, first: str | None = None) -> dict[str, Any]:
        result: dict[str, Any] = {}
        while True:
            if first is not None:
                content, first = first, None
            else:
                line = self.peek()
                if line is None or line[0] < indent:
                    return result
                if line[0] > indent:
                    raise YamlSubsetError(f"{self.where()}: unexpected indentation")
                content = line[1]
                if _LIST_ITEM.match(content):
                    return result
            match = _MAP_KEY.match(content)
            if match is None:
                raise YamlSubsetError(f"{self.where()}: expected 'key: value'")
            key, value = match.group(1), (match.group(2) or "").strip()
            if key in result:
                raise YamlSubsetError(f"{self.where()}: duplicate key {key!r}")
            line_index = self.index
            self.index += 1
            if value in {">", ">-", ">+", "|", "|-", "|+"}:
                result[key] = self.block(indent, value)
            elif value and not value.startswith("#"):
                result[key] = _scalar(value, self.where(line_index))
            else:
                child = self.peek()
                if child is not None and child[0] == indent and _LIST_ITEM.match(child[1]):
                    result[key] = self.sequence(indent)
                elif child is not None and child[0] > indent:
                    result[key] = self.node(child[0])
                else:
                    result[key] = None

    def sequence(self, indent: int) -> list[Any]:
        result: list[Any] = []
        while True:
            line = self.peek()
            if line is None or line[0] != indent or not _LIST_ITEM.match(line[1]):
                if line is not None and line[0] > indent:
                    raise YamlSubsetError(f"{self.where()}: unexpected indentation")
                return result
            rest = line[1][1:].lstrip()
            if not rest or rest.startswith("#"):
                self.index += 1
                result.append(self.node(indent + 1))
            elif _MAP_KEY.match(rest) and not rest.startswith('"'):
                offset = len(line[1]) - len(rest)
                result.append(self.mapping(indent + offset, first=rest))
            else:
                result.append(_scalar(rest, self.where()))
                self.index += 1

    def block(self, parent_indent: int, style: str) -> str:
        lines: list[str] = []
        block_indent: int | None = None
        while self.index < len(self.raw):
            line = self.raw[self.index]
            if line.strip():
                indent = len(line) - len(line.lstrip(" "))
                if indent <= parent_indent:
                    break
                if block_indent is None:
                    block_indent = indent
                if indent < block_indent:
                    raise YamlSubsetError(f"{self.where()}: block text is less indented than its first line")
                lines.append(line[block_indent:])
            else:
                lines.append("")
            self.index += 1
        while lines and lines[-1] == "":
            lines.pop()
        if style.startswith("|"):
            text = "\n".join(lines)
        else:
            text = ""
            for position, line in enumerate(lines):
                if position == 0:
                    text = line
                elif line == "":
                    text += "\n"
                elif lines[position - 1] == "":
                    text += line
                else:
                    text += " " + line
        if not text:
            return ""
        return text if style.endswith("-") else text + "\n"


def load_yaml_subset(path: str) -> Any:
    with open(path, encoding="utf-8") as handle:
        text = handle.read()
    return _Parser(text, os.path.relpath(path, SKILL_ROOT)).parse()


def _text(value: Any) -> str:
    return value.strip() if isinstance(value, str) else ""


def _texts(value: Any) -> list[str]:
    return [item.strip() for item in value if isinstance(item, str)] if isinstance(value, list) else []


def _variants(value: Any) -> list[dict[str, str]]:
    variants: list[dict[str, str]] = []
    if not isinstance(value, list):
        return variants
    for item in value:
        if not isinstance(item, dict):
            continue
        image = _text(item.get("image"))
        if not image:
            continue
        variant = {
            "image": image,
            "config": _text(item.get("config")),
            "shows": _text(item.get("shows")),
        }
        example = _text(item.get("example"))
        if example:
            variant["example"] = example
        variants.append(variant)
    return variants


def load_families() -> list[dict[str, Any]]:
    document = load_yaml_subset(os.path.join(REFERENCES, "plots.yaml"))
    groups = document.get("groups") if isinstance(document, dict) else None
    if not isinstance(groups, dict):
        raise YamlSubsetError("references/plots.yaml has no groups")
    families: list[dict[str, Any]] = []
    for family, group in groups.items():
        if not isinstance(group, dict):
            continue
        families.append(
            {
                "family": family,
                "title": _text(group.get("title")),
                "use_when": _text(group.get("use_when")),
                "catalog": _text(group.get("catalog")),
            }
        )
    return families


def load_contracts() -> dict[str, Any]:
    with open(CONTRACTS, encoding="utf-8") as handle:
        contracts = json.load(handle)
    if not isinstance(contracts, dict):
        raise ValueError("references/template_contracts.json must be an object")
    return contracts


_CATALOG: list[dict[str, Any]] | None = None


def load_catalog() -> list[dict[str, Any]]:
    """Every template in family order, merged with its input roles."""
    global _CATALOG
    if _CATALOG is not None:
        return _CATALOG
    contracts = load_contracts()
    by_id = {
        str(item.get("id")): item
        for item in contracts.get("templates", [])
        if isinstance(item, dict) and item.get("id")
    }
    templates: list[dict[str, Any]] = []
    for family in load_families():
        document = load_yaml_subset(os.path.join(SKILL_ROOT, family["catalog"]))
        records = document.get("templates") if isinstance(document, dict) else None
        for record in records if isinstance(records, list) else []:
            if not isinstance(record, dict) or not _text(record.get("id")):
                continue
            template_id = _text(record.get("id"))
            contract = by_id.get(template_id, {})
            roles = contract.get("roles") if isinstance(contract.get("roles"), dict) else {}
            templates.append(
                {
                    "template_id": template_id,
                    "family": family["family"],
                    "title": _text(record.get("title")),
                    "aliases": _texts(record.get("aliases")),
                    "use_when": _text(record.get("use_when")),
                    "avoid_when": _text(record.get("avoid_when")),
                    "input_shape": _text(record.get("input_shape")),
                    "examples": _texts(record.get("examples")),
                    "source": _text(record.get("source")),
                    "preview": _text(record.get("preview")),
                    "variants": _variants(record.get("variants")),
                    "status": _text(record.get("status")),
                    "roles": {
                        "required": _texts(roles.get("required")),
                        "optional": _texts(roles.get("optional")),
                    },
                    "risks": _texts(contract.get("risks")),
                }
            )
    _CATALOG = templates
    return templates


SUMMARY_KEYS = ("template_id", "title", "aliases", "use_when", "avoid_when", "input_shape")


def summary(template: dict[str, Any]) -> dict[str, Any]:
    """What choosing between templates needs, without paths."""
    return {key: template[key] for key in SUMMARY_KEYS}


def entry(template: dict[str, Any]) -> dict[str, Any]:
    """The full public form of one template: roles, example data, and shipped preview."""
    result = summary(template)
    result["family"] = template["family"]
    result["roles"] = template["roles"]
    result["examples"] = [os.path.join(SKILL_ROOT, path) for path in template["examples"]]
    preview = installed_preview({"preview": template["preview"]})
    if preview is not None:
        result["preview"] = preview
        result["preview_markdown"] = f"![{template['template_id']}]({preview})"
    variants = []
    for item in template.get("variants") or []:
        image = installed_preview({"preview": item["image"]})
        if image is None:
            continue
        stem = os.path.splitext(os.path.basename(image))[0]
        name = stem[len("preview_") :] if stem.startswith("preview_") else stem
        variant = {
            "config": item["config"],
            "shows": item["shows"],
            "preview": image,
            "preview_markdown": f"![{template['template_id']}:{name}]({image})",
        }
        example = item.get("example")
        if example:
            variant["example"] = os.path.join(SKILL_ROOT, example)
        variants.append(variant)
    if variants:
        result["variants"] = variants
    if template["risks"]:
        result["risks"] = template["risks"]
    return result


def normalize_name(value: str) -> str:
    return re.sub(r"\s+", " ", value).strip().casefold()


def lookup_alias(name: str) -> list[dict[str, Any]]:
    wanted = normalize_name(name)
    return [
        template
        for template in load_catalog()
        if wanted
        in {normalize_name(template["template_id"]), normalize_name(template["title"])}
        | {normalize_name(alias) for alias in template["aliases"]}
    ]
