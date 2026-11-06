"""Facts about input tables for the omics visualization CLI.

Everything here is counted, not guessed: column types, missing values, ranges,
distinct values, identifier overlap between tables, and whether an
agent-supplied role=column mapping fits a template's declared roles.
Missing values follow readr's default: an empty cell or "NA".
"""

from __future__ import annotations

import csv
import math
import os
import sys
from typing import Any, Iterator

sys.dont_write_bytecode = True

TABLE_EXTENSIONS = {".tsv", ".tab", ".csv", ".txt"}
MISSING = {"", "NA"}
MAX_DETAILED_COLUMNS = 200
MAX_DISTINCT = 10_000
MAX_EXAMPLES = 5

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))


def is_table(path: str) -> bool:
    return os.path.splitext(path)[1].lower() in TABLE_EXTENSIONS


def delimiter_of(path: str) -> str:
    suffix = os.path.splitext(path)[1].lower()
    if suffix in {".tsv", ".tab"}:
        return "\t"
    if suffix == ".csv":
        return ","
    with open(path, encoding="utf-8-sig", errors="replace", newline="") as handle:
        sample = handle.read(8192)
    lines = sample.splitlines()
    first = lines[0] if lines else ""
    if first.count("\t") > first.count(","):
        return "\t"
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,;").delimiter
    except csv.Error:
        return ","


def read_table(path: str) -> tuple[list[str], Iterator[list[str]]]:
    handle = open(path, encoding="utf-8-sig", errors="replace", newline="")
    reader = csv.reader(handle, delimiter=delimiter_of(path))
    try:
        header = [cell.strip() for cell in next(reader)]
    except StopIteration:
        handle.close()
        return [], iter(())

    def rows() -> Iterator[list[str]]:
        try:
            for row in reader:
                if row:
                    yield row
        finally:
            handle.close()

    return header, rows()


def read_header(path: str) -> list[str]:
    if not is_table(path):
        return []
    with open(path, encoding="utf-8-sig", errors="replace", newline="") as handle:
        try:
            return [cell.strip() for cell in next(csv.reader(handle, delimiter=delimiter_of(path)))]
        except StopIteration:
            return []


def parse_number(value: str) -> float | None:
    try:
        return float(value)
    except ValueError:
        return None


def json_float(value: float) -> float | int:
    return int(value) if value.is_integer() else value


class ColumnStats:
    def __init__(self, name: str, detailed: bool) -> None:
        self.name = name
        self.detailed = detailed
        self.missing = 0
        self.numeric = 0
        self.text = 0
        self.non_finite = 0
        self.minimum: float | None = None
        self.maximum: float | None = None
        self.distinct: set[str] = set()
        self.distinct_capped = False
        self.examples: list[str] = []
        self.text_examples: list[str] = []

    def add(self, value: str) -> None:
        value = value.strip()
        if value in MISSING:
            self.missing += 1
            return
        number = parse_number(value)
        if number is None:
            self.text += 1
            if len(self.text_examples) < MAX_EXAMPLES and value not in self.text_examples:
                self.text_examples.append(value)
        else:
            self.numeric += 1
            if math.isfinite(number):
                self.minimum = number if self.minimum is None else min(self.minimum, number)
                self.maximum = number if self.maximum is None else max(self.maximum, number)
            else:
                self.non_finite += 1
        if not self.detailed:
            return
        if value not in self.distinct:
            if len(self.distinct) >= MAX_DISTINCT:
                self.distinct_capped = True
            else:
                self.distinct.add(value)
                if len(self.examples) < MAX_EXAMPLES:
                    self.examples.append(value)

    @property
    def kind(self) -> str:
        if self.numeric == 0 and self.text == 0:
            return "empty"
        return "number" if self.text == 0 else "text"

    def summary(self) -> dict[str, Any]:
        present = self.numeric + self.text
        result: dict[str, Any] = {"name": self.name, "type": self.kind, "missing": self.missing}
        if self.kind == "text" and self.numeric:
            result["numeric_values"] = self.numeric
        if self.minimum is not None and self.maximum is not None:
            result["min"] = json_float(self.minimum)
            result["max"] = json_float(self.maximum)
        if self.non_finite:
            result["non_finite"] = self.non_finite
        if self.detailed:
            if self.distinct_capped:
                result["distinct"] = f">{MAX_DISTINCT}"
            else:
                result["distinct"] = len(self.distinct)
                result["unique"] = present > 0 and len(self.distinct) == present
            result["examples"] = self.examples
        return result


def profile_table(path: str) -> dict[str, Any]:
    name = os.path.basename(path)
    if not is_table(path):
        return {"name": name, "path": path, "kind": "file", "bytes": os.path.getsize(path)}
    header, rows = read_table(path)
    columns = [ColumnStats(column, index < MAX_DETAILED_COLUMNS) for index, column in enumerate(header)]
    count = 0
    ragged = 0
    for row in rows:
        count += 1
        if len(row) != len(header):
            ragged += 1
        for index, stats in enumerate(columns):
            stats.add(row[index] if index < len(row) else "")
    result: dict[str, Any] = {
        "name": name,
        "path": path,
        "kind": "table",
        "rows": count,
        "column_count": len(header),
        "columns": [stats.summary() for stats in columns[:MAX_DETAILED_COLUMNS]],
    }
    if len(columns) > MAX_DETAILED_COLUMNS:
        rest = columns[MAX_DETAILED_COLUMNS:]
        result["more_columns"] = {
            "count": len(rest),
            "number": sum(1 for stats in rest if stats.kind == "number"),
            "text": sum(1 for stats in rest if stats.kind == "text"),
            "empty": sum(1 for stats in rest if stats.kind == "empty"),
        }
    duplicated = sorted({column for column in header if header.count(column) > 1})
    if duplicated:
        result["duplicate_headers"] = duplicated
    if ragged:
        result["ragged_rows"] = ragged
    return result


def sibling_files(paths: list[str]) -> list[str]:
    given = {os.path.realpath(path) for path in paths}
    found: set[str] = set()
    for directory in {os.path.dirname(os.path.realpath(path)) for path in paths}:
        try:
            names = os.listdir(directory)
        except OSError:
            continue
        for name in names:
            full = os.path.join(directory, name)
            if name.startswith(".") or not os.path.isfile(full) or os.path.realpath(full) in given:
                continue
            found.add(full)
    return sorted(found)


def column_values(path: str, column: str) -> list[str]:
    header, rows = read_table(path)
    index = header.index(column)
    return [row[index].strip() if index < len(row) else "" for row in rows]


def tables_by_name(paths: list[str]) -> dict[str, str]:
    names: dict[str, str] = {}
    for path in paths:
        name = os.path.basename(path)
        if name in names and os.path.realpath(names[name]) != os.path.realpath(path):
            raise ValueError(f"Two input tables are both named {name}; give them distinct file names.")
        names[name] = path
    return names


def parse_join(spec: str, names: dict[str, str]) -> tuple[tuple[str, str], tuple[str, str]]:
    left, separator, right = spec.partition("=")
    if not separator:
        raise ValueError(f'A join must look like "links.tsv:from=nodes.tsv:id"; got "{spec}".')
    sides: list[tuple[str, str]] = []
    for side in (left, right):
        table, colon, column = side.strip().partition(":")
        if not colon or not table or not column:
            raise ValueError(f'A join must look like "links.tsv:from=nodes.tsv:id"; got "{spec}".')
        if table not in names:
            raise ValueError(f'The join names {table}, which is not one of the given tables: {", ".join(names)}.')
        header = read_header(names[table])
        if column not in header:
            raise ValueError(f"{table} has no column {column}. Its columns are: {', '.join(header)}.")
        sides.append((table, column))
    return sides[0], sides[1]


def join_overlap(spec: str, names: dict[str, str]) -> dict[str, Any]:
    (left_table, left_column), (right_table, right_column) = parse_join(spec, names)
    left = {value for value in column_values(names[left_table], left_column) if value not in MISSING}
    right = {value for value in column_values(names[right_table], right_column) if value not in MISSING}
    left_only = sorted(left - right)
    right_only = sorted(right - left)
    return {
        "join": spec,
        "left_distinct": len(left),
        "right_distinct": len(right),
        "shared": len(left & right),
        "left_only": len(left_only),
        "right_only": len(right_only),
        "left_only_examples": left_only[:MAX_EXAMPLES],
        "right_only_examples": right_only[:MAX_EXAMPLES],
    }


def profile(paths: list[str], joins: list[str]) -> dict[str, Any]:
    names = tables_by_name(paths)
    result: dict[str, Any] = {"tables": [profile_table(path) for path in paths]}
    siblings = sibling_files(paths)
    if siblings:
        result["other_files"] = [os.path.basename(path) for path in siblings]
    if joins:
        result["joins"] = [join_overlap(spec, names) for spec in joins]
    return result


def parse_mapping(entries: list[str]) -> dict[str, list[str]]:
    mapping: dict[str, list[str]] = {}
    for item in entries:
        role, separator, column = item.partition("=")
        role, column = role.strip(), column.strip()
        if not separator or not role or not column:
            raise ValueError(f'Each column mapping must look like "role=column"; got "{item}".')
        mapping.setdefault(role, [])
        if column not in mapping[role]:
            mapping[role].append(column)
    return mapping


def check_columns(
    template_roles: dict[str, list[str]],
    role_types: dict[str, str],
    tables: list[str],
    entries: list[str],
) -> dict[str, Any]:
    """Check an agent's role=column mapping. Problems are reported, never fixed."""
    mapping = parse_mapping(entries)
    required = template_roles.get("required", [])
    optional = template_roles.get("optional", [])
    headers = {path: read_header(path) for path in tables}
    problems: list[str] = []
    mapped: list[dict[str, Any]] = []
    for role, columns in mapping.items():
        if role not in required and role not in optional:
            problems.append(f"{role} is not a role of this template (roles: {', '.join(required + optional)}).")
            continue
        kind = role_types.get(role, "any")
        if len(columns) > 1 and kind != "numbers":
            problems.append(f"{role} takes one column; got {', '.join(columns)}.")
        for column in columns:
            found = [path for path in tables if column in headers[path]]
            record: dict[str, Any] = {"role": role, "column": column, "expects": kind}
            if not found:
                problems.append(f"{role}: no input table has a column named {column}.")
                record["tables"] = []
                mapped.append(record)
                continue
            record["tables"] = [os.path.basename(path) for path in found]
            if kind in {"number", "numbers"}:
                for path in found:
                    stats = ColumnStats(column, detailed=False)
                    for value in column_values(path, column):
                        stats.add(value)
                    if stats.text:
                        examples = ", ".join(stats.text_examples)
                        problems.append(
                            f"{role} expects numbers; {os.path.basename(path)}:{column} has "
                            f"{stats.text} non-numeric value(s), e.g. {examples}."
                        )
            mapped.append(record)
    not_checked = [role for role in required if role_types.get(role) == "table"]
    missing = [
        role
        for role in required
        if role not in mapping and role_types.get(role) not in {"table", "numbers"}
    ]
    for role in missing:
        problems.append(f"Required role {role} has no column.")
    unmapped_numbers = [role for role in required if role_types.get(role) == "numbers" and role not in mapping]
    result: dict[str, Any] = {"ok": not problems, "mapped": mapped, "problems": problems}
    if not_checked or unmapped_numbers:
        result["not_checked"] = not_checked + unmapped_numbers
    return result
