#!/usr/bin/env python3
"""Merge multi-step protein design metrics by lineage ID, filter with a
Chaperone-compatible boolean expression DSL, and rank/select diverse candidates.

Pure-Python implementation (no third-party dependencies required):
  - Multi-table left/outer join on design_id (supports TSV, CSV, JSON arrays)
  - Safe recursive-descent boolean DSL (`and`, `or`, `not`, parentheses,
    `>= <= > < == !=`, dotted keys, backtick-quoted keys, negative numbers)
  - Three-valued null semantics matching Chaperone/Polars FilterBlock:
      True or None  -> True
      False and None -> False
      overall None  -> dropped
  - Ranking (`min`/`max`), per-group diversity cap (`--group-by` + `--max-per-group`),
    and `--top-n` truncation
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


# --- Table loading ---

def _load_table(path: Path) -> list[dict[str, Any]]:
    suffix = path.suffix.lower()
    if suffix == ".json":
        data = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(data, list):
            raise ValueError(f"JSON table must be a list of objects: {path}")
        return [dict(row) for row in data]
    delimiter = "\t" if suffix in (".tsv", ".tab") else ","
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        return [dict(row) for row in reader]


def _parse_table_spec(spec: str) -> tuple[str | None, Path]:
    if ":" in spec:
        maybe_prefix, maybe_path = spec.split(":", 1)
        if maybe_prefix and not ("/" in maybe_prefix or "\\" in maybe_prefix):
            return maybe_prefix.strip(), Path(maybe_path)
    return None, Path(spec)


def _to_number(value: Any) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return None if math.isnan(value) or math.isinf(value) else float(value)
    s = str(value).strip()
    if not s or s.lower() in ("na", "nan", "null", "none"):
        return None
    try:
        v = float(s)
        return None if math.isnan(v) or math.isinf(v) else v
    except ValueError:
        return None


def merge_tables(
    table_specs: list[str],
    id_column: str,
    join_mode: str = "outer",
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    merged: dict[str, dict[str, Any]] = {}
    order: list[str] = []
    source_counts: dict[str, int] = {}

    for idx, spec in enumerate(table_specs):
        prefix, path = _parse_table_spec(spec)
        rows = _load_table(path)
        source_counts[str(path)] = len(rows)

        for row in rows:
            raw_id = row.get(id_column)
            if raw_id is None or str(raw_id).strip() == "":
                continue
            design_id = str(raw_id).strip()
            if design_id not in merged:
                if idx > 0 and join_mode == "left":
                    continue
                merged[design_id] = {id_column: design_id}
                order.append(design_id)

            target = merged[design_id]
            for k, v in row.items():
                if k == id_column:
                    continue
                if prefix:
                    out_key = f"{prefix}.{k}"
                else:
                    out_key = k
                if v is not None and str(v).strip() != "":
                    target[out_key] = v
                elif out_key not in target:
                    target[out_key] = ""

    return [merged[k] for k in order], source_counts


# --- Safe recursive-descent parser for Chaperone FilterBlock expr DSL ---

class ExprSyntaxError(ValueError):
    """Raised when a filter expression cannot be parsed."""


Token = tuple[str, str]  # (kind, value)


def tokenize_expr(expr: str) -> list[Token]:
    tokens: list[Token] = []
    i = 0
    n = len(expr)
    while i < n:
        ch = expr[i]
        if ch.isspace():
            i += 1
            continue
        if ch == "(":
            tokens.append(("LPAREN", "("))
            i += 1
            continue
        if ch == ")":
            tokens.append(("RPAREN", ")"))
            i += 1
            continue
        if expr[i : i + 2] in (">=", "<=", "==", "!="):
            tokens.append(("OP", expr[i : i + 2]))
            i += 2
            continue
        if ch in (">", "<"):
            tokens.append(("OP", ch))
            i += 1
            continue
        if ch == "`":
            end = expr.find("`", i + 1)
            if end == -1:
                raise ExprSyntaxError("Unterminated backtick identifier")
            tokens.append(("IDENT", expr[i + 1 : end]))
            i = end + 1
            continue
        if ch in ("'", '"'):
            quote = ch
            j = i + 1
            buf: list[str] = []
            while j < n and expr[j] != quote:
                if expr[j] == "\\" and j + 1 < n:
                    j += 1
                buf.append(expr[j])
                j += 1
            if j >= n:
                raise ExprSyntaxError("Unterminated string literal")
            tokens.append(("STRING", "".join(buf)))
            i = j + 1
            continue
        # Signed or unsigned number
        if ch.isdigit() or (ch in ("-", "+") and i + 1 < n and (expr[i + 1].isdigit() or expr[i + 1] == ".")):
            j = i + 1 if ch in ("-", "+") else i
            has_dot = False
            while j < n and (expr[j].isdigit() or (expr[j] == "." and not has_dot)):
                if expr[j] == ".":
                    has_dot = True
                j += 1
            if j < n and expr[j] in ("e", "E"):
                j += 1
                if j < n and expr[j] in ("+", "-"):
                    j += 1
                while j < n and expr[j].isdigit():
                    j += 1
            tokens.append(("NUMBER", expr[i:j]))
            i = j
            continue
        # Bare identifier or boolean keyword (supports dots, slashes, hyphens inside names)
        j = i
        while j < n and (expr[j].isalnum() or expr[j] in ("_", ".", "/", "-")):
            j += 1
        if j == i:
            raise ExprSyntaxError(f"Unexpected character {ch!r} at position {i}")
        word = expr[i:j]
        lower = word.lower()
        if lower in ("and", "or", "not"):
            tokens.append((lower.upper(), lower))
        else:
            tokens.append(("IDENT", word))
        i = j
    return tokens


class _ExprParser:
    def __init__(self, tokens: list[Token]) -> None:
        self.tokens = tokens
        self.pos = 0

    def peek(self) -> Token | None:
        return self.tokens[self.pos] if self.pos < len(self.tokens) else None

    def consume(self, expected_kind: str | None = None) -> Token:
        tok = self.peek()
        if tok is None:
            raise ExprSyntaxError("Unexpected end of expression")
        if expected_kind and tok[0] != expected_kind:
            raise ExprSyntaxError(f"Expected {expected_kind}, got {tok[0]} ({tok[1]!r})")
        self.pos += 1
        return tok

    def parse(self) -> Any:
        node = self.parse_or()
        if self.peek() is not None:
            raise ExprSyntaxError(f"Unexpected trailing token: {self.peek()}")
        return node

    def parse_or(self) -> Any:
        node = self.parse_and()
        while self.peek() and self.peek()[0] == "OR":  # type: ignore[index]
            self.consume("OR")
            rhs = self.parse_and()
            node = ("or", node, rhs)
        return node

    def parse_and(self) -> Any:
        node = self.parse_not()
        while self.peek() and self.peek()[0] == "AND":  # type: ignore[index]
            self.consume("AND")
            rhs = self.parse_not()
            node = ("and", node, rhs)
        return node

    def parse_not(self) -> Any:
        if self.peek() and self.peek()[0] == "NOT":  # type: ignore[index]
            self.consume("NOT")
            return ("not", self.parse_not())
        return self.parse_primary()

    def parse_primary(self) -> Any:
        tok = self.peek()
        if tok is None:
            raise ExprSyntaxError("Expected condition or '('")
        if tok[0] == "LPAREN":
            self.consume("LPAREN")
            node = self.parse_or()
            self.consume("RPAREN")
            return node
        # Comparison: IDENT OP (NUMBER | STRING | IDENT)
        lhs = self.consume("IDENT")[1]
        op = self.consume("OP")[1]
        rhs_tok = self.consume()
        if rhs_tok[0] == "NUMBER":
            rhs: Any = ("num", float(rhs_tok[1]))
        elif rhs_tok[0] == "STRING":
            rhs = ("str", rhs_tok[1])
        elif rhs_tok[0] == "IDENT":
            rhs = ("ident", rhs_tok[1])
        else:
            raise ExprSyntaxError(f"Invalid right-hand operand: {rhs_tok}")
        return ("cmp", op, lhs, rhs)


def compile_expr(expr: str) -> Any:
    expr = expr.strip()
    if not expr:
        return None
    tokens = tokenize_expr(expr)
    return _ExprParser(tokens).parse()


def _three_val_and(a: bool | None, b: bool | None) -> bool | None:
    if a is False or b is False:
        return False
    if a is True and b is True:
        return True
    return None


def _three_val_or(a: bool | None, b: bool | None) -> bool | None:
    if a is True or b is True:
        return True
    if a is False and b is False:
        return False
    return None


def eval_ast(node: Any, row: dict[str, Any]) -> bool | None:
    if node is None:
        return True
    tag = node[0]
    if tag == "and":
        return _three_val_and(eval_ast(node[1], row), eval_ast(node[2], row))
    if tag == "or":
        return _three_val_or(eval_ast(node[1], row), eval_ast(node[2], row))
    if tag == "not":
        v = eval_ast(node[1], row)
        return None if v is None else (not v)
    if tag == "cmp":
        _, op, lhs_key, (rhs_kind, rhs_val) = node
        raw_lhs = row.get(lhs_key)
        if raw_lhs is None or str(raw_lhs).strip() == "":
            return None
        if rhs_kind == "str":
            l_str = str(raw_lhs).strip()
            if op == "==":
                return l_str == rhs_val
            if op == "!=":
                return l_str != rhs_val
            return None
        l_num = _to_number(raw_lhs)
        r_num = rhs_val if rhs_kind == "num" else _to_number(row.get(rhs_val))
        if l_num is None or r_num is None:
            return None
        if op == ">=":
            return l_num >= r_num
        if op == "<=":
            return l_num <= r_num
        if op == ">":
            return l_num > r_num
        if op == "<":
            return l_num < r_num
        if op == "==":
            return l_num == r_num
        if op == "!=":
            return l_num != r_num
    raise ExprSyntaxError(f"Unknown AST node: {node}")


def rank_and_select(
    rows: list[dict[str, Any]],
    rank_by: str,
    rank_direction: str,
    group_by: str,
    max_per_group: int | None,
    top_n: int | None,
) -> list[dict[str, Any]]:
    selected = list(rows)
    if rank_by:
        reverse = rank_direction.lower() == "max"

        def sort_key(r: dict[str, Any]) -> tuple[int, float]:
            val = _to_number(r.get(rank_by))
            if val is None:
                return (1, 0.0)
            return (0, -val if reverse else val)

        selected.sort(key=sort_key)

    if group_by and max_per_group is not None and max_per_group > 0:
        counts: dict[str, int] = {}
        diverse: list[dict[str, Any]] = []
        for r in selected:
            grp = str(r.get(group_by, "")).strip() or "__ungrouped__"
            c = counts.get(grp, 0)
            if c < max_per_group:
                counts[grp] = c + 1
                diverse.append(r)
        selected = diverse

    if top_n is not None and top_n > 0:
        selected = selected[:top_n]

    return selected


def _write_tsv(path: Path, rows: list[dict[str, Any]], fallback_columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    seen: set[str] = set()
    for col in fallback_columns:
        if col not in seen:
            fieldnames.append(col)
            seen.add(col)
    for row in rows:
        for col in row.keys():
            if col not in seen:
                fieldnames.append(col)
                seen.add(col)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if row.get(k) is None else row.get(k)) for k in fieldnames})


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge multi-step protein design metrics, filter via boolean DSL, and select top diverse candidates."
    )
    parser.add_argument(
        "--tables",
        nargs="+",
        required=True,
        help="Input metric tables (TSV/CSV/JSON). Optional 'prefix:path' adds step prefix to column names.",
    )
    parser.add_argument("--id-column", default="design_id", help="Lineage join key column (default: design_id).")
    parser.add_argument("--join-mode", choices=("outer", "left"), default="outer", help="Table join mode.")
    parser.add_argument("--expr", default="", help="Chaperone-compatible boolean filter expression.")
    parser.add_argument("--rank-by", default="", help="Optional numeric column to rank passing designs.")
    parser.add_argument("--rank-direction", choices=("max", "min"), default="max", help="Ranking direction.")
    parser.add_argument("--group-by", default="", help="Optional grouping column for diversity cap (e.g. parent_design_id).")
    parser.add_argument("--max-per-group", type=int, default=0, help="Max selected designs per group (0 = unlimited).")
    parser.add_argument("--top-n", type=int, default=0, help="Keep top N selected designs (0 = all).")
    parser.add_argument("--out-merged", required=True, help="Output path for merged lineage TSV.")
    parser.add_argument("--out-filtered", required=True, help="Output path for filtered/selected TSV.")
    parser.add_argument("--out-summary", required=True, help="Output path for funnel summary JSON.")
    args = parser.parse_args()

    merged_rows, source_counts = merge_tables(args.tables, id_column=args.id_column, join_mode=args.join_mode)
    ast = compile_expr(args.expr)

    passing_rows: list[dict[str, Any]] = []
    for row in merged_rows:
        verdict = eval_ast(ast, row) is True
        row["passed_filter"] = "true" if verdict else "false"
        if verdict:
            passing_rows.append(dict(row))

    selected_rows = rank_and_select(
        rows=passing_rows,
        rank_by=args.rank_by,
        rank_direction=args.rank_direction,
        group_by=args.group_by,
        max_per_group=args.max_per_group if args.max_per_group > 0 else None,
        top_n=args.top_n if args.top_n > 0 else None,
    )

    _write_tsv(Path(args.out_merged), merged_rows, [args.id_column, "passed_filter"])
    _write_tsv(Path(args.out_filtered), selected_rows, [args.id_column, "passed_filter"])

    summary = {
        "id_column": args.id_column,
        "join_mode": args.join_mode,
        "expr": args.expr,
        "rank_by": args.rank_by,
        "rank_direction": args.rank_direction,
        "group_by": args.group_by,
        "max_per_group": args.max_per_group,
        "top_n": args.top_n,
        "source_row_counts": source_counts,
        "total_designs": len(merged_rows),
        "passed_filter": len(passing_rows),
        "selected_designs": len(selected_rows),
    }
    out_summary = Path(args.out_summary)
    out_summary.parent.mkdir(parents=True, exist_ok=True)
    out_summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
