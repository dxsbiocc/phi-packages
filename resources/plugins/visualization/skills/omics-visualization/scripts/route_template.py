#!/usr/bin/env python3
"""Recommend omics visualization templates from lightweight data contracts.

This script is intentionally dependency-free so it can run before R packages or
plotting libraries are installed. It shortlists templates; it does not replace
the full catalog review required for low-confidence or publication decisions.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path
from typing import Any

DEFAULT_CONTRACTS = Path(__file__).resolve().parents[1] / "references" / "template_contracts.json"
SIDECAR_STEMS = {
    "annot",
    "colinfo",
    "cytoband",
    "domains",
    "edges",
    "enrichment",
    "features",
    "genes",
    "groups",
    "karyotype",
    "links",
    "loops",
    "nodes",
    "points",
    "rowinfo",
    "windows",
}
SIDECAR_SUFFIXES = {".csv", ".tab", ".tsv"}
SIDECAR_ROLE_EXPORTS = {
    "alteration_type",
    "axis_entity",
    "category",
    "chromosome",
    "count",
    "cytoband",
    "direction",
    "feature_type",
    "flow_value",
    "genomic_end",
    "genomic_start",
    "mutation_class",
    "numeric_x",
    "numeric_y",
    "pair_id",
    "residue_position",
    "significance",
    "signed_association",
    "site_position",
    "source_entity",
    "target_entity",
    "track",
    "weight",
    "zoom_end",
    "zoom_start",
}
SIDECAR_DEPENDENT_SHAPES = {
    "annotation_sidecars",
    "enrichment_sidecar",
    "enrichment_zoom_sidecars",
    "genomic_annotation_sidecars",
    "node_link_sidecars",
    "protein_domain_sidecar",
    "synteny_sidecars",
}
SIDECAR_SHAPE_INTENT_KEYWORDS = {
    "annotation_sidecars": ("annotation", "annotated", "grouped", "注释", "分组"),
    "enrichment_sidecar": ("enrichment", "pathway", "go", "kegg", "富集", "通路"),
    "enrichment_zoom_sidecars": ("enrichment", "zoom", "aligned", "pathway", "富集", "通路", "对齐"),
    "genomic_annotation_sidecars": ("genomic", "locus", "coverage", "track", "基因组", "位点"),
    "node_link_sidecars": ("node", "link", "edge", "network", "graph", "节点", "连边", "网络"),
    "protein_domain_sidecar": ("protein", "domain", "residue", "lollipop", "蛋白", "结构域"),
    "synteny_sidecars": ("synteny", "homology", "assembly", "共线性", "同源"),
}

ROLE_HINTS: dict[str, tuple[str, ...]] = {
    "focal_entity": (
        "focal_entity",
        "focus_entity",
        "central_entity",
        "focal",
        "focus",
        "center",
        "centre",
        "central",
        "hub",
        "anchor",
        "query",
        "seed",
        "primary",
        "index",
        "reference",
    ),
    "source_entity": (
        "source_entity",
        "source_node",
        "source_feature",
        "source",
        "from",
        "node1",
        "node_a",
        "entity1",
        "feature1",
    ),
    "parent_entity": (
        "parent_entity",
        "parent_node",
        "parent",
        "ancestor",
        "root",
        "path",
        "level",
        "level1",
        "level_1",
        "from",
    ),
    "axis_entity": (
        "axis",
        "axis_entity",
        "dimension",
        "metric",
        "variable",
        "name",
        "label",
    ),
    "target_entity": (
        "target_entity",
        "related_entity",
        "related_item",
        "target_item",
        "target",
        "related",
        "item",
        "label",
        "feature",
        "entity",
        "name",
        "symbol",
        "term",
        "pathway",
        "protein",
        "metabolite",
        "transcript",
        "gene",
        "node",
        "to",
        "description",
    ),
    "x_category": (
        "x",
        "x_axis",
        "axis_x",
        "var1",
        "variable1",
        "feature1",
        "entity1",
        "item1",
        "category_x",
        "class_x",
        "group_x",
        "column",
        "variable_x",
        "set_x",
    ),
    "y_category": (
        "y",
        "y_axis",
        "axis_y",
        "var2",
        "variable2",
        "feature2",
        "entity2",
        "item2",
        "category_y",
        "class_y",
        "group_y",
        "row",
        "variable_y",
        "set_y",
    ),
    "category": (
        "category",
        "class",
        "group",
        "type",
        "set",
        "cluster",
        "collection",
        "family",
        "ontology",
        "signature",
        "condition",
        "cohort",
    ),
    "secondary_category": (
        "secondary_category",
        "sub_category",
        "subcategory",
        "group",
        "subgroup",
        "series",
        "component",
        "condition",
        "cohort",
        "class",
        "type",
        "facet",
        "stratum",
        "level",
        "arm",
    ),
    "component": (
        "component",
        "segment",
        "part",
        "ring",
        "subtype",
        "state",
        "type",
        "signature",
        "set",
        "class",
        "gender",
    ),
    "svg_marker": (
        "svg",
        "svg_file",
        "svg_basename",
        "icon",
        "icon_file",
        "glyph",
        "glyph_file",
        "marker",
        "marker_file",
        "image",
        "image_file",
        "basename",
        "file",
    ),
    "node_color": (
        "color",
        "colour",
        "itemstyle_color",
        "itemstyle_colour",
        "fill",
        "hex",
    ),
    "signed_association": (
        "correlation",
        "corr",
        "rho",
        "spearman",
        "pearson",
        "coefficient",
        "coef",
        "association",
        "effect",
        "estimate",
        "score",
        "beta",
        "delta",
        "r",
    ),
    "significance": (
        "p",
        "pvalue",
        "p_value",
        "pval",
        "p_adj",
        "padj",
        "adjusted_p",
        "adjusted_p_value",
        "q",
        "qval",
        "qvalue",
        "q_value",
        "fdr",
        "false_discovery_rate",
        "neglogp",
        "minus_log10_p",
    ),
    "effect_size": (
        "logfc",
        "log2fc",
        "log2foldchange",
        "fold_change",
        "effect_size",
        "effect",
        "estimate",
        "beta",
        "delta",
    ),
    "magnitude": (
        "magnitude",
        "size",
        "strength",
        "weight",
        "importance",
        "score_abs",
        "abs_score",
        "neglogp",
        "minus_log10_p",
    ),
    "average_abundance": (
        "average_abundance",
        "average_expression",
        "average",
        "avg",
        "mean",
        "mean_expression",
        "mean_count",
        "mean_intensity",
        "basemean",
        "base_mean",
        "abundance",
    ),
    "ratio": (
        "ratio",
        "gene_ratio",
        "rich_factor",
        "enrichment_ratio",
        "fraction",
        "proportion",
    ),
    "count": (
        "count",
        "counts",
        "n",
        "number",
        "size",
        "hits",
        "overlap",
        "gene_count",
        "frequency",
        "freq",
        "observed",
        "expected",
    ),
    "uncertainty": (
        "uncertainty",
        "error",
        "errorbar",
        "se",
        "sem",
        "stderr",
        "std_error",
        "standard_error",
        "sd",
        "std",
        "ci",
        "ci95",
        "confidence_interval",
        "lower",
        "upper",
        "ymin",
        "ymax",
        "lcl",
        "ucl",
    ),
    "pair_id": (
        "pair",
        "paired",
        "pair_id",
        "subject",
        "patient",
        "donor",
        "individual",
        "cell_line",
        "line",
        "unit",
        "sample_id",
        "sample",
        "id",
        "replicate",
    ),
    "rank": (
        "rank",
        "order",
        "position",
        "index",
        "ordinal",
    ),
    "direction": (
        "direction",
        "sign",
        "regulation",
        "regulation_direction",
        "change",
        "up_down",
        "updown",
    ),
    "numeric_x": (
        "x",
        "x1",
        "value_x",
        "measure_x",
        "expression_x",
        "score_x",
        "variable_x",
        "umap_1",
        "umap1",
        "tsne_1",
        "tsne1",
        "pc1",
        "pca_1",
        "pca1",
    ),
    "numeric_y": (
        "y",
        "x2",
        "value_y",
        "measure_y",
        "expression_y",
        "score_y",
        "variable_y",
        "umap_2",
        "umap2",
        "tsne_2",
        "tsne2",
        "pc2",
        "pca_2",
        "pca2",
        "value",
    ),
    "time": (
        "time",
        "date",
        "year",
        "timepoint",
        "time_point",
        "stage",
        "phase",
        "period",
        "week",
        "days",
        "months",
        "duration",
        "followup",
        "follow_up",
        "os_time",
        "pfs_time",
    ),
    "interval_start": (
        "start",
        "start_time",
        "start_day",
        "begin",
        "from_time",
        "xstart",
        "xmin",
        "t_start",
    ),
    "interval_end": (
        "end",
        "end_time",
        "end_day",
        "finish",
        "stop",
        "to_time",
        "xend",
        "xmax",
        "t_end",
    ),
    "chromosome": (
        "chromosome",
        "chrom",
        "chr",
        "seqname",
        "seqnames",
        "contig",
        "scaffold",
    ),
    "genomic_start": (
        "start",
        "chrom_start",
        "chromstart",
        "genomic_start",
        "genome_start",
        "bp_start",
        "begin",
    ),
    "genomic_end": (
        "end",
        "chrom_end",
        "chromend",
        "genomic_end",
        "genome_end",
        "bp_end",
        "stop",
    ),
    "site_position": (
        "pos",
        "position",
        "site",
        "locus",
        "coordinate",
        "bp",
    ),
    "zoom_start": (
        "zstart",
        "zoom_start",
        "inner_start",
        "nested_start",
    ),
    "zoom_end": (
        "zend",
        "zoom_end",
        "inner_end",
        "nested_end",
    ),
    "track": (
        "track",
        "assay",
        "layer",
        "signal",
        "mark",
        "histone_mark",
        "lane",
    ),
    "cytoband": (
        "band",
        "stain",
        "gie_stain",
        "g_stain",
        "arm",
        "cytoband",
    ),
    "transcript_id": (
        "transcript",
        "transcript_id",
        "isoform",
        "mrna",
        "model",
    ),
    "feature_type": (
        "feature_type",
        "feature",
        "type",
        "region",
        "element",
        "exon",
        "cds",
        "utr",
    ),
    "residue_position": (
        "aa",
        "residue",
        "residue_position",
        "protein_position",
        "aa_position",
        "aapos",
        "codon",
    ),
    "mutation_site": (
        "site",
        "mutation_site",
        "residue",
        "residue_label",
        "position_label",
    ),
    "mutation_class": (
        "class",
        "mutation_class",
        "variant_class",
        "variant_classification",
        "alteration_class",
        "type",
    ),
    "alteration_type": (
        "alteration",
        "alteration_type",
        "event_type",
        "event",
        "mutation",
        "variant",
        "cnv",
        "copy_number",
        "type",
    ),
    "status": (
        "status",
        "event",
        "censor",
        "censored",
        "death",
        "os_status",
        "pfs_status",
    ),
    "weight": (
        "weight",
        "score",
        "strength",
        "combined_score",
        "confidence",
    ),
}


def normalize(value: str) -> str:
    value = value.strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")


def looks_non_ascii(value: str) -> bool:
    return any(ord(char) > 127 for char in value)


def contains_keyword(query: str, keyword: str) -> bool:
    query_lower = query.lower()
    keyword_lower = keyword.lower()
    if looks_non_ascii(keyword_lower):
        return keyword_lower in query_lower
    return normalize(keyword_lower) in normalize(query_lower)


def infer_delimiter(path: Path, sample: str) -> str:
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".tab"}:
        return "\t"
    if suffix == ".csv":
        return ","
    first_line = sample.splitlines()[0] if sample.splitlines() else ""
    if first_line.count("\t") > first_line.count(","):
        return "\t"
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,;").delimiter
    except csv.Error:
        return ","


def read_rows(path: Path, limit: int = 20000) -> tuple[list[str], list[dict[str, str]], int]:
    sample = path.read_text(encoding="utf-8-sig", errors="replace")
    delimiter = infer_delimiter(path, sample[:8192])
    reader = csv.DictReader(sample.splitlines(), delimiter=delimiter)
    headers = list(reader.fieldnames or [])
    rows: list[dict[str, str]] = []
    total = 0
    for row in reader:
        total += 1
        if len(rows) < limit:
            rows.append({header: (row.get(header) or "").strip() for header in headers})
    return headers, rows, total


def parse_number(value: str) -> float | None:
    if value is None:
        return None
    stripped = value.strip()
    if not stripped or stripped.lower() in {"na", "nan", "null", "none", "inf", "-inf"}:
        return None
    try:
        number = float(stripped.replace(",", ""))
    except ValueError:
        return None
    if not math.isfinite(number):
        return None
    return number


def column_values(rows: list[dict[str, str]], column: str) -> list[str]:
    return [row.get(column, "").strip() for row in rows if row.get(column, "").strip()]


def numeric_summary(rows: list[dict[str, str]], column: str) -> dict[str, Any]:
    values = column_values(rows, column)
    numbers = [number for value in values if (number := parse_number(value)) is not None]
    summary: dict[str, Any] = {
        "non_missing": len(values),
        "numeric": len(numbers),
        "numeric_fraction": round(len(numbers) / len(values), 3) if values else 0.0,
        "is_numeric": bool(values) and len(numbers) / len(values) >= 0.8,
    }
    if numbers:
        summary.update(
            {
                "min": min(numbers),
                "max": max(numbers),
                "positive": sum(1 for number in numbers if number > 0),
                "negative": sum(1 for number in numbers if number < 0),
                "zero": sum(1 for number in numbers if number == 0),
            }
        )
    return summary


def unique_count(rows: list[dict[str, str]], column: str) -> int:
    return len(set(column_values(rows, column)))


def is_bounded_category_column(rows: list[dict[str, str]], column: str) -> bool:
    count = unique_count(rows, column)
    return 1 < count <= max(20, len(rows) // 2)


def is_classification_column_name(column: str) -> bool:
    norm = normalize(column)
    tokens = set(norm.split("_"))
    return bool(tokens & {"class", "condition", "cohort", "group", "status", "type"})


def is_explicit_component_column(column: str) -> bool:
    norm = normalize(column)
    tokens = set(norm.split("_"))
    return bool(tokens & {"component", "segment", "part", "ring", "subtype", "state"})


def score_column_name(column: str, role: str) -> int:
    norm = normalize(column)
    if not norm:
        return 0
    tokens = set(norm.split("_"))
    score = 0
    for hint in ROLE_HINTS.get(role, ()):
        hint_norm = normalize(hint)
        if not hint_norm:
            continue
        if norm == hint_norm:
            score = max(score, 8)
        elif hint_norm in tokens:
            score = max(score, 6)
        elif norm.startswith(f"{hint_norm}_") or norm.endswith(f"_{hint_norm}"):
            score = max(score, 5)
        elif len(hint_norm) >= 3 and hint_norm in norm:
            score = max(score, 3)
    if role == "focal_entity" and norm.startswith(("focal_", "focus_", "center_", "centre_", "central_", "hub_")):
        score = max(score, 9)
    if role == "target_entity" and norm.startswith(("target_", "related_", "item_", "feature_", "to_")):
        score = max(score, 9)
    if role == "source_entity" and norm.startswith(("source_", "from_")):
        score = max(score, 9)
    if role == "parent_entity" and norm.startswith(("parent_", "ancestor_", "root_", "path_", "level_")):
        score = max(score, 9)
    if role == "pair_id" and is_classification_column_name(column) and not norm.endswith(("_id", "id")):
        score = min(score, 2)
    return score


def infer_role_mapping(headers: list[str], rows: list[dict[str, str]]) -> dict[str, str]:
    summaries = {header: numeric_summary(rows, header) for header in headers}
    mapping: dict[str, str] = {}

    for role in ROLE_HINTS:
        best_column = ""
        best_score = 0
        for header in headers:
            score = score_column_name(header, role)
            if role in {
                "numeric_x",
                "numeric_y",
                "signed_association",
                "effect_size",
                "significance",
                "magnitude",
                "average_abundance",
                "ratio",
                "count",
                "weight",
                "uncertainty",
                "rank",
                "genomic_start",
                "genomic_end",
                "site_position",
                "zoom_start",
                "zoom_end",
                "residue_position",
            }:
                if not summaries[header]["is_numeric"]:
                    score -= 4
            if role in {
                "category",
                "secondary_category",
                "component",
                "x_category",
                "y_category",
                "target_entity",
                "focal_entity",
                "source_entity",
                "parent_entity",
                "axis_entity",
                "pair_id",
                "direction",
                "status",
                "time",
                "svg_marker",
                "node_color",
                "interval_start",
                "interval_end",
                "chromosome",
                "track",
                "cytoband",
                "transcript_id",
                "feature_type",
                "mutation_site",
                "mutation_class",
                "alteration_type",
            }:
                if summaries[header]["is_numeric"] and score < 8:
                    score -= 2
            if score > best_score:
                best_score = score
                best_column = header
        if best_score >= 3 and best_column:
            mapping[role] = best_column

    numeric_columns = [header for header in headers if summaries[header]["is_numeric"]]
    used = {
        mapping.get(role)
        for role in (
            "significance",
            "effect_size",
            "signed_association",
            "magnitude",
            "average_abundance",
            "rank",
            "uncertainty",
            "genomic_start",
            "genomic_end",
            "site_position",
            "zoom_start",
            "zoom_end",
            "residue_position",
            "count",
        )
    }
    available_numeric = [header for header in numeric_columns if header not in used]
    if "numeric_x" not in mapping and available_numeric:
        mapping["numeric_x"] = available_numeric[0]
    if "numeric_y" not in mapping and available_numeric:
        mapping["numeric_y"] = available_numeric[1] if len(available_numeric) > 1 else available_numeric[0]
    protected_numeric_columns = {
        mapping.get(role)
        for role in (
            "significance",
            "effect_size",
            "signed_association",
            "magnitude",
            "average_abundance",
            "rank",
            "uncertainty",
            "genomic_start",
            "genomic_end",
            "site_position",
            "zoom_start",
            "zoom_end",
            "residue_position",
            "count",
        )
    }
    for role in ("numeric_x", "numeric_y"):
        if mapping.get(role) in protected_numeric_columns:
            mapping.pop(role, None)

    if "target_entity" not in mapping:
        text_columns = [header for header in headers if not summaries[header]["is_numeric"]]
        if text_columns:
            mapping["target_entity"] = max(text_columns, key=lambda column: unique_count(rows, column))

    if "category" not in mapping:
        text_columns = [header for header in headers if not summaries[header]["is_numeric"]]
        bounded = [column for column in text_columns if is_bounded_category_column(rows, column)]
        if bounded:
            mapping["category"] = min(bounded, key=lambda column: unique_count(rows, column))

    text_columns = [header for header in headers if not summaries[header]["is_numeric"]]
    bounded_text = [column for column in text_columns if is_bounded_category_column(rows, column)]
    if len(numeric_columns) == 1 and len(bounded_text) >= 2:
        cardinalities = {column: unique_count(rows, column) for column in bounded_text}
        high_cardinality = max(bounded_text, key=lambda column: cardinalities[column])
        low_cardinality = min(bounded_text, key=lambda column: cardinalities[column])
        if cardinalities[high_cardinality] > cardinalities[low_cardinality] and cardinalities[low_cardinality] <= 12:
            if mapping.get("pair_id") not in {high_cardinality, low_cardinality}:
                mapping["category"] = high_cardinality
                mapping["secondary_category"] = low_cardinality

    used_entity_columns = {
        mapping.get("category"),
        mapping.get("target_entity"),
        mapping.get("focal_entity"),
        mapping.get("source_entity"),
        mapping.get("parent_entity"),
        mapping.get("pair_id"),
    }
    if mapping.get("secondary_category") == mapping.get("category") or "secondary_category" not in mapping:
        candidates = [column for column in bounded_text if column not in used_entity_columns]
        if candidates:
            mapping["secondary_category"] = min(candidates, key=lambda column: unique_count(rows, column))
        elif mapping.get("secondary_category") == mapping.get("category"):
            mapping.pop("secondary_category", None)

    if mapping.get("component") == mapping.get("category") or "component" not in mapping:
        candidates = [
            column
            for column in bounded_text
            if column not in {mapping.get("category"), mapping.get("secondary_category"), mapping.get("target_entity")}
        ]
        if candidates:
            mapping["component"] = min(candidates, key=lambda column: unique_count(rows, column))
        elif mapping.get("component") == mapping.get("category"):
            mapping.pop("component", None)

    if (
        "pair_id" in mapping
        and mapping["pair_id"] in {mapping.get("category"), mapping.get("secondary_category"), mapping.get("component")}
        and is_classification_column_name(mapping["pair_id"])
    ):
        mapping.pop("pair_id", None)

    if "signed_association" in mapping and "target_entity" in mapping and "category" in mapping:
        if "x_category" not in mapping:
            mapping["x_category"] = mapping["target_entity"]
        if "y_category" not in mapping:
            mapping["y_category"] = mapping["category"]

    return mapping


def role_present(role: str, mapping: dict[str, str], shapes: set[str]) -> bool:
    if role in mapping:
        return True
    if role == "numeric_matrix":
        return (
            "numeric_matrix" in shapes
            or "wide_numeric_table" in shapes
            or "radar_wide_table" in shapes
            or "dendrogram_matrix" in shapes
        )
    if role == "flow_value":
        return "weight" in mapping or "count" in mapping or "numeric_y" in mapping
    if role == "ternary_components":
        return "ternary_components" in shapes
    if role in {
        "cytoband_table",
        "genomic_density_windows",
        "genomic_heatmap_table",
        "genomic_interval_table",
        "genomic_locus_table",
        "genomic_track_table",
        "group_split_matrix",
        "mutation_energy_matrix",
        "oncoprint_events",
        "protein_lollipop_table",
        "site_score_table",
        "synteny_blocks",
        "transcript_feature_table",
        "two_set_sample_matrix",
    }:
        return role in shapes
    return False


def looks_count_like_column(column: str) -> bool:
    norm = normalize(column)
    tokens = set(norm.split("_"))
    count_tokens = {"count", "counts", "total", "unique", "core", "size", "overlap", "hits", "n"}
    return bool(tokens & count_tokens)


def has_explicit_pair_axis_columns(mapping: dict[str, str]) -> bool:
    x_column = normalize(mapping.get("x_category", ""))
    y_column = normalize(mapping.get("y_category", ""))
    x_names = {"x", "x_axis", "var1", "variable1", "feature1", "entity1", "item1"}
    y_names = {"y", "y_axis", "var2", "variable2", "feature2", "entity2", "item2"}
    return x_column in x_names and y_column in y_names


def has_embedding_axis_columns(mapping: dict[str, str]) -> bool:
    x_column = normalize(mapping.get("numeric_x", ""))
    y_column = normalize(mapping.get("numeric_y", ""))
    x_names = {"umap_1", "umap1", "tsne_1", "tsne1", "pc1", "pca_1", "pca1"}
    y_names = {"umap_2", "umap2", "tsne_2", "tsne2", "pc2", "pca_2", "pca2"}
    return x_column in x_names and y_column in y_names


def has_position_axis_columns(mapping: dict[str, str]) -> bool:
    x_column = normalize(mapping.get("numeric_x", ""))
    y_column = normalize(mapping.get("numeric_y", ""))
    x_names = {"x", "x1", "umap_1", "umap1", "tsne_1", "tsne1", "pc1", "pca_1", "pca1"}
    y_names = {"y", "x2", "umap_2", "umap2", "tsne_2", "tsne2", "pc2", "pca_2", "pca2"}
    return x_column in x_names and y_column in y_names


def has_synteny_block_columns(headers: list[str]) -> bool:
    norms = {normalize(header) for header in headers}
    numbered = {"chr1", "start1", "end1", "chr2", "start2", "end2"} <= norms
    prefixed = {
        "source_chr",
        "source_start",
        "source_end",
        "target_chr",
        "target_start",
        "target_end",
    } <= norms
    query_target = {
        "query_chr",
        "query_start",
        "query_end",
        "target_chr",
        "target_start",
        "target_end",
    } <= norms
    return numbered or prefixed or query_target


def amino_acid_numeric_columns(headers: list[str], numeric_columns: list[str]) -> list[str]:
    amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    return [
        column
        for column in numeric_columns
        if normalize(column).upper() in amino_acids and len(normalize(column)) == 1
    ]


def discover_sidecars(input_path: Path, sidecar_dir: Path | None = None) -> dict[str, dict[str, Any]]:
    directory = sidecar_dir or input_path.parent
    if not directory.is_dir():
        return {}
    input_resolved = input_path.resolve(strict=False)
    sidecars: dict[str, dict[str, Any]] = {}
    for path in sorted(directory.iterdir(), key=lambda item: item.name.lower()):
        if not path.is_file() or path.suffix.lower() not in SIDECAR_SUFFIXES:
            continue
        if path.resolve(strict=False) == input_resolved:
            continue
        if normalize(path.stem) not in SIDECAR_STEMS:
            continue
        headers, rows, total_rows = read_rows(path, limit=5000)
        if not headers:
            continue
        mapping = infer_role_mapping(headers, rows)
        shapes, _details = detect_shapes(headers, rows, mapping)
        sidecars[path.name] = {
            "path": str(path.resolve(strict=False)),
            "row_count": total_rows,
            "sampled_rows": len(rows),
            "column_count": len(headers),
            "columns": headers,
            "role_mapping": mapping,
            "shapes": sorted(shapes),
        }
    return sidecars


def collect_sidecar_role_mapping(sidecars: dict[str, dict[str, Any]]) -> dict[str, str]:
    role_mapping: dict[str, str] = {}
    for filename, profile in sidecars.items():
        mapping = profile.get("role_mapping", {})
        if not isinstance(mapping, dict):
            continue
        for role, column in mapping.items():
            if role not in SIDECAR_ROLE_EXPORTS or role in role_mapping:
                continue
            role_mapping[role] = f"{filename}:{column}"
    return role_mapping


def unique_values(rows: list[dict[str, str]], column: str) -> set[str]:
    return {value for value in column_values(rows, column)}


def first_identifier_column(
    headers: list[str],
    rows: list[dict[str, str]],
    mapping: dict[str, str],
    preferred_roles: tuple[str, ...],
) -> str:
    for role in preferred_roles:
        column = mapping.get(role)
        if column and unique_values(rows, column):
            return column
    summaries = {header: numeric_summary(rows, header) for header in headers}
    for header in headers:
        if not summaries[header]["is_numeric"] and unique_values(rows, header):
            return header
    return headers[0] if headers else ""


def alignment_status(
    relationship: str,
    overlap_count: int,
    sidecar_count: int,
    reference_count: int,
    missing_count: int,
    reference_coverage: float,
) -> str:
    if sidecar_count == 0 or reference_count == 0:
        return "skipped"
    if overlap_count == 0:
        return "error"
    if missing_count > 0:
        return "warning"
    if relationship == "matrix_columns" and reference_coverage < 0.8:
        return "warning"
    return "ok"


def alignment_check(
    *,
    sidecar: str,
    relationship: str,
    sidecar_column: str,
    sidecar_values: set[str],
    reference: str,
    reference_values: set[str],
) -> dict[str, Any]:
    overlap = sidecar_values & reference_values
    missing = sidecar_values - reference_values
    extra = reference_values - sidecar_values
    sidecar_count = len(sidecar_values)
    reference_count = len(reference_values)
    overlap_count = len(overlap)
    sidecar_coverage = round(overlap_count / sidecar_count, 3) if sidecar_count else 0.0
    reference_coverage = round(overlap_count / reference_count, 3) if reference_count else 0.0
    status = alignment_status(
        relationship,
        overlap_count,
        sidecar_count,
        reference_count,
        len(missing),
        reference_coverage,
    )
    message = (
        f"{sidecar} {sidecar_column} overlaps {overlap_count}/{sidecar_count} "
        f"sidecar IDs and {overlap_count}/{reference_count} reference IDs"
    )
    if status == "error":
        message = f"{sidecar} {sidecar_column} has no overlap with {reference}"
    elif status == "warning":
        message = f"{sidecar} {sidecar_column} is only partially aligned with {reference}"
    return {
        "sidecar": sidecar,
        "relationship": relationship,
        "sidecar_column": sidecar_column,
        "reference": reference,
        "sidecar_count": sidecar_count,
        "reference_count": reference_count,
        "overlap_count": overlap_count,
        "sidecar_coverage": sidecar_coverage,
        "reference_coverage": reference_coverage,
        "missing_count": len(missing),
        "extra_count": len(extra),
        "missing_examples": sorted(missing)[:5],
        "extra_examples": sorted(extra)[:5],
        "status": status,
        "message": message,
    }


def sidecar_profile_rows(profile: dict[str, Any]) -> tuple[list[str], list[dict[str, str]]]:
    path = Path(str(profile.get("path", "")))
    if not path:
        return [], []
    headers, rows, _total = read_rows(path, limit=5000)
    return headers, rows


def sidecar_by_stem(sidecars: dict[str, dict[str, Any]]) -> dict[str, tuple[str, dict[str, Any]]]:
    return {normalize(Path(filename).stem): (filename, profile) for filename, profile in sidecars.items()}


def evaluate_sidecar_alignment(
    input_path: Path,
    headers: list[str],
    rows: list[dict[str, str]],
    mapping: dict[str, str],
    details: dict[str, Any],
    sidecars: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    checks: list[dict[str, Any]] = []
    if not sidecars:
        return {"status": "not_checked", "checks": []}

    main_row_column = first_identifier_column(
        headers,
        rows,
        mapping,
        ("target_entity", "axis_entity", "focal_entity", "pair_id", "category", "source_entity", "parent_entity"),
    )
    main_row_ids = unique_values(rows, main_row_column) if main_row_column else set()
    matrix_column_ids = {column for column in details.get("numeric_columns", []) if column in headers}
    by_stem = sidecar_by_stem(sidecars)

    for stem in ("rowinfo", "annot"):
        if stem not in by_stem or not main_row_ids:
            continue
        filename, profile = by_stem[stem]
        sidecar_headers, sidecar_rows = sidecar_profile_rows(profile)
        sidecar_mapping = profile.get("role_mapping", {}) if isinstance(profile.get("role_mapping"), dict) else {}
        sidecar_column = first_identifier_column(
            sidecar_headers,
            sidecar_rows,
            sidecar_mapping,
            ("target_entity", "axis_entity", "focal_entity", "pair_id", "category"),
        )
        if sidecar_column:
            checks.append(
                alignment_check(
                    sidecar=filename,
                    relationship="matrix_rows",
                    sidecar_column=sidecar_column,
                    sidecar_values=unique_values(sidecar_rows, sidecar_column),
                    reference=f"{input_path.name}:{main_row_column}",
                    reference_values=main_row_ids,
                )
            )

    if "colinfo" in by_stem and matrix_column_ids:
        filename, profile = by_stem["colinfo"]
        sidecar_headers, sidecar_rows = sidecar_profile_rows(profile)
        sidecar_mapping = profile.get("role_mapping", {}) if isinstance(profile.get("role_mapping"), dict) else {}
        sidecar_column = first_identifier_column(
            sidecar_headers,
            sidecar_rows,
            sidecar_mapping,
            ("pair_id", "target_entity", "axis_entity", "category"),
        )
        if sidecar_column:
            checks.append(
                alignment_check(
                    sidecar=filename,
                    relationship="matrix_columns",
                    sidecar_column=sidecar_column,
                    sidecar_values=unique_values(sidecar_rows, sidecar_column),
                    reference=f"{input_path.name}:numeric columns",
                    reference_values=matrix_column_ids,
                )
            )

    node_ids: set[str] = set()
    node_reference = ""
    input_stem = normalize(input_path.stem)
    if input_stem == "nodes" and main_row_ids:
        node_ids = main_row_ids
        node_reference = f"{input_path.name}:{main_row_column}"
    elif "nodes" in by_stem:
        filename, profile = by_stem["nodes"]
        sidecar_headers, sidecar_rows = sidecar_profile_rows(profile)
        sidecar_mapping = profile.get("role_mapping", {}) if isinstance(profile.get("role_mapping"), dict) else {}
        sidecar_column = first_identifier_column(
            sidecar_headers,
            sidecar_rows,
            sidecar_mapping,
            ("target_entity", "axis_entity", "focal_entity", "pair_id", "category", "source_entity"),
        )
        if sidecar_column:
            node_ids = unique_values(sidecar_rows, sidecar_column)
            node_reference = f"{filename}:{sidecar_column}"

    for stem in ("links", "edges"):
        if stem not in by_stem or not node_ids:
            continue
        filename, profile = by_stem[stem]
        sidecar_headers, sidecar_rows = sidecar_profile_rows(profile)
        sidecar_mapping = profile.get("role_mapping", {}) if isinstance(profile.get("role_mapping"), dict) else {}
        source_column = sidecar_mapping.get("source_entity")
        target_column = sidecar_mapping.get("target_entity")
        if source_column and target_column:
            endpoint_values = unique_values(sidecar_rows, source_column) | unique_values(sidecar_rows, target_column)
            checks.append(
                alignment_check(
                    sidecar=filename,
                    relationship="node_link_endpoints",
                    sidecar_column=f"{source_column},{target_column}",
                    sidecar_values=endpoint_values,
                    reference=node_reference,
                    reference_values=node_ids,
                )
            )

    statuses = {check["status"] for check in checks}
    if "error" in statuses:
        status = "error"
    elif "warning" in statuses:
        status = "warning"
    elif checks:
        status = "ok"
    else:
        status = "not_checked"
    return {"status": status, "checks": checks}


def sidecar_alignment_messages(profile: dict[str, Any], limit: int = 3) -> list[str]:
    alignment = profile.get("sidecar_alignment", {})
    checks = alignment.get("checks", []) if isinstance(alignment, dict) else []
    messages = [check.get("message", "") for check in checks if check.get("status") in {"warning", "error"}]
    return [message for message in messages if message][:limit]


def sidecar_shapes(sidecars: dict[str, dict[str, Any]], input_stem: str = "") -> set[str]:
    shapes: set[str] = set()
    stems = {normalize(Path(filename).stem) for filename in sidecars}
    stems_with_input = set(stems)
    if input_stem:
        stems_with_input.add(input_stem)
    if {"nodes", "links"} <= stems_with_input or {"nodes", "edges"} <= stems_with_input:
        shapes.add("node_link_sidecars")
        shapes.add("network_edges")
    if {"rowinfo", "colinfo", "enrichment"} <= stems:
        shapes.add("annotation_sidecars")
        shapes.add("enrichment_zoom_sidecars")
    elif stems & {"annot", "colinfo", "groups", "rowinfo"}:
        shapes.add("annotation_sidecars")
    if "enrichment" in stems:
        shapes.add("enrichment_sidecar")
    if stems & {"cytoband", "features", "genes", "loops", "points", "windows"}:
        shapes.add("genomic_annotation_sidecars")
    if "domains" in stems:
        shapes.add("protein_domain_sidecar")
    if "karyotype" in stems:
        shapes.add("synteny_sidecars")

    for filename, profile in sidecars.items():
        mapping = profile.get("role_mapping", {})
        profile_shapes = set(profile.get("shapes", []))
        stem = normalize(Path(filename).stem)
        if not isinstance(mapping, dict):
            continue
        if stem in {"edges", "links"} or {"source_entity", "target_entity"} <= set(mapping):
            shapes.add("network_edges")
            if "signed_association" in mapping:
                shapes.add("signed_network_edges")
            if "weight" in mapping or "count" in mapping or "numeric_y" in mapping:
                shapes.add("flow_table")
        if "cytoband_table" in profile_shapes:
            shapes.add("cytoband_table")
        if "synteny_blocks" in profile_shapes:
            shapes.add("synteny_blocks")
            shapes.add("synteny_sidecars")
    return shapes


def detect_shapes(headers: list[str], rows: list[dict[str, str]], mapping: dict[str, str]) -> tuple[set[str], dict[str, Any]]:
    summaries = {header: numeric_summary(rows, header) for header in headers}
    numeric_columns = [header for header, summary in summaries.items() if summary["is_numeric"]]
    text_columns = [header for header in headers if header not in numeric_columns]
    position_columns = set()
    if has_position_axis_columns(mapping):
        position_columns = {mapping.get("numeric_x"), mapping.get("numeric_y")}
    position_columns.update(
        column
        for role in (
            "genomic_start",
            "genomic_end",
            "site_position",
            "zoom_start",
            "zoom_end",
            "residue_position",
        )
        if (column := mapping.get(role))
    )
    measurement_numeric_columns = [column for column in numeric_columns if column not in position_columns]
    summary_numeric_columns = {
        mapping.get(role)
        for role in (
            "significance",
            "effect_size",
            "signed_association",
            "magnitude",
            "average_abundance",
            "rank",
            "uncertainty",
            "count",
        )
    }
    matrix_layer_columns = [
        column for column in measurement_numeric_columns if column not in summary_numeric_columns
    ]
    has_numeric_measure = bool(measurement_numeric_columns) or (
        "numeric_y" in mapping and mapping["numeric_y"] not in position_columns
    )
    shapes: set[str] = set()
    if rows:
        shapes.add("long_table")
    if len(numeric_columns) >= 2 and len(rows) >= 3:
        shapes.add("two_numeric_observation")
    if len(numeric_columns) >= 3 and len(text_columns) >= 1:
        shapes.add("wide_numeric_table")
    if len(numeric_columns) >= max(3, len(headers) - 1) and len(rows) >= 2:
        shapes.add("numeric_matrix")
    if "numeric_matrix" in shapes or "wide_numeric_table" in shapes:
        shapes.add("dendrogram_matrix")
    if len(numeric_columns) >= 3 and len(rows) >= 3:
        shapes.add("wide_many_variables")
    nonnegative_numeric = [
        column
        for column in numeric_columns
        if summaries[column].get("numeric", 0) > 0 and summaries[column].get("negative", 0) == 0
    ]
    if len(nonnegative_numeric) >= 3 and len(rows) >= 1:
        shapes.add("ternary_components")
    if has_embedding_axis_columns(mapping) and "category" in mapping:
        shapes.add("embedding_with_tracks")
    if len(numeric_columns) >= 4 and len(rows) >= 3 and len(text_columns) <= 1:
        shapes.add("two_set_sample_matrix")

    has_target = "target_entity" in mapping
    has_signed = "signed_association" in mapping and summaries[mapping["signed_association"]]["is_numeric"]
    if has_target and has_signed:
        shapes.add("long_association")
    if "x_category" in mapping and "y_category" in mapping and has_signed:
        shapes.add("association_grid")
    elif has_signed:
        categorical = [column for column in text_columns if unique_count(rows, column) > 1]
        if len(categorical) >= 2:
            shapes.add("association_grid")
    if ((has_numeric_measure or "magnitude" in mapping or "count" in mapping) and len(text_columns) >= 2):
        shapes.add("category_grid")

    if "focal_entity" in mapping and has_target and has_signed and unique_count(rows, mapping["focal_entity"]) <= 1:
        shapes.add("one_to_many_association")
    elif has_target and has_signed:
        if "association_grid" not in shapes and len(rows) > 1:
            shapes.add("one_to_many_association")

    if has_target and "effect_size" in mapping and "significance" in mapping:
        shapes.add("feature_level_testing")
    if has_target and "effect_size" in mapping and "average_abundance" in mapping:
        shapes.add("mean_difference_table")
    if has_target and ("rank" in mapping or "effect_size" in mapping or "signed_association" in mapping):
        shapes.add("ranked_value")
    if "category" in mapping and has_numeric_measure:
        shapes.add("category_distribution")
        shapes.add("category_magnitude")
        has_second_category = (
            "secondary_category" in mapping
            and mapping["secondary_category"] != mapping.get("category")
        )
        has_explicit_component = (
            "component" in mapping
            and mapping["component"] != mapping.get("category")
            and is_explicit_component_column(mapping["component"])
        )
        if len(text_columns) == 1 and len(numeric_columns) == 1:
            shapes.add("part_to_whole")
        if has_second_category or has_explicit_component:
            shapes.add("grouped_category_value")
            shapes.add("grouped_distribution")
            if "count" in mapping or "magnitude" in mapping or "numeric_y" in mapping:
                shapes.add("single_axis_series")
        if has_explicit_component:
            shapes.add("stacked_composition")
            shapes.add("nested_composition")
            shapes.add("part_to_whole")
        if "uncertainty" in mapping:
            shapes.add("category_uncertainty")
        if "pair_id" in mapping:
            pair_values = column_values(rows, mapping["pair_id"])
            if len(set(pair_values)) < len(pair_values):
                shapes.add("paired_distribution")
        series_columns = [
            column
            for column in numeric_columns
            if column
            not in {
                mapping.get("uncertainty"),
                mapping.get("significance"),
                mapping.get("count"),
                mapping.get("rank"),
            }
        ]
        has_row_level_label = any(
            column and unique_count(rows, column) == len(rows)
            for column in (mapping.get("target_entity"), mapping.get("category"))
        )
        has_xy_position_names = (
            normalize(mapping.get("numeric_x", "")) in {"x", "umap_1", "umap1", "tsne_1", "tsne1"}
            and normalize(mapping.get("numeric_y", "")) in {"y", "umap_2", "umap2", "tsne_2", "tsne2"}
        )
        if len(series_columns) == 2 and has_row_level_label and not has_xy_position_names:
            shapes.add("category_two_series")
        count_like_columns = [column for column in numeric_columns if looks_count_like_column(column)]
        if count_like_columns:
            shapes.add("set_membership_summary")
            if len(rows) >= 5 or len(count_like_columns) >= 2:
                shapes.add("multi_set_overlap")
    if has_target and "category" in mapping and len(text_columns) >= 2 and not numeric_columns:
        shapes.add("set_membership_table")
    if has_target and "ratio" in mapping and "significance" in mapping:
        shapes.add("enrichment_terms")
    if has_target and "significance" in mapping and "effect_size" not in mapping:
        shapes.add("enrichment_ranked_terms")
        if "category" in mapping:
            shapes.add("enrichment_grouped_terms")
    if has_target and "category" in mapping and "significance" in mapping and "direction" in mapping:
        shapes.add("enrichment_term_grid")
    if "time" in mapping and "status" in mapping and "category" in mapping:
        shapes.add("survival_table")
    if "time" in mapping and has_numeric_measure:
        shapes.add("time_series")
        if "category" in mapping or "secondary_category" in mapping:
            shapes.add("grouped_time_series")
        series_numeric_columns = [
            column for column in measurement_numeric_columns if column != mapping.get("time")
        ]
        if len(series_numeric_columns) >= 2:
            shapes.add("wide_two_series")
        if "rank" in mapping or has_target or len(series_numeric_columns) >= 2:
            shapes.add("rank_time_series")
    if "source_entity" in mapping and has_target:
        source_col = mapping["source_entity"]
        target_col = mapping["target_entity"]
        if source_col != target_col:
            shapes.add("network_edges")
            if has_signed:
                shapes.add("signed_network_edges")
            if "weight" in mapping or "count" in mapping or "numeric_y" in mapping:
                shapes.add("flow_table")
            if "count" in mapping and "significance" in mapping:
                shapes.add("classification_enrichment_tree")
    if len(text_columns) >= 3 and "count" in mapping:
        shapes.add("parallel_sets_table")
    if "parent_entity" in mapping and has_target and mapping["parent_entity"] != mapping["target_entity"]:
        shapes.add("hierarchy_edges")
        if "numeric_y" in mapping or "count" in mapping:
            shapes.add("hierarchy_area")
            if "node_color" in mapping:
                shapes.add("colored_hierarchy_area")
        empty_parent_rows = sum(1 for row in rows if not row.get(mapping["parent_entity"], "").strip())
        if empty_parent_rows >= 2:
            shapes.add("multi_root_hierarchy")
    if ("axis_entity" in mapping or has_target or "category" in mapping) and len(numeric_columns) >= 2:
        shapes.add("radar_wide_table")
    if has_signed and has_explicit_pair_axis_columns(mapping):
        shapes.add("triangular_correlation")
        if "category" in mapping:
            shapes.add("grouped_correlation")
    if "pair_id" in mapping and "interval_start" in mapping and "interval_end" in mapping:
        shapes.add("interval_timeline")
    if "svg_marker" in mapping and ("numeric_y" in mapping or {"numeric_x", "numeric_y"} <= set(mapping)):
        shapes.add("svg_marked_observation")
    if has_target and "pair_id" in mapping and "alteration_type" in mapping and not numeric_columns:
        shapes.add("oncoprint_events")
    if "mutation_site" in mapping and len(amino_acid_numeric_columns(headers, numeric_columns)) >= 8:
        shapes.add("mutation_energy_matrix")
    if "site_position" in mapping and has_numeric_measure and "chromosome" not in mapping:
        shapes.add("site_score_table")
    if "residue_position" in mapping and ("mutation_class" in mapping or "category" in mapping):
        if "count" in mapping or has_numeric_measure:
            shapes.add("protein_lollipop_table")
    if "transcript_id" in mapping and "feature_type" in mapping and "genomic_start" in mapping and "genomic_end" in mapping:
        shapes.add("transcript_feature_table")

    has_synteny_blocks = has_synteny_block_columns(headers)
    has_genomic_interval = "chromosome" in mapping and "genomic_start" in mapping and "genomic_end" in mapping
    if has_genomic_interval and not has_synteny_blocks:
        shapes.add("genomic_interval_table")
        if "cytoband" in mapping:
            shapes.add("cytoband_table")
        if has_target and "category" in mapping:
            shapes.add("genomic_locus_table")
        if len(measurement_numeric_columns) == 1 and "track" not in mapping and "cytoband" not in mapping:
            shapes.add("genomic_density_windows")
        if "track" in mapping and has_numeric_measure:
            shapes.add("genomic_track_table")
            shapes.add("genomic_coverage_tracks")
        if len(measurement_numeric_columns) >= 2:
            shapes.add("genomic_heatmap_table")
        if "zoom_start" in mapping and "zoom_end" in mapping:
            shapes.add("nested_genomic_windows")
    if has_target and "category" in mapping and "significance" in mapping and len(matrix_layer_columns) >= 3:
        shapes.add("group_split_matrix")
    if has_synteny_blocks:
        shapes.add("synteny_blocks")

    profile_details = {
        "numeric_columns": numeric_columns,
        "text_columns": text_columns,
        "numeric_summaries": summaries,
        "unique_counts": {header: unique_count(rows, header) for header in headers},
    }
    return shapes, profile_details


def build_profile(path: Path, sidecar_dir: Path | None = None) -> dict[str, Any]:
    headers, rows, total_rows = read_rows(path)
    mapping = infer_role_mapping(headers, rows)
    shapes, details = detect_shapes(headers, rows, mapping)
    sidecars = discover_sidecars(path, sidecar_dir)
    detected_sidecar_shapes = sidecar_shapes(sidecars, normalize(path.stem))
    shapes.update(detected_sidecar_shapes)
    sidecar_role_mapping = collect_sidecar_role_mapping(sidecars)
    sidecar_alignment = evaluate_sidecar_alignment(path, headers, rows, mapping, details, sidecars)
    return {
        "path": str(path),
        "row_count": total_rows,
        "sampled_rows": len(rows),
        "column_count": len(headers),
        "columns": headers,
        "role_mapping": mapping,
        "sidecar_dir": str((sidecar_dir or path.parent).resolve(strict=False)),
        "sidecars": sidecars,
        "sidecar_shapes": sorted(detected_sidecar_shapes),
        "sidecar_role_mapping": sidecar_role_mapping,
        "sidecar_alignment": sidecar_alignment,
        "shapes": sorted(shapes),
        **details,
    }


def load_contracts(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def score_template(template: dict[str, Any], profile: dict[str, Any], query: str, mode: str) -> dict[str, Any]:
    shapes = set(profile["shapes"])
    mapping = {**profile.get("sidecar_role_mapping", {}), **profile["role_mapping"]}
    preferred_shapes = set(template.get("preferred_shapes", []))
    matched_shapes = sorted(preferred_shapes & shapes)
    required_roles = template.get("roles", {}).get("required", [])
    optional_roles = template.get("roles", {}).get("optional", [])
    score = 0
    rationale: list[str] = []
    risks: list[str] = []
    confidence_cap = ""

    template_id = template["id"]
    exact_id = template_id.lower() in query.lower()
    if exact_id:
        score += 100
        rationale.append(f"query names exact template id {template_id}")

    if matched_shapes:
        shape_score = 20 * len(matched_shapes)
        score += shape_score
        rationale.append(f"data shape matches {', '.join(matched_shapes)}")
        decisive_shapes = sorted(
            set(matched_shapes)
            & {
                "classification_enrichment_tree",
                "colored_hierarchy_area",
                "embedding_with_tracks",
                "enrichment_zoom_sidecars",
                "enrichment_terms",
                "genomic_density_windows",
                "genomic_heatmap_table",
                "genomic_locus_table",
                "genomic_annotation_sidecars",
                "genomic_track_table",
                "group_split_matrix",
                "interval_timeline",
                "mean_difference_table",
                "multi_root_hierarchy",
                "mutation_energy_matrix",
                "node_link_sidecars",
                "one_to_many_association",
                "oncoprint_events",
                "protein_domain_sidecar",
                "protein_lollipop_table",
                "rank_time_series",
                "site_score_table",
                "synteny_blocks",
                "synteny_sidecars",
                "ternary_components",
                "transcript_feature_table",
            }
        )
        if decisive_shapes:
            decisive_shape_weights = {
                "cytoband_table": 25,
                "genomic_density_windows": 25,
                "genomic_heatmap_table": 25,
                "genomic_annotation_sidecars": 15,
                "genomic_locus_table": 25,
                "genomic_track_table": 25,
                "group_split_matrix": 25,
                "mutation_energy_matrix": 25,
                "node_link_sidecars": 20,
                "oncoprint_events": 25,
                "protein_domain_sidecar": 15,
                "protein_lollipop_table": 25,
                "site_score_table": 25,
                "synteny_blocks": 35,
                "synteny_sidecars": 20,
                "transcript_feature_table": 25,
            }
            score += sum(decisive_shape_weights.get(shape, 10) for shape in decisive_shapes)
            rationale.append(f"decisive shape match: {', '.join(decisive_shapes)}")

        sidecar_intent_hits: list[str] = []
        for shape in sorted(set(matched_shapes) & SIDECAR_DEPENDENT_SHAPES):
            for keyword in SIDECAR_SHAPE_INTENT_KEYWORDS.get(shape, ()):
                if contains_keyword(query, keyword):
                    sidecar_intent_hits.append(f"{shape}:{keyword}")
        if sidecar_intent_hits:
            score += min(24, 8 * len(sidecar_intent_hits))
            rationale.append("sidecar intent match: " + ", ".join(sidecar_intent_hits[:5]))

        sidecar_status = profile.get("sidecar_alignment", {}).get("status")
        if set(matched_shapes) & SIDECAR_DEPENDENT_SHAPES and sidecar_status in {"warning", "error"}:
            if sidecar_status == "error":
                score -= 28
                confidence_cap = "low"
            else:
                score -= 12
                confidence_cap = "medium"
            messages = sidecar_alignment_messages(profile)
            if messages:
                risks.append("sidecar alignment issue: " + "; ".join(messages))
            else:
                risks.append(f"sidecar alignment status is {sidecar_status}")
    else:
        score -= 18
        risks.append("no preferred input shape matched the table profile")

    for role in required_roles:
        if role_present(role, mapping, shapes):
            score += 10
            column = mapping.get(role, role)
            rationale.append(f"required role {role} is present ({column})")
        else:
            score -= 16
            risks.append(f"missing required role {role}")

    for role in optional_roles:
        if role_present(role, mapping, shapes):
            score += 4
            column = mapping.get(role, role)
            rationale.append(f"optional role {role} is available ({column})")

    keyword_hits = [keyword for keyword in template.get("intent_keywords", []) if contains_keyword(query, keyword)]
    if keyword_hits:
        score += 8 * len(keyword_hits)
        rationale.append("intent keyword match: " + ", ".join(keyword_hits[:5]))

    if mode == "preview" and template.get("preview"):
        score += 5
        rationale.append("preview asset is available for fast comparison")
    elif mode == "publication":
        score -= 3
        risks.append("publication mode still requires catalog and final-size visual QA")

    if template_id == "graph-focus":
        query_mentions_network = any(contains_keyword(query, keyword) for keyword in ("network", "graph", "edge", "网络", "互作"))
        if "network_edges" not in shapes and not query_mentions_network:
            score -= 12
            risks.append("graph-focus needs an explicit network/edge interpretation")

    if template_id == "scatter-one2many" and "association_grid" in shapes and "one_to_many_association" not in shapes:
        score -= 8
        risks.append("association grid may read better as heatmap-corr-bubble than a focal layout")

    if template_id == "heatmap-corr-bubble" and "one_to_many_association" in shapes and not exact_id:
        score -= 6
        risks.append("one focal entity may read more directly as scatter-one2many")

    if score >= 80 or exact_id:
        confidence = "high"
    elif score >= 45:
        confidence = "medium"
    else:
        confidence = "low"
    if confidence_cap:
        confidence_rank = {"low": 0, "medium": 1, "high": 2}
        if confidence_rank[confidence] > confidence_rank[confidence_cap]:
            confidence = confidence_cap

    return {
        "id": template_id,
        "family": template.get("family"),
        "title": template.get("title"),
        "source": template.get("source"),
        "preview": template.get("preview"),
        "score": score,
        "confidence": confidence,
        "matched_shapes": matched_shapes,
        "required_roles": required_roles,
        "optional_roles": optional_roles,
        "role_mapping": {role: mapping[role] for role in sorted(mapping) if role in set(required_roles + optional_roles)},
        "rationale": rationale[:12],
        "risks": risks[:8],
        "use_when": template.get("use_when"),
        "avoid_when": template.get("avoid_when"),
    }


def recommend(profile: dict[str, Any], contracts: dict[str, Any], query: str, mode: str, top: int) -> list[dict[str, Any]]:
    scored = [score_template(template, profile, query, mode) for template in contracts.get("templates", [])]
    scored.sort(key=lambda item: (item["score"], item["confidence"] == "high"), reverse=True)
    return scored[:top]


def render_text(payload: dict[str, Any]) -> str:
    lines = [
        f"Input: {payload['input_profile']['path']}",
        f"Rows x columns: {payload['input_profile']['row_count']} x {payload['input_profile']['column_count']}",
        "Detected shapes: " + ", ".join(payload["input_profile"]["shapes"]),
        "Sidecars: " + ", ".join(payload["input_profile"].get("sidecars", {}).keys()),
        "Sidecar alignment: " + payload["input_profile"].get("sidecar_alignment", {}).get("status", "not_checked"),
        "Top recommendations:",
    ]
    for index, item in enumerate(payload["recommendations"], start=1):
        lines.append(f"{index}. {item['id']} ({item['family']}) score={item['score']} confidence={item['confidence']}")
        if item["matched_shapes"]:
            lines.append("   shapes: " + ", ".join(item["matched_shapes"]))
        if item["role_mapping"]:
            role_bits = [f"{role}={column}" for role, column in item["role_mapping"].items()]
            lines.append("   roles: " + ", ".join(role_bits))
        if item["rationale"]:
            lines.append("   why: " + "; ".join(item["rationale"][:3]))
        if item["risks"]:
            lines.append("   check: " + "; ".join(item["risks"][:2]))
    return "\n".join(lines)


def build_payload(args: argparse.Namespace) -> dict[str, Any]:
    input_path = Path(args.input).expanduser().resolve()
    sidecar_dir = Path(args.sidecar_dir).expanduser().resolve() if args.sidecar_dir else None
    contracts_path = Path(args.contracts).expanduser().resolve() if args.contracts else DEFAULT_CONTRACTS
    profile = build_profile(input_path, sidecar_dir)
    contracts = load_contracts(contracts_path)
    recommendations = recommend(profile, contracts, args.query, args.mode, args.top)
    return {
        "mode": args.mode,
        "query": args.query,
        "contracts": str(contracts_path),
        "input_profile": profile,
        "recommendations": recommendations,
    }


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Recommend an omics visualization template from a result table.")
    parser.add_argument("--input", required=True, help="CSV/TSV result table to profile.")
    parser.add_argument("--query", required=True, help="User request or intended figure purpose.")
    parser.add_argument(
        "--mode",
        choices=("preview", "publication", "template-dev"),
        default="preview",
        help="Routing strictness. Preview is fast; publication requires manual catalog confirmation.",
    )
    parser.add_argument("--top", type=int, default=4, help="Number of recommendations to return.")
    parser.add_argument("--contracts", help="Override template_contracts.json path.")
    parser.add_argument(
        "--sidecar-dir",
        help="Optional directory containing companion files such as nodes.tsv, links.tsv, rowInfo.tsv, or cytoband.tsv.",
    )
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    payload = build_payload(args)
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    else:
        print(render_text(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
