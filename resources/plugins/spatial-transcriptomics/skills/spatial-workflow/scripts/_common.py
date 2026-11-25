#!/usr/bin/env python3
"""Shared helpers for the spatial-workflow script toolkit.

Rules for this toolkit (mirrors scanpy-workflow/_common.py):
  * stdout carries exactly one JSON object; every progress line goes to stderr.
  * a missing optional package is reported in the JSON result, never fatal.
  * arguments arrive as `--flag value` pairs built from the declared schema;
    a flag with no value is boolean True, a repeated flag collects a list.

Spatial inputs are dispatched to the format readers in ``readers/`` (Visium,
Visium HD, Stereo-seq GEF/GEM). The notebook load fragment carries the same
logic inline so a rendered notebook stays self-contained in the user's project.
"""

import datetime
import hashlib
import importlib.util
import json
import os
import struct
import sys

import numpy as np

TOOL_MAP = {
    "project_inspect.py": "st_py_inspect",
    "qc_filter.py": "st_py_qc",
    "preprocess_cluster.py": "st_py_preprocess",
    "spatial_statistics.py": "st_py_spatial",
    "niche.py": "st_py_niche",
    "trajectory.py": "st_py_trajectory",
    "find_markers.py": "st_py_markers",
    "annotate_clusters.py": "st_py_annotate",
    "deconvolve.py": "st_py_deconvolve",
    "communication.py": "st_py_communication",
    "build_notebook.py": "st_py_notebook",
    "env_probe.py": "st_py_env_probe",
}


def script_dir():
    return os.path.dirname(os.path.abspath(__file__))


def info(*parts):
    print("[spatial] " + "".join(str(p) for p in parts), file=sys.stderr)


def fail(*parts):
    message = "".join(str(p) for p in parts)
    print("[spatial] error: " + message, file=sys.stderr)
    print(json.dumps({"error": message}, indent=2, ensure_ascii=False))
    sys.exit(1)


def _clean(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return [_clean(v) for v in value.tolist()]
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_clean(v) for v in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def _write_artifact_sidecar(path, tool, skill, input_path=None):
    ext = os.path.splitext(path)[1].lower()
    figures = {".png": "image/png", ".pdf": "application/pdf", ".svg": "image/svg+xml"}
    tables = {".tsv": "text/tab-separated-values", ".csv": "text/csv"}
    if ext in figures:
        kind, media = "figure", figures[ext]
    elif ext in tables:
        kind, media = "table", tables[ext]
    elif ext == ".html":
        kind, media = "report", "text/html"
    else:
        return None
    stem = os.path.splitext(os.path.basename(path))[0].replace("-", " ").replace("_", " ").strip()
    title = (f"{tool}: {stem}" if stem else tool)[:200]
    desc = {
        "contractVersion": "1.0.0",
        "kind": kind,
        "file": os.path.basename(path),
        "mediaType": media,
        "title": title,
        "provenance": {
            "createdAt": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "tool": tool,
            "skill": skill,
        },
    }
    if kind == "figure":
        fig_block = {"format": ext[1:]}
        if ext == ".png":
            try:
                with open(path, "rb") as fh:
                    hdr = fh.read(24)
                if hdr[:8] == b"\x89PNG\r\n\x1a\n":
                    w, h = struct.unpack(">II", hdr[16:24])
                    fig_block["widthPx"], fig_block["heightPx"] = int(w), int(h)
            except OSError:
                pass
        desc["figure"] = fig_block
    elif kind == "table":
        try:
            sep = "," if ext == ".csv" else "\t"
            with open(path, "r", encoding="utf-8", errors="replace") as fh:
                first = fh.readline()
                n_cols = len(first.rstrip("\r\n").split(sep)) if first else 0
                n_rows = sum(1 for _ in fh)
            desc["table"] = {"rows": max(0, n_rows), "columns": n_cols}
        except OSError:
            pass
    if input_path and isinstance(input_path, str) and os.path.isfile(input_path):
        try:
            rel = os.path.relpath(os.path.realpath(input_path), os.path.realpath(os.getcwd()))
            if not rel.startswith("..") and not os.path.isabs(rel):
                with open(input_path, "rb") as fh:
                    sha = hashlib.sha256(fh.read()).hexdigest()
                desc["provenance"]["inputs"] = [{"path": rel, "sha256": sha}]
        except OSError:
            pass
    sidecar = path + ".phi-artifact.json"
    tmp = sidecar + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(desc, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    os.replace(tmp, sidecar)
    return path


def _attach_artifacts(result, skill="spatial-workflow"):
    if not isinstance(result, dict) or "artifacts" in result or result.get("error"):
        return result
    candidates = []

    def _add(val):
        if isinstance(val, str) and val:
            candidates.append(val)
        elif isinstance(val, (list, tuple)):
            for item in val:
                if isinstance(item, str):
                    candidates.append(item)
                elif isinstance(item, dict):
                    for k in ("path", "table", "file"):
                        if isinstance(item.get(k), str):
                            candidates.append(item[k])
        elif isinstance(val, dict):
            for k in ("path", "table", "file"):
                if isinstance(val.get(k), str):
                    candidates.append(val[k])

    for key in (
        "figures",
        "table",
        "tables",
        "outputs",
        "written",
        "svg",
        "niche",
        "communication",
        "deconvolution",
        "trajectory",
        "gene_scores",
        "spatial_diff",
    ):
        if key in result:
            _add(result[key])
    script_name = os.path.basename(sys.argv[0]) if sys.argv else ""
    tool = TOOL_MAP.get(script_name, "st_py_" + os.path.splitext(script_name)[0])
    inp = result.get("input") if isinstance(result.get("input"), str) else None
    seen = set()
    written = []
    for cand in candidates:
        if cand in seen or not os.path.isfile(cand):
            continue
        seen.add(cand)
        art = _write_artifact_sidecar(cand, tool=tool, skill=skill, input_path=inp)
        if art:
            written.append(art)
    if written:
        result = {**result, "artifacts": written[:8]}
    return result


def emit(result):
    print(json.dumps(_clean(_attach_artifacts(result)), indent=2, ensure_ascii=False))


# --- argument parsing ---------------------------------------------------------


def parse_flags(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    out = {}
    index = 0
    while index < len(argv):
        token = argv[index]
        if not token.startswith("--"):
            index += 1
            continue
        key = token[2:].replace("-", "_")
        has_value = index + 1 < len(argv) and not argv[index + 1].startswith("--")
        if not has_value:
            out[key] = True
            index += 1
        else:
            value = argv[index + 1]
            if key not in out:
                out[key] = value
            elif isinstance(out[key], list):
                out[key].append(value)
            else:
                out[key] = [out[key], value]
            index += 2
    return out


def flag_str(flags, name, default=None):
    value = flags.get(name)
    if value is None or value is True:
        return default
    if isinstance(value, list):
        value = value[0]
    return str(value)


def flag_req(flags, name):
    value = flag_str(flags, name)
    if not value:
        fail("missing required argument --" + name.replace("_", "-"))
    return value


def flag_int(flags, name, default=None):
    value = flag_str(flags, name)
    if value is None:
        return default
    try:
        return int(value)
    except ValueError:
        fail("--" + name.replace("_", "-") + " must be an integer, got " + value)


def flag_num(flags, name, default=None):
    value = flag_str(flags, name)
    if value is None:
        return default
    try:
        return float(value)
    except ValueError:
        fail("--" + name.replace("_", "-") + " must be a number, got " + value)


def flag_bool(flags, name, default=False):
    value = flags.get(name)
    if value is None:
        return default
    if value is True:
        return True
    if isinstance(value, list):
        value = value[0]
    return str(value).lower() in ("true", "1", "yes", "on")


def flag_vec(flags, name, default=()):
    value = flags.get(name)
    if value is None or value is True:
        return list(default)
    if isinstance(value, list):
        return [str(v) for v in value]
    return [str(value)]


def flag_nums(flags, name, default=()):
    values = flag_vec(flags, name, default)
    try:
        return [float(v) for v in values]
    except ValueError:
        fail("--" + name.replace("_", "-") + " takes numbers only")


# --- packages -----------------------------------------------------------------


def optional_module(name):
    if importlib.util.find_spec(name) is not None:
        return True
    info("optional package not installed: " + name)
    return False


def require_scanpy():
    try:
        import scanpy as sc
    except ImportError:
        fail("scanpy is not installed in this environment; use the plugin:spatial environment")
    return sc


def require_squidpy():
    try:
        import squidpy as sq
    except ImportError:
        fail("squidpy is not installed in this environment; use the plugin:spatial environment")
    return sq


# --- objects ------------------------------------------------------------------


def detect_kind(path):
    """Classify a spatial input. Sniffs HDF5 keys for .h5/.gef to tell platforms apart."""
    from readers import sniff_h5

    if os.path.isdir(path):
        if os.path.isdir(os.path.join(path, "spatial")):
            return "visium-dir"
        nested = os.path.join(path, "filtered_feature_bc_matrix")
        if os.path.isdir(nested) and os.path.isdir(os.path.join(path, "spatial")):
            return "visium-dir"
        markers = ("features.tsv.gz", "genes.tsv.gz", "features.tsv", "genes.tsv")
        has_markers = any(os.path.exists(os.path.join(path, m)) for m in markers)
        has_matrix = os.path.exists(os.path.join(path, "matrix.mtx.gz")) or os.path.exists(
            os.path.join(path, "matrix.mtx")
        )
        if has_markers and has_matrix:
            return "10x-matrix-dir"
        return "unknown-directory"
    lower = path.lower()
    if lower.endswith(".h5ad"):
        return "anndata"
    if lower.endswith(".gef"):
        return "stereoseq-gef"
    if lower.endswith((".gem", ".gem.gz")):
        return "stereoseq-gem"
    if lower.endswith(".h5"):
        return sniff_h5(path)
    if lower.endswith(".loom"):
        return "loom"
    if lower.endswith((".csv", ".tsv", ".txt", ".csv.gz", ".tsv.gz")):
        return "dense-table"
    return "unknown-file"


def read_object(path, platform="auto", bin_size=None, cellbin=False, var_names="gene_symbols"):
    """Load any supported spatial input as AnnData with a counts layer + obsm['spatial']."""
    from readers import read_spatial

    if not os.path.exists(path):
        fail("input not found: " + path)
    kind = detect_kind(path)
    info("input kind: " + kind)
    adata, meta = read_spatial(
        path, kind, platform=platform, bin_size=bin_size, cellbin=cellbin, var_names=var_names
    )
    if "counts" not in adata.layers:
        info("no counts layer; snapshotting X as counts")
        adata.layers["counts"] = adata.X.copy()
    adata.uns["spatial_platform"] = meta
    return adata


def write_object(adata, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    adata.write_h5ad(path)
    info("wrote " + path)
    return path


def ensure_dir(path, fallback):
    target = path if path else fallback
    os.makedirs(target, exist_ok=True)
    return target


# --- cluster levels -----------------------------------------------------------


def sorted_levels(values):
    uniques = sorted({str(v) for v in values})
    try:
        return sorted(uniques, key=float)
    except ValueError:
        return uniques