#!/usr/bin/env python3
"""Shared helpers for the scatac-workflow script toolkit.

Rules for this toolkit (mirrors spatial-workflow/_common.py):
  * stdout carries exactly one JSON object; every progress line goes to stderr.
  * a missing optional package or genome resource is reported in the JSON
    result, never fatal.
  * arguments arrive as `--flag value` pairs built from the declared schema;
    a flag with no value is boolean True, a repeated flag collects a list.

scATAC inputs are fragments (.fragments.gz / .tsv.gz), a position-sorted BAM,
a saved snapatac2 .h5ad, or a 10x peak/tile matrix directory/.h5. The notebook
load fragment carries the same import logic inline so a rendered notebook stays
self-contained in the user's project.

snapatac2 AnnData is file-backed: .obs is a lazy element that supports
``obs["col"].to_numpy()`` and ``"col" in adata.obs`` but NOT iteration, .columns
or .keys() (those panic in the Rust layer). Use obs_col / obs_has below.
"""

import importlib.util
import json
import os
import sys

import numpy as np


def script_dir():
    return os.path.dirname(os.path.abspath(__file__))


def info(*parts):
    print("[scatac] " + "".join(str(p) for p in parts), file=sys.stderr)


def fail(*parts):
    message = "".join(str(p) for p in parts)
    print("[scatac] error: " + message, file=sys.stderr)
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


def emit(result):
    print(json.dumps(_clean(result), indent=2, ensure_ascii=False))


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


def require_snapatac2():
    try:
        import snapatac2 as snap
    except ImportError:
        fail("snapatac2 is not installed in this environment; use the plugin:scatac environment")
    return snap


def require_scanpy():
    try:
        import scanpy as sc
    except ImportError:
        fail("scanpy is not installed in this environment; use the plugin:scatac environment")
    return sc


# --- backed AnnData obs/var access --------------------------------------------


def obs_has(adata, name):
    """Membership test on a snapatac2 backed obs element (safe; no iteration)."""
    try:
        return name in adata.obs
    except Exception:
        return False


def obs_col(adata, name):
    """Return an obs column as a numpy array, or None when absent."""
    if not obs_has(adata, name):
        return None
    try:
        return np.asarray(adata.obs[name].to_numpy())
    except Exception:
        try:
            return np.asarray(adata.obs[name])
        except Exception:
            return None


def var_col(adata, name):
    """Return a var column as numpy, or None. var membership is unreliable on
    the backed element, so probe by item access instead."""
    if adata.n_vars == 0:
        return None
    try:
        return np.asarray(adata.var[name].to_numpy())
    except Exception:
        return None


def obs_columns(adata, candidates):
    """Enumerate obs columns by probing a candidate list (the backed element
    exposes no .columns). Returns the candidates that are present."""
    return [c for c in candidates if obs_has(adata, c)]


# --- genome resources ---------------------------------------------------------


BUNDLED_GENOMES = ("hg38", "mm10", "hg19", "mm39")


def resolve_chrom_sizes(snap, genome, chrom_sizes_path):
    """chrom_sizes for import: a bundled genome dict (offline) or a user file.

    snapatac2 ships chrom_sizes for hg38/mm10/hg19/mm39 as a hardcoded dict; the
    annotation and FASTA are lazy downloads, so passing the dict keeps import
    offline. A custom genome needs a two-column chrom<TAB>size file.
    """
    if chrom_sizes_path:
        if not os.path.exists(chrom_sizes_path):
            fail("chrom_sizes file not found: " + chrom_sizes_path)
        sizes = {}
        with open(chrom_sizes_path, encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        sizes[parts[0]] = int(parts[1])
                    except ValueError:
                        continue
        if not sizes:
            fail("no chromosome sizes parsed from " + chrom_sizes_path)
        return sizes, {"source": "file", "path": chrom_sizes_path, "n_chrom": len(sizes)}
    name = (genome or "hg38").lower()
    if name not in BUNDLED_GENOMES:
        fail(
            "genome '" + str(genome) + "' has no bundled chrom_sizes; "
            "pass --chrom-sizes with a two-column chrom<TAB>size file"
        )
    g = getattr(snap.genome, name)
    return dict(g.chrom_sizes), {"source": "bundled", "genome": name, "n_chrom": len(g.chrom_sizes)}


# --- objects ------------------------------------------------------------------


def detect_kind(path):
    """Classify a scATAC input."""
    if os.path.isdir(path):
        for sub in ("filtered_peak_bc_matrix", "filtered_feature_bc_matrix", "raw_peak_bc_matrix"):
            if os.path.isdir(os.path.join(path, sub)):
                return "10x-matrix-dir"
        markers = ("features.tsv.gz", "genes.tsv.gz", "peaks.bed.gz", "features.tsv", "genes.tsv")
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
    if lower.endswith(".bam"):
        return "bam"
    if lower.endswith((".bed", ".bed.gz", ".narrowpeak", ".broadpeak")):
        return "peak-file"
    if lower.endswith(".h5"):
        return "10x-h5"
    if "fragment" in lower and lower.endswith((".gz", ".tsv", ".txt", ".bgz")):
        return "fragment"
    if lower.endswith((".fragments.gz", ".fragments.tsv.gz")):
        return "fragment"
    if lower.endswith((".tsv.gz", ".txt.gz", ".tsv", ".csv")):
        return "fragment-like-table"
    return "unknown-file"


def read_object(path, snap, genome="hg38", chrom_sizes_path=None, min_num_fragments=200,
                backing=None, is_paired=True, n_jobs=1):
    """Load any supported scATAC input as a snapatac2 AnnData (backed)."""
    if not os.path.exists(path):
        fail("input not found: " + path)
    kind = detect_kind(path)
    info("input kind: " + kind)
    if kind == "anndata":
        return snap.read(path, backed="r+"), kind
    if kind == "10x-matrix-dir":
        return snap.read_10x_mtx(path), kind
    if kind == "10x-h5":
        return snap.read_10x_mtx(path), kind
    if kind in ("fragment", "fragment-like-table"):
        sizes, _ = resolve_chrom_sizes(snap, genome, chrom_sizes_path)
        adata = snap.pp.import_fragments(
            path, chrom_sizes=sizes, sorted_by_barcode=False,
            min_num_fragments=min_num_fragments, file=backing, is_paired=is_paired, n_jobs=n_jobs,
        )
        return adata, kind
    if kind == "bam":
        if not backing:
            fail("BAM import needs a backing .h5ad path (pass --output)")
        frag = backing + ".fragments.tsv.gz"
        snap.pp.make_fragment_file(path, frag, is_paired=is_paired, source="10x")
        sizes, _ = resolve_chrom_sizes(snap, genome, chrom_sizes_path)
        adata = snap.pp.import_fragments(
            frag, chrom_sizes=sizes, sorted_by_barcode=False,
            min_num_fragments=min_num_fragments, file=backing, is_paired=is_paired, n_jobs=n_jobs,
        )
        return adata, kind
    fail("unsupported scATAC input kind: " + kind + " (" + path + ")")


def write_object(adata, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    adata.write(path)
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