import importlib.util
import os
import sys
from importlib.metadata import version

import numpy as np
import pandas as pd
import polars as pl
import snapatac2 as snap
import scanpy as sc

sc.settings.verbosity = 2
sc.settings.autoshow = False
np.random.seed(seed)


def info(*parts):
    print("[scatac] " + "".join(str(p) for p in parts), file=sys.stderr)


def optional_module(name):
    if importlib.util.find_spec(name) is not None:
        return True
    info("optional package not installed: " + name)
    return False


def sorted_levels(values):
    uniques = sorted({str(v) for v in values})
    try:
        return sorted(uniques, key=float)
    except ValueError:
        return uniques


def detect_kind(path):
    """Classify a scATAC input (mirrors scripts/_common.py for notebook self-containment)."""
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


# snapatac2 keeps a file-backed AnnData. Its .obs is a lazy Rust element that
# supports obs["col"].to_numpy() and "col" in adata.obs, but NOT .columns,
# .keys() or iteration (those panic). These helpers are the only safe access.
def obs_has(adata, name):
    try:
        return name in adata.obs
    except Exception:
        return False


def obs_col(adata, name):
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
    if adata.n_vars == 0:
        return None
    try:
        return np.asarray(adata.var[name].to_numpy())
    except Exception:
        return None


BUNDLED_GENOMES = ("hg38", "mm10", "hg19", "mm39")


def resolve_chrom_sizes():
    """chrom_sizes for import: a bundled genome dict (offline) or a user file.
    snapatac2 ships chrom_sizes for hg38/mm10/hg19/mm39; annotation and FASTA
    are lazy downloads, so the dict keeps import offline."""
    if chrom_sizes:
        if not os.path.exists(chrom_sizes):
            raise ValueError("chrom_sizes file not found: " + chrom_sizes)
        sizes = {}
        with open(chrom_sizes, encoding="utf-8") as handle:
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
            raise ValueError("no chromosome sizes parsed from " + chrom_sizes)
        return sizes, {"source": "file", "path": chrom_sizes, "n_chrom": len(sizes)}
    name = (genome or "hg38").lower()
    if name not in BUNDLED_GENOMES:
        raise ValueError(
            "genome '" + str(genome) + "' has no bundled chrom_sizes; pass chrom_sizes"
        )
    g = getattr(snap.genome, name)
    return dict(g.chrom_sizes), {"source": "bundled", "genome": name, "n_chrom": len(g.chrom_sizes)}


# Bookkeeping shared with the script tools: the same fragments run there, and
# the report layer reads these variables to build the JSON result.
figure_paths = {}
analysis_steps = []
analysis_notes = []
qc_scored = []
qc_filters = []
obj_before = None
adata = None
peak_adata = None
gene_adata = None
loaded_kind = None
genome_info = None
cluster_columns = []
embeddings = []
spectral_variance = None
peak_tables = None
merged_peaks = None
markers = None
marker_tables = []
diff_results = None
annotation_panel = None
annotation_scores = None
annotation_assignments = None
annotation_automated = None
motif_results = None
integration_info = None


# snapatac2 plotting returns a plotly Figure ONLY with show=False; otherwise it
# renders inline and returns None. show_plot takes that Figure, writes an
# interactive HTML copy when figures_dir is set, optionally a static PNG, and
# displays it inline in a notebook. Returns the HTML path (or None).
def show_plot(fig, name):
    if fig is None:
        return None
    path = None
    if figures_dir:
        os.makedirs(figures_dir, exist_ok=True)
        path = os.path.join(figures_dir, name + ".html")
        try:
            fig.write_html(path, include_plotlyjs="cdn")
            figure_paths[name] = path
        except Exception as error:
            info("could not write HTML for " + name + ": " + str(error))
        if static_png:
            try:
                png = os.path.join(figures_dir, name + ".png")
                fig.write_image(png)
                figure_paths[name + "-png"] = png
            except Exception as error:
                info("kaleido PNG export skipped for " + name + ": " + str(error))
    try:
        from IPython.display import display

        display(fig)
    except Exception:
        pass
    return path


print(f"snapatac2 {version('snapatac2')}, scanpy {version('scanpy')}, MACS3 {version('macs3')}, seed {seed}")