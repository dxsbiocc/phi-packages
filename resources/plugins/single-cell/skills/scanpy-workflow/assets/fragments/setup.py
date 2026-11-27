import importlib.util
import os
import sys
from importlib.metadata import version

import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 2
sc.set_figure_params(dpi=120, facecolor="white")
sc.settings.autosave = False
np.random.seed(seed)


def info(*parts):
    print("[scanpy] " + "".join(str(p) for p in parts), file=sys.stderr)


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


# Bookkeeping shared with the script tools: the same fragments run there, and
# the report layer reads these variables to build the JSON result.
figure_paths = {}
analysis_steps = []
analysis_notes = []
qc_scored = []
qc_filters = []
obj_before = None
doublet_result = None
pca_variance = None
cluster_columns = []
embeddings = []
markers = None
marker_tables = []
annotation_panel = None
annotation_scores = None
annotation_assignments = None
annotation_automated = None
pb_info = None
integration_info = None
traj_info = None
ccc_info = None
reg_info = None
loaded_kind = None


def _figure_of(result):
    items = result
    if isinstance(items, dict):
        items = list(items.values())
    if not isinstance(items, (list, tuple)):
        items = [items]
    flat = []
    for item in items:
        if isinstance(item, np.ndarray):
            flat.extend(item.ravel().tolist())
        elif isinstance(item, (list, tuple)):
            flat.extend(item)
        else:
            flat.append(item)
    return next((item.figure for item in flat if hasattr(item, "figure")), None)


# Show a plot inline and, when figures_dir is set, also write a PNG so a run
# leaves reviewable figures behind. Returns the path (or None).
def show_plot(result, name):
    import matplotlib.pyplot as plt

    figure = _figure_of(result)
    if figure is None:
        return None
    path = None
    if figures_dir:
        os.makedirs(figures_dir, exist_ok=True)
        path = os.path.join(figures_dir, name + ".png")
        figure.savefig(path, dpi=150, bbox_inches="tight")
        figure_paths[name] = path
    if sc.settings.autoshow:
        try:
            from IPython.display import display

            display(figure)
        except ImportError:
            pass
    plt.close(figure)
    return path


print(f"scanpy {version('scanpy')}, seed {seed}")