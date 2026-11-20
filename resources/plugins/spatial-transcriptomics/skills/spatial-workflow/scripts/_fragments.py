"""Fragment runner: the assets/fragments library is the single analysis
source for both the notebook renderer and the script tools.

A fragment is plain teaching code executed in a shared namespace that already
holds the config globals (same names as the notebook config cell). Fragments
leave bookkeeping behind (figure_paths, analysis_steps, qc_filters, ...); the
reports/ layer serializes it into the tool JSON contract.
"""

import contextlib
import os
import sys

from _common import info, script_dir


def fragment_path(step):
    return os.path.join(script_dir(), "..", "assets", "fragments", step + ".py")


def run_fragments(steps, namespace):
    # stdout carries exactly one JSON object: fragment progress (print, scanpy
    # logs) goes to stderr, the same discipline as the R adapters' sink(stderr()).
    with contextlib.redirect_stdout(sys.stderr):
        for step in steps:
            path = fragment_path(step)
            if not os.path.exists(path):
                raise ValueError("unknown fragment: " + step)
            info("fragment " + step)
            with open(path, encoding="utf-8") as handle:
                source = handle.read()
            exec(compile(source, path, "exec"), namespace)


def run_report(name, namespace):
    path = os.path.join(script_dir(), "reports", name + ".py")
    with open(path, encoding="utf-8") as handle:
        source = handle.read()
    exec(compile(source, path, "exec"), namespace)