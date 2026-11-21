"""Build the fragment namespace from steps.json defaults + CLI flags.

The notebook renderer and the script tools share assets/fragments as the single
analysis source. This helper gives the script path the same config globals the
notebook config cell defines, so a fragment behaves identically in both. Flag
names match the parameter names in steps.json; structural globals (data_path,
output_object, ...) and adapter-only globals (emit, tool_input, tool_output) are
passed as overrides and win over the library defaults.
"""

import json
import os

from _common import flag_bool, flag_int, flag_num, flag_nums, flag_str, script_dir


def registry():
    path = os.path.join(script_dir(), "..", "assets", "steps.json")
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def config_namespace(flags, overrides=None):
    ns = {}
    for name, spec in registry()["parameters"].items():
        kind = spec["type"]
        default = spec["default"]
        if kind == "int":
            ns[name] = flag_int(flags, name, default)
        elif kind == "number":
            ns[name] = flag_num(flags, name, default)
        elif kind == "bool":
            ns[name] = flag_bool(flags, name, bool(default))
        elif kind == "list":
            ns[name] = flag_nums(flags, name, default or [])
        else:
            ns[name] = flag_str(flags, name, default)
    if overrides:
        ns.update(overrides)
    return ns


def structural(flags, input_path, output):
    """The globals every adapter shares: paths, figures and the report hooks."""
    from _common import emit

    return {
        "data_path": input_path,
        "project_name": flag_str(flags, "project", "phi"),
        "output_object": output,
        "sample_column": flag_str(flags, "sample_column", "sample"),
        "figures_dir": flag_str(flags, "figures_dir"),
        "emit": emit,
        "tool_input": input_path,
        "tool_output": None,
    }