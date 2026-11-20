# JSON contract for sc_py_trajectory.
emit(
    {
        "input": tool_input,
        "output": tool_output,
        **traj_info,
        "figures": [
            f
            for f in (figure_paths.get("trajectory-paga"), figure_paths.get("trajectory-pseudotime"))
            if f
        ],
        "notes": [
            "pseudotime is relative to the root cluster, not absolute ordering",
            "unreachable cells (inf pseudotime) sit in disconnected components",
        ],
    }
)