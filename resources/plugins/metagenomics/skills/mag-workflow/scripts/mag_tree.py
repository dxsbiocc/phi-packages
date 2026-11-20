"""Patristic distances between MAGs and to nearest GTDB references from a Newick tree."""
import argparse
import re
from pathlib import Path

import pandas as pd

from _common import bin_stem, emit, fail, read_tsv, write_artifact

REF_PATTERN = re.compile(r"^(RS_|GB_)?GC[FA]_")
DELIMS = ":,();[] \t\n\r"


def read_label_length(text, i):
    n = len(text)
    start = i
    if i < n and text[i] == "'":
        i += 1
        while i < n and text[i] != "'":
            i += 1
        label = text[start + 1:i]
        i += 1
        if i < n and text[i] == "'":
            i += 1
    else:
        while i < n and text[i] not in DELIMS:
            i += 1
        label = text[start:i]
    length = 0.0
    if i < n and text[i] == ":":
        i += 1
        start = i
        while i < n and text[i] not in DELIMS:
            i += 1
        try:
            length = float(text[start:i])
        except ValueError:
            fail(f"bad branch length near: {text[start:i][:40]}")
    return label, length, i


def parse_newick(text):
    """Return (names, lengths, parent) arrays; children are not retained."""
    names, lengths, parent = [], [], []
    stack = []

    def new_node(par):
        idx = len(names)
        names.append(None)
        lengths.append(0.0)
        parent.append(par)
        return idx

    cur = new_node(None)
    i, n = 0, len(text)
    while i < n:
        ch = text[i]
        if ch == "(":
            stack.append(cur)
            cur = new_node(cur)
            i += 1
        elif ch == ",":
            cur = new_node(stack[-1])
            i += 1
        elif ch == ")":
            if not stack:
                fail("unbalanced Newick: ')' without '('")
            cur = stack.pop()
            i += 1
            label, length, i = read_label_length(text, i)
            lengths[cur] = length
        elif ch == ";":
            break
        elif ch in " \t\n\r":
            i += 1
        elif ch == "[":
            end = text.find("]", i)
            if end < 0:
                fail("unterminated Newick comment")
            i = end + 1
        else:
            label, length, i = read_label_length(text, i)
            names[cur] = label
            lengths[cur] = length
    if stack:
        fail("unbalanced Newick: unclosed '('")
    return names, lengths, parent


def ancestor_path(tip, parent, lengths):
    path = {}
    cum = 0.0
    node = tip
    while node is not None:
        path[node] = cum
        cum += lengths[node]
        node = parent[node]
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True)
    ap.add_argument("--summary")
    ap.add_argument("--output_dir", required=True)
    a = ap.parse_args()

    tree_path = Path(a.tree)
    if not tree_path.is_file():
        fail(f"tree file not found: {tree_path}")
    text = tree_path.read_text(encoding="utf-8", errors="replace").strip()
    if not text.startswith("(") and ";" not in text:
        fail(f"{tree_path} does not look like Newick")
    names, lengths, parent = parse_newick(text)
    tips = [i for i, name in enumerate(names) if name]
    if not tips:
        fail("tree has no labeled tips")

    if a.summary:
        s = read_tsv(a.summary)
        if "user_genome" not in s.columns:
            fail("summary lacks a user_genome column")
        wanted = {bin_stem(v) for v in s["user_genome"].astype(str)}
        identification = "summary"
    else:
        wanted = None
        identification = "heuristic"
    user_tips, ref_tips = [], []
    for t in tips:
        is_user = (bin_stem(names[t]) in wanted) if wanted is not None \
            else not REF_PATTERN.match(names[t])
        (user_tips if is_user else ref_tips).append(t)
    if not user_tips:
        fail("no user MAG tips found in the tree"
             + (" (check --summary user_genome names)" if a.summary else ""))

    user_names = [bin_stem(names[t]) for t in user_tips]
    paths = [ancestor_path(t, parent, lengths) for t in user_tips]

    n = len(user_tips)
    dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            node, cum = user_tips[j], 0.0
            while node not in paths[i]:
                cum += lengths[node]
                node = parent[node]
            d = cum + paths[i][node]
            dist[i][j] = dist[j][i] = d

    node_users = {}
    for ui, path in enumerate(paths):
        for node, cum in path.items():
            node_users.setdefault(node, []).append((ui, cum))
    best = [None] * n
    if ref_tips:
        for r in ref_tips:
            node, cum = r, 0.0
            while node is not None:
                entries = node_users.get(node)
                if entries:
                    for ui, cu in entries:
                        d = cum + cu
                        if best[ui] is None or d < best[ui][0]:
                            best[ui] = (d, names[r])
                cum += lengths[node]
                node = parent[node]

    outdir = Path(a.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    inputs = [a.tree] + ([a.summary] if a.summary else [])
    params = {"user_identification": identification}
    artifacts = []
    outputs = {}

    dm = pd.DataFrame(dist, index=user_names, columns=user_names)
    dm.index.name = "bin"
    dm_path = outdir / "mag-patristic-distances.tsv"
    dm.to_csv(dm_path, sep="\t")
    outputs["distance_matrix"] = str(dm_path)
    artifacts.append(write_artifact(dm_path, "mg_mag_tree",
                                    "Patristic distances between MAGs (GTDB tree)",
                                    skill="mag-workflow", inputs=inputs,
                                    parameters=params, shape=dm.shape))

    if ref_tips:
        rows = [{"bin": user_names[i],
                 "nearest_reference": best[i][1] if best[i] else None,
                 "patristic_distance": best[i][0] if best[i] else None}
                for i in range(n)]
        nr = pd.DataFrame(rows, columns=["bin", "nearest_reference", "patristic_distance"])
        nr_path = outdir / "mag-nearest-reference.tsv"
        nr.to_csv(nr_path, sep="\t", index=False)
        outputs["nearest_reference"] = str(nr_path)
        artifacts.append(write_artifact(nr_path, "mg_mag_tree",
                                        "Nearest GTDB reference per MAG (patristic)",
                                        skill="mag-workflow", inputs=inputs,
                                        parameters=params, shape=nr.shape))

    emit({
        "tool": "mg_mag_tree",
        "tree": str(tree_path),
        "n_tips": len(tips),
        "n_user_tips": n,
        "n_reference_tips": len(ref_tips),
        "user_identification": identification,
        "outputs": outputs,
        "artifacts": artifacts,
        "caveat": "Patristic distances come from the GTDB backbone placement tree: "
                  "they reflect marker-gene divergence, not ANI, and placement of "
                  "incomplete bins is uncertain. No tree figure is rendered here.",
    })


if __name__ == "__main__":
    main()