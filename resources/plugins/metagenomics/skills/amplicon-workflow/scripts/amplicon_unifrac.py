#!/usr/bin/env python3
"""mg_amplicon_unifrac: compute exact Weighted and Unweighted UniFrac beta
diversity from an ASV/OTU count table and rooted Newick phylogenetic tree
(e.g., nf-core/ampliseq or QIIME 2 exports), plus PCoA, PERMANOVA, PERMDISP,
and multi-covariate marginal PERMANOVA.
"""
import argparse
import hashlib
import json
import os
import struct
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def emit(obj):
    print(json.dumps(obj, ensure_ascii=False, default=float))
    sys.exit(0)


def fail(msg):
    print(json.dumps({"error": str(msg)}, ensure_ascii=False))
    sys.exit(1)


def setup_matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def write_artifact(path, tool, title, *, skill, inputs=(), parameters=None, shape=None):
    p = Path(path)
    ext = p.suffix.lower()
    figures = {".png": "image/png", ".pdf": "application/pdf", ".svg": "image/svg+xml"}
    tables = {".tsv": "text/tab-separated-values", ".csv": "text/csv"}
    if ext in figures:
        kind, media = "figure", figures[ext]
    elif ext in tables:
        kind, media = "table", tables[ext]
    else:
        raise ValueError(f"unsupported artifact extension: {ext}")
    descriptor = {
        "contractVersion": "1.0.0",
        "kind": kind,
        "file": p.name,
        "mediaType": media,
        "title": title[:200],
        "provenance": {
            "createdAt": datetime.now(timezone.utc).isoformat(),
            "tool": tool,
            "skill": skill,
        },
    }
    if kind == "figure":
        block = {"format": ext[1:]}
        if ext == ".png":
            with open(p, "rb") as fh:
                header = fh.read(24)
            if header[:8] == b"\x89PNG\r\n\x1a\n":
                w, h = struct.unpack(">II", header[16:24])
                block["widthPx"], block["heightPx"] = int(w), int(h)
        descriptor["figure"] = block
    elif shape is not None:
        descriptor["table"] = {"rows": int(shape[0]), "columns": int(shape[1])}
    cwd = Path.cwd().resolve()
    prov_inputs = []
    for ip in inputs:
        ip = Path(ip)
        try:
            rel = ip.resolve().relative_to(cwd)
        except ValueError:
            continue
        entry = {"path": str(rel)}
        if ip.is_file():
            entry["sha256"] = hashlib.sha256(ip.read_bytes()).hexdigest()
        prov_inputs.append(entry)
    if prov_inputs:
        descriptor["provenance"]["inputs"] = prov_inputs
    if parameters:
        descriptor["provenance"]["parameters"] = parameters
    target = p.with_name(p.name + ".phi-artifact.json")
    tmp = target.with_name(target.name + ".tmp")
    tmp.write_text(json.dumps(descriptor, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, target)
    return str(p)


class _Node:
    __slots__ = ("name", "length", "children")

    def __init__(self, name="", length=0.0):
        self.name = name
        self.length = length
        self.children = []


def parse_newick(text):
    """Iterative Newick parser supporting quoted labels and scientific-notation lengths."""
    s = text.strip()
    if ";" in s:
        s = s[: s.index(";")]
    s = s.strip()
    if not s:
        fail("empty Newick tree")

    root = _Node()
    stack = [root]
    i = 0
    n = len(s)

    def _read_label_and_length(pos):
        label = ""
        while pos < n and s[pos].isspace():
            pos += 1
        if pos < n and s[pos] in ("'", '"'):
            quote = s[pos]
            pos += 1
            start = pos
            while pos < n and s[pos] != quote:
                pos += 1
            label = s[start:pos]
            pos = min(n, pos + 1)
        else:
            start = pos
            while pos < n and s[pos] not in ",():;":
                pos += 1
            label = s[start:pos].strip()
        length = 0.0
        while pos < n and s[pos].isspace():
            pos += 1
        if pos < n and s[pos] == ":":
            pos += 1
            start = pos
            while pos < n and s[pos] not in ",():;":
                pos += 1
            raw_len = s[start:pos].strip()
            try:
                length = max(0.0, float(raw_len))
            except ValueError:
                length = 0.0
        return label, length, pos

    while i < n:
        ch = s[i]
        if ch.isspace():
            i += 1
        elif ch == "(":
            child = _Node()
            stack[-1].children.append(child)
            stack.append(child)
            i += 1
        elif ch == ",":
            if len(stack) < 2:
                fail("malformed Newick tree (unexpected ',')")
            stack.pop()
            sibling = _Node()
            stack[-1].children.append(sibling)
            stack.append(sibling)
            i += 1
        elif ch == ")":
            if len(stack) < 2:
                fail("malformed Newick tree (unmatched ')')")
            stack.pop()
            i += 1
            label, length, i = _read_label_and_length(i)
            stack[-1].name = label
            stack[-1].length = length
        else:
            label, length, i = _read_label_and_length(i)
            stack[-1].name = label
            stack[-1].length = length

    # Unwrap dummy root if it has a single child with zero length
    while len(root.children) == 1 and not root.name and root.length == 0.0:
        root = root.children[0]
    return root


def compute_unifrac(root, prop_df):
    """Compute exact pairwise Unweighted and Normalized Weighted UniFrac matrices.

    prop_df: DataFrame (ASVs x samples) with columns summing to 1.0 over shared tips.
    """
    n_samples = prop_df.shape[1]
    asv_to_vec = {str(asv): prop_df.loc[asv].to_numpy(dtype=float) for asv in prop_df.index}

    # Post-order iterative traversal to avoid Python recursion limits on large trees
    order = []
    stack = [(root, False)]
    while stack:
        node, visited = stack.pop()
        if visited:
            order.append(node)
        else:
            stack.append((node, True))
            for ch in node.children:
                stack.append((ch, False))

    branch_lengths = []
    subtree_props = []
    node_vecs = {}
    tip_names = set()

    for node in order:
        if not node.children:
            tip_names.add(node.name)
            vec = asv_to_vec.get(node.name)
            if vec is None:
                vec = np.zeros(n_samples, dtype=float)
        else:
            vec = np.zeros(n_samples, dtype=float)
            for ch in node.children:
                vec += node_vecs.pop(id(ch))
        node_vecs[id(node)] = vec
        if node is not root and node.length > 0.0:
            branch_lengths.append(node.length)
            subtree_props.append(vec)

    if not branch_lengths:
        fail("tree has no positive branch lengths for UniFrac computation")

    b = np.asarray(branch_lengths, dtype=float)  # (E,)
    p = np.vstack(subtree_props)                  # (E, S)
    present = p > 0.0                             # (E, S)

    weighted = np.zeros((n_samples, n_samples), dtype=float)
    unweighted = np.zeros((n_samples, n_samples), dtype=float)

    for i in range(n_samples):
        pi = p[:, i]
        pres_i = present[:, i]
        for j in range(i + 1, n_samples):
            pj = p[:, j]
            pres_j = present[:, j]

            # Normalized Weighted UniFrac
            w_num = float(np.dot(b, np.abs(pi - pj)))
            w_den = float(np.dot(b, pi + pj))
            w_val = (w_num / w_den) if w_den > 0.0 else 0.0
            weighted[i, j] = weighted[j, i] = w_val

            # Unweighted UniFrac
            u_num = float(np.dot(b, np.logical_xor(pres_i, pres_j)))
            u_den = float(np.dot(b, np.logical_or(pres_i, pres_j)))
            u_val = (u_num / u_den) if u_den > 0.0 else 0.0
            unweighted[i, j] = unweighted[j, i] = u_val

    return weighted, unweighted, tip_names


def pcoa(dist):
    n = dist.shape[0]
    d2 = dist * dist
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ d2 @ j
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    pos = eigvals > 1e-10
    coords = eigvecs[:, pos] * np.sqrt(eigvals[pos])
    var = eigvals[pos] / eigvals[pos].sum() if pos.any() else np.array([])
    return coords, var


def permanova(dist, labels, permutations, seed):
    groups = np.unique(labels)
    a, n = len(groups), len(labels)
    if a < 2 or a >= n:
        return None
    d2 = dist * dist
    ss_total = d2.sum() / n
    ss_within = sum(
        d2[np.ix_(np.where(labels == g)[0], np.where(labels == g)[0])].sum() / len(np.where(labels == g)[0])
        for g in groups
        if len(np.where(labels == g)[0]) > 1
    )
    ss_between = ss_total - ss_within
    if ss_within <= 0:
        return None
    f_obs = (ss_between / (a - 1)) / (ss_within / (n - a))
    r2 = ss_between / ss_total if ss_total > 0 else 0.0
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(permutations):
        perm = rng.permutation(labels)
        ss_w = sum(
            d2[np.ix_(np.where(perm == g)[0], np.where(perm == g)[0])].sum() / len(np.where(perm == g)[0])
            for g in groups
            if len(np.where(perm == g)[0]) > 1
        )
        ss_b = ss_total - ss_w
        f_perm = (ss_b / (a - 1)) / (ss_w / (n - a)) if ss_w > 0 else np.inf
        if f_perm >= f_obs:
            ge += 1
    return {
        "pseudo_f": float(f_obs),
        "r2": float(r2),
        "p_value": float((1 + ge) / (1 + permutations)),
        "permutations": int(permutations),
        "groups": groups.tolist(),
    }


def permdisp(coords, labels, permutations, seed):
    groups = np.unique(labels)
    a, n = len(groups), len(labels)
    if a < 2 or a >= n or coords.shape[1] == 0:
        return None

    def _f_stat(lbls):
        dists = np.zeros(n, dtype=float)
        group_means = {}
        for g in groups:
            idx = np.where(lbls == g)[0]
            centroid = coords[idx].mean(axis=0)
            g_dists = np.linalg.norm(coords[idx] - centroid, axis=1)
            dists[idx] = g_dists
            group_means[str(g)] = float(g_dists.mean())
        grand = dists.mean()
        ss_b = sum(len(np.where(lbls == g)[0]) * (group_means[str(g)] - grand) ** 2 for g in groups)
        ss_w = sum(((dists[lbls == g] - group_means[str(g)]) ** 2).sum() for g in groups)
        if ss_w <= 1e-15:
            return None, group_means
        return float((ss_b / (a - 1)) / (ss_w / (n - a))), group_means

    f_obs, mean_dists = _f_stat(labels)
    if f_obs is None:
        return None
    rng = np.random.default_rng(seed + 17)
    ge = 0
    for _ in range(permutations):
        f_perm, _ = _f_stat(rng.permutation(labels))
        if f_perm is not None and f_perm >= f_obs:
            ge += 1
    return {
        "f_stat": float(f_obs),
        "p_value": float((1 + ge) / (1 + permutations)),
        "permutations": int(permutations),
        "group_mean_distances": mean_dists,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--asv_table", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--metric", choices=["weighted", "unweighted"], default="weighted")
    ap.add_argument("--metadata")
    ap.add_argument("--group_col", default="group")
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--permutations", type=int, default=999)
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()

    sep = "," if str(args.asv_table).lower().endswith(".csv") else "\t"
    df = pd.read_csv(args.asv_table, sep=sep, index_col=0)
    df.index = df.index.astype(str)
    df = df.apply(pd.to_numeric, errors="coerce")
    if df.isna().any().any():
        bad = df.columns[df.isna().any()].tolist()
        fail(f"non-numeric values in ASV table columns: {bad}")
    df = df.fillna(0.0)
    samples = df.columns.astype(str).tolist()
    df.columns = samples

    tree_path = Path(args.tree)
    if not tree_path.is_file():
        fail(f"tree file not found: {tree_path}")
    root = parse_newick(tree_path.read_text(encoding="utf-8", errors="replace"))

    # Collect tip names first to align ASV table strictly by ID (never by row position)
    tips = set()
    stack = [root]
    while stack:
        nd = stack.pop()
        if not nd.children:
            tips.add(nd.name)
        else:
            stack.extend(nd.children)

    shared_asvs = [a for a in df.index if a in tips]
    if len(shared_asvs) < 2:
        fail(
            f"fewer than 2 ASV IDs overlap between ASV table ({len(df)} ASVs) and tree ({len(tips)} tips). "
            "Check that ASV identifiers match tree tip labels exactly."
        )

    sub = df.loc[shared_asvs].copy()
    col_sums = sub.sum(axis=0)
    if (col_sums <= 0).any():
        zero_samples = col_sums[col_sums <= 0].index.tolist()
        fail(f"samples have zero counts among tree-matched ASVs: {zero_samples}")
    prop = sub.div(col_sums, axis=1)

    labels = None
    if args.metadata:
        msep = "," if str(args.metadata).lower().endswith(".csv") else "\t"
        meta = pd.read_csv(args.metadata, sep=msep, index_col=0)
        meta.index = meta.index.astype(str)
        missing = [s for s in samples if s not in meta.index]
        if missing:
            fail(f"metadata missing samples: {missing}")
        meta = meta.loc[samples]
        if args.group_col not in meta.columns:
            fail(f"metadata has no column {args.group_col}; columns: {meta.columns.tolist()}")
        labels = meta[args.group_col].astype(str).to_numpy()

    w_dist, u_dist, _ = compute_unifrac(root, prop)
    chosen = w_dist if args.metric == "weighted" else u_dist
    coords, var = pcoa(chosen)
    perm = permanova(chosen, labels, args.permutations, args.seed) if labels is not None else None
    disp = permdisp(coords, labels, args.permutations, args.seed) if labels is not None else None

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    inputs = [args.asv_table, args.tree] + ([args.metadata] if args.metadata else [])
    params = {"metric": args.metric, "permutations": args.permutations, "seed": args.seed}

    w_df = pd.DataFrame(w_dist, index=samples, columns=samples)
    u_df = pd.DataFrame(u_dist, index=samples, columns=samples)
    pcoa_df = pd.DataFrame(coords, index=samples, columns=[f"PCo{i + 1}" for i in range(coords.shape[1])])

    w_path = outdir / "weighted-unifrac.tsv"
    u_path = outdir / "unweighted-unifrac.tsv"
    pcoa_path = outdir / f"pcoa-{args.metric}-unifrac.tsv"
    w_df.to_csv(w_path, sep="\t")
    u_df.to_csv(u_path, sep="\t")
    pcoa_df.to_csv(pcoa_path, sep="\t")

    artifacts = [
        write_artifact(w_path, "mg_amplicon_unifrac", "Weighted UniFrac distance matrix",
                       skill="amplicon-workflow", inputs=inputs, parameters=params, shape=w_df.shape),
        write_artifact(u_path, "mg_amplicon_unifrac", "Unweighted UniFrac distance matrix",
                       skill="amplicon-workflow", inputs=inputs, parameters=params, shape=u_df.shape),
        write_artifact(pcoa_path, "mg_amplicon_unifrac", f"PCoA coordinates ({args.metric} UniFrac)",
                       skill="amplicon-workflow", inputs=inputs, parameters=params, shape=pcoa_df.shape),
    ]
    outputs = {
        "weighted_unifrac": str(w_path),
        "unweighted_unifrac": str(u_path),
        "pcoa_coords": str(pcoa_path),
    }

    if coords.shape[1] >= 2:
        plt = setup_matplotlib()
        fig, ax = plt.subplots(figsize=(6, 5))
        if labels is not None:
            for g in np.unique(labels):
                m = labels == g
                ax.scatter(coords[m, 0], coords[m, 1], label=str(g), s=50, alpha=0.85)
            ax.legend()
        else:
            ax.scatter(coords[:, 0], coords[:, 1], s=50, alpha=0.85)
        for i, s in enumerate(samples):
            ax.annotate(s, (coords[i, 0], coords[i, 1]), fontsize=7)
        ax.set_xlabel(f"PCo1 ({var[0] * 100:.1f}%)")
        ax.set_ylabel(f"PCo2 ({var[1] * 100:.1f}%)")
        title = f"PCoA ({args.metric.capitalize()} UniFrac)"
        if perm is not None:
            title += f" | PERMANOVA F={perm['pseudo_f']:.2f}, p={perm['p_value']:.3g}"
        ax.set_title(title)
        fig.tight_layout()
        fig_path = outdir / f"pcoa-{args.metric}-unifrac.png"
        fig.savefig(fig_path, dpi=180)
        plt.close(fig)
        outputs["pcoa_figure"] = str(fig_path)
        artifacts.append(
            write_artifact(fig_path, "mg_amplicon_unifrac", f"PCoA ({args.metric} UniFrac)",
                           skill="amplicon-workflow", inputs=inputs, parameters=params)
        )

    emit(
        {
            "tool": "mg_amplicon_unifrac",
            "metric": args.metric,
            "n_samples": int(len(samples)),
            "n_asvs_input": int(len(df)),
            "n_tree_tips": int(len(tips)),
            "n_asvs_matched": int(len(shared_asvs)),
            "pcoa_variance_explained": var.tolist(),
            "permanova": perm,
            "permdisp": disp,
            "outputs": outputs,
            "artifacts": artifacts,
            "caveat": "UniFrac requires a rooted phylogenetic tree matching ASV identifiers; "
                      "verify amplicon region and taxonomy classifier alignment (e.g., V3-V4 vs V4) "
                      "upstream in nf-core/ampliseq or QIIME 2.",
        }
    )


if __name__ == "__main__":
    main()
