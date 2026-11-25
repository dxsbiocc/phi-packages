"""Shared utilities and 4PL Hill curve fitting for pharmacogenomics (pgx-workflow)."""
from __future__ import annotations

import hashlib
import json
import math
import os
import struct
import sys
from datetime import datetime, timezone
from pathlib import Path


def emit(payload: dict) -> None:
    sys.stdout.write(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    sys.stdout.flush()


def fail(message: str, code: int = 1) -> "NoReturn":
    sys.stderr.write(f"[pgx-workflow error] {message}\n")
    sys.stderr.flush()
    raise SystemExit(code)


def bool_arg(val) -> bool:
    return str(val).lower() not in ("false", "0", "no")


def read_table(path: str | Path):
    import pandas as pd

    p = Path(path)
    if not p.is_file():
        fail(f"Input file not found: {path}")
    name = p.name.lower()
    if name.endswith(".csv") or name.endswith(".csv.gz"):
        sep = ","
    elif name.endswith(".tsv") or name.endswith(".tsv.gz") or name.endswith(".txt") or name.endswith(".maf"):
        sep = "\t"
    else:
        with open(p, "rt", encoding="utf-8", errors="replace") as fh:
            first = ""
            for line in fh:
                if not line.startswith("#") and line.strip():
                    first = line
                    break
        sep = "\t" if "\t" in first else ","
    return pd.read_csv(p, sep=sep, comment="#")


def bh_fdr(pvals) -> list[float]:
    import numpy as np

    arr = np.asarray(pvals, dtype=float)
    n = arr.size
    if n == 0:
        return []
    safe = np.where(np.isfinite(arr), np.clip(arr, 0.0, 1.0), 1.0)
    order = np.argsort(safe)
    ranked = safe[order] * n / np.arange(1, n + 1, dtype=float)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(ranked, 0.0, 1.0)
    return out.tolist()


def setup_matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def detect_viability_scale(values, requested: str = "auto") -> tuple[str, float]:
    import numpy as np

    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        fail("No finite viability values found.")
    if requested == "percent":
        return "percent", 0.01
    if requested == "fraction":
        return "fraction", 1.0
    p95 = float(np.percentile(arr, 95))
    if p95 > 2.5:
        return "percent", 0.01
    return "fraction", 1.0


def hill_4pl_log10(log10_d, e0: float, einf: float, log10_ec50: float, hill: float):
    import numpy as np

    exponent = np.clip(hill * (np.asarray(log10_d, dtype=float) - log10_ec50), -50.0, 50.0)
    return einf + (e0 - einf) / (1.0 + np.power(10.0, exponent))


def fit_hill_4pl(doses, viabilities, *, fix_e0: float | None = None) -> dict:
    """Fit 4-parameter log-logistic Hill curve to positive doses and fractional viabilities (~0..1)."""
    import numpy as np
    from scipy import optimize

    d = np.asarray(doses, dtype=float)
    y = np.asarray(viabilities, dtype=float)
    mask = np.isfinite(d) & np.isfinite(y) & (d > 0)
    d = d[mask]
    y = y[mask]
    if d.size < 4:
        fail(f"At least 4 positive dose points are required to fit a 4PL Hill curve (found {d.size}).")

    order = np.argsort(d)
    d = d[order]
    y = y[order]
    log_d = np.log10(d)
    min_d, max_d = float(d[0]), float(d[-1])
    min_ld, max_ld = float(log_d[0]), float(log_d[-1])

    # Collapse replicate doses for trapezoidal AUC
    uniq_ld, inv = np.unique(log_d, return_inverse=True)
    if uniq_ld.size < 3:
        fail(f"At least 3 distinct positive concentrations are required (found {uniq_ld.size}).")
    mean_y = np.array([float(y[inv == i].mean()) for i in range(uniq_ld.size)], dtype=float)
    span_ld = max(1e-9, float(uniq_ld[-1] - uniq_ld[0]))
    trapz_fn = getattr(np, "trapezoid", None) or getattr(np, "trapz")
    auc_emp = float(np.clip(trapz_fn(np.clip(mean_y, 0.0, 1.2), uniq_ld) / span_ld, 0.0, 1.0))

    e0_init = float(np.clip(fix_e0 if fix_e0 is not None else np.median(mean_y[: max(1, min(2, len(mean_y)))]), 0.65, 1.25))
    einf_init = float(np.clip(np.min(mean_y), 0.0, min(1.0, e0_init - 0.01)))
    mid_target = 0.5 * (e0_init + einf_init)
    mid_idx = int(np.argmin(np.abs(mean_y - mid_target)))
    ec50_ld_init = float(uniq_ld[mid_idx])

    if fix_e0 is not None:
        e0_lo, e0_hi = max(0.05, fix_e0 - 0.02), fix_e0 + 0.02
    else:
        e0_lo, e0_hi = 0.55, 1.35

    lo = [e0_lo, 0.0, min_ld - 2.0, 0.15]
    hi = [e0_hi, 1.05, max_ld + 2.0, 8.0]

    best_popt = None
    best_sse = float("inf")
    starts = [
        [np.clip(e0_init, e0_lo + 1e-4, e0_hi - 1e-4), np.clip(einf_init, 0.0, 1.0), ec50_ld_init, 1.2],
        [np.clip(1.0, e0_lo + 1e-4, e0_hi - 1e-4), 0.05, 0.5 * (min_ld + max_ld), 1.0],
        [np.clip(1.0, e0_lo + 1e-4, e0_hi - 1e-4), float(np.clip(mean_y[-1], 0.0, 1.0)), max_ld, 2.0],
    ]
    for p0 in starts:
        p0_c = [float(np.clip(v, lo[i] + 1e-5, hi[i] - 1e-5)) for i, v in enumerate(p0)]
        try:
            popt, _ = optimize.curve_fit(
                hill_4pl_log10,
                log_d,
                y,
                p0=p0_c,
                bounds=(lo, hi),
                method="trf",
                maxfev=6000,
            )
            pred = hill_4pl_log10(log_d, *popt)
            sse = float(np.sum((y - pred) ** 2))
            if sse < best_sse:
                best_sse = sse
                best_popt = popt
        except Exception:
            continue

    if best_popt is None:
        # Flat fallback
        e0_f = float(np.clip(np.mean(y), 0.0, 1.2))
        best_popt = np.array([e0_f, e0_f, max_ld + 1.0, 1.0], dtype=float)

    e0, einf, log10_ec50, hill = [float(v) for v in best_popt]
    # Enforce monotonic inhibition (e0 >= einf); if inverted on noisy flat curve, flatten
    if e0 < einf:
        avg_v = float(np.clip(np.mean(y), 0.0, 1.1))
        e0, einf, log10_ec50, hill = avg_v, avg_v, max_ld + 1.0, 1.0

    y_hat = hill_4pl_log10(log_d, e0, einf, log10_ec50, hill)
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    ss_res = float(np.sum((y - y_hat) ** 2))
    r2 = float(max(0.0, 1.0 - ss_res / ss_tot)) if ss_tot > 1e-10 else 0.0

    ec50 = float(10.0 ** log10_ec50)
    # Absolute IC50 where predicted viability crosses 0.50
    ic50_raw: float | None = None
    if e0 > 0.50 and einf < 0.50 and (e0 - einf) > 1e-6:
        ratio = (e0 - 0.50) / (0.50 - einf)
        if ratio > 0:
            log10_ic50 = log10_ec50 + math.log10(ratio) / hill
            if log10_ic50 < 20.0:
                ic50_raw = float(10.0 ** log10_ic50)

    beyond_range = bool(ic50_raw is None or ic50_raw > max_d)
    ic50_censored = float(max_d if beyond_range else ic50_raw)

    # Model-integrated AUC on dense log-dose grid clipped to [0, 1]
    grid_ld = np.linspace(min_ld, max_ld, 200)
    grid_y = np.clip(hill_4pl_log10(grid_ld, e0, einf, log10_ec50, hill), 0.0, 1.0)
    auc_fit = float(np.clip(trapz_fn(grid_y, grid_ld) / max(1e-9, max_ld - min_ld), 0.0, 1.0))

    return {
        "e0": e0,
        "einf": einf,
        "emax_inhibition": float(max(0.0, 1.0 - einf)),
        "ec50": ec50,
        "log10_ec50": log10_ec50,
        "hill": hill,
        "ic50_raw": ic50_raw,
        "ic50_censored": ic50_censored,
        "ln_ic50_censored": float(math.log(max(ic50_censored, 1e-15))),
        "beyond_tested_range": beyond_range,
        "min_dose": min_d,
        "max_dose": max_d,
        "n_points": int(d.size),
        "n_doses": int(uniq_ld.size),
        "auc": auc_fit,
        "auc_empirical": auc_emp,
        "aac": float(1.0 - auc_fit),
        "r2": r2,
    }


def write_artifact(
    path: str | Path,
    tool: str,
    title: str,
    *,
    skill: str = "pgx-workflow",
    inputs=(),
    parameters: dict | None = None,
    shape: tuple[int, int] | dict | None = None,
    scientific_context: dict | None = None,
) -> str:
    p = Path(path)
    ext = p.suffix.lower()
    figures = {".png": "image/png", ".pdf": "application/pdf", ".svg": "image/svg+xml"}
    tables = {".tsv": "text/tab-separated-values", ".csv": "text/csv"}
    if ext in figures:
        kind, media = "figure", figures[ext]
    elif ext in tables:
        kind, media = "table", tables[ext]
    else:
        raise ValueError(f"Unsupported artifact extension: {ext}")

    descriptor: dict = {
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
        block: dict = {"format": ext[1:]}
        if ext == ".png":
            with open(p, "rb") as fh:
                header = fh.read(24)
            if header[:8] == b"\x89PNG\r\n\x1a\n":
                w, h = struct.unpack(">II", header[16:24])
                block["widthPx"], block["heightPx"] = int(w), int(h)
        descriptor["figure"] = block
    elif shape is not None:
        if isinstance(shape, dict):
            rows_n = int(shape.get("rows", 0))
            cols_n = int(shape.get("columns", shape.get("cols", 0)))
        else:
            rows_n, cols_n = int(shape[0]), int(shape[1])
        descriptor["table"] = {"rows": rows_n, "columns": cols_n}

    if scientific_context:
        descriptor["scientificContext"] = scientific_context

    cwd = Path.cwd().resolve()
    prov_inputs = []
    for ip in inputs:
        if not ip:
            continue
        ip_p = Path(ip)
        try:
            rel = ip_p.resolve().relative_to(cwd)
        except ValueError:
            continue
        entry = {"path": str(rel)}
        if ip_p.is_file():
            entry["sha256"] = hashlib.sha256(ip_p.read_bytes()).hexdigest()
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
