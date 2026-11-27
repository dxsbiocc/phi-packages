"""Centralized path resolution keeping code and mutable state/data strictly separated.

Code files (`schema_v1.sql`, `pipeline/sync_policy.json`, `scanners/*.py`) always resolve
relative to `CODE_ROOT` (`data-pipelines/omics-sniffer/`).

Mutable runtime state (`data/`, `serving/`, `logs/`) resolves via environment variables so
production deployments (e.g. `phi-site` on `/var/lib/phi-omics-sniffer`) never write into
the git checkout:
- `PHI_OMICS_SNIFFER_STATE_DIR` / `OMICS_SNIFFER_STATE_DIR` / `OMICS_SNIFFER_ROOT`:
  base state directory containing `data/`, `serving/`, and `logs/` (defaults to `CODE_ROOT`
  for standalone local development).
- `OMICS_SNIFFER_DATA_DIR`: explicit override for SQLite staging/main `data/` directory.
- `PHI_OMICS_SNIFFER_DATA_DIR` / `OMICS_SNIFFER_SERVING_DIR`: explicit override for the
  Parquet `serving/` directory (`releases/v_*` + atomic `current` symlink) shared with
  `phi-omics-gateway`.
- `OMICS_SNIFFER_LOGS_DIR`: explicit override for `logs/` directory.
"""

from __future__ import annotations

import os
from pathlib import Path

CODE_ROOT = Path(__file__).resolve().parent.parent
SCHEMA_PATH = CODE_ROOT / "schema_v1.sql"
POLICY_PATH = CODE_ROOT / "pipeline" / "sync_policy.json"


def state_root() -> Path:
    for var in ("PHI_OMICS_SNIFFER_STATE_DIR", "OMICS_SNIFFER_STATE_DIR", "OMICS_SNIFFER_ROOT"):
        val = os.environ.get(var, "").strip()
        if val:
            return Path(val).resolve()
    return CODE_ROOT


def data_dir() -> Path:
    val = os.environ.get("OMICS_SNIFFER_DATA_DIR", "").strip()
    return Path(val).resolve() if val else (state_root() / "data")


def serving_dir() -> Path:
    for var in ("OMICS_SNIFFER_SERVING_DIR", "PHI_OMICS_SNIFFER_DATA_DIR"):
        val = os.environ.get(var, "").strip()
        if val:
            return Path(val).resolve()
    return state_root() / "serving"


def logs_dir() -> Path:
    val = os.environ.get("OMICS_SNIFFER_LOGS_DIR", "").strip()
    return Path(val).resolve() if val else (state_root() / "logs")
