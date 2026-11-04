#!/usr/bin/env python3
"""Emit a grounded wiki compile request for the invoking agent (no API keys)."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main() -> int:
    # Delegate to wiki_ingest --request-only
    args = [sys.executable, str(HERE / "wiki_ingest.py"), "--request-only", *sys.argv[1:]]
    # Ensure --compile agent if not provided
    if not any(a == "--compile" or a.startswith("--compile=") for a in sys.argv[1:]):
        args[2:2] = ["--compile", "agent"]
    # Fix: insert after script name
    cmd = [sys.executable, str(HERE / "wiki_ingest.py")]
    rest = list(sys.argv[1:])
    if "--request-only" not in rest:
        rest = ["--request-only", *rest]
    if not any(a == "--compile" or a.startswith("--compile=") for a in rest):
        rest = ["--compile", "agent", *rest]
    return subprocess.call([sys.executable, str(HERE / "wiki_ingest.py"), *rest])


if __name__ == "__main__":
    sys.exit(main())
