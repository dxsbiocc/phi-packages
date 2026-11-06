#!/usr/bin/env python3
"""CLI: init library.sqlite schema."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from db import init_db


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--db_path", required=True)
    args = p.parse_args()
    result = init_db(args.db_path)
    result["mode"] = "live"
    print(json.dumps(result, ensure_ascii=False))
    return 0 if result.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
