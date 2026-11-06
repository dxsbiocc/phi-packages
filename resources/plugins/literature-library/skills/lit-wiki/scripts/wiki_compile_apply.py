#!/usr/bin/env python3
"""Apply agent-produced compile JSON into the wiki vault (validate + write)."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main() -> int:
    rest = list(sys.argv[1:])
    # Map --response / positional convenience
    if "--llm-json" not in rest and "--apply-compile" not in rest:
        # look for --response PATH
        if "--response" in rest:
            i = rest.index("--response")
            if i + 1 < len(rest):
                path = rest[i + 1]
                rest = rest[:i] + rest[i + 2 :]
                rest = ["--llm-json", path, *rest]
        else:
            print('{"ok": false, "error": "provide --llm-json / --apply-compile / --response PATH"}')
            return 1
    if not any(a == "--compile" or a.startswith("--compile=") for a in rest):
        rest = ["--compile", "agent", *rest]
    return subprocess.call([sys.executable, str(HERE / "wiki_ingest.py"), *rest])


if __name__ == "__main__":
    sys.exit(main())
