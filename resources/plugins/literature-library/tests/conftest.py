"""Make the tool scripts and their internal modules importable for every test."""
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / "skills" / "manage-library" / "scripts"
for path in (SCRIPTS / "lib", SCRIPTS):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))
