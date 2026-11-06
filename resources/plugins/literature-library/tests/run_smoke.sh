#!/usr/bin/env bash
# Smoke: init demo collection, ingest bib + ids, mock search, query, export.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEMO="$ROOT/demo-project"
FIX="$ROOT/tests/fixtures"
OUT="$ROOT/tests/out"
PY="${PYTHON:-python3}"
mkdir -p "$OUT"
rm -rf "$DEMO"
mkdir -p "$DEMO"

echo "== unit test_dedup =="
"$PY" "$ROOT/tests/test_dedup.py"

echo "== create collection =="
COL=$("$PY" "$ROOT/skills/literature-library/scripts/collection.py" \
  --action create --project_root "$DEMO" --name "Demo HCC" --description "smoke")
echo "$COL" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o["ok"] and o["mode"]=="live", o; print(o["collection_id"])' > "$OUT/collection_id.txt"
CID=$(cat "$OUT/collection_id.txt")
echo "collection_id=$CID"

echo "== ingest bib =="
ING=$("$PY" "$ROOT/skills/literature-library/scripts/ingest.py" \
  --project_root "$DEMO" --collection_id "$CID" \
  --input_path "$FIX/sample.bib" --format bib)
echo "$ING" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o["ok"], o; c=o["counts"]; assert c["inserted"]>=1 and c["merged"]>=1, o; print("counts", c)'

echo "== unit: test_ids_coidentity =="
"$PY" "$ROOT/tests/test_ids_coidentity.py" -q
echo "== ingest id_list =="
ING2=$("$PY" "$ROOT/skills/literature-library/scripts/ingest.py" \
  --project_root "$DEMO" --collection_id "$CID" \
  --input_path "$FIX/id_list.txt" --format ids)
echo "$ING2" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o["ok"], o; assert o.get("n_parsed")==2, o; print("id_list", o["counts"], "n_parsed", o.get("n_parsed"))'

echo "== ingest ris =="
ING3=$("$PY" "$ROOT/skills/literature-library/scripts/ingest.py" \
  --project_root "$DEMO" --collection_id "$CID" \
  --input_path "$FIX/sample.ris" --format ris)
echo "$ING3" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o["ok"] and o["counts"]["inserted"]>=1, o'

echo "== run-config + mock search =="
RC=$("$PY" "$ROOT/skills/literature-library/scripts/run_config.py" \
  --project_root "$DEMO" --collection_id "$CID" \
  --query "T cell exhaustion liver" --sources pubmed,europepmc --limit 5)
RID=$(echo "$RC" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o["ok"], o; print(o["run_id"])')
SR=$("$PY" "$ROOT/skills/literature-library/scripts/search.py" \
  --project_root "$DEMO" --collection_id "$CID" --run_id "$RID" --mock)
echo "$SR" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o["ok"] and o["mode"]=="mock" and o["n_hits"]==2, o; assert o.get("ingest",{}).get("ok"), o; print("mock hits", o["n_hits"])'

echo "== query =="
Q=$("$PY" "$ROOT/skills/literature-library/scripts/query.py" \
  --project_root "$DEMO" --collection_id "$CID" --q "T Cells" --limit 20)
echo "$Q" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o["ok"] and o["n"]>=1, o; print("query n", o["n"], "backend", o["backend"])'

echo "== export =="
EX=$("$PY" "$ROOT/skills/literature-library/scripts/export.py" \
  --project_root "$DEMO" --collection_id "$CID" \
  --out_path "$DEMO/exports/search_log.jsonl")
echo "$EX" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o["ok"] and o["n"]>=1, o; print("exported", o["n"], o["out_path"])'

echo "== expand preview (mock, no write) =="
ST=$("$PY" "$ROOT/skills/literature-library/scripts/expand.py" \
  --project_root "$DEMO" --collection_id "$CID" --seed_id "10.1016/j.cell.2017.05.035" --mock)
echo "$ST" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o["ok"] and o.get("stub") is False and o.get("wrote") is False, o; print("expand preview", o.get("n_hits"))'

echo "== merge_ids dry-run =="
MR=$("$PY" "$ROOT/skills/literature-library/scripts/merge_ids_cli.py" --project_root "$DEMO")
echo "$MR" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o["ok"] and o["mode"]=="dry-run", o; print("actionable", o.get("n_actionable"))'

echo "SMOKE PASSED — demo DB at $DEMO/library.sqlite"
