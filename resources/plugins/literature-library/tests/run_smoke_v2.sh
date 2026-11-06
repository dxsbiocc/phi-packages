#!/usr/bin/env bash
# V2/V4 smoke: wiki vault init, mock LLM compile (no API key), query, suggest/apply tree, lint.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEMO="$ROOT/demo-project-v2"
FIX="$ROOT/tests/fixtures"
OUT="$ROOT/tests/out-v2"
PY="${PYTHON:-python3}"
mkdir -p "$OUT"
rm -rf "$DEMO"
mkdir -p "$DEMO"

echo "== create collection (default wiki-init) =="
COL=$("$PY" "$ROOT/skills/lit-collection/scripts/collection.py" \
  --action create --project_root "$DEMO" --name "Wiki Demo" --description "v2 smoke")
echo "$COL" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o["ok"] and o.get("wiki_vault"), o; print(o["collection_id"]); open("'"$OUT"'/vault.txt","w").write(o["wiki_vault"])'
CID=$(echo "$COL" | "$PY" -c 'import json,sys; print(json.load(sys.stdin)["collection_id"])')
VAULT=$(cat "$OUT/vault.txt")
echo "collection_id=$CID"
echo "vault=$VAULT"
test -f "$VAULT/AGENTS.md"
test -f "$VAULT/wiki/index.md"
test -f "$VAULT/wiki/log.md"

echo "== unit: test_wiki_llm =="
"$PY" "$ROOT/tests/test_wiki_llm.py"

echo "== unit: compile_extract heuristics (fallback path) =="
"$PY" - << PY
import sys
sys.path.insert(0, "$ROOT/skills/lit-wiki-lib/scripts")
from wiki_lib import compile_extract, condense_summary
title = "Landscape of Infiltrating T Cells in Liver Cancer Revealed by Single-Cell Sequencing"
abstract = (
    "Here we report a single-cell RNA sequencing analysis of T cells from HCC patients. "
    "We identified exhausted CD8+ T cells and Tregs preferentially enriched in HCC, "
    "and LAYN is linked to immunosuppressive T cells. Using TraCeR we assembled TCR sequences."
)
ex = compile_extract(title, abstract)
assert "placeholder" not in ex["summary"].lower(), ex["summary"]
assert ex["claims"], ex
assert any("scRNA" in m or "single-cell" in m.lower() or "TraCeR" in m for m in ex["methods"]), ex["methods"]
assert any("LAYN" == e or e.upper()=="LAYN" for e in ex["entities"]), ex["entities"]
assert any("exhaust" in c.lower() or "hepatocellular" in c.lower() or "T cell" in c for c in ex["concepts"]), ex["concepts"]
assert ex["tags"], ex["tags"]
noisy = compile_extract(
    "PD-L1 and IFN-γ in HCC",
    "PD-L1 blockade restored IFN-γ production; NKG2A and HLA-E were upregulated.",
)
elow = [e.lower() for e in noisy["entities"]]
assert "l1" not in elow and "ifn" not in elow, noisy["entities"]
assert any(e in ("PD-L1", "IFNG", "NKG2A", "HLA-E") for e in noisy["entities"]), noisy["entities"]
s = condense_summary(title, abstract)
assert "Here we report" in s
print("compile_extract OK", {"n_claims": len(ex["claims"]), "methods": ex["methods"][:4], "entities": ex["entities"][:6], "tags": ex["tags"][:6]})
PY

echo "== complementary DOI+PMID merge (shells then joint) =="
"$PY" "$ROOT/tests/test_merge_smoke.py" --project_root "$DEMO" --collection_id "$CID"

echo "== ingest bib into sqlite =="
ING=$("$PY" "$ROOT/skills/lit-ingest/scripts/ingest.py" \
  --project_root "$DEMO" --collection_id "$CID" \
  --input_path "$FIX/sample.bib" --format bib)
echo "$ING" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o["ok"], o; print(o["counts"])'

echo "== resolve Zheng publication id =="
PID=$("$PY" "$ROOT/skills/lit-query/scripts/query.py" \
  --project_root "$DEMO" --collection_id "$CID" --q "Landscape" --limit 5 \
  | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o["n"]>=1, o; print(o["publications"][0]["id"])')
echo "publication_id=$PID"

echo "== wiki-ingest (mock LLM compile; no API key) =="
WI=$("$PY" "$ROOT/skills/lit-wiki-ingest/scripts/wiki_ingest.py" \
  --project_root "$DEMO" --collection_id "$CID" \
  --publication_id "$PID" --compile mock)
echo "$WI" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o["ok"] and o["source_page"].startswith("wiki/sources/"), o; assert o.get("n_claims",0)>=1, o; assert o.get("concept_pages_created") or o.get("entity_pages_created"), o; assert o.get("tags_persisted") is True, o; assert "placeholder" not in (o.get("summary_preview") or "").lower(), o; assert o.get("compile_backend")=="mock" or str(o.get("compile","")).startswith("mock"), o; print(o["safe_id"], "claims", o["n_claims"], "compile", o.get("compile"), "concepts+", len(o.get("concept_pages_created") or []), "entities+", len(o.get("entity_pages_created") or []), "tags", o.get("effective_tags"))'
SAFE=$(echo "$WI" | "$PY" -c 'import json,sys; print(json.load(sys.stdin)["safe_id"])')
SRC="$VAULT/wiki/sources/${SAFE}.md"
test -f "$SRC"
test -f "$VAULT/raw/publications/${SAFE}/meta.json"
grep -q "## Summary" "$SRC"
! grep -q "Summary (placeholder)" "$SRC"
! grep -q "_Stub summary generated" "$SRC"
grep -q "## Key claims" "$SRC"
grep -qE "mock-v0\.[34]" "$SRC"
# concepts or entities dirs should have pages
COMP_N=$(find "$VAULT/wiki/concepts" -name '*.md' | wc -l)
ENT_N=$(find "$VAULT/wiki/entities" -name '*.md' | wc -l)
echo "concept_pages=$COMP_N entity_pages=$ENT_N"
test "$COMP_N" -ge 1
test "$ENT_N" -ge 1
! test -f "$VAULT/wiki/entities/l1.md"
! test -f "$VAULT/wiki/entities/ifn.md"
test -f "$VAULT/wiki/entities/layn.md"

echo "== agent compile request → apply (fixture; no API key) =="
set +e
REQ_OUT=$("$PY" "$ROOT/skills/lit-wiki-ingest/scripts/wiki_compile_request.py" \
  --project_root "$DEMO" --collection_id "$CID" --publication_id "$PID" --ensure_vault)
REQ_RC=$?
set -e
# exit 2 = waiting_for_agent (expected)
test "$REQ_RC" -eq 2 -o "$REQ_RC" -eq 0
echo "$REQ_OUT" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o.get("waiting_for_agent") and o.get("request_path"), o; print("request", o["request_path"])'
REQ_PATH=$(echo "$REQ_OUT" | "$PY" -c 'import json,sys; print(json.load(sys.stdin)["request_path"])')
RESP_PATH=$(echo "$REQ_OUT" | "$PY" -c 'import json,sys; print(json.load(sys.stdin)["response_path"])')
"$PY" - << PY
import json
from pathlib import Path
req = json.loads(Path("$REQ_PATH").read_text())
assert req["publication_id"] == "$PID"
# agent-filled fixture (same shape the host agent would write)
resp = {
  "publication_id": "$PID",
  "summary": "Agent-model synthesis of Zheng 2017: a scRNA-seq atlas of infiltrating T cells in HCC highlighting exhaustion programs and LAYN on immunosuppressive subsets.",
  "key_claims": [
    "Deep scRNA-seq profiled thousands of T cells from HCC patients.",
    "Exhausted CD8+ T cells and Tregs are enriched in tumor.",
    "LAYN (layilin) is upregulated on immunosuppressive T cells."
  ],
  "methods": ["scRNA-seq"],
  "concepts": [{"name": "T cell exhaustion"}, {"name": "hepatocellular carcinoma"}],
  "entities": [{"symbol": "LAYN", "name": "layilin", "kind": "gene"}, {"symbol": "CD8", "kind": "protein"}],
  "tags": ["HCC", "exhaustion", "LAYN", "scRNA-seq"],
  "confidence_notes": "smoke agent fixture"
}
Path("$RESP_PATH").write_text(json.dumps(resp, ensure_ascii=False, indent=2), encoding="utf-8")
print("wrote", "$RESP_PATH")
PY
APPLY=$("$PY" "$ROOT/skills/lit-wiki-ingest/scripts/wiki_compile_apply.py" \
  --project_root "$DEMO" --collection_id "$CID" --publication_id "$PID" \
  --response "$RESP_PATH")
echo "$APPLY" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o.get("ok"), o; assert o.get("compile_backend")=="agent" or str(o.get("compile","")).startswith("agent"), o; print("applied", o.get("compile"), "backend", o.get("compile_backend"))'

echo "== tags persisted in sqlite =="
"$PY" - << PY
import sqlite3, json
db = sqlite3.connect("$DEMO/library.sqlite")
row = db.execute("SELECT tags_json FROM publications WHERE id=?", ("$PID",)).fetchone()
tags = json.loads(row[0] or "[]")
assert tags, row
print("sqlite tags", tags[:8])
db.close()
PY

echo "== wiki-query exhaustion / LAYN =="
WQ=$("$PY" "$ROOT/skills/lit-wiki-query/scripts/wiki_query.py" \
  --project_root "$DEMO" --collection_id "$CID" --q "exhaustion" --limit 20)
echo "$WQ" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o["ok"] and o["embeddings"] is False, o; assert o["n_wiki"]>=1 or o["n_publications"]>=1, o; print("exhaustion wiki", o["n_wiki"], "pubs", o["n_publications"])'
WQ2=$("$PY" "$ROOT/skills/lit-wiki-query/scripts/wiki_query.py" \
  --project_root "$DEMO" --collection_id "$CID" --q "LAYN" --limit 20)
echo "$WQ2" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o["ok"], o; assert o["n_wiki"]>=1 or o["n_publications"]>=1, o; print("LAYN wiki", o["n_wiki"], "pubs", o["n_publications"])'

echo "== suggest-tree (no manual tags — should use auto tags) =="
SUG=$("$PY" "$ROOT/skills/lit-wiki-suggest-tree/scripts/wiki_suggest_tree.py" \
  --project_root "$DEMO" --collection_id "$CID" \
  --out_path "$OUT/proposals.json" --max_branches 8)
echo "$SUG" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o["ok"] and o["wrote"] is False and o["n_proposals"]>=1, o; print("proposals", o["n_proposals"])'
test -f "$OUT/proposals.json"

echo "== apply-tree dry-run (no --confirm) =="
DRY=$("$PY" "$ROOT/skills/lit-wiki-apply-tree/scripts/wiki_apply_tree.py" \
  --project_root "$DEMO" --collection_id "$CID" --proposals_path "$OUT/proposals.json")
echo "$DRY" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o["ok"] and o["wrote"] is False, o'

echo "== apply-tree --confirm =="
APP=$("$PY" "$ROOT/skills/lit-wiki-apply-tree/scripts/wiki_apply_tree.py" \
  --project_root "$DEMO" --collection_id "$CID" --proposals_path "$OUT/proposals.json" --confirm)
echo "$APP" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o["ok"] and o["wrote"] is True and o["n_applied"]>=1, o; print("applied", o["n_applied"])'
grep -q "### " "$VAULT/wiki/index.md"

echo "== wiki-lint =="
LINT=$("$PY" "$ROOT/skills/lit-wiki-lint/scripts/wiki_lint.py" \
  --project_root "$DEMO" --collection_id "$CID")
echo "$LINT" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o["ok"], o; assert o["n_errors"]==0, o; ph=[i for i in o["issues"] if i.get("code")=="placeholder_summary"]; assert not ph, ph; print("issues", o["n_issues"], "errors", o["n_errors"], "sources", o["n_sources"])'

echo "== idempotent wiki-init =="
WI2=$("$PY" "$ROOT/skills/lit-wiki-init/scripts/wiki_init.py" \
  --project_root "$DEMO" --collection_id "$CID")
echo "$WI2" | "$PY" -c 'import json,sys; o=json.load(sys.stdin); assert o["ok"], o'

echo "SMOKE V2 PASSED — vault at $VAULT"
