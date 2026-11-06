# UX notes — literature-library 0.5.3-draft

Problem → fix mapping from real acceptance friction (Zheng / live v052 run).

| # | Problem | Fix |
|--|--|--|
| P0.1 | Concept merge wrote `[nsclc](nsclc.md)` while page was `non-small-cell-lung-cancer.md` → lint `broken_link` | `resolve_related_concepts` (exact / slugify / aliases / light fuzzy); only link existing `wiki/concepts/*.md`; unresolved → `related_unresolved` in apply JSON |
| P0.2 | Agents/users typed `--query` but CLI only had `--q` | `lit-query` + `lit-wiki-query` accept `--query` alias; SKILL one-liners show `--query` |
| P0.3 | Local PDF `file_path` drifts after moves; identity is hash not path | `skills/literature-library/scripts/rebind_files.py` scan_dir → match `content_hash` → update path; dry-run default, `--confirm` writes |
| P1.4 | Hard to see “what wiki work is left” | `wiki_agent_status.py`: missing enriched sources + pending merges + suggested commands; `--emit_requests` batches compile_request.json |
| P1.5 | Gold recall 2/7 on search-only looked like a silent failure | coverage `hints[]` when recall &lt; 1.0 / shells / uncontacted / OOW leak — honest, no fake recall |
| P1.6 | README buried agent-compile; demos leaned on mock | Top 「最短真实工作流」; explicit: mock = offline CI only |
| P2.7 | Agent entities with “implied / not claimed” noise | Soft-strip in `validate_compile_payload` + `entity_warnings` |
| P2.8 | Wiki query AND semantics undocumented | argparse help + SKILL note |

## Deferred (not in 0.5.3) → landed in 0.6.0

See `NOTES-0.6.0.md`. This package is the first landing in phi-packages. Remaining: live SSO soak and a host-provided `LIT_MCP_CMD`.

## Tests

`tests/test_ux_v053.py` — related resolve, `--query` alias, rebind dry+confirm, coverage hints, agent_status basic.
