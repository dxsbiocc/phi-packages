# Notes — literature-library 0.6.0

Implements the five deferred items from UX-NOTES-0.5.3.

| # | Item | Landed |
|--|--|--|
| 1 | Background path-watch auto-rebind | `skills/lit-ingest/scripts/watch_rebind.py` (stdlib poll; `--once` / `--interval`; dry-run / `--confirm`; reuses `rebind_pub_files`) |
| 2 | PDF binary text + optional OCR | `pdf_extract.py` + ingest path: **pdftotext** default; tesseract / `LIT_PDF_OCR_CMD` soft-fail; TextIn **not** required / not KB |
| 3 | Campus / SSO authorized full-text | `skills/lit-institutional/` — `--institutional_access` gate; methods `vpn_proxy` \| `session_reuse` \| `saml_browser` \| `user_provided_file`; `LIT_INSTITUTIONAL_PROXY`; credentials stay in user vault |
| 4 | Optional MCP client + reev consumer | `providers_mcp.py`; `--backend mcp` / `LIT_SEARCH_BACKEND`; `--mcp_mock_fixture`; HTTP remains default. research-evidence `literature-search` can `--delegate_lit` / `--hits_jsonl` |
| 5 | Packaging | version **0.6.0** landed in phi-packages (first package landing; `-draft` dropped); schema **v6**; `demo-*` trees excluded |

## Schema

`SCHEMA_VERSION = 6`: `pub_files.access_origin`, `pub_files.access_method`.

## Tests

`tests/test_deferred_v060.py` (+ existing suite). Offline only.

## Still open (honest)

- Live SAML/CAS Playwright = host/agent concern (plugin records method only)
- Production `LIT_MCP_CMD` bridge must be provided by Phi host (not shipped)
- Real campus proxy soak-test (not in CI)
- Host `validatePlugin` is not in this repository; local contract check covers the manifest
