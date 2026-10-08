# Initial import validation

Validated with the Phi package builder and installer validators from the source working tree on 2026-10-08.

- Imported source and reference files: 7885; all copied files matched the captured SHA-256 values and sizes.
- High-confidence credential-pattern scan: no findings in the exported files.
- Package build: 206 packages (14 skills, 18 MCP connectors, 1 plugin, 173 wrapper packages).
- Archive staging checks: every generated package passed size, SHA-256, safe extraction, per-file allowlist/hash and manifest validation.
- Authored repository metadata passed whitespace checks. Existing upstream content is preserved byte-for-byte, including TSV empty columns and test-snapshot whitespace.

## Existing build diagnostics

The vendored RNA-seq workflow has four `plugin/nf-schema` include references reported as `not-found` by the local include scanner. They are Nextflow plugin references rather than relative source files; the existing workflow's plugin configuration is preserved. `.nf-core.yml` is reported as unattributed support metadata. These diagnostics did not prevent package generation or archive validation.

## Verification scope

No managed environments were built, no Nextflow workflows or plotting scripts were executed, and no remote registry was signed or connected to Phi. Generated package archives were used for local validation and are not committed in this source repository.
