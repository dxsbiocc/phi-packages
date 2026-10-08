# BioMCP connector validation

Added the standalone `mcp:biomcp@1.0.0` adapter on 2026-10-08 without changes to Phi application source or the user's Phi configuration.

- Phi's canonical package validator accepted the manifest with zero errors and warnings.
- All 23 offline launcher tests passed, including pinned download checks, safe extraction, concurrent startup, offline reuse, executable repair and stdout discipline.
- The official macOS ARM64 BioMCP 0.9.1 wheel passed its pinned size and SHA-256 checks.
- The native program, adapter launcher and installed managed stdio entry completed actual MCP initialization and advertised seven tools.
- The existing Phi installer installed the package and generated an enabled managed stdio entry in an isolated temporary agent directory. The Python execution environment in this check used fixture metadata and an existing Python executable; no micromamba build was performed.
- Existing distribution boundary tests remained green (11 tests).
- The distribution ZIP contains one local registry and the connector archive; the package itself contains no native binaries or copied core environment files.

## Verification limits

macOS ARM64 is the runtime-tested platform. Intel macOS and Linux x86-64 downloads are pinned but their native programs were not executed here. No upstream biomedical queries or credentialed API calls were made. The local registry is unsigned/imported; automatic remote registry support is unchanged.

## Earlier content verification

# Content boundary correction validation

Corrected the public source repository on 2026-10-08:

- Removed the application palette tree, core runtime environment tree and local-only differential-expression image recipes.
- Preserved all 652 callable wrapper adapters, including the eight differential-expression wrappers and the RNA-seq QC include.
- Boundary and process-dependency tests: 11 passed. JavaScript syntax and changed-file whitespace checks passed.
- Package builder and installer staging checks: all 206 packages passed (14 skills, 18 MCP connectors, 1 plugin, 173 wrapper packages).
- Thirty-two R script, parameter and wrapper-entrypoint files remained byte-identical to the previous source.
- Public OCI and native SIF dependency builds succeeded. The final layered OCI image build loaded DESeq2, edgeR, limma, BiocParallel, variancePartition, ggplot2, pheatmap, ashr and lme4 in its runtime stage.
- The original local runtime reports ashr 2.2.63 and lme4 2.0.6; the new image resolves the same versions, now also pinned for local Conda execution.
- OCI manifest digest and Linux AMD64 platform were verified. Maximum compressed layer is approximately 205 MB.
- A focused independent source/dependency review found no material issue.

## Verification limits

No complete Nextflow analysis was run with real inputs, no native SIF execution was tested,
and no native Conda environments were built on other platforms. The local OCI pull/check was
cancelled after Docker's credentials helper blocked; library loading was verified in the
successful public image build instead. The Phi application repository remains unchanged.

## Historical initial import checks

The following record describes the original broader import before this boundary correction.

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
