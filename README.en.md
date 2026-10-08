# Phi Packages

[中文](README.md) | **English**

Phi's independent public content repository, containing source files for installable domain capabilities.

## Current contents

| Content | Count |
| --- | ---: |
| Standalone skills | 14 |
| MCP connector definitions | 19 |
| Composite Phi plugins | 1 |
| Skills within plugins | 1 |
| Wrapper adapters | 652 |

```text
resources/
  skills/                Standalone skills, references, and scripts
  connectors/            MCP connector declarations
  plugins/               Visualization plugin and its components
  wrappers/              Nextflow modules, subworkflows, workflows, and adapters
docs/
  contracts/             Reference snapshots of Phi content contracts
  design/                Content distribution, runtime, and wrapper designs
  decisions/             Recorded distribution decisions
  roadmap/               Reference content distribution implementation plans
SOURCE.json              Import provenance, version, scope, and exclusions
```

Each installable package has its own version. A repository commit or release does not require every package to be upgraded together.

## Relationship to the Phi application

Application repository: [dxsbiocc/phi](https://github.com/dxsbiocc/phi).

The initial import was a snapshot of Git-tracked working-tree files. See [SOURCE.json](SOURCE.json) for the source commit and included uncommitted content. Phi's installable domain content now comes from this repository's signed online catalog by default.

The official source uses the [catalog-v1 release](https://github.com/dxsbiocc/phi-packages/releases/tag/catalog-v1). Catalog browsing fetches the signed index, icons, and connector metadata; installation downloads the selected archives and necessary dependencies into managed paths under `~/.phi`. Offline browsing uses a verified cache. See [publishing instructions](docs/publishing.md) for the script, signing procedure, and installed paths. Core runtime code, the installer, validators, and application UI remain in the Phi repository.

Phi supplies the official shared environments and application-level palettes; their core definitions are not copied here. Resources required by a plugin remain with that plugin. Contracts and designs under `docs/` are reference snapshots, with Phi remaining the authoritative source. Runtime binaries, installed environments, account credentials, and run outputs were not imported.

## Distribution boundaries

This repository contains skills, connectors, plugins, wrappers, and the resources those components need. The application-level `resources/palettes/` and core `resources/runtime/` trees are outside its distribution scope.

The old local image directory for the differential-expression wrappers has been removed. See [differential-expression dependencies](resources/wrappers/modules/local/differential-expression/README.md) for the dependency migration.

Run the boundary checks with `node --test tests/content-boundaries.test.mjs`.

## Import exclusions and third-party content

- `resources/plugins/office/` was not imported: the existing Office skill licenses explicitly prohibit copying and third-party distribution. The entire Office plugin was excluded to preserve component integrity. It can be included after obtaining the necessary authorization or replacing the restricted components.
- `resources/skills/create-wrapper/` remains in Phi as a built-in engine authoring skill.
- Imported files retain their original licenses, citations, authors, and upstream provenance. This repository includes upstream content such as nf-core; importing it does not change its existing license terms.
- This repository currently has no single repository-wide open-source license.

See [SOURCE.json](SOURCE.json) for the full scope and [VALIDATION.md](VALIDATION.md) for validation results.

## Installable connector

[BioMCP](resources/connectors/biomcp/README.md) provides a standalone
[local registry distribution](https://github.com/dxsbiocc/phi-packages/releases/tag/biomcp-v1.0.0)
that can be added through Phi's existing directory installation flow without application source changes.

## Content conventions

- Follow the existing directory layout and content formats where possible. Package manifests follow the [Package contract](docs/contracts/package.md).
- Keep wrapper include targets, support files, and test data together so copying an adapter alone does not leave missing dependencies.
- Check component references when modifying plugins. Local MCP services continue to run in managed environments on the user's machine.
- Do not commit secrets, account configuration, installed environments, or Nextflow run outputs.
- Default package-local `icon.svg`, `icon.png`, `icon.webp`, `icon.jpg`, or `icon.jpeg` files travel with content; see [content icons](docs/content-icons.md) for the convention and existing brand attribution.
