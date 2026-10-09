# BioMCP connector for Phi

**English** | [中文](README.zh-CN.md)

This `phi-packages` connector runs the official [BioMCP](https://github.com/genomoncology/biomcp)
native executable in its own Phi-managed environment.

## Versions and platforms

- Connector package: **1.1.0**; requires Phi **1.0.1** or newer.
- Upstream package: **biomcp-cli 0.9.1**, distributed as native wheels on [PyPI](https://pypi.org/project/biomcp-cli/0.9.1/).
- Platforms: macOS Apple Silicon, macOS Intel, and Linux x86-64. Linux requires glibc 2.28 or newer.
- The [official website icon](../../../docs/content-icons.md#biomcp-website-icon) stays with the package.

## Installation and startup

Install or update BioMCP from Phi's **Connectors** catalog. Phi verifies the package,
then prepares its declared native environment before enabling the connector. The
installation downloads approximately 15–17 MB for the selected platform, verifies
the pinned size and SHA-256, and extracts only the declared native executable.
Installation progress and failures belong to Phi's managed installation flow; a
failed preparation leaves the connector unavailable until installation is retried.

`environment.yml` is the sole artifact-pin document. It declares the native backend,
the executable `biomcp`, and exact official wheel URLs, sizes, hashes, and ZIP members
for each supported platform. `locks/<platform>.txt` contains `@EXPLICIT` with zero
Conda dependencies. The binary is installed under Phi's managed runtime prefix;
the package launches it directly with `biomcp serve`.

This package has no Python launcher, Python runtime dependency, `pip install`,
Conda package installation, Git clone, or first-start software download. Once
prepared, startup reuses the installed native executable. Online biomedical
queries still require their upstream services. Older 1.0.x connector launcher
caches are left untouched and are no longer used by this package.

The 1.1.0 source must be rebuilt into the [signed catalog](../../../docs/publishing.md)
before the app can install it. The earlier standalone
[1.0.0 release](https://github.com/dxsbiocc/phi-packages/releases/tag/biomcp-v1.0.0)
contains the retired launcher and does not provide this managed installation flow.

## Credentials and provenance

The connector exposes BioMCP's read-only biomedical tools. Many queries need no
key; credentialed upstream services still require their own API keys. See the
[upstream API-key guide](https://biomcp.org/getting-started/api-keys/).

The wheel version, URLs, SHA-256 values, and sizes match the former 1.0.x pins.
Each verified wheel contains `biomcp_cli-0.9.1.data/scripts/biomcp`; only that member
is installed. `SOURCE.json` records the official repository, PyPI distribution,
version, and icon provenance. The upstream MIT copyright and license remain in
[LICENSE.upstream](LICENSE.upstream). Data-source terms remain separate; see
[upstream source licensing](https://biomcp.org/reference/source-licensing/).

## Verification

With Phi and `phi-packages` side by side, use Phi's existing parser and builder:

```sh
cd ../Phi
node --import ./scripts/test-loader.mjs --test ../phi-packages/tests/biomcp-native.test.mjs
node --import ./scripts/test-loader.mjs --test tests/envs-application-native.test.ts tests/envs-application-artifacts.test.ts tests/envs-ensure.test.ts
```

The content suite verifies the native declarations, original pins, every platform
lock, environment ownership and identity, license, and archive allowlist. Shared
Phi tests cover download integrity, cache reuse, concurrency, bounded ZIP-member
extraction, cancellation, and runtime readiness. The former launcher tests moved
to these managed-runtime responsibilities instead of retaining a second installer.
To check actual previously verified wheels without network access, set
`PHI_BIOMCP_WHEEL_DIR` to a disposable directory containing the three official
wheel filenames; the content suite then exercises Phi's native installer on them.
