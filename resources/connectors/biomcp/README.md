# BioMCP connector for Phi

**English** | [中文](README.zh-CN.md)

A Phi MCP connector for [BioMCP](https://github.com/genomoncology/biomcp), maintained
as an adapter in `phi-packages`. It uses Phi's existing package and stdio mechanisms;
no Phi application source changes are required.

## Versions and supported platforms

- Connector package: **1.0.0**.
- Upstream BioMCP: **0.9.1**, pinned in `upstream.json`.
- macOS Apple Silicon, macOS Intel, and Linux x86-64 are supported. The Linux binary
  requires glibc 2.28 or newer, matching its upstream wheel tag.
- The connector references Phi's managed `phi:python@1` environment. Its definitions
  and locks remain in Phi; this package adds no Python library dependencies.

## Install without changing Phi

Download [biomcp-registry-1.0.0.zip](https://github.com/dxsbiocc/phi-packages/releases/download/biomcp-v1.0.0/biomcp-registry-1.0.0.zip)
from the [connector release](https://github.com/dxsbiocc/phi-packages/releases/tag/biomcp-v1.0.0),
extract it, open Phi's **Skills** catalog, click the folder-icon **Add** button, and
select the extracted directory. This existing picker registers sources for all
package types. A source containing only this connector has no skill rows; that is
expected. Close the dialog, reopen the **Connectors** catalog, and install BioMCP.
The directory contains `index.json` and `mcp-biomcp-1.0.0.tar.gz`.
The current connector dialog itself has no local-directory picker. If Phi reports that its Python environment is not
ready, build the referenced environment using Phi's existing environment controls.

The local registry is unsigned and is treated as an imported source. This release
uses the existing local-directory installation route; Phi's automatic remote registry
fetching is not yet connected.

## Startup and caching

The launcher uses Python's standard library to select the official platform wheel
from PyPI, check its pinned size and SHA-256, and extract only the native `biomcp`
executable. It does not run `pip install` or modify a managed environment.

The executable and verified wheel are stored outside the package, in:

- macOS: `~/Library/Caches/Phi/connectors/biomcp/<version>/<platform>/`.
- Linux: `~/.cache/Phi/connectors/biomcp/<version>/<platform>/`.

The Linux default intentionally stays independent of Phi's managed XDG cache so
prewarming works across execution environments. A file lock serializes concurrent starts. Each start verifies the cached wheel and
executable; corrupted binaries are repaired from the verified wheel. The launcher
then replaces itself with `biomcp serve`, preserving the stdio MCP stream. Bootstrap
messages go only to stderr.

The first start downloads approximately 15–17 MB. For a slow connection or offline
use, warm the cache first with the same user account that runs Phi:

```bash
python /path/to/biomcp/server.py --prepare
```

After a Phi package install, the script is under
`~/.phi/packages/mcp/biomcp/1.0.0/server.py`. A dedicated cache can be selected with
`--cache-dir /path/to/biomcp-cache` when testing or preparing an offline environment.
The initial download requires internet access; online biomedical queries continue
to require access to their upstream data sources.

## Capabilities and credentials

This package enables BioMCP's standard read-only biomedical MCP surface, including
literature, gene, variant, drug and trial retrieval. Many upstream queries work
without keys. BioMCP's optional credentialed sources still require their own keys;
this initial package does not add credential fields to Phi's stdio package contract
or put credentials in source files. See the [upstream API-key guide](https://biomcp.org/getting-started/api-keys/).

## Provenance and verification

`upstream.json` records the exact official `biomcp-cli` wheel URLs, sizes and SHA-256
checksums. BioMCP is MIT-licensed; its copyright and license notice are retained in
[LICENSE.upstream](LICENSE.upstream). This notice applies to upstream BioMCP and does
not establish a repository-wide license for the adapter. Data-source terms remain
separate: see [upstream source licensing](https://biomcp.org/reference/source-licensing/).

Run the adapter's offline tests with:

```bash
python3 -B -m unittest discover -s tests -p test_biomcp_connector.py
```

The macOS ARM64 native binary and launcher are also checked with actual MCP
initialization and tool discovery. No real biomedical query or credentialed endpoint
is required for that protocol check.
