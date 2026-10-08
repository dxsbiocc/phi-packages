# Connector contract

contractVersion: 1.0.0

A **connector** gives Phi an MCP server: a remote HTTP server (most catalog entries) or
a local stdio server that runs in a managed environment. Connectors are distributed as
packages (`type: mcp`, added here as package contract 1.1.0), so the catalog can grow
without an app release. Public REST databases are API skills, not connectors (content
distribution design §8.3). See the [package](package.md),
[environment](environment.md), and [enablement](enablement.md) contracts.

## 1. Package layout

```text
<connector-id>/
  phi-package.yaml          # type: mcp, with the `connector` block (§ 2)
  icon.svg | icon.png       # optional, shown in the catalog
  environment.yml           # stdio only, optional (§ 3)
  locks/<platform>.txt      # with environment.yml
  <server files>            # stdio only: what `command` runs, when not provided by the environment
```

## 2. The `connector` block

Shared package fields (package contract § 1) plus:

| Field                   | Required   | Rule                                                                                                                                                                      |
| ----------------------- | ---------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `connector.transport`   | yes        | `http` or `stdio`                                                                                                                                                         |
| `connector.publisher`   | yes        | 1–80 characters                                                                                                                                                           |
| `connector.category`    | yes        | one of the catalog categories: `生产力`, `沟通协作`, `设计创作`, `健康与生命科学`, `科研数据`                                                                             |
| `connector.homepage`    | no         | `https://` URL                                                                                                                                                            |
| `connector.url`         | http only  | `https://` URL of the MCP endpoint                                                                                                                                        |
| `connector.auth`        | http only  | `none`, `oauth` (the server's own sign-in flow, as today), or `header` (the user enters a token once; stored like other provider secrets)                                 |
| `connector.environment` | stdio only | `phi:<name>@<major>` or `./environment.yml` (environment contract; a package-local environment uses scope `mcp` and the package id as owner — environment contract 1.3.0) |
| `connector.command`     | stdio only | the program to run, resolved **inside** the environment, or `./<path>` relative to the package                                                                            |
| `connector.args`        | stdio only | list of strings; `${package}` expands to the installed package directory                                                                                                  |
| `connector.secrets`     | no         | **reserved** for stdio servers that need credentials; v1 rejects it                                                                                                       |

Unknown keys under `connector` are errors.

## 3. Install and run

- **HTTP**: installing writes the server into `mcp.json` exactly as adding it from the
  catalog does today (name = package id; `type: http`; `url`); sign-in follows `auth`.
  The `mcp.json` entry records the package id so updates and uninstall find it.
  Uninstall removes that entry; a server the user added by hand is never touched.
- **stdio**: installing places the package under `~/.phi/packages/mcp/<id>/<version>/`
  and writes an `mcp.json` entry built by the managed stdio mechanism: the command is
  resolved inside the environment and started with exactly the execution contract's
  variables, so host variables never reach the server. The environment is built on
  first use with the usual prompt; until then the connector shows as "环境未构建". When
  the environment's envId changes, the entry is regenerated.
- Enablement (`mcp:<package id>`) controls whether the server's `mcp.json` entry is
  enabled. Defaults: catalog connectors the user added are enabled; nothing is added
  without the user's action.

## 4. Catalog

The connector catalog shown to the user is the list of `type: mcp` packages from the
registries Phi knows: the bundled registry built from `resources/connectors/` and any
local or remote registry the user added. An entry already present in `mcp.json` shows
as added. Entries that need an app newer than this one show but cannot be added.

## 5. Versioning

`contractVersion` follows the content distribution design §4.4: minor versions are
additive only (for example allowing `secrets`); anything else needs a decision record
and a deprecation window.
