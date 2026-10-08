# Plugin contract

contractVersion: 1.1.0

A plugin is a **composite capability**: one or more specialist agents, the skills
they use, and the environments they run in, installed and removed as one unit.
Single tools are not plugins; they ship as standalone skills, wrappers, or
connectors. This contract defines the plugin directory, its `phi-package.yaml`,
and what the engine does when a plugin is installed, upgraded, or removed. It
builds on the [agent](agent.md), [skill](skill.md), and
[environment](environment.md) contracts. See
[runtime foundation](../design/phi-runtime-foundation.md) §7 and §10, and
[content distribution design](../design/phi-content-distribution-design.md) §5.

## 1. Layout

```text
<plugin-id>/
  phi-package.yaml          # required (§ 2)
  README.md                 # optional, shown in the plugins page
  agents/<Name>.md          # agent definitions (agent contract)
  skills/<name>/            # skills (skill contract), validated as "inside a plugin"
  environments/<name>/      # environment specs and locks (environment contract)
    environment.yml
    locks/<platform>.txt
  assets/                   # optional data the plugin's skills read
```

`mcp/`, `wrappers/`, and `orchestrator/` are **reserved**: a v1 plugin must not
contain them. Nothing outside the listed components is executed. Files a plugin
needs at run time live under `skills/<name>/` (scripts, references, assets) or
`assets/`.

## 2. `phi-package.yaml`

UTF-8 YAML. Unknown top-level keys are errors. Since 1.1.0 the shared distribution fields of the [package contract](package.md) § 1 (`minAppVersion`, `requires`, `dependsOn`, `files`) are also accepted, with the meaning that contract gives them.

| Field           | Required | Rule                                                                                                                                                                                                                                |
| --------------- | -------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `schemaVersion` | yes      | `1`                                                                                                                                                                                                                                 |
| `id`            | yes      | `^[a-z][a-z0-9-]{1,63}$`; unique among installed plugins; equals the directory name of a bundled plugin                                                                                                                             |
| `type`          | yes      | `plugin`                                                                                                                                                                                                                            |
| `version`       | yes      | semantic version `MAJOR.MINOR.PATCH` (optional `-prerelease`)                                                                                                                                                                       |
| `title`         | yes      | 1–80 characters, shown to the user                                                                                                                                                                                                  |
| `summary`       | yes      | 1–300 characters                                                                                                                                                                                                                    |
| `toolPrefix`    | yes      | skill contract pattern `^[a-z][a-z0-9]{1,11}$`, not an engine-reserved prefix; the prefix of every script tool the plugin's skills declare. Unique among installed plugins and standalone skills                                    |
| `components`    | yes      | `{ agents: [<path>…], skills: [<path>…] }`: paths relative to the plugin root, `agents/<Name>.md` and `skills/<name>`; at least one agent or skill; every listed path exists; no file under `agents/` or `skills/` is left unlisted |
| `environments`  | no       | map `<name>` → `{ spec: environments/<name>/environment.yml }`. `<name>` uses the environment name pattern. The spec's `name` equals `<name>`, and `environments/<name>/locks/<platform>.txt` exists for every platform Phi ships   |

Example:

```yaml
schemaVersion: 1
id: visualization
type: plugin
version: 1.0.2
title: 科研绘图
summary: 用模板生成可直接发表的组学图表，配有专属绘图智能体。
toolPrefix: viz
components:
  agents: [agents/Visualization.md]
  skills: [skills/omics-visualization]
```

The bundled visualization plugin deliberately omits `environments`: its agent
and `omics-visualization` skill both declare `phi:r@1`. The shared `phi-r`
environment contains R, the visualization packages, and Python 3.12 for the
plugin's standard-library-only command-line scripts. Other plugins may still
declare private environments through the optional field above.

## 3. Rules across components

1. **Environments are private.** `plugin:<name>` in a component resolves only to an
   environment declared by **the same plugin**. A reference to another plugin's
   environment, or to an undeclared name, is an error at validation and at run time.
2. **Agents** follow the agent contract. Their names are unique across all agents
   Phi loads; a clash with an installed agent makes the plugin uninstallable until
   resolved. `visibility: internal` stays reserved (agent contract § 2.2).
3. **Skills** are validated with "inside a plugin" (skill contract § 2.2): they must
   not declare `toolPrefix`; their script tools are named `<plugin toolPrefix>_<name>`.
   `attachTo` may name `main` and the plugin's own agents only. Skill names are
   unique across all loaded skills.
4. A plugin's agents may list its skills in `skills:`; they may also use standalone
   skills that are installed separately.

## 4. Lifecycle

Installed plugins live in `~/.phi/packages/plugin/<id>/<version>/`, one active
version per id, recorded in `~/.phi/packages/plugins.json`. The engine owns every
step below; a plugin contains no install code.

1. **Validate** the directory against this contract and every component against its
   own contract; collect all problems (`errors` and `warnings`, like the other
   validators). Errors stop the install.
2. **Install** copies the plugin's allowlisted files (§ 1; never caches, run
   leftovers, or `__pycache__`) into the version directory, then registers its
   agents, skills, and script tools. Environments are not built at install: they
   are built on first use with the usual prompt and progress, or from the
   environment panel. Each environment is recorded as referenced by the plugin.
3. **Enable / disable** keeps the files and hides or restores the plugin's agents,
   skills, and tools.
4. **Upgrade** installs the new version beside the old one, builds the new
   version's environments whose envIds changed **before** switching, then switches
   the active version atomically and drops the old version's references; old
   environments are collected once unreferenced. A failed build leaves the old
   version active.
5. **Uninstall** unregisters the components, removes the version directories and
   the plugin's environment references; garbage collection removes environments no
   longer referenced.

**Bundled plugins** ship with the app under `resources/plugins/<id>/` and are
installed through the same loader: on start, a bundled plugin that is not
installed, or whose version is newer than the installed one, is installed or
upgraded. A user who uninstalled a bundled plugin keeps it uninstalled.

For the visualization upgrade that removes its former private environment, the
atomic switch drops the installed version's `plugin:viz` reference. The normal
environment GC then removes that old environment once it has no other
referrers; there is no visualization-specific migration or deletion path.

## 5. Validation

`validatePlugin(dir)` checks everything above and returns all problems; a plugin is
valid when `errors` is empty. It runs in `bun run lint` for every plugin in the
repository and before every install or upgrade.

## 6. Versioning

`contractVersion` follows the content distribution design §4.4: minor versions are
additive only (new optional fields, allowing a reserved directory or value); anything
else needs a decision record and a deprecation window.

## Changes

- **1.1.0** (2026-10-02): accepts the package contract's shared distribution fields (`minAppVersion`, `requires`, `dependsOn`, `files`). Additive.
