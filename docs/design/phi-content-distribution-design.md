# Phi Content Distribution: Installable Units, Plugins, and Environments

Date: 2026-09-29
Status: **Design draft — needs owner sign-off before implementation.** Phases 0–5
are proposed for the internal beta; phases 6–9 are post-beta (see §15).

Chinese version: [phi-content-distribution-design.zh-CN.md](phi-content-distribution-design.zh-CN.md).
Decision record: [content-distribution.md](../decisions/content-distribution.md).

## 1. Summary

Phi's domain capabilities (skills, wrappers, database connectors, agent
definitions, visualization templates) are shipped inside the app bundle under
`resources/` and are all loaded for every user. This design turns them into
**content that is distributed and loaded on demand**:

- The app keeps a small **engine** (kernel): agent runtime, generic core
  tools, installer, environment manager, approvals, remote execution, and a
  handful of built-in skills. The engine contains no domain-specific logic.
- The engine exposes a fixed set of **versioned contracts** for plugins,
  skills, wrappers, connectors, agent definitions, environments, and core
  services (§4.4). Contracts are validated with reference implementations,
  then frozen; afterwards new capability arrives as content, not as engine
  changes.
- Everything domain-specific becomes an installable **unit** (skill, wrapper,
  MCP connector) or, for composite capabilities only, a **plugin**.
- Units and plugins share one package format, one installer, and one registry.
  The registry is local (generated from `resources/`) first and a signed remote
  registry later, so content can be updated without shipping a new app.
- Script-bearing content runs in **managed, locked micromamba environments**
  chosen automatically by Phi.
- The `db_*` connector toolchain is replaced by **API skills plus one core
  fetch tool**, based on the evaluation in Appendix A.

Goals, in order: (1) content is not locked into the installer, (2) the main
agent's routing load scales with what the user installed, not with everything
Phi knows, (3) scripts run reproducibly without the user fixing environments.

## 2. Non-goals

- A public marketplace with accounts, ratings, payments, or third-party
  publishing. The first remote registry is first-party only.
- Loading arbitrary executable TypeScript from content packages. Content is
  declarative plus sandboxable scripts and MCP servers (§11.2).
- Windows support for managed environments (bioconda is effectively
  unavailable there).
- Rewriting how Nextflow wrappers execute. Nextflow keeps owning per-process
  environments (`conda` / `docker` / `singularity` profiles).

## 3. Current state

| Content | Location | How it reaches the agent | Problem |
|---|---|---|---|
| 21 skills | `resources/skills/` | whole directory appended to the runtime skill paths (`appendExistingBundledSkillPaths`) | every skill description is in the main prompt; bundle carries 65M (57M is `omics-visualization`) |
| ~150 module / subworkflow wrappers | `resources/wrappers/` | `wrapper_agent` + `wrapper_search`; overlay packs under `~/.phi/wrappers/packs/<version>/` replace the whole tree | low routing cost, but one monolithic pack; install copies whole directories, including dev leftovers |
| 32 db connectors | `resources/db-connectors/` | Database agent with 7 `db_*` tools over ~20k lines in `src/main/agent/db/` | slower and less reliable than plain fetch (Appendix A) |
| 3 agent definitions | `resources/agents/` | scanned at session creation | Visualization is coupled to bundled tools and skill |
| pi plugins | npm / git via runtime | extensions, prompts, themes | developer extensions, not research capabilities |
| Remote MCP connectors | `~/.phi/mcp.json` via catalog dialog | MCP tools (deferrable) | already the "add what you need" model this design generalises |

Script dependencies are unmanaged: `viz_render` calls `Rscript` / `python3`
from `PATH`, and skill scripts run with whatever Python the agent's `bash`
finds.

## 4. Concepts

### 4.1 Engine, units, plugins

- **Engine** — code that ships with the app and is versioned with it. Includes
  generic core tools (`wrapper_*`, the core fetch tool, `skill_run`, artifact
  presentation, file and shell tools), the installer, the environment manager,
  and **built-in skills**
  that extend Phi itself (`create-wrapper`). Built-in skills are not packages
  and never appear in catalogs.
- **Unit** — one piece of content of one type: a skill, a wrapper tool family,
  or an MCP connector. Units are what the agent actually uses at run time.
- **Plugin** — a package that bundles several units of different types that
  must be installed, upgraded, and removed together, optionally with its own
  specialist agent and its own environment.

Plugins answer *how capability is acquired*; units answer *what the agent can
use*. The two layers coexist: official and shared content is distributed, user
and project content is authored directly as units.

### 4.2 Unit or plugin?

A package is a plugin only if at least one holds:

1. it bundles two or more component types that depend on each other
   (agent + skill + templates, skill + wrappers);
2. it needs a dedicated environment that conflicts with `phi-python` (§9);
3. it needs a dedicated specialist agent (own workflow, own instructions).

Everything else is a unit. A single tool never becomes a plugin. Having
scripts alone does not make a skill a plugin (§9.4). Today only
`omics-visualization` (with the Visualization agent and its templates)
qualifies.

### 4.3 Where content goes at run time

Contents are routed by type, regardless of whether they came from a plugin or
were added standalone:

| Content | Routed to |
|---|---|
| specialist agent | main agent, as one delegation entry |
| wrapper | `wrapper_agent`'s search pool |
| API / knowledge skill | the specialist it declares (e.g. Database), else the main agent |
| script skill needing a dedicated env | its plugin's agent (env is bound to that session) |
| MCP connector | its plugin's agent, else the main agent (deferred tool loading) |

Consequences: the main agent sees specialist entries plus a small number of
knowledge skills, so its prompt grows with *installed specialists*, not with
tool count. Outside an orchestration run, specialists never call each other;
only the main agent orchestrates, and hand-off between specialists is by
project files plus a short description. Inside a plugin's orchestration run
(§11.4), the plugin's orchestrator coordinates the plugin's own agents under
engine-enforced limits.

Names are namespaced by package: `visualization/omics-visualization`,
`gatk4/haplotypecaller`.

### 4.4 Kernel contracts and stability

The engine is fixed first. Everything content can rely on is a written
contract with a JSON schema, a validator, and conformance tests; content may
use nothing else.

| Contract | Defines | Starting point |
|---|---|---|
| Package | `phi-package.yaml`, file layout, `files.json`, versioning, `dependsOn`, `requires` (§5) | new |
| Skill | `SKILL.md` frontmatter (name, description, `attachTo` specialist, environment, declared scripts), `references/` / `scripts/` / `assets/` layout, `skill_run` semantics | runtime skill format + Phi additions |
| Wrapper | `wrapper.yaml`, `params.json`, `main.nf` entrypoint, outputs, profiles, family packaging and `dependsOn` | existing wrapper technical design |
| Connector | MCP (standard protocol) + Phi metadata: transport, auth kind, category, environment for stdio servers | existing MCP catalog |
| Agent definition | frontmatter (`name`, `description`, `tools`, `spawns`, `model`, `thinkingLevel` as in omp; plus Phi's `visibility`, `skills`, `delegationMode`, `delegation`, `fallback`, environment binding), structured result via `outputSchema` | existing `resources/agents` format aligned with omp agent definitions |
| Plugin | components, environment bindings, entry agent, `requires.coreTools` (§11) | new |
| Environment | `environment.yml` + per-platform lock, variables injected at run time (§9.4) | new |
| Core services | the engine tools content may call: fetch, `skill_run`, `wrapper_*`, artifact presentation, file / shell / approval behaviour | partly existing |
| Orchestration | multi-agent run model on top of omp's `task` / registry / `hub`: internal agents, generated `spawns`, run-scoped store, budgets, checkpoints, optional orchestrator component (§11.4) | extension points reserved in v1; the contract is added as a minor, additive version and frozen with its reference implementation (phase 9) |
| Artifact | typed outputs content hands to the UI: `table`, `figure`, `structure`, `molecule`, `network`, `report`, with a JSON sidecar describing the file | generalises today's db result viewers |

Rules:

1. **No domain logic in the engine.** Domain behaviour lives in content. The
   current `viz_*` tools and the `db_*` toolchain are the two existing
   violations; they move out (§10, §11.3). The UI renders typed artifacts; it
   never needs code from a plugin.
2. **Validate before freezing.** Each contract gets one reference
   implementation before v1: an API skill (`protein-apis`), a wrapper family
   (`samtools`), an MCP connector (one remote catalog entry plus one stdio
   server), and the visualization plugin. Gaps found here change the draft,
   not v1. v1 also reserves the orchestration extension points (plugin
   `orchestrator` slot, agent `visibility`, the `orchestration.*` service
   namespace) so that adding multi-agent orchestration later is additive.
3. **Contract versions.** Every contract carries `contractVersion` (semver)
   and packages declare the versions they target. Minor versions are
   additive only (new optional fields, new artifact kinds, new core services);
   unknown optional fields are ignored.
4. **Breaking changes are exceptional.** A major version needs an ADR, a
   deprecation window of at least two app releases in which both versions are
   accepted, and a migration note or automatic shim.
5. **One validator for everyone.** The same validator (`phi validate`) runs in
   the package builder, in CI for the registry, and on user-authored content,
   so authors see the same errors Phi enforces.
6. **Conformance suite.** Each contract has fixture packages that must load,
   route, and run; the suite is part of `npm test` and gates engine changes.

## 5. Package model

### 5.1 Manifest

Every package has `phi-package.yaml` at its root:

```yaml
schemaVersion: 1
id: protein-apis                 # unique within the registry, [a-z0-9-]
type: skill                      # skill | wrapper | mcp | plugin
version: 1.2.0                   # semver
title: 蛋白质数据库 API
summary: UniProt, PDB, AlphaFold and InterPro lookups
minAppVersion: 0.9.0
requires:
  coreTools: [fetch, skill_run]  # engine capabilities this package relies on
dependsOn:                       # other packages, resolved at install time
  - id: samtools
    version: ">=1.0.0"
environment:                     # optional, see §9
  spec: environment/environment.yml
  lock: environment/conda-lock.yml
files: files.json                # generated allowlist with sha256 per file
```

A plugin manifest additionally lists its components:

```yaml
type: plugin
components:
  agents: [agents/Visualization.md]   # visibility: entry | internal (frontmatter)
  skills: [skills/omics-visualization]
  wrappers: []
  mcp: []
  orchestrator: null                  # reserved: workflow.yaml or a script (§11.4)
```

Visualization deliberately declares no private environment: its agent and
skill frontmatter both bind to the shared official `phi:r@1` environment.

### 5.2 Allowlisted contents

Packages are built from an allowlist (git-tracked files, or the existing
`index.json` digest list for wrappers), never by copying a directory. The
builder rejects `.nextflow/`, `.nextflow.log*`, `work/`, `results/`,
`.DS_Store`, `__pycache__`, and anything untracked. CI fails a package that
exceeds a per-type size budget unless the manifest carries an explicit
justification. This replaces `copyWrapperSourceTree`'s whole-directory copy.

### 5.3 Versioning and compatibility

- Semver per package. `minAppVersion` and `requires.coreTools` are checked
  before install; an incompatible package is shown but not installable.
- Core tools carry an API version; removing or changing a tool signature is a
  breaking engine change and bumps that version.

## 6. Registry and installer

### 6.1 Registry

A registry is an `index.json` listing, per package: `id`, `type`, `version`,
`summary`, `sha256`, `size`, `url`, `dependsOn`, `minAppVersion`, and catalog
metadata (category, preview image). Two sources implement the same interface:

- **Local registry** — generated at build time from `resources/`. Used in
  phases 2–6; lets every install/upgrade path be tested without network.
- **Remote registry** — phase 7. Static files on object storage + CDN; the
  index is signed with an ed25519 key held in CI; the app embeds the public
  key and rejects unsigned or mis-signed indexes. Mirrors are allowed only if
  they serve content signed with the same key.

### 6.2 Install flow

1. Resolve dependencies; show the plan (packages, sizes, environment cost).
2. Download (remote) or read (local) each archive; verify `sha256`.
3. Extract into a temporary directory under `~/.phi/.staging/`, rejecting
   absolute paths, `..` segments, and symlinks leaving the root.
4. Verify every file against `files.json`.
5. Atomically rename into place; write `.source.json`
   (`{ registry, id, version, sha256, installedAt, installedBy: user|dependency }`).
6. Build or reuse the environment if declared (§9), possibly deferred until
   first use with explicit progress.

Uninstall removes a package only when no installed package depends on it;
dependency-only packages are garbage-collected when their last dependent goes.
Upgrade installs the new version side by side and switches atomically.

### 6.3 On-disk layout

```text
~/.phi/
  packages/<type>/<id>/<version>/  # installed packages; type is skill, wrapper, mcp, or plugin
  wrappers/tree/                   # assembled wrapper tree (§8.2)
  skills/<name>/                   # user-authored skills (not packages)
  wrappers/custom/<id>/            # user-authored wrappers
  runtime/                         # micromamba root: envs/, pkgs/, mambarc (runtime foundation §2)
  state/enabled.json               # enablement, global + per-project overrides
  .staging/
<project>/.phi/skills/...          # project-scoped units
```

## 7. Loading, enablement, and catalog UI

- One **content loader** merges four sources — built-in, installed packages
  (units and plugin components), user-authored, project — and records
  provenance for every item.
- Installation is global; enablement is global with per-project overrides
  (a single-cell project can enable only the relevant capabilities).
- Only enabled items reach the runtime. The main prompt may carry one short
  line listing *installable but not installed* capabilities so the agent can
  suggest "install X" instead of assuming a capability does not exist.
- UI:
  - **Skills / Wrappers / Connectors pages** list every unit with its source
    (`plugin: visualization`, `mine`, `project`, `imported`), allow enable /
    disable, and offer an "Add from catalog" dialog modelled on
    `McpConnectorCatalogDialog`. Items that belong to a plugin can be disabled
    but not uninstalled individually.
  - **Plugins page** lists composite capabilities only, with the example
    outputs and environment cost on the card.
  - The existing pi plugin page is renamed **Developer extensions** and moved
    to advanced settings. The word "plugin" in the UI means Phi plugins only.
- A user can later "package these units as a plugin" to share them; the
  two-layer model is what makes that path possible.

## 8. Unit types

### 8.1 Skills

- `SKILL.md` + optional `references/`, `scripts/`, `assets/`.
- Script-bearing skills declare their environment (§9) and document scripts
  as `skill_run` calls. A third-party skill that still says `python x.py`
  keeps working through `bash`, without the environment guarantee.
- Built-in (not packaged): `create-wrapper`.

### 8.2 Wrappers

- Package granularity is the **tool family** as it exists under
  `modules/<provider>/<tool>/` (e.g. `samtools`, `gatk4`, `bcftools`); each
  subworkflow is its own package that `dependsOn` the modules it includes.
  Scenario bundles ("variant calling") are meta-packages with only
  dependencies.
- Subworkflows include modules by relative path
  (`../../../../modules/nf-core/...`). Installed wrapper packages are
  therefore assembled into one tree at `~/.phi/wrappers/tree/` that mirrors
  the source layout, rather than rewriting include paths.
- This revises `composition/packs.ts`: the tree is merged per package and
  verified per package (each package's `files.json`) instead of being
  replaced and verified as one pack. The overlay-pack mechanism is retired
  once per-package install ships.
- `wrapper_search` only returns wrappers from installed, enabled packages.

### 8.3 MCP connectors

- Remote HTTP connectors keep the current catalog flow; the catalog becomes
  registry-backed so entries can be added without an app release.
- Local stdio connectors are packages with an environment (§9) and a launch
  command that the environment manager resolves.
- MCP is used for sources that need identity or credentials (lab LIMS,
  institutional or commercial databases), for protocols that are genuinely
  hard to call directly, and for existing official servers (PubMed, bioRxiv,
  ClinicalTrials already in the catalog). Public REST databases are API skills
  (§10), not MCP servers.

## 9. Execution environments

> **Superseded in detail by [phi-runtime-foundation.md](phi-runtime-foundation.md)**
> (runtime layout, explicit locks, read-only prefixes, activation snapshot,
> `runInEnvironment`, binding rules). Where this chapter differs, the
> foundation document wins.

### 9.1 Manager

- Phi bundles **micromamba** (single static binary, BSD license). No user
  conda installation is required or used.
- Environments live at `~/.phi/runtime/envs/<envId>/`, where `envId` is derived from the
  solved lock. Identical requirements share one environment; micromamba
  installs from its package cache with hard links, so each extra environment
  usually costs tens to hundreds of MB, not a full copy.
- Environments are **immutable**. Changed requirements produce a new
  environment; the old one is collected when no package references it.
- Channel mirrors (tuna / ustc) and proxies are configurable. Doctor compares
  an environment with its lock and offers rebuild on drift.

### 9.2 Shared base environments

- `phi-python`: the locked scientific Python stack used by official skills
  (numpy, pandas, matplotlib, scikit-learn, scanpy, scvelo, rdkit, openpyxl,
  defusedxml, pillow, markitdown, …), derived from the actual imports of the
  current skills. Maintained by Phi developers, versioned, locked per
  platform (osx-arm64, osx-64, linux-64).
- `phi-r`: optional, built on demand; shared by R notebooks, scanpy R
  interoperability, and the visualization plugin.

### 9.3 Resolving an environment for a package

- **Official packages** must be satisfied by the current `phi-python` / `phi-r` or
  bring their own lock (plugins). CI runs
  `micromamba install -n phi-python --file environment.yml --dry-run --json` and
  fails the package if the result is not "no changes". Clients never solve
  environments for official skills.
- **User-added skills** are resolved on add:
  1. no environment file → `phi-python`, flagged "dependencies undeclared";
  2. dry-run against `phi-python` reports no changes → use it;
  3. dry-run against each existing managed environment → reuse the first
     with no changes;
  4. solve `phi-python` spec + the skill's spec together → new environment;
  5. if that conflicts, solve the skill's spec alone → new environment.
  Steps 4–5 show the size and time estimate and need user confirmation.
- Environments cannot be layered: Python `site-packages` and R libraries do
  not stack across conda environments, so "base + extras" is always a new,
  fully solved environment.
- `pip:` sections are allowed for user skills but are outside lock
  guarantees; official packages use conda-forge / bioconda only. Pixi may be
  evaluated later if PyPI-heavy packages become common.

### 9.4 Using the environment

| Path | Mechanism |
|---|---|
| skill scripts | `skill_run(skill, script, args)` resolves the skill's environment and runs the script by absolute path |
| plugin agent's `bash` | the session is created with the plugin environment on `PATH`, `CONDA_PREFIX` set, `PYTHONNOUSERSITE=1`, `R_LIBS_USER` cleared |
| core services that execute content (`skill_run`, `wrapper_*`) | receive a run context with the environment's absolute interpreter paths instead of looking up `PATH` |
| stdio MCP servers | launched with the environment's interpreter |
| wrappers | unchanged; Nextflow manages per-process environments |

Package environments are read-only. Extra packages requested during a project
go into a project environment (itself a fully solved environment, §9.3).
Run records (`wrappers/reproducibility.ts` and equivalents) store the
environment hash.

### 9.5 Remote and offline hosts

Phase 8. The same lock is rebuilt on the remote host with micromamba; for
offline clusters the environment is packed (`conda-pack`) locally or on a login
node and uploaded through the existing remote bundle path. linux-64 locks
must stay compatible with glibc 2.17 where such clusters are targeted.

## 10. Data access redesign

### 10.1 Evidence

Appendix A: on 14 structured lookups × 2 repetitions, plain URL fetching was
as accurate as the `db_*` toolchain (27/28 and 26/28 vs 25/28) at roughly half
the latency and a third of the tool calls; `db_query` failed on 40% of calls.
The toolchain's failures were connector coverage gaps; fetch's failures were
upstream errors without retry.

### 10.2 Target shape

- **Core fetch tool (engine)** — HTTP GET/POST with retry and backoff,
  per-host rate limits, domain allowlist and egress audit (reusing
  `policy-*.ts`), pagination and bulk download (from `db_download`), optional
  credential injection from the credential store, and registration of results
  with the UI result viewers (`DbQueryResultPreview` and viewers).
- **API skills (units)** — about six domain skills replacing 32 connectors:
  `protein-apis`, `genomics-apis`, `chemistry-apis`, `pathway-network-apis`,
  `clinical-cancer-apis`, `ontology-apis`. `SKILL.md` is short; per-database
  detail lives in `references/<db>.md` (endpoints, query syntax, fields,
  examples, pitfalls). Existing connector manifests are the source material
  for these references.
- **Helper scripts** (standard library only, run via `skill_run`) for GraphQL
  (gnomAD), SPARQL (UniProt, WikiPathways), and Entrez batch/rate-limited
  access.
- **Database agent** stays as the single main-agent entry; its tools become
  the core fetch tool, `skill_run`, and the API skills.
- **MCP** for authenticated or private sources (§8.3).

### 10.3 Retirement gate

The `db_*` tools and the connector manifest format are retired only after an `api-skill` arm is at least as good as the
`db` arm on the extended evaluation (§14), including batch, pagination,
GraphQL, SPARQL, and rate-limited tasks. Until then both paths coexist.

## 11. Plugins

> Plugin directory layout and how plugin agents use environments are defined
> in [phi-runtime-foundation.md](phi-runtime-foundation.md) §6–§7.

### 11.1 Shape

A plugin is a package of `type: plugin` whose components are installed under
`~/.phi/packages/plugin/<id>/<version>/` and routed by type (§4.3). A plugin has at most one
**entry** agent, which is what the main agent sees; it may ship any number of
**internal** agents, which are visible only to the plugin's orchestrator
(§11.4). A private environment is declared in the manifest; component bindings
live in agent or skill frontmatter. Visualization instead binds to `phi:r@1`.

### 11.2 Code boundary

Plugins do not ship TypeScript that runs inside Phi. New programmatic
capability comes from (a) core tools referenced in `requires.coreTools`,
(b) scripts run in the component's resolved managed environment, or (c) MCP servers. This keeps
review, security, and compatibility tractable.

### 11.3 First plugin: visualization

`agents/Visualization.md` + `skills/omics-visualization` (scripts, templates,
previews), both bound to the shared `phi:r@1` environment. Today's `viz_*` tools
(`src/main/agent/visualization/`, ~1.1k lines) are domain logic in the engine;
they become the command-line program `scripts/viz.py` declared as script tools
(runtime foundation §5.1), keeping their names, and figures are returned as
`figure` artifacts. This happens in step 4 of the implementation plan, before
plugins exist; step 6 only packages it. The visualization plugin is the
reference implementation for the Plugin, official Environment reference, script-tool, and
Artifact contracts.

Candidate later plugins: single-cell analysis (agent + scanpy/scvi skills +
starsolo wrappers + torch environment), bulk RNA-seq (agent + skills + wrapper
meta-package).

### 11.4 Multi-agent orchestration (reserved)

Some capabilities are not one specialist but a team working over many rounds,
for example an AI Co-Scientist–style system: a supervisor drives generation,
reflection, ranking (pairwise debates in a tournament), evolution, proximity,
and meta-review agents for hours, keeping a growing set of hypotheses and
scores. Phi should run such systems, and users should be able to define their
own, **without engine changes per system**.

#### 11.4.1 What the omp runtime already provides

Phi runs on `@oh-my-pi/pi-coding-agent` (18.1.10), which already ships a
multi-agent runtime. Phi should build on it rather than re-implement it:

| omp capability | Where | What it gives Phi |
|---|---|---|
| `task` tool | `src/task/` | spawn named agents; batch mode runs many tasks in parallel with a shared `context`; per-spawn `outputSchema` with `strict` / `permissive` validation for structured results; optional worktree isolation |
| Agent definitions | `src/task/agents.ts`, discovery | Markdown + frontmatter (`name`, `description`, `tools`, `spawns`, `model`, `thinkingLevel`, `blocking`); sources `bundled` / `user` / `project` |
| Spawn policy | `task/spawn-policy.ts` | `spawns` frontmatter restricts which agents an agent may spawn; `task.maxRecursionDepth` (default 2) bounds nesting |
| Agent registry and lifecycle | `src/registry/` | process-global agents with stable ids; `running` / `idle` / `parked` / `aborted`; parked agents keep their transcript and are **revived on demand** |
| `hub` tool + IRC bus | `tools/hub/`, `irc/bus.ts` | agent-to-agent messaging (`send`, `wait`, `inbox`, `list`) that wakes idle agents or revives parked ones; background jobs (`start`, `ps`, `logs`, `stop`) |
| Limits | `task/executor.ts`, `task/provider-concurrency.ts` | soft per-agent request budgets with wrap-up steering; per-provider concurrency semaphores |
| Long-running loops | `src/goals/`, `src/autoresearch/` | goal mode with persisted state and token usage; an autoresearch loop with dashboard and resume |

Phi currently bypasses most of this: specialists are created with
`restrictToolNames` and Phi's own `agents/registry.ts` re-implements parallel,
background, steer, and stop. Phase 1 must decide whether Phi specialists move
onto omp's `task` / registry (proposed, §11.4.4).

#### 11.4.2 Phi owns the contract, omp is the implementation

The kernel-stability rule (§4.4) still holds: plugins depend only on Phi
contracts, never on omp internals, because omp moves fast. The Orchestration
and Agent-definition contracts are designed as a **Phi-owned subset that maps
thinly onto omp**, and one adapter layer in the engine does the mapping:

- Phi agent frontmatter keeps omp's field names where they exist (`name`,
  `description`, `tools`, `spawns`, `model`, `thinkingLevel`) and adds Phi
  fields (`visibility`, `skills`, `environment`, `delegation*`).
- Structured results use `outputSchema` (JSON Schema) instead of a
  Phi-specific parsing protocol where possible.
- The adapter is covered by the conformance suite; an omp upgrade that breaks
  it fails CI before release. omp stays pinned.

Phi adds what omp does not have, as engine services:

- **Budgets as hard limits** per orchestration run (tokens, cost, wall time,
  total agent runs), on top of omp's soft request budgets.
- **Run-scoped typed store**: collections with JSON schemas declared in the
  manifest (hypotheses, reviews, matches), append-only history, browsable in
  the UI. Agents access it through a store tool, not ad-hoc files.
- **Approvals, project boundaries, remote-project guards** for every spawned
  agent (existing Phi extensions).
- **Human checkpoints** surfaced in the UI (`ask` / `checkpoint`).
- **Run view**: agent tree, messages, store, budget, with pause / steer / stop,
  and main-agent wake-up on completion (existing `run-host.ts`).
- **Isolation by package**: a plugin's agents may spawn only agents from the
  same plugin plus public specialists listed in `uses`; enforced by
  generating `spawns` from the manifest, not trusting agent files.

#### 11.4.3 Three authoring levels, one runtime

1. **Agent-driven (default, available first).** The entry agent is a
   supervisor whose `spawns` lists the plugin's internal agents. It uses
   `task` (batch, `outputSchema`) and `hub` messaging, and the store tool.
   No new runtime concept beyond the Phi services above.
2. **Declarative** `workflow.yaml` (later): sequence, parallel fan-out /
   fan-in, map over store records, loop until condition or budget, human
   checkpoint. For users who want repeatable pipelines without writing code.
3. **Programmatic** orchestrator (later): a script in the plugin environment
   that drives the same primitives over a local JSON-RPC channel, for
   deterministic procedures (Elo tournaments, evolutionary selection).

All three spawn agents through the same adapter, so limits, approvals,
store, and the run view behave identically.

#### 11.4.4 Decision for phase 1

Before freezing v1, prototype Phi's existing specialists (Database, Wrapper,
Visualization) on omp's `task` + registry through the adapter. If the
prototype keeps today's behaviour (approval routing, remote guards, run
cards, steer / stop, wake-up), retire the duplicated parts of
`agents/registry.ts` and base both specialist delegation and orchestration
on omp. If not, record the gaps in the ADR and keep Phi's registry for v1.

#### 11.4.5 Sketch of a Co-Scientist–style plugin

```yaml
id: co-scientist
type: plugin
components:
  agents:
    - agents/CoScientist.md        # visibility: entry; supervisor
    - agents/Generation.md         # visibility: internal
    - agents/Reflection.md
    - agents/Ranking.md
    - agents/Evolution.md
    - agents/MetaReview.md
  skills: [skills/elo-tournament]  # stdlib script run via skill_run
  orchestrator: null               # level 1: the supervisor agent orchestrates
uses: [Database]
store:
  hypotheses: schemas/hypothesis.json
  reviews: schemas/review.json
  matches: schemas/match.json
limits:
  concurrentAgents: 4
  totalAgentRuns: 200
  wallTime: 6h
  tokens: 5000000
```

The supervisor spawns Generation tasks in a batch (`outputSchema` =
hypothesis) → Reflection reviews each hypothesis → Ranking agents debate pairs
and record matches; the Elo update runs as a skill script over the store →
Evolution refines the top-k → MetaReview writes feedback that seeds the next
round → stop on budget or convergence → emit a `report` artifact and ranked
hypotheses. A later version could move the loop into a level-3 orchestrator
without changing the agents.

## 12. Security and trust

- Trust tiers, aligned with wrapper `trustTier`: `builtin`, `official`
  (signed registry), `user` (authored locally), `imported` (third-party
  files).
- `skill_run` participates in approvals: in `ask` mode, the first run of an
  `imported` skill's script needs confirmation; `official` content is trusted.
- The core fetch tool keeps the current URL policy and egress audit.
- Package extraction is path-contained (§6.2); signatures are mandatory for
  the remote registry.

## 13. Migration

- On first launch after the change, skills that appear in the user's session
  history are installed from the local registry; everything else moves to the
  catalogs. Beta data may be reset instead if simpler (the beta roadmap
  already allows clean resets).
- `~/.phi/wrappers/installed/` bundled entries are replaced by per-package
  installs; `custom` wrappers move to `~/.phi/wrappers/custom/`.
- Overlay packs under `~/.phi/wrappers/packs/` become inert and are removed by
  GC.
- `resources/` becomes the source tree of the local registry; the installer
  keeps only engine content (built-in skills, agent definitions needed by the
  engine, micromamba) once the remote registry ships.

## 14. Evaluation as a gate

`scripts/eval/` becomes a standing regression suite:

- **data access** — `db-vs-fetch.ts` extended with an `api-skill` arm and
  tasks for batch, pagination, GraphQL, SPARQL, rate limits;
- **routing** — does the main agent pick the right specialist / skill with a
  given installed set, and does it suggest installing a missing capability;
- **context cost** — first-turn prompt size and per-run tokens via
  `npm run usage:report`.

Each phase's exit criterion (§15) cites these numbers. Runs use at least two
models where credentials allow.

## 15. Phases

> **Sequencing is superseded by the prioritised plan in
> [content-distribution-implementation.md](../roadmap/content-distribution-implementation.md)
> (P0–P9, contracts frozen just before their first consumer).** The table
> below remains as the logical grouping of work.

| Phase | Scope | Exit criterion |
|---|---|---|
| 0 | This document, ADR in `docs/decisions/`, roadmap update | owner sign-off |
| 1 | **Kernel contracts v1** (§4.4): specs, JSON schemas, `phi validate`, conformance fixtures, reference implementations of each contract; prototype specialists on omp `task` / registry (§11.4.4); then freeze | every reference implementation validates and runs; omp decision recorded in the ADR; contracts tagged v1 |
| 2 | Content hygiene (allowlist builder, size budget), package format, local registry, installer, content loader with provenance, eval suite | install / upgrade / uninstall of one skill end to end; eval suite runs with one command |
| 3 | Catalog UI for skills, wrappers, connectors; wrapper per-family packages and tree assembly; enablement; migration | routing eval not worse; main first-turn prompt measurably smaller |
| 4 | micromamba, `phi-python`, environment resolution, `skill_run`, environment page | script skills run on a clean machine without manual setup |
| 5 | Core fetch tool, API skills, helper scripts, extended eval, db toolchain retirement | `api-skill` ≥ `db` on the extended eval; `src/main/agent/db/` substantially reduced |
| 6 | Plugin loader, visualization plugin (viz logic moved out of the engine), `phi-r`, developer-extensions rename | visualization installs, runs in shared `phi-r`, uninstalls cleanly; no `viz_*` code left in the engine |
| 7 | Remote signed registry, updates, mirrors, offline import, installer slimming | fresh install obtains content from the server; offline first launch works |
| 8 | Remote / HPC parity for skills and environments | skill scripts run on an SSH project and on an offline cluster |
| 9 | Orchestration on omp (§11.4): adapter-generated `spawns`, hard budgets, store tool, human checkpoints, run view; Co-Scientist–style reference plugin at level 1; freeze the Orchestration contract. Levels 2–3 follow when needed | the reference plugin completes a multi-round run within budget, survives an app restart (parked agents revived), and can be steered and stopped |

Dependencies: 0 → 1 → 2 → 3; 2 → 4 and 2 → 5 (4 and 5 can run in parallel);
4 → 6; 3, 5, 6 → 7 → 8; 6 → 9. Phases 2–8 implement frozen contracts; if one of them
needs a contract change, it goes through the rules in §4.4 instead of being
patched into the engine. Proposed: phases 0–5 in the internal beta, 6–9 after.

## 16. Open questions (with proposed defaults)

1. Beta scope — phases 0–5 in beta? *Proposed: yes.*
2. Run phases 4 and 5 in parallel? *Proposed: yes.*
3. Keep bundled db connectors as offline fallback during phase 5?
   *Proposed: yes, until the retirement gate passes.*
4. Update policy — auto or prompted? *Proposed: prompted; no silent upgrades
   for projects that recorded versions.*
5. Project lock of used packages/environments for collaborators?
   *Proposed: record in run records in phases 2–5; project-level lock later.*
6. Registry hosting — static object storage + CDN, no backend?
   *Proposed: yes; a domestic CDN for users in China.*
7. Signing key custody — CI secret vs offline signing? *Proposed: CI secret
   for the first-party registry.*
8. UI naming — "插件" for Phi plugins; namespace separator `/`.
   *Proposed: yes.*
9. Move `viz_*` out of the engine into the visualization plugin?
   *Confirmed: yes, as a CLI plus declared script tools in implementation step
   4; the engine keeps only artifact presentation.*
10. Contract format — JSON Schema files under `docs/contracts/` plus a
    validator in the engine? *Proposed: yes; schemas are the source of truth
    and the prose specs link to them.*
11. Orchestration authoring — start with agent-driven orchestration on omp
    (level 1) and add declarative / programmatic levels later?
    *Proposed: yes. Programmatic orchestrators, when added, are limited to
    `official` and `user` trust tiers until sandboxing is reviewed.*
12. Move Phi's specialist delegation onto omp's `task` / registry and retire
    the duplicated parts of `agents/registry.ts`? *Proposed: decide in phase 1
    from the prototype (§11.4.4).*

## Appendix A: db toolchain vs fetch evaluation (2026-09-29)

Harness: `scripts/eval/db-vs-fetch.ts` (tasks in `db-vs-fetch-tasks.ts`,
scoring in `db-vs-fetch-score.py`). Each arm runs as an isolated in-memory
agent session with only its own tools. Model: `cursor/cursor-grok-4.6-fast`
(through `scripts/eval/cursor-bridge.mjs`); Claude via Cursor was region
blocked. 14 tasks × 2 repetitions; ground truth verified against live APIs.

| Arm | Correct | Tool calls / task | Failed calls / task | Mean s | p90 s | Output tokens / task |
|---|---|---|---|---|---|---|
| `db` (Database prompt + 7 `db_*` tools) | 25/28 | 8.0 | 1.2 | 45.7 | 89.3 | 6891 |
| `fetch` (URL read only) | 27/28 | 3.0 | 0.2 | 23.9 | 44.4 | 2941 |
| `fetch-hints` (URL read + ~15 lines of API base URLs) | 26/28 | 2.5 | 0.1 | 20.9 | 37.8 | 2338 |

- `db_query` failed on 32 of 79 calls (filter/field format). Discovery tools
  (`db_search`, `db_routes`, `db_domain`, `db_resolve`) made 124 calls versus
  79 data calls.
- `db` failures: PDBe `entry_summary` lacks resolution (×2); AlphaFold domain
  lacks global pLDDT (×1).
- `fetch` failures: AlphaFold 403 (×2) and Ensembl 500 (×1) without retry.
- Not covered: batch downloads, pagination, credentialed APIs, obscure
  databases, rate limiting, UI result viewers. Cursor reports input tokens
  as 0, so only output tokens are comparable. Single fast-tier model.
