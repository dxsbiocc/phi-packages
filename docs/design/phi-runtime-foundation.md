# Phi Runtime Foundation: Environments, Execution, Binding, and Plugin Structure

Date: 2026-09-29
Updated: 2026-10-07 — visualization uses the shared `phi-r` environment.
Chinese version: [phi-runtime-foundation.zh-CN.md](phi-runtime-foundation.zh-CN.md).
Status: **Confirmed (2026-09-29; all six decisions in §9 accepted).** This
document precedes the [content distribution design](phi-content-distribution-design.md)
and the [implementation plan](../roadmap/content-distribution-implementation.md);
wherever those touch environments, execution, or plugin structure, this
document wins.

## 0. Why this document comes first

Content distribution (installing skills, wrappers, connectors, plugins) is the
outer layer. It depends on questions that were not settled:

1. Where does the bundled micromamba live, how is it invoked, and how is it
   isolated from the user's conda?
2. What is an environment, how is it identified, created, and collected?
3. What exactly does "run in an environment" mean, and how do we guarantee it
   does not leak into the host environment?
4. How do skills, specialist agents, core tools, notebooks, wrappers, and MCP
   servers each use environments?
5. Who declares an environment, who binds it, and who wins on conflict?
6. What is the plugin layout, and how does a plugin's agent manage its
   environment?

The answers go **inside-out**: L0 runtime layout → L1 environment model → L2
execution primitive → L3 consumers → L4 binding rules → L5 plugin structure →
L6 distribution. Each layer depends only on inner layers, and **building
follows the same order, small to large** (§8).

## 1. Current state

| Consumer            | How it gets a runtime today                                                                     | Problem                                              |
| ------------------- | ----------------------------------------------------------------------------------------------- | ---------------------------------------------------- |
| skill scripts       | agent `bash` runs host `python` / `Rscript`                                                     | depends on what the user installed; not reproducible |
| `viz_render`        | finds `Rscript`, `python3` on `PATH`                                                            | same; can only report "missing dependency"           |
| notebook kernels    | host `jupyter kernelspec`                                                                       | same                                                 |
| wrappers (Nextflow) | host `nextflow`; `-profile conda` uses host conda                                               | users must install nextflow and conda                |
| environment panel   | `environment/detect.ts` probes host micromamba, nextflow, docker, singularity, jupyter, Rscript | can detect, cannot provide                           |

Phi has no runtime of its own; everything depends on the host. That is the
root of the problems above it.

## 2. L0 Runtime layout

### 2.1 The micromamba binary

- **Bundled with the app** at `resources/runtime/micromamba/<platform>-<arch>/micromamba`
  (`darwin-arm64`, `darwin-x64`, `linux-x64`), shipped via `extraResources`
  outside the asar: `process.resourcesPath/runtime/micromamba/...` when
  packaged, the same repository path in development.
- Version and sha256 live in `resources/runtime/manifest.json`. A build script
  downloads and verifies the binary; it is not committed to git.
- Only the build platform's binary is bundled. The linux binary for remote
  hosts is downloaded and uploaded on demand in step 8 (§7).
- The engine locates it with one function, `getMicromambaPath()`; nothing
  looks it up on `PATH`.

### 2.2 Full isolation from the user's conda

Every micromamba call sets:

- `MAMBA_ROOT_PREFIX=~/.phi/runtime`;
- `--rc-file ~/.phi/runtime/mambarc`, with `CONDARC`, `MAMBARC`, and `CONDA_*`
  cleared, so the user's `~/.condarc` is never read;
- no `shell init`, no changes to the user's shell configuration, nothing added
  to the user's `PATH`.

### 2.3 Directory layout

```text
~/.phi/runtime/
  mambarc                 # channels (conda-forge, bioconda), channel_priority: strict, mirrors, proxy
  pkgs/                   # package cache; every environment hard-links from here
  envs/
    <envId>/              # environment prefix, read-only after creation
      .phi/env.json       # environment metadata (§3.3)
  state/environments.json # index: envId → status, referrers
  logs/                   # build logs
```

`envId = <name>-<hash12>`, where `hash` covers platform plus lock content
(§3.2).

## 3. L1 Environment model

### 3.1 Kinds

| Kind    | Examples                                             | Defined by                                              | Mutability                                               |
| ------- | ---------------------------------------------------- | ------------------------------------------------------- | -------------------------------------------------------- |
| base    | `phi-python`, `phi-r`, `phi-nextflow`, `phi-jupyter` | Phi developers, shipped with the app or registry        | immutable; a new version is a new environment            |
| package | a plugin-local or skill-owned environment            | plugin or skill author                                  | immutable                                                |
| project | extra packages a project needs                       | the user (requested by an agent, confirmed by the user) | immutable; adding dependencies creates a new environment |

`phi-r` is the shared base for R notebooks, scanpy R interoperability, and
the bundled visualization plugin. Visualization does not declare a private
package environment.

### 3.2 Specs and locks

- **Spec**: `environment.yml` (name, channels, dependencies; optional `pip:`).
- **Lock**: one **explicit lock** per platform (`@EXPLICIT`, one URL + md5 per
  line), produced in CI. Clients install with `micromamba create -p <prefix>
-f <lock>` and **never solve**: results are deterministic and installable
  offline from the cache.
- Official base and package environments must ship locks. A user project
  environment without a prebuilt lock is solved locally once, exported as an
  explicit lock, and installed from that lock afterwards.
- **Host dependencies** conda cannot provide (LibreOffice on macOS, Docker,
  Singularity) are declared in a `host:` section; they are checked and their
  paths recorded, never installed.
- **Source packages**: R packages without a conda build are declared in a
  `sourcePackages:` section, installed from source when the environment is
  built, and frozen together with the conda part (§3.5).

### 3.3 Metadata `.phi/env.json`

```json
{
  "envId": "phi-r-3f9a1c2b7d10",
  "name": "phi-r",
  "kind": "base",
  "platform": "darwin-arm64",
  "lockSha256": "…",
  "createdAt": "…",
  "micromambaVersion": "2.x",
  "activation": { "set": { "CONDA_PREFIX": "…", "JAVA_HOME": "…" }, "pathPrepend": ["…/bin"] },
  "host": {},
  "sourcePackages": [
    { "language": "r", "name": "ggsankey", "source": "github", "ref": "…", "sha256": "…" }
  ],
  "status": "ready"
}
```

- `activation` is an **activation snapshot**: after building, run
  `micromamba run -p <prefix> env` once and diff it against the sanitised base
  environment to capture what activation sets, including `activate.d` scripts
  (`JAVA_HOME`, `GDAL_DATA`, `PROJ_LIB`, …). Every later execution applies the
  snapshot directly: no shell activation, fast, injectable into any process.
- States: `absent` → `building` → `ready`; plus `failed` and `drifted` (doctor
  found a mismatch with the lock).

### 3.4 Lifecycle

- **ensure**: spec + lock → `envId` → return if ready; otherwise build (with
  progress and logs): install the conda part from the lock → install source
  packages (§3.5) → capture activation → make the prefix read-only → mark
  `ready`. If any step fails, the whole prefix is removed and the state is
  `failed`; no half-built environment is left behind.
- **Read-only**: after creation the prefix is read-only, so an agent's
  `pip install` into it fails. Immutability is enforced by the file system,
  not by prompts.
- **Reference counting and GC**: `environments.json` records which packages,
  plugins, and projects reference each environment; unreferenced environments
  are deleted. The package cache is trimmed by size separately.
- **Repair**: doctor compares an environment with its lock and rebuilds on
  drift.

### 3.5 Source packages

For packages without a conda build; currently R only, from CRAN or GitHub.

```yaml
# Phi extension section in environment.yml (stripped before locking, never passed to micromamba)
sourcePackages:
  - language: r
    name: ggsankey
    source: github # cran | github
    repo: davidsjoberg/ggsankey # github only
    ref: 5a3b1c… # github: full commit sha; cran: exact version such as 1.2.3
    sha256: 9f2e… # sha256 of the source archive
```

Rules:

- **Fully pinned**: GitHub accepts only full commit shas (no branches or
  tags); CRAN accepts only exact versions (fetched from the CRAN archive);
  every archive needs a sha256. The whole `sourcePackages` section is part of
  the envId hash, so any change yields a new environment.
- **No automatic dependencies**: installation uses `R CMD INSTALL`
  (equivalent to `install.packages(..., repos = NULL, dependencies = FALSE)`).
  Dependencies must come from the conda part or from earlier entries in the
  same section. A missing dependency is an error, never fetched from the
  network, so the environment is fully determined by its spec.
- **Compilation**: packages with native code require the matching compiler
  toolchain in the conda part (e.g. `compilers` alongside `r-base`); builds use
  only the environment's compilers, never the host's.
- **Cache and offline**: archives are cached at
  `~/.phi/runtime/sources/<sha256>` and used only after verification; a cached
  archive needs no network. Once the remote registry exists, archives can ship
  with packages.
- **Frozen**: activation is captured and the prefix made read-only only after
  source packages are installed, so they cannot change afterwards.
- **Recorded**: installed source packages (name, source, ref, sha256) are
  written to `env.json` and checked by doctor.

## 4. L2 Execution primitive

The engine has **one** implementation of "run in an environment", used by
every consumer:

```ts
environmentVariables(envRef, options): Record<string, string>   // for consumers that spawn processes themselves
runInEnvironment(envRef, argv, { cwd, stdin, timeout, signal }) // run directly
```

### 4.1 Sanitised base variables

Child processes do **not** inherit the Phi process environment. Instead:

- **keep**: `HOME`, `USER`, `LOGNAME`, `TMPDIR`, `TERM`, `LANG`, `LC_*` (UTF-8
  default, reusing `utf8LocaleEnv` from `render.ts`), `HTTP(S)_PROXY`,
  `NO_PROXY`, `SSH_AUTH_SOCK`;
- **clear**: `PYTHONPATH`, `PYTHONHOME`, `VIRTUAL_ENV`, `CONDA_*`, `MAMBA_*`,
  `R_LIBS`, `R_LIBS_USER`, `R_LIBS_SITE`, `R_HOME`, `JAVA_HOME`, `PERL5LIB`,
  `LD_LIBRARY_PATH`, `DYLD_*`, `NODE_PATH`;
- **set**: `PYTHONNOUSERSITE=1`, `R_LIBS_USER=` (empty), `MPLBACKEND=Agg`;
- **PATH**: `<env>/bin` + other activation paths + directories of declared
  host dependencies + the minimal system path (`/usr/bin:/bin:/usr/sbin:/sbin`).
  **Never** the user's conda, Homebrew, pyenv, or `~/.local/bin`.

### 4.2 Verifiable

- Every execution returns `envId` and the resolved interpreter path and writes
  them to the run record (reproducibility).
- Conformance tests (§8 step 2) prove the interpreter is inside the prefix, a
  host-only package fails to import, and execution works on a clean machine
  without any system scientific stack.

## 5. L3 Consumers

Every consumer executes through L2, with two exceptions:

- **The main agent's bash** represents the user working on their own machine
  and keeps the host environment.
- **User-facing interactive tools where the user explicitly picks the host
  version**: notebook kernels and Nextflow default to managed environments,
  but the user may explicitly choose a host kernel or host nextflow (§5.2).

**Content execution** (skill scripts, script tools, plugin and specialist
agents) **always** uses managed environments and never reuses the host, even
when the host has the same tool.

| Consumer                         | Today                                                     | Target                                                                                                                                                                                                                                                                               |
| -------------------------------- | --------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| skill scripts                    | agent `bash` with host interpreters                       | declared script tools (§5.1) or `skill_run(skill, script, args)` → `runInEnvironment(skill environment)`                                                                                                                                                                             |
| specialist / plugin agent `bash` | host shell                                                | the session is bound to an environment: a Phi extension rewrites bash calls in the `tool_call` event and merges `environmentVariables(agent environment)` into the bash `env` parameter (omp supports returning replacement input). `python x.py` then uses the environment's python |
| domain tools (e.g. `viz_render`) | TypeScript tools inside the engine, `Rscript` from `PATH` | leave the engine; become command-line programs in the plugin, registered as script tools (§5.1); tool names unchanged                                                                                                                                                                |
| MCP stdio servers                | —                                                         | started with `environmentVariables`                                                                                                                                                                                                                                                  |
| notebook kernels                 | host kernelspecs                                          | Jupyter server in `phi-jupyter`; default kernels from managed analysis environments (`ipykernel` in `phi-python`, `irkernel` in `phi-r`); existing host kernels listed for explicit selection (§5.2)                                                                                 |
| wrappers (Nextflow)              | host `nextflow`, host conda                               | default `nextflow` and Java from `phi-nextflow`; `-profile conda` sets `conda.useMicromamba = true` with the bundled micromamba and `conda.cacheDir` under `~/.phi/runtime`; the user may explicitly point at a host nextflow (§5.2); Docker / Singularity remain host dependencies  |
| environment panel                | probes host tools                                         | shows managed environments (status, size, referrers, host dependencies); host probing kept only for host dependencies (Docker, Singularity, LibreOffice) and tools the user may explicitly choose (nextflow, Jupyter kernels)                                                        |

### 5.1 Script tools

The engine holds no domain tools. Skills and plugins get typed tools by
**declaring script tools** in SKILL.md; the engine registers them from the
declaration, with no per-domain code.

```yaml
# SKILL.md frontmatter (excerpt). Phi fields live under the phi block; see docs/contracts/skill.md
name: omics-visualization
description: …
phi:
  environment: phi:r@1
  scripts:
    - name: route
      run: [python, ./scripts/viz.py, route]
      description: Profile a result table and shortlist fitting figure templates
      args: # JSON Schema
        type: object
        required: [data, purpose]
        additionalProperties: false
        properties:
          data: { type: string, format: input-path } # read-only, must be inside the project
          purpose: { type: string }
      approval: read # read | write
      output: ./schemas/route-result.json # JSON Schema
    - name: render
      run: [python, ./scripts/viz.py, render]
      args:
        type: object
        required: [script, output]
        additionalProperties: false
        properties:
          script: { type: string, format: input-path }
          output: { type: string, format: project-path } # writable, must be inside the project
      approval: write
```

Engine responsibilities (all generic):

- **Registration**: tool name `<toolPrefix>_<name>` (§10), registered for the
  agents that declare the skill; for the main agent when the skill's
  `attachTo` points there.
- **Arguments**: tool parameters generated from and validated against `args`.
- **Paths**: `input-path` is read-only, `project-path` writable; both must
  resolve inside the project or the call is rejected.
- **Approvals**: `approval` maps onto the existing approval flow.
- **Execution**: `run` executes through `runInEnvironment` in the skill's
  environment, with arguments passed as `--<name> <value>`.
- **Output**: the program writes one JSON object to stdout, validated against
  `output`; on failure it exits non-zero with `{"error": "..."}`.
- **Artifacts**: the program writes `<file>.phi-artifact.json` (kind, title,
  provenance); the engine presents the artifact in the UI.

An agent in an environment-bound session can still run these programs with
`bash`; script tools are the more reliable path: typed arguments, constrained
paths, and graded approvals.

### 5.2 Explicitly choosing host versions

- **Notebook kernels**: the kernel list shows managed kernels and existing host
  kernels (from host kernelspecs). Host kernels are labelled "host
  (unmanaged)" and must be chosen explicitly; they are never the default.
- **Nextflow**: reuse the existing `customPaths` in `~/.phi/environment.json`;
  the user sets a host nextflow path in the environment panel. The version is
  checked against wrapper minimums and rejected if too old; otherwise it is
  labelled "host (unmanaged)". This mainly serves HPC sites with their own
  nextflow.
- Host versions carry no reproducibility guarantee; run records note the host
  version and path.

## 6. L4 Binding and resolution

### 6.1 Who declares environments

| Where                         | Syntax                                                                                          | Meaning                                                                       |
| ----------------------------- | ----------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------- |
| skill frontmatter `phi` block | `environment: phi:python@1` / `environment: ./environment.yml` / `environment: plugin:statistics` | where the skill's scripts and script tools run                                |
| agent frontmatter             | `environment: phi:r@1` / `environment: plugin:statistics`                                      | which environment the agent session's bash and core tools use                 |
| plugin manifest               | `environments: { statistics: {...} }`                                                          | named environments shipped by the plugin, referenced by its skills and agents |
| project                       | `.phi/environment.yml` (optional), referenced as `project:default`                              | extra dependencies for the project                                            |

`plugin:<name>` resolves only within the same plugin; cross-plugin references
are not allowed.

The bundled visualization plugin uses the official `phi:r@1` reference in
both its skill and agent frontmatter and therefore omits the manifest's
optional `environments` field.

### 6.2 Resolution order

- **`skill_run(skill, …)`**: the skill's own environment → the environment
  bound to the calling agent session → `phi-python` (with an "undeclared
  dependencies" warning).
- **agent session**: the agent's declared environment → unbound (host, like
  the main agent).
- Two environments in one session never mix: the agent's bash uses the
  agent's environment; a skill that declares a different one runs separately
  in that environment through `skill_run`.

### 6.3 Extra packages

Environments are read-only; agents cannot install into them. For extra
dependencies:

1. the agent calls `env_request(packages, reason)` (a reserved core
   tool);
2. the user confirms in the UI;
3. the engine solves "original spec + extras" into a new **project**
   environment, exports a lock, and records it in project state;
4. from then on, that agent / skill in that project binds to the new
   environment. The original environment and other projects are untouched.

## 7. L5 Plugin structure

### 7.1 Layout

```text
<plugin-id>/
  phi-package.yaml
  agents/
    Visualization.md            # frontmatter: visibility, environment, skills, tools, spawns
  skills/
    omics-visualization/
      SKILL.md                  # frontmatter: environment: phi:r@1, scripts: [...]
      scripts/  references/  assets/
  mcp/                          # optional: stdio server definitions
  wrappers/                     # optional
  orchestrator/                 # reserved
```

### 7.2 Visualization manifest

```yaml
id: visualization
type: plugin
version: 1.0.2
toolPrefix: viz # script tool prefix, see §10
components:
  agents: [agents/Visualization.md]
  skills: [skills/omics-visualization]
```

The visualization manifest has no `environments` field. Its agent and skill
frontmatter bind both execution paths to `phi:r@1`.

### 7.3 How the visualization plugin uses the shared environment

1. **On plugin install** the engine registers the agent, skill, and script
   tools. There is no plugin-private environment to build.
2. **On agent session creation** the engine reads `environment: phi:r@1` →
   ensures the official environment is ready (asking the user to build it
   rather than silently falling back to the host) → attaches the bash
   injection extension → passes the same environment to core tools through
   the run context → the plugin's skill resolves to the same environment via
   `skill_run`.
3. **While running** the environment is read-only; extra packages go through
   §6.3.
4. **On the upgrade that removes the legacy private environment** the active
   version switches atomically, its old `plugin:viz` reference is dropped,
   and the existing GC collects that environment once it is unreferenced.
5. **On uninstall** plugin components are removed. The plugin does not own the
   shared official `phi-r` environment.

The agent does not "manage" its environment; it only **declares** which one
it uses. Creation, activation, isolation, and collection belong to the engine,
so plugin authors write no environment code and the engine needs no
per-plugin special cases.

## 8. Build order (inside-out, small to large)

Each step starts only after the previous one is complete and tested.

| Step | Layer   | Scope                                                                                                                                                                                                                                | Done when                                                                                     |
| ---- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------- |
| 1    | L0 + L1 | bundle micromamba; `~/.phi/runtime` layout and `mambarc`; create from explicit locks, activation snapshot, read-only, index, GC                                                                                                      | a python-only lock builds an environment on a clean account; repeated ensure does not rebuild |
| 2    | L2      | `environmentVariables` / `runInEnvironment`; sanitisation; canary tests (interpreter location, negative import, clean machine)                                                                                                       | all three isolation tests pass                                                                |
| 3    | L3      | `skill_run`; first `phi-python`; migrate one skill (e.g. scanpy)                                                                                                                                                                     | the skill works on a machine without a host Python scientific stack                           |
| 4    | L4      | agent frontmatter `environment`; bash injection extension; visualization rewritten as the command-line program `scripts/viz.py`, declared as script tools, running in `phi-r`; `src/main/agent/visualization/` deleted from the engine | visualization renders on a machine without host R; no visualization code in the engine        |
| 5    | L3      | remaining consumers: notebook kernels (`phi-jupyter`), Nextflow (`phi-nextflow` + `conda.useMicromamba`), MCP stdio; explicit host versions (§5.2)                                                                                   | notebooks and wrappers work on a machine without host jupyter / nextflow / conda              |
| 6    | L5      | plugin layout and manifest; package the already rewritten visualization locally and install / run / uninstall it per §7.3                                                                                                            | the full plugin lifecycle works                                                               |
| 7    | L6      | content distribution: units, catalogs, installer, registry (see the content distribution design)                                                                                                                                     | per that design                                                                               |
| 8    | remote  | remote / HPC: upload linux micromamba, build from locks remotely; offline clusters via package cache or conda-pack                                                                                                                   | step-3 skills run in an SSH project and on an offline cluster                                 |

Steps 1–6 are the kernel. After them, how skills, agents, plugins, wrappers,
and notebooks use environments is fixed; outer distribution just moves files
on top of it.

## 9. Confirmed foundational decisions (2026-09-29)

1. **Only the main agent's bash uses the host environment**; all other
   execution runs in managed environments.
2. **Official environments use explicit locks; clients never solve**; only
   user project environments are solved locally, once.
3. **Environment prefixes are read-only after creation**; extra dependencies
   go through `env_request` into a new project environment.
4. **LibreOffice, Docker, and Singularity are host dependencies**: checked,
   not installed.
5. **Nextflow and Jupyter default to managed environments**; notebook kernels
   and Nextflow may explicitly use a host version (Nextflow after a version
   check), labelled unmanaged. Content execution (skills, script tools,
   plugin and specialist agents) always uses managed environments.
   (Revised 2026-09-29.)
6. micromamba is bundled with the app, **current platform only**; remote
   hosts get theirs on demand.
7. **No domain tools in the engine**: domain logic lives in command-line
   programs inside plugins or skills, exposed as typed tools by declaring
   script tools in SKILL.md (§5.1). `viz_*` migrates first, keeping its tool
   names.
8. **Unified naming conventions**: see §10.

## 10. Naming conventions

| Object                               | Rule                                                                                                                                                                                                                                                                                                  | Examples                                                                                        |
| ------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------- |
| core tools (engine)                  | snake_case, `<domain>_<verb>`                                                                                                                                                                                                                                                                         | `skill_run`, `env_request`, `http_fetch`, `wrapper_search`, `agent_status`                      |
| engine-reserved tool prefixes        | not usable by content                                                                                                                                                                                                                                                                                 | `skill_`, `env_`, `http_`, `wrapper_`, `agent_`, `db_` (retired but still reserved), `mcp__` (omp MCP tools) |
| script tools                         | `<toolPrefix>_<script name>`; `toolPrefix` declared in the package manifest, 2–12 lowercase letters or digits, unique in the registry, no clash with reserved prefixes                                                                                                                                | `viz_route`, `viz_render`                                                                       |
| package ids, skill names, plugin ids | kebab-case                                                                                                                                                                                                                                                                                            | `omics-visualization`, `visualization`, `protein-apis`                                          |
| agent names                          | PascalCase                                                                                                                                                                                                                                                                                            | `Visualization`, `Database`, `Wrapper`                                                          |
| environment names                    | kebab-case; Phi-maintained environments use a `phi-` prefix and are named by purpose; plugin-local environments have no prefix                                                                                                                                                                        | `phi-python`, `phi-r`, `phi-nextflow`, `phi-jupyter`; `statistics`                              |
| environment references               | `<scope>:<name>[@<major>]`, or a relative path                                                                                                                                                                                                                                                        | `phi:python@1`, `phi:r@1`, `plugin:statistics`, `project:default`, `./environment.yml`          |
| envId                                | `<scope>-<owner>-<name>-<hash12>`; owner omitted for scope `phi`                                                                                                                                                                                                                                      | `phi-python-3f9a1c2b7d10`, `plugin-reports-statistics-…`, `project-<short project id>-default-…` |
| package manifest                     | `phi-package.yaml` for every type, distinguished by `type`                                                                                                                                                                                                                                            | `type: skill` / `wrapper` / `mcp` / `plugin`                                                    |
| frontmatter and manifest fields      | camelCase, aligned with omp; exception: SKILL.md standard fields keep the Agent Skills spelling (`allowed-tools`, `disable-model-invocation`) and all Phi fields go under the `phi` block; legacy snake_case fields (e.g. `delegation_mode`) are read as aliases and the validator suggests migration | `thinkingLevel`, `outputSchema`, `attachTo`, `toolPrefix`, `delegationMode`                     |
| script entry points                  | `scripts/<toolPrefix>.py`, subcommands named after the script tools                                                                                                                                                                                                                                   | `scripts/viz.py route`                                                                          |
| artifact descriptors                 | `<file>.phi-artifact.json`                                                                                                                                                                                                                                                                            | `fig.png.phi-artifact.json`                                                                     |
| contract files                       | `docs/contracts/<name>.schema.json` (kebab-case); prose specs `<name>.md`                                                                                                                                                                                                                             | `environment.schema.json`, `env-metadata.schema.json`, `execution.md`                           |
| user directories                     | `~/.phi/runtime/` (micromamba root, environments, package cache), `~/.phi/packages/<type>/<id>/<version>/` (all packages, including plugins), `~/.phi/state/` (enablement and other state)                                                                                                            |                                                                                                 |
| repository directories               | `resources/runtime/environments/<env name>/`, `resources/plugins/<plugin id>/`                                                                                                                                                                                                                        | `resources/runtime/environments/phi-python/`                                                    |

The implementation plan follows §8:
[content-distribution-implementation.md](../roadmap/content-distribution-implementation.md).
The environment (§9) and plugin (§11) chapters of the content distribution
design defer to this document.
