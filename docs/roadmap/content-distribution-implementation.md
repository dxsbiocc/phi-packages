# Runtime and Content Distribution Implementation Plan

Date: 2026-09-29 (rewritten per [runtime foundation](../design/phi-runtime-foundation.md) §8; same day: R source packages, docx and scvi-tools skills removed, nf-core and nf-test in `phi-nextflow`)
Updated: 2026-10-07 — visualization dependencies consolidated into `phi-r`.
Based on: [runtime foundation](../design/phi-runtime-foundation.md) (confirmed), [content distribution design](../design/phi-content-distribution-design.md), [decisions](../decisions/content-distribution.md).
Chinese version: [content-distribution-implementation.zh-CN.md](content-distribution-implementation.zh-CN.md).

## 1. Principles

1. **Inside-out, small to large.** Step order follows the runtime foundation's
   layers: runtime → execution primitive → consumers → binding → plugins →
   distribution → remote. Each step depends only on completed inner layers.
2. **No outer work before the inner step passes.** Every step has explicit
   acceptance criteria.
3. **Contracts freeze with their layer.** When a layer is done, its contract is
   frozen as v1 and changes only under content distribution design §4.4.
4. **Work outside the main line is listed separately.** Work that neither
   belongs to the layered line nor changes the kernel (data access redesign,
   multi-agent orchestration) is a side track and does not reorder the main
   line.
5. Each task is 0.5–3 person-days, merges and reverts on its own, and leaves
   the app working.

## 2. Overview

| Step | Layer | Goal | Contracts frozen | Done when | Est. | Beta |
|---|---|---|---|---|---|---|
| 0 | prep | cleanup and baselines, no new architecture | — | clean installer; dependency inventory and context baseline exist | 3 d | yes |
| 1 | L0 + L1 | bundled micromamba, `~/.phi/runtime`, environments from locks | Environment | an environment built from a lock on a clean account; idempotent, read-only, collectable | 2 w | yes |
| 2 | L2 | single execution primitive and isolation | Execution (`runInEnvironment`) | three isolation tests pass locally and in CI | 1 w | yes |
| 3 | L3 | `skill_run` and script tools; skill scripts run in managed environments | Skill, `skill_run`, script tools | official skills run on a machine without a host Python scientific stack | 2.5 w | yes |
| 4 | L4 | agent sessions bound to environments; visualization rewritten as a CLI plus script tools and removed from the engine | Agent definition, Artifact | visualization renders on a machine without host R; no visualization code in the engine | 3.5 w | yes |
| 5 | L3 | remaining consumers: Nextflow, Jupyter, MCP stdio; explicit host versions; environment panel | — | wrappers and notebooks work without host nextflow, conda, or jupyter | 2 w | yes |
| 6 | L5 | plugin structure and loader; visualization installed as a plugin | Plugin | plugin install / run / upgrade / uninstall | 1.5 w | no |
| 7 | L6 | content distribution: packages, installer, catalogs, enablement, remote registry | Package, enablement state, Wrapper, Connector | per content distribution design | batches | partly |
| 8 | remote | runtime and environments on remote / HPC hosts | — | step-3 skills run in an SSH project and on an offline cluster | 2 w | no |

After steps 1–6 the kernel is fixed (naming conventions: [runtime foundation §10](../design/phi-runtime-foundation.md)): how skills, agents, plugins, wrappers,
and notebooks use environments no longer changes; steps 7 and 8 only
distribute content on top of it and extend it to remote hosts.

Estimates are single-person. Proposed beta scope: steps 0–5.

## 3. Details

Format: scope / main files / acceptance / estimate.

### Step 0 Preparation

| ID | Scope | Main files | Acceptance | Est. |
|---|---|---|---|---|
| 0.1 | Remove run leftovers from `resources/`; add `scripts/check-resources.mjs` (fails on untracked files) to `npm run lint` | `resources/`, `scripts/` | `git status --ignored resources` clean | 0.5 d |
| 0.2 | Exclude leftovers in `electron-builder.yml` `files` | `electron-builder.yml` | no `.nextflow` in `build:unpack` output | 0.5 d |
| 0.3 | Delete `resources/skills/create-database-connector` and its test references | `resources/skills/`, `tests/resources.test.ts`, `tests/phi-agents.test.ts` | `npm test` passes | 0.5 d |
| 0.4 | Commit `scripts/eval/` with a README and `bun run eval:fetch` | `scripts/eval/` | one command reproduces the URL-reading evaluation | 0.5 d |
| 0.5 | **Dependency inventory** per skill and visualization script: Python packages, R packages, external commands (known: `soffice`, pandoc, poppler, tesseract, Node `docx-js`, R) | `docs/runtime/dependency-inventory.md` | inventory for all script content; input to steps 3 and 4 | 1 d |

### Step 1 Runtime and environment model (L0 + L1)

New module: `src/main/agent/envs/`. The existing `src/main/agent/environment/`
(host tool detection) shrinks in step 5.

| ID | Scope | Main files | Acceptance | Est. |
|---|---|---|---|---|
| 1.1 | **Freeze Environment contract v1**: supported `environment.yml` fields, per-platform explicit lock format, `host:` section, `sourcePackages:` section (R source packages, foundation §3.5), `env.json` schema, `envId` rule, state machine | `docs/contracts/environment.schema.json`, `docs/contracts/env-metadata.schema.json` | schemas and fixtures validate | 1 d |
| 1.2 | Bundle micromamba: `resources/runtime/manifest.json` (version, per-platform URL and sha256); `scripts/runtime/fetch-micromamba.mjs` downloads and verifies for build and dev; `extraResources`; `getMicromambaPath()` | `resources/runtime/`, `scripts/runtime/`, `electron-builder.yml`, `envs/paths.ts` | `micromamba --version` runs in dev and packaged builds | 2 d |
| 1.3 | Runtime directory and invocation wrapper: create `~/.phi/runtime/{envs,pkgs,logs,state}`; generate `mambarc` (conda-forge, bioconda, strict, mirrors and proxy from settings); `micromamba(args)` pins `MAMBA_ROOT_PREFIX` and `--rc-file` and clears `CONDA_*` / `MAMBA_*` | `envs/runtime.ts` | on a machine with user conda, wrapped calls ignore `~/.condarc` (asserted) | 1.5 d |
| 1.4 | Lock tooling: `scripts/runtime/lock-env.mjs` uses `micromamba --platform <p> --dry-run` to produce explicit locks for three platforms | `scripts/runtime/` | a python-only fixture spec yields three locks | 1.5 d |
| 1.5 | `ensureEnvironment(spec, lock)`: compute `envId` → return if ready → file lock against concurrent builds → `micromamba create -p … -f lock` → install source packages (1.7) → capture activation snapshot → make prefix read-only → write `env.json` → update `state/environments.json`; progress events and logs | `envs/ensure.ts`, `envs/activation.ts`, `envs/index-store.ts` | fixture environment built on a clean account; repeated calls do not rebuild; writes into the prefix fail; any failed step removes the prefix and records `failed` | 3 d |
| 1.6 | Reference counting, GC, doctor: record referrers; delete unreferenced environments; compare `micromamba list --json` with the lock and rebuild on drift | `envs/gc.ts`, `envs/doctor.ts` | tests for GC, drift detection, rebuild | 1.5 d |
| 1.7 | Source packages (R): fetch the CRAN exact version or GitHub full-commit archive per `sourcePackages` → verify sha256 → cache at `~/.phi/runtime/sources/<sha256>` → `R CMD INSTALL` in declared order (no automatic dependencies) → record in `env.json`; only the environment's compilers are used | `envs/source-packages.ts` | a real install of a small pure-R package; tests for sha256 mismatch, missing dependency, cache hit without network | 2 d |

### Step 2 Execution primitive (L2)

| ID | Scope | Main files | Acceptance | Est. |
|---|---|---|---|---|
| 2.1 | **Freeze Execution contract v1**: sanitisation rules (kept / cleared / set variables, `PATH` composition); `runInEnvironment` inputs and outputs (`envId`, interpreter path, exit code, output truncation) | `docs/contracts/execution.md` | spec reviewed | 0.5 d |
| 2.2 | Sanitisation module with unit tests | `envs/sanitize.ts` | a test per rule | 1 d |
| 2.3 | `environmentVariables(envRef)`, `runInEnvironment(envRef, argv, opts)`: timeout, abort, streaming; host-dependency directories on `PATH` | `envs/run.ts` | unit tests | 1.5 d |
| 2.4 | Canary environment and isolation tests: (a) `sys.executable` and `sys.prefix` under the prefix; (b) a host-only package fails to import; (c) user site disabled; R `.libPaths()` once an R environment exists | `tests/runtime/`, `npm run test:runtime` | all pass locally | 1.5 d |
| 2.5 | CI on clean macOS and Ubuntu runners for `test:runtime` | CI config | CI green; gates step 3 | 1 d |

### Step 3 Skills run in environments (L3)

| ID | Scope | Main files | Acceptance | Est. |
|---|---|---|---|---|
| 3.1 | **Freeze Skill, `skill_run`, and script-tool declaration contracts v1**: frontmatter `environment` (`phi:python@1` / `./environment.yml` / `plugin:<name>` / `project:default`), `attachTo`, `scripts` (`name`, `run`, `description`, `args`, `approval`, `output`; path `format` `input-path` / `project-path`); script skills must declare an environment; `skill_run` parameters, result (with `envId` and interpreter path), errors, approvals; `validateSkill()` | `docs/contracts/skill.schema.json`, `src/main/agent/content/validate-skill.ts` | fixtures pass; `npm run lint` validates all skills | 2.5 d |
| 3.2 | `phi-python` v1 from the 0.5 inventory (including pandoc, poppler, tesseract, nodejs, ipykernel); locks in `resources/runtime/environments/phi-python/` | `resources/runtime/environments/` | locks build on all three platforms | 2 d |
| 3.3 | `skill_run` core tool: resolution (skill → session binding → `phi-python` with a warning); if the environment is not ready, ask the user to build it with size and progress, never fall back to the host; approvals in ask mode; registered for the main agent and specialists | `src/main/agent/content/skill-run.ts` | unit tests plus one real run | 3 d |
| 3.4 | Script-tool registration: generate `<toolPrefix>_<name>` tools from `scripts`; argument validation; `input-path` / `project-path` constraints; approval levels; execution through `runInEnvironment`; JSON output validation; read `<file>.phi-artifact.json` | `src/main/agent/content/script-tools.ts` | a fixture skill covers arguments, out-of-project paths, approvals, output validation, error exits | 3 d |
| 3.5 | Migrate scanpy: declare its environment; SKILL.md script calls become `skill_run` | `resources/skills/scanpy/` | example runs without a host Python scientific stack | 1 d |
| 3.6 | Migrate the remaining Python skills (pptx, xlsx, pdf, markitdown, matplotlib, scikit-learn, scvelo, rdkit); LibreOffice as a host dependency | `resources/skills/` | every skill passes `validateSkill()` and runs one example via `skill_run` | 3 d |

### Step 4 Agents bound to environments; visualization leaves the engine (L4)

| ID | Scope | Main files | Acceptance | Est. |
|---|---|---|---|---|
| 4.1 | **Spike**: prototype Wrapper and Visualization on omp `task` / registry; decide whether delegation moves to omp and whether duplicated parts of `agents/registry.ts` retire; record in the decision record | `src/main/agent/agents/` | decision record updated | 1 w |
| 4.2 | **Freeze Agent-definition contract v1**: omp fields (`name`, `description`, `tools`, `spawns`, `model`, `thinkingLevel`) plus Phi `environment`, `visibility`, `skills`, `delegationMode`, `delegation`, `fallback` (legacy `delegation_mode` read as an alias); structured results via `outputSchema` | `docs/contracts/agent.schema.json`, `agents/definition.ts` | the three existing agents validate | 2 d |
| 4.3 | bash injection extension: in bound sessions, merge `environmentVariables` into the bash call's `env` in the `tool_call` event | `src/main/agent/agents/`, `omp/omp-sdk-worker.ts` | `which python` points into the environment in bound sessions; the main agent's bash is unaffected | 2 d |
| 4.4 | Create the plugin-shaped directory `resources/plugins/visualization/`: move in `resources/agents/Visualization.md` and `resources/skills/omics-visualization`; merge the complete visualization dependency set into the official `phi-r` spec and locks (159 R scripts plus Python 3.12 for standard-library scripts). Keep the seven direct R packages without conda builds — gground, ggideogram, ggcor, linkET, ggsankey, ggsvg, ggmagnify — plus ggmagnify's gridGeometry dependency as eight ordered, pinned `sourcePackages`; preserve `phi-r`'s existing GenomeInfoDbData source package and settle the ggideogram / ggplot2 4.x incompatibility here. The plugin declares no private environment | `resources/plugins/visualization/`, `resources/runtime/environments/phi-r/` | `phi-r` builds on all three platforms; all templates pass a smoke render; IRkernel, Seurat, and SingleCellExperiment still load | 4 d |
| 4.5 | **Freeze Artifact contract v1**: fields of `<file>.phi-artifact.json` (kind `figure` / `table` / `structure` / `molecule` / `network` / `report`, title, provenance); the engine presents artifacts | `docs/contracts/artifact.schema.json`, presentation layer | fixtures pass; generic artifact viewers read artifacts | 2 d |
| 4.6 | Rewrite visualization as the command-line program `scripts/viz.py` (subcommands `examples`, `route`, `prepare`, `render`): `route` calls the existing `route_template.py`; `prepare` and `examples` are ported from TypeScript to Python; `render` runs `Rscript` in the same environment plus QA and writes a `figure` artifact | `resources/plugins/visualization/skills/omics-visualization/scripts/` | tests per subcommand; outputs validate against their JSON Schemas | 4 d |
| 4.7 | Declare four script tools in SKILL.md (`toolPrefix: viz`, names kept as `viz_examples`, `viz_route`, `viz_prepare`, `viz_render`); the Visualization agent and `omics-visualization` skill both declare `environment: phi:r@1`; workflow-based tool filtering moves into the agent's instructions; session creation ensures the environment or asks to build it | SKILL.md, `Visualization.md`, session creation | renders on a machine without host R; no regression on the visualization eval | 2 d |
| 4.8 | Delete `src/main/agent/visualization/` and visualization-specific logic such as `visualizationToolNamesForWorkflow` from the engine | `src/main/agent/` | no visualization code in the engine; `npm test` passes | 1 d |
| 4.9 | `env_request`: agent asks for extra packages → user confirms → solve "original spec + extras" into a new project environment with a lock → the project rebinds | `src/main/agent/envs/`, project state | tests for request, confirmation, rebinding | 2.5 d |

Until step 6, `resources/plugins/visualization/` is mounted by temporary built-in loading (registering its agent and skill directly); step 6 replaces it with the real plugin loader.

### Step 5 Remaining consumers (L3)

| ID | Scope | Main files | Acceptance | Est. |
|---|---|---|---|---|
| 5.1 | `phi-nextflow` environment (nextflow + openjdk + nf-core + nf-test; nf-core and nf-test are used to author and test wrapper modules, about 380 MB) used by default by the wrapper executor; explicit host nextflow kept via `customPaths` in `environment.json`, accepted only if it meets wrapper minimum versions and labelled "host (unmanaged)" | `wrappers/composition/executor.ts`, `environment/store.ts` | wrappers run without host nextflow; an outdated host nextflow is rejected | 2.5 d |
| 5.2 | `-profile conda` sets `conda.useMicromamba = true`, the micromamba path, and `conda.cacheDir` under `~/.phi/runtime` | wrapper `nextflow.config` generation | `npm run smoke:wrappers` passes without host conda | 2 d |
| 5.3 | Jupyter server in `phi-jupyter`; default kernels from `phi-python` (ipykernel) and `phi-r` (irkernel); existing host kernels listed as "host (unmanaged)", selectable only explicitly | `notebook/analysis-kernels.ts` | notebooks work without host jupyter; host kernels selectable but never default | 3 d |
| 5.4 | MCP stdio servers start through `environmentVariables` (mechanism only; user-configured servers in `mcp.json` stay as they are) | MCP code | unit tests | 1 d |
| 5.5 | Environment panel: managed environments (status, size, referrers, rebuild, clean); host dependencies (Docker, Singularity, LibreOffice); host tools that may be chosen explicitly (nextflow, Jupyter kernels) | `features/environment/`, `environment/detect.ts` | usable UI; old detection fields migrated | 2.5 d |

Wrappers in remote projects keep using the remote host's nextflow until step 8.

### Step 6 Plugins (L5)

| ID | Scope | Main files | Acceptance | Est. |
|---|---|---|---|---|
| 6.1 | **Freeze Plugin contract v1**: layout, `phi-package.yaml` (`type: plugin`, `toolPrefix`), `environments`, `components`, agent `visibility`, reserved `orchestrator` | `docs/contracts/plugin.schema.json` | fixtures pass | 2 d |
| 6.2 | Plugin loader: install from a local directory into `~/.phi/packages/plugin/<id>/<version>/`; ensure environments; register agents, skills, and script tools; upgrade builds the new environment before switching; uninstall removes references and collects environments | `src/main/agent/plugins/` (new) | tests for install, upgrade, uninstall | 3 d |
| 6.3 | Add `phi-package.yaml` without `environments` to `resources/plugins/visualization/` and install it through the plugin loader; remove step 4's temporary loading. Bump the bundled plugin version so upgrade drops the legacy `plugin:viz` reference; the existing GC collects that environment once unreferenced | `resources/`, `src/main/agent/` | works after installing or upgrading through the loader; no visualization-private environment remains referenced | 1 d |
| 6.4 | Plugins page (list, install from local, uninstall); rename the pi plugin page to Developer extensions under advanced settings | `features/plugin/` | usable UI | 2 d |

### Step 7 Content distribution (L6)

Moves content on top of the kernel; details in the content distribution
design. Batches:

1. **Package contract and local registry**: freeze Package v1; allowlist
   packaging from git-tracked files; installer (verify, staging, atomic
   install, `.source.json`, upgrade, uninstall, GC) calling step 1's
   `ensureEnvironment`.
2. **Skill catalog and enablement**: freeze the enablement-state contract v1;
   filter bundled skills by enablement; "Add from catalog" and toggles on the
   Skills page without rewriting skill files; migrate from session history.
3. **Wrapper and connector packaging** (post-beta): freeze Wrapper and
   Connector v1; per-family wrapper packages assembled into one tree;
   registry-backed MCP catalog.
4. **Remote registry** (post-beta): signed index, static storage + CDN,
   update prompts, mirrors, offline import, slim installer.

### Step 8 Remote / HPC

- Download the linux micromamba on demand and upload it; the remote runtime
  directory mirrors the local layout.
- Build environments from locks remotely; offline clusters use a prefetched
  package cache or conda-pack; locks compatible with glibc 2.17.
- Remote wrappers use the remote `phi-nextflow` environment.

## 4. Side tracks (no kernel changes, no reordering)

| Track | Needs | Scope |
|---|---|---|
| Data access redesign | independent side track | The old database toolchain is removed in this release; see [data-access-implementation.md](data-access-implementation.md) for public URL reading and the next-release source-family MCP connector plan. |
| Multi-agent orchestration | step 6 | freeze the Orchestration contract; hard budgets, store tool, human checkpoints, and run view on omp; Co-Scientist–style reference plugin |

## 5. Contract freeze schedule

| Contract | Frozen at | Layer |
|---|---|---|
| Environment | 1.1 | L1 |
| Execution | 2.1 | L2 |
| Skill / `skill_run` / script tools | 3.1 | L3 |
| Agent definition | 4.2 | L4 |
| Artifact | 4.5 | L4 |
| Plugin | 6.1 | L5 |
| Package, enablement state | step 7 batches 1–2 | L6 |
| Wrapper, Connector | step 7 batch 3 | L6 |
| Core fetch tool | data access track | side |
| Orchestration | orchestration track | side |

## 6. First iteration

All of step 0, then 1.1 → 1.2 → 1.3 → 1.4 → 1.5 → 1.6 → 1.7 (about two and a half weeks).

Outcome: a clean installer and a dependency inventory; Phi has its own
runtime and can build read-only, reproducible, collectable environments from
locks on a clean account without reading the user's conda configuration.
