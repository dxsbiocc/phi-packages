# Phi Content Distribution Decisions

Date: 2026-09-29
Status: Accepted direction. Open questions are tracked in the design's §16.
Design: [phi-content-distribution-design.md](../design/phi-content-distribution-design.md)
(Chinese: [phi-content-distribution-design.zh-CN.md](../design/phi-content-distribution-design.zh-CN.md)).
Chinese version of this record: [content-distribution.zh-CN.md](content-distribution.zh-CN.md).

## Context

Phi ships all domain capabilities inside the app bundle (`resources/skills`,
`resources/wrappers`, `resources/db-connectors`, `resources/agents`) and loads
them for every user. Content cannot be updated without an app release, every
skill description sits in the main agent's prompt, and script dependencies
are whatever the user's `PATH` happens to provide. Domain logic has also
grown inside the engine (`src/main/agent/db/`, `src/main/agent/visualization/`).
An evaluation on 2026-09-29 showed the `db_*` toolchain was no more accurate
than plain URL fetching and about twice as slow (design Appendix A).

## Decisions

1. **Engine and content are separate.** The engine (kernel) holds the agent
   runtime, generic core tools, installer, environment manager, approvals,
   remote execution, and built-in skills such as `create-wrapper`. It holds no
   domain logic. Everything domain-specific is content.
2. **Kernel contracts first, then freeze.** Package, Skill, Wrapper,
   Connector, Agent definition, Plugin, Environment, Core services, Artifact,
   and (reserved) Orchestration are written contracts with JSON schemas, one
   validator (`phi validate`), and conformance tests. Each is validated by a
   reference implementation and frozen as v1 before its first consumer is
   built. Afterwards minor versions are
   additive only; a breaking change needs an ADR and a deprecation window of
   at least two app releases.
3. **Units and plugins.** Single tools are standalone units: skill, wrapper
   tool family, or MCP connector. A plugin exists only for composite
   capability: several interdependent component types, a dedicated
   environment, or a dedicated specialist agent. Units and plugins share one
   package format, one installer, and one registry.
4. **Local registry first, remote later.** The registry is generated from
   `resources/` until the installer, loader, and catalogs are stable; a signed
   remote registry on static storage follows after the beta.
5. **Allowlist packaging.** Packages are built from git-tracked files or
   manifests, never by copying directories, with a size budget in CI.
6. **Routing by type.** The main agent sees specialist entries and a few
   knowledge skills; wrappers go to `wrapper_agent`, API skills to their
   specialist. Only installed, enabled content reaches the runtime.
7. **Managed environments.** Phi bundles micromamba. Official script skills
   share a locked, developer-maintained `phi-python`, checked in CI; user skills
   get an environment chosen by dry-run, or a new content-addressed,
   immutable one. Scripts run through `skill_run`.
8. **Data access via API skills.** The `db_*` toolchain is replaced by a core
   fetch tool (retry, rate limits, allowlist and audit, pagination and bulk
   download, result viewers) plus about six domain API skills. MCP is used for
   authenticated or private sources and existing official servers. The old
   toolchain retires only after an `api-skill` arm is at least as good on the
   extended evaluation. The `create-database-connector` skill is deleted
   immediately.
9. **No domain tools in the engine.** Domain logic becomes command-line
   programs inside plugins or skills, exposed as typed tools by declaring
   script tools in SKILL.md; the engine only registers, validates, runs them in
   their environment, and presents typed artifacts. Visualization migrates
   first as `scripts/viz.py`, keeping the tool names `viz_examples`,
   `viz_route`, `viz_prepare`, `viz_render`.
10. **Multi-agent orchestration builds on omp.** Orchestration is reserved in
    v1 (plugin `orchestrator` slot, agent `visibility`, `orchestration.*`
    namespace). It uses omp's `task`, `spawns`, agent registry, and `hub`
    through a thin Phi adapter; Phi adds hard budgets, a typed run store,
    approvals, human checkpoints, isolation by package, and the run view.
    Agent-driven orchestration comes first; declarative and programmatic
    levels later.
11. **Pi plugins are developer extensions.** The existing pi plugin page is
    renamed and moved to advanced settings; "plugin" in the UI means Phi
    plugins.
12. **Runtime foundation first** ([phi-runtime-foundation.md](../design/phi-runtime-foundation.md), confirmed 2026-09-29):
    only the main agent's `bash` uses the host environment; all other
    execution runs in managed environments. Official environments install from
    explicit locks without client-side solving. Environment prefixes are
    read-only; extra dependencies create a new project environment. LibreOffice,
    Docker, and Singularity are host dependencies (checked, not installed).
    Nextflow and Jupyter default to managed environments; notebook kernels and
    Nextflow may explicitly use a host version (Nextflow after a version
    check), labelled unmanaged, while content execution always stays managed. Only the current
    platform's micromamba is bundled; remote hosts get theirs on demand. Work is
    built inside-out: runtime → execution → consumers → binding → plugins →
    distribution.
13. **One naming convention** for tools, environments, references, manifests,
    fields, and directories, defined in the runtime foundation §10: core tools
    `<domain>_<verb>` with reserved prefixes; script tools
    `<toolPrefix>_<name>`; environments `phi-<purpose>` referenced as
    `phi:<name>@<major>`, `plugin:<name>`, `project:default`; one
    `phi-package.yaml` for every package type; camelCase fields aligned with
    omp; all packages under `~/.phi/packages/<type>/<id>/<version>/`.

14. **Specialist delegation stays on Phi's runtime; omp `task` is for
    orchestration** (step 4.1 spike, 2026-09-30, omp 18.1.10). A spike ran a
    real omp `task` subagent from a Phi-shaped parent session. Findings:
    - A parent's in-process `customTools` (how Phi builds `viz_*`, wrapper,
      database, and script tools, each closing over `requestHost`) are **not**
      visible to the subagent ("tool not found").
    - A tool registered by an **extension** (`pi.registerTool`) **is** visible
      and runs in the parent process, closures intact. A parent extension's
      `tool_call` hook also sees every subagent call. Both hold only while the
      parent session does not set `restrictToolNames`, because omp forwards
      prepared extensions only then.
    - Subagents run with `tools.approvalMode: "yolo"`. Phi's ask-mode approvals
      still apply, but only through that forwarded extension hook, and the hook
      loses Phi's run context (`agentRunId`, the chat card).
    - Agents are discovered only from `.omp/agents` directories and omp
      plugins, not from Phi's scanned definitions.
    - omp has no per-agent environment binding.

    Moving Phi's three specialists onto `task` would therefore mean
    re-registering every Phi tool through an extension and scoping it per
    agent. It would also mean rebuilding what `agents/registry.ts` and the UI
    already provide: per-conversation limits, background runs with chat cards
    and the jobs panel, steer, the report protocol, fallback policy, usage
    logs, remote-project guards, and parent-bound approvals. That is a
    rewrite with beta risk and little user-visible gain.

    Decision: the beta keeps Phi's delegation runtime, and nothing in
    `agents/registry.ts` retires now. The Agent-definition contract (4.2) uses
    omp's field names (`name`, `description`, `tools`, `spawns`, `model`,
    `thinkingLevel`) so definitions stay loadable by omp. Orchestration
    (decision 10) uses omp `task` / `hub` with Phi tools provided through a
    tool-registering extension and an unrestricted orchestrator session; its
    adapter must restore run context for approvals.

## Consequences

- The internal beta roadmap changes: catalog install / enable / disable for
  skills, wrappers, and connectors replaces the "read-only Skills and MCP
  pages" rule; the DB connector prototype is frozen pending retirement.
- Work follows the layered plan
  ([content-distribution-implementation.md](../roadmap/content-distribution-implementation.md)),
  derived from the runtime foundation: runtime → execution → consumers →
  binding → plugins → distribution → remote. Steps 0–5 are proposed for the
  beta; plugins, content distribution, and remote / HPC (steps 6–8) follow.
  Data access and orchestration are side tracks. Each contract is frozen when
  its layer is complete.
- The omp spike (step 4.1) is done: specialist delegation stays on Phi's
  runtime (decision 14); `agents/registry.ts` is kept.
- Content authors get one validator and one package format; engine changes
  are gated by the conformance suite.

## Alternatives considered

- **Keep everything bundled** — simplest, but content stays tied to releases
  and routing load grows with every addition.
- **Everything is a plugin** — forces packaging on single tools and user
  content; rejected in favour of units plus composite plugins.
- **Rewrite db connectors as MCP servers** — MCP is a transport; the
  evaluation showed the cost comes from the query abstraction, so MCP is kept
  for authenticated / complex sources only.
- **Phi-native orchestration runtime** — duplicates what omp already
  provides; rejected in favour of an adapter over omp.
- **One shared environment for all content** — cannot satisfy conflicting
  user skills; rejected in favour of `phi-python` plus content-addressed
  environments.
