# Agent-definition contract

contractVersion: 1.0.0

This contract defines a Phi agent: a specialist the main agent delegates to. It
covers the file, its frontmatter, and what the runtime does with each field. It
builds on the [environment contract](environment.md), the
[execution contract](execution.md), and the [skill contract](skill.md). See
[runtime foundation](../design/phi-runtime-foundation.md) §7 and §10, and
decisions 10 and 14 in the
[content distribution decision record](../decisions/content-distribution.md).

## 1. File

One Markdown file per agent: YAML frontmatter, then the system prompt as the body.
The file name without `.md` equals `name`. Phi scans its own agent directories
(bundled `resources/agents/`, and plugin `agents/` directories from step 6); files
from other tools' agent directories are read in **compatibility mode** (§ 6).

## 2. Fields

Field names follow omp where omp has the field, and camelCase everywhere
(runtime foundation §10), so a Phi agent file stays loadable by omp's `task`.

### 2.1 Fields shared with omp

| Field           | Required | Rule                                                                                                                                                                                                                |
| --------------- | -------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `name`          | yes      | `^[A-Z][A-Za-z0-9]*$`, must not end in `Agent`; equals the file name. It is also the name of the main agent's delegation tool                                                                                       |
| `description`   | yes      | 1–1024 characters; what the specialist does and when the main agent should use it                                                                                                                                   |
| `tools`         | yes      | non-empty list of tool names: SDK built-ins (`read`, `bash`, …), Phi core tools (`skill_run`, …), and Phi tool functions. Names the runtime does not provide are dropped with a warning when the session is created |
| `spawns`        | no       | **reserved** for orchestration (decision 10). v1 rejects the field                                                                                                                                                  |
| `model`         | no       | a model selector as in omp (`provider/id`, or a list tried in order). Absent: the delegating session's model. Unresolvable: the delegating session's model, with a warning                                          |
| `thinkingLevel` | no       | `off`, `minimal`, `low`, `medium`, `high`, or `xhigh`. Absent: the delegating session's level                                                                                                                       |

### 2.2 Phi fields

| Field            | Required | Rule                                                                                                                                                                                                                                                       |
| ---------------- | -------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `skills`         | no       | list of skill names the agent may use; default none. Skills that are not installed and enabled are dropped with a warning                                                                                                                                  |
| `environment`    | no       | an environment reference (environment contract): `phi:<name>@<major>`, `plugin:<name>`, or `project:<name>`. See § 3                                                                                                                                       |
| `visibility`     | no       | `entry` (default): the main agent gets a delegation tool for it. **`internal` is reserved** for orchestration and rejected in v1                                                                                                                           |
| `delegationMode` | no       | `required-first` (matching work goes to this agent before generic tools), `preferred`, or `optional` (default)                                                                                                                                             |
| `delegation`     | no       | guidance for the main agent on when and how to hand work over; added to the main agent's instructions                                                                                                                                                      |
| `fallback`       | no       | `{ afterFailures: <positive integer>, tools: [<main-agent tool>…], match: [<token>…] }`: after that many consecutive failed, blocked, or not-found runs, the main agent may use `tools` for inputs containing one of the `match` tokens (case-insensitive) |
| `outputSchema`   | no       | **reserved** for structured results. v1 rejects the field; a later minor version defines it                                                                                                                                                                |

### 2.3 Other keys

Legacy snake_case spellings are read as aliases and reported as a **warning** naming
the new spelling: `delegation_mode` → `delegationMode`, `fallback.after_failures` →
`fallback.afterFailures`. Any other unknown key is a **warning** and ignored.
Anything in §§ 2.1–2.2 that breaks its rule is an **error**.

## 3. Environment binding

When `environment` is set, the agent's sessions are **bound** to that environment:

1. **Session creation.** Before the session starts, the environment is resolved and
   must be ready. If it is not built, the user is asked to build it, with the same
   prompt and progress as `skill_run` (skill contract § 3.3). If the user declines
   or the build fails, the delegation fails with `environment <ref> is not ready`.
   The session never falls back to the host.
2. **`bash`.** Every `bash` call of the session runs with the environment's
   variables as defined by the execution contract: host variables outside the
   execution contract's allowlist are removed, and `PATH` and the isolation
   variables come from the environment. Variables the call sets explicitly are
   kept, except the names the execution contract reserves (`PATH`, `PYTHON*`,
   `R_*`, `CONDA_*`, `MAMBA_*`, `LD_*`, `DYLD_*`, `PHI_ENV_*`), which the binding
   always sets.
3. **Skills.** `skill_run` and script tools called by this session use the agent's
   environment as the session environment (skill contract § 3.3 step 2). A skill's
   own `phi.environment` still wins.
4. **Scope.** The binding applies only to this agent's sessions. The main agent's
   `bash` always uses the host.

Without `environment`, the agent's `bash` uses the host, like the main agent.

## 4. Results

The agent's final message is its report to the main agent. The runtime appends
Phi's report protocol to the system prompt: the report ends with a status
(`done`, `blocked`, `not-found`, or `failed`) and may name missing inputs or a
next agent. The protocol text belongs to the runtime, not to this contract.

## 5. Validation

`validateAgent(filePath)` checks everything above and returns all problems, split
into `errors` and `warnings`; an agent is valid when `errors` is empty. It runs in
`npm run lint` for every agent in the repository, and when agents are scanned: an
invalid agent is skipped with its errors in the log, and the others still load.

## 6. Compatibility mode

Agent files from other tools (for example `~/.claude/agents`) are read leniently:
the name is derived from `name` or the file name (`code-reviewer` → `CodeReviewer`);
`tools` is mapped to Phi built-ins where an equivalent exists (unknown names are
dropped) and defaults to `read`, `glob`, `grep`, `bash`, `write`, `edit`; every Phi
field is ignored. Compatibility agents are never bound to an environment.

## 7. Versioning

`contractVersion` follows the content distribution design §4.4: minor versions are
additive only (new optional fields, or allowing a reserved field or value); anything
else needs a decision record and a deprecation window.
