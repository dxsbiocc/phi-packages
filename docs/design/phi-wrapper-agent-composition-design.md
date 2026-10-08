# Agent-Facing Nextflow Wrappers

Date: 2026-09-15
Status: **Decision draft — align implementation to this shape before adding composition code.**

## 1. Decision

Phi should not compose wrappers by chaining multiple Phi runs outside
Nextflow. It should expose only agent-facing **wrappers** that are runnable
Nextflow entrypoints, and it should let Nextflow own execution, monitoring,
retry, caching, logs, resources, and scheduler integration.

The core shape is:

```text
modules/<provider>/<tool>/
  main.nf
  meta.yml
  environment.yml
  tests/
  wrapper/
    main.nf
    params.json
    wrapper.yaml

subworkflows/<provider>/<workflow>/
  main.nf
  meta.yml
  tests/
  wrapper/
    main.nf
    params.json
    wrapper.yaml
```

`modules/` and `subworkflows/` stay ordinary Nextflow/nf-core assets.
`wrapper/` is the thin adapter Phi exposes to the agent. Phi does not expose
raw modules as tools, because modules are `include` targets, not standalone
pipelines.

## 2. Why This Shape

The agent needs a small, reliable calling surface. A wrapper provides that:

- `wrapper/main.nf` turns agent-friendly params into the channel/module or
  subworkflow call.
- `wrapper/params.json` is a real default parameter file that runs on tiny
  test data with no edits.
- `wrapper/wrapper.yaml` tells the agent which params to replace and where
  outputs land.

This keeps the agent in its comfort zone: discover a wrapper, inspect the
minimal parameter contract, override input/output params, and ask Phi to
plan the run. If the user wants parameter tuning, debugging, or composition,
the agent can load source and tests progressively.

## 3. Wrapper Contract

`wrapper.yaml` is intentionally small. It is not a full workflow manifest,
not a replacement for `meta.yml`, and not a Nextflow schema. It only
describes the agent-facing parameters and declared outputs.

```yaml
id: nf-core/modules/fastqc
name: FastQC
summary: Run FastQC quality checks on FASTQ files.

params:
  reads:
    kind: input
    type: path_glob
    required: true
    description: FASTQ files to analyze.
  outdir:
    kind: output
    type: path
    required: true
    description: Output directory.
  threads:
    kind: option
    type: integer
    required: false
    minimum: 1
    maximum: 16
    description: CPU threads.

outputs:
  reports:
    type: directory
    path: '${outdir}'
    primary: true
```

Rules:

- Use fixed filenames: `wrapper/main.nf`, `wrapper/params.json`,
  `wrapper/wrapper.yaml`.
- Do not put `entrypoint`, `paramsFile`, defaults, environment, or agent
  hints in `wrapper.yaml`.
- Default values live only in `wrapper/params.json`.
- `params` uses one table for inputs, outputs, and ordinary options.
- `kind` is only `input`, `output`, or `option`.
- Requiredness is `required: true`; it is not a separate kind.
- Outputs do not need `required`; `primary: true` outputs must exist after a
  successful run.
- Environment and detailed tool metadata stay in `environment.yml`,
  `meta.yml`, and Nextflow config.

## 4. Discovery And Loading

Phi should scan only:

```text
modules/**/wrapper/wrapper.yaml
subworkflows/**/wrapper/wrapper.yaml
```

If a module or subworkflow has no `wrapper/`, it is not directly callable by
the agent.

Agent tools should be generic and progressively loaded:

- `wrapper_search`: returns compact results only: id, name, summary, key
  input/output params, and primary output.
- `wrapper_inspect`: returns `wrapper.yaml` plus default
  `wrapper/params.json`.
- `wrapper.plan_run`: accepts param overrides, merges them with
  `wrapper/params.json`, validates the merged params, and creates a plan.
- `wrapper.load_source`: optional, used only for composition or debugging;
  can load `wrapper/main.nf`, component `main.nf`, `tests/`, `meta.yml`, or
  `environment.yml`.

Do not register one LLM tool per wrapper. The agent should use the generic
tools, like skills: search first, inspect only selected wrappers, and load
source only when needed.

## 5. Execution

Every wrapper must be smoke-testable with the fixed command:

```bash
nextflow run wrapper/main.nf -params-file wrapper/params.json
```

`wrapper/params.json` must use real tiny test data, not placeholders. That
file is the default run configuration. For real data, the agent usually
overrides only `kind: input` params and commonly overrides `kind: output`
params; `kind: option` params stay at their defaults unless the user asks
for tuning.

Phi validation stays light:

- required params are present after merging defaults and overrides
- simple type checks: string, integer, boolean, path, path_glob, file,
  directory, html, json, csv, fastq_glob, bam, bai, fasta, gtf
- enum/range checks when declared
- local path/glob existence where applicable
- output path can be created or written
- primary declared outputs exist after a successful run

Phi should treat `wrapper/main.nf` as a black-box Nextflow entrypoint in the
first implementation. Correct channel semantics are proven by the wrapper
smoke test and by Nextflow at runtime, not by a Phi-side Nextflow parser.

## 6. Composition

When the agent needs multiple tools, it should generate a temporary
subworkflow component with the same layout:

```text
plans/<planId>/subworkflows/agent-generated/
  main.nf
  wrapper/
    main.nf
    params.json
    wrapper.yaml
```

The temporary component does not write back into the global workspace unless
the user explicitly saves or curates it later.

Composition rules:

- The agent first searches and inspects only the selected wrappers.
- If needed, it uses `wrapper.load_source` for the selected components'
  source/tests.
- The generated subworkflow should include the original modules or
  subworkflows, not other wrappers.
- Wrapper adapters should not be nested inside other wrapper adapters.
- The temporary wrapper still runs as a single Nextflow run via
  `-params-file`, so Nextflow owns monitoring and orchestration.

This preserves a single user approval surface for the generated run plan,
instead of asking the user to approve a chain of separate Phi runs.

## 7. Tests As Usage Evidence

Component `tests/` remain the source of truth for module/subworkflow
correctness and calling patterns. The wrapper adds one extra agent-facing
smoke path:

```bash
nextflow run wrapper/main.nf -params-file wrapper/params.json
```

That proves the default params, wrapper adapter, and declared outputs work
together. The agent does not need to read tests for normal execution; tests
are loaded only for composition, debugging, or deeper parameter work.

## 8. Implementation Implications

This replaces the earlier plan to orchestrate multiple independent wrapper
runs from Phi. The implementation should move toward:

- generic wrapper tools instead of eager `wrapper_<id>` execute tools
- scanning `wrapper/` directories under `modules/` and `subworkflows/`
- plan creation from overrides merged into `wrapper/params.json`
- a minimal `wrapper.yaml` parser/validator separate from the existing full
  `wrapper.yaml` package manifest, or a migration of the existing manifest
  toward this smaller per-wrapper contract
- wrapper smoke validation before a wrapper is made available to the agent
- temporary composition plans that create a `subworkflows/agent-generated`
  component with its own `wrapper/`

Existing bundled wrappers can be migrated incrementally. The next slice
should be small: one module wrapper using this layout, indexed through the
generic search/inspect/plan path, with a smoke test proving the fixed
Nextflow command works.

## 9. Phi Agents (the `Wrapper` agent)

The generic tools in section 4 are **not** handed to the main agent. Wrapper work
belongs to a specialist agent, `Wrapper`, and the main agent is the leader that
delegates to it. This is a general Phi mechanism, not a wrapper special case.

### Naming

- **Agents** are capitalised and never carry an `Agent` suffix: `Wrapper`,
  `CodeReviewer`. The name is enforced (`^[A-Z][A-Za-z0-9]*$`, not ending in
  `Agent`) and must equal the file name.
- **Tool functions** stay snake_case (`wrapper_search`, `read`, `bash`). The
  two are never confused: an agent is invoked by its name, a function by its.

### Definition and scanning (`src/main/agent/agents/`)

An agent is a Markdown file: frontmatter plus the system prompt as the body.

```markdown
---
name: Wrapper
description: One sentence the main agent sees.
tools: [read, glob, grep, bash, write, edit, wrapper_search, wrapper_inspect, wrapper_run]
skills: [create-wrapper, nextflow]
delegation: |
  Guidance for the main agent: when to hand work over, what not to do itself.
---
You are Wrapper, ...
```

Phi scans **its own** locations first; the SDK's `.omp/agents` mechanism is not
used. Legacy layouts are then read for compatibility only, the same way Phi
skills also honour legacy config directories. The first definition of a name
wins, so a Phi agent always beats a compat one.

| Order | Location | Kind |
| --- | --- | --- |
| 1 | `<project>/.phi/agents/` | Phi |
| 2 | `~/.phi/agents/` | Phi |
| 3 | bundled `resources/agents/` (`Wrapper.md`) | Phi |
| 4 | `<project>/{.omp,.pi,.claude}/agents/` | compat |
| 5 | `~/.omp/agent/agents/`, `~/.pi/agent/agents/`, `~/.claude/agents/` | compat |

Compat definitions are normalised rather than rejected: the name becomes
`CodeReviewer` (`code-reviewer`, `planner-agent` → `Planner`), a comma-separated
`tools:` string is accepted, foreign tool names map onto Phi built-ins
(`Read`→`read`, `MultiEdit`→`edit`, `WebSearch`→`web_search`; the rest are
dropped), and a missing `tools:` defaults to the standard file/shell toolbox.
An invalid file becomes a diagnostic (logged as `agent_definition_invalid`) and
never blocks the others. Codex `.toml` agents are not read (different schema).

### How the main agent leads

The **main process scans once per session** and uses that single result twice,
so the prompt and the tools cannot disagree:

1. `buildAgentLeaderPrompt` appends an `<phi_agents>` block to the main system
   prompt (beside the notebook and DB-connector runtime prompts): delegate work
   in a specialist's remit, write self-contained tasks (the specialist cannot see
   the conversation or ask questions), relay reports faithfully, then continue
   downstream. Each agent's `description` and `delegation` are listed. It never
   names a specialist's own tool functions.
2. The definitions are sent to the worker (`session.create` → `phiAgents`), which
   exposes **one delegation tool per agent, named after it** (`Wrapper`), with
   `loadMode: 'essential'` (custom tools default to `discoverable`, hidden behind
   `read xd://`) and `approval: 'read'`.

### How an agent runs

- **Own session.** `createPhiAgentSession` starts an in-memory session with the
  parent's model and thinking level. The definition drives it: the body is the
  system prompt, `tools` the restricted toolbox (`restrictToolNames` +
  `allowRestrictedCustomTools`; no MCP, LSP, web or task), `skills` the only
  skills exposed. Phi tool functions in `tools:` are resolved in the worker
  (`phiToolFunctions`) and are never given to the main agent, so wrapper
  catalogs and Nextflow logs stay out of the main conversation; only the short
  final report returns.
- **Approvals unchanged.** The session's approval extension is bound to the
  *parent's* `sessionId`, so its shell and file writes go through the same
  approval flow, in the chat the user is watching.
- **Cancel / time limit / progress.** The parent tool call's abort signal aborts
  the session; a 3 hour wall-clock backstop applies; each tool call the agent makes
  is streamed to the tool card as a one-line update.

The SDK's built-in `task` sub-agents cannot host this: custom tools registered
through `customTools` are not inherited by them, and they are discovered from
`.omp/agents`. They are left untouched.

Creating or changing a wrapper edits `resources/wrappers/` and `tests/` in the
Phi source tree, so it only works from a Phi checkout, not from a packaged app.

### Runs are background jobs

A wrapper run can take hours, so `wrapper_run` does not block. The agent starts a
run and gets a run id back at once; the main conversation is never held up.

- **Owner.** `WrapperJobManager` (`composition/job-manager.ts`) is a singleton in
  the **main process**. It cannot live in the agent worker: the bridge stops the
  worker whenever it is idle, which would kill the run with it. The worker's tools
  reach the manager through the existing host-request channel
  (`createHostJobClient` → `wrapperJob.*` → `wrapperJobHostHandlers`, which
  validates what arrives from the worker).
- **Tools** (the `Wrapper` agent's, never the main agent's):

  | Tool | Behaviour |
  | --- | --- |
  | `wrapper_run` | Validates, records, starts Nextflow, returns the run id at once. |
  | `wrapper_status` | State, progress, outputs and log tail of one run; without an id, the recent runs. |
  | `wrapper_wait` | Blocks until the run ends or `timeout_seconds` (max 600); aborting the call stops the wait, never the run. |
  | `wrapper_cancel` | Stops the run. |

  By default `Wrapper` does not wait: it reports the run id, output directory and
  "running in the background". It waits (repeated `wrapper_wait`) only when the task
  needs the result ("run it, then summarise…"). The leader learns this from the
  agent's `delegation` text and follows up with "report the status of run <id>".
- **Progress.** Parsed from Nextflow's own console output (Nextflow 26
  `[PROCESS ab/123456] NAME` and the classic `Submitted process > NAME`):
  `started` = distinct processes begun, `total` = distinct process nodes in
  `wrapper/dag.mmd`. It is an estimate of how far the run has got, not a completion
  percentage. Persisted in `run.json` (throttled) and shown in the Wrappers view.
- **Log.** The full output is appended to `runs/<runId>/nextflow.log`, so status
  is answerable from disk after a restart.
- **Cancelling stops Nextflow.** The process is started in its own process group;
  cancel sends SIGTERM to the whole group and SIGKILL after a grace period.
  Verified against real Nextflow + Docker: the run is `cancelled` within
  milliseconds and no container is left behind. Quitting the app stops every live
  run and records it `cancelled`.
- **Recorded in the run store**, the one the Wrappers view reads
  (`composition/run-record.ts`): `origin: 'composition'`, `actor: 'agent'`,
  `planId: ''`, `wrapper.canonicalId` = the composition id. `params.json`,
  `outputs.json` and `summary.json` (with the manifest) land in the run
  directory, so reproducibility export works.
- **States:** `running` → `cancelling` → `cancelled`, or `completed` | `failed`
  (non-zero exit, or a primary output is missing). A run still non-terminal at the
  next app start belonged to a process that is gone and is marked `lost` (outcome
  unknown, so not `failed`).
- **Limits:** at most 3 runs at once.
- **UI.** The main process broadcasts `wrappers:runsChanged`; the run history
  reloads on it (coalesced), shows `2/6 步 · HISAT2_ALIGN`, and offers 取消运行
  for a running background run.

### When a run ends

`WrapperJobManager.onFinish` fires once for a run that ends on its own
(completed, failed, cancelled) — not on app shutdown, when nobody is left to tell.
The main process (`composition/job-notify.ts`) then does two independent things:

1. **A notice in the conversation that started the run.** Every run is stamped with
   the runtime session whose `Wrapper` agent started it (`originSessionId`, set by
   the worker's job client). The main process registers each runtime session's
   conversation when the session is created, appends a `wrapper_run_finished`
   event to that conversation's timeline (persisted, so it is still there when the
   user switches back) and pushes it to the window. The renderer shows it as an
   info banner: outcome, wrapper, elapsed time, output directory, run id, and for a
   failure the cause and "ask Wrapper to look at the log". Live and restored
   history share one function (`chatItems.ts`), and the event does not touch the
   conversation's run status.
2. **An OS notification, only when Phi is not in the foreground.** A run the user
   just cancelled produces no pop-up (they know), but is still recorded in the
   conversation. A run whose conversation is unknown still gets the notification.

A failure in one channel never blocks the other.

#### Waking the conversation

Telling the user is not enough for a chain like "run fastqc, then plot the results":
the agent has to carry on when the run ends. So, by default, Phi also **wakes the
conversation that started the run**: it submits a message to the agent, wrapped in
`<phi_wrapper_run_finished>` and labelled as coming from Phi, not the user. It carries
the outcome, the exit code, the output directory and, for a failure, the missing
outputs. The leader learns this from the `Wrapper` agent's `delegation` text: after
delegating it should say the run is going and end its turn, with no waiting or
polling.

- **Core.** `agent:prompt` was split: `submitPromptRun` runs a prompt in *any*
  conversation, by session key, without touching which conversation is on screen; the
  IPC handler is a thin wrapper around it. Waking uses the same path, so the run gets
  the usual run tracking, approvals, model settings and events, and a conversation in
  the background is woken in the background (it is never switched to).
- **What the user sees.** No user bubble: the message is not the user's words, so no
  `user_message` event is recorded (when the timeline holds text, history restore skips
  the runtime user message that duplicates it). The `wrapper_run_finished` banner is
  already there, followed by the agent's reply.
- **Busy conversation.** A run that ends while the conversation is mid-turn is queued,
  never injected. When that turn ends the queue is flushed as **one** message covering
  every run that ended meanwhile. If the user pressed stop, the queue is dropped.
- **Already reported.** A run whose outcome `wrapper_wait` or `wrapper_status` already
  handed to the agent (`WrapperJobManager.hasBeenReported`) is left out, so a result the
  agent just reported is not announced twice.
- **Who is not woken:** a cancelled run (the user did it), a run started with
  `continue_when_done: false` (the `Wrapper` agent uses this when the task says nothing
  should happen afterwards), and a conversation Phi no longer knows.
- **Loop guard.** At most `MAX_AUTOMATIC_CONTINUATIONS` (5) wake-ups in a row per
  conversation; a real user message starts the count over. Without it an agent that
  starts a run every time it is woken would run unattended indefinitely.
- **Policy and message** are in `composition/job-continue.ts`; the queue and the call
  into `submitPromptRun` are in `index.ts`.

Not built yet: the sidebar does not mark the conversation unread. Local runs do not
survive quitting the app; remote runs do (next section).

### Running on an HPC cluster (`target: "remote"`)

`wrapper_run` takes `target: "local" | "remote"`. A remote run goes to the **HPC
connection saved on the project** the chat belongs to (Wrappers page, remote settings)
and uses the same job manager, run records, progress, notifications and wake-ups as a
local run: `remote-job.ts` returns the same `WrapperProcess` the local runner does.

- **Controller** (`composition/remote-controller.ts`; the connection's "Nextflow 主进程运行位置"). Two ways to host the Nextflow head process; everything after the launch (log streaming, outputs, reconnecting, cancel, reattach) is shared.
  - `login` (default): started on the login node under `setsid` (`detached_ssh`), identified by its pid. With the `slurm` scheduler it then submits every task to Slurm itself.
  - `sbatch` (`slurm-controller`): the head process is a Slurm job of its own (`job.sbatch`), identified by its job id, for sites that forbid long-lived processes on login nodes. It costs a queue wait before anything starts, and the cluster must allow a job to submit jobs (nested `sbatch`).
  Phi only polls over SSH, so a dropped link or a closed app does not touch the run.
- **Where things live** (`composition/remote-config.ts`): under the project's remote
  workspace root, `wrappers/bundles/<hash>/` holds the wrapper source tree and
  `wrappers/runs/<runId>/` holds `params.json`, `phi_remote.config`, `launch.sh`, `logs/`,
  `work/` and `results/`.
- **Source bundle** (`remote-bundle.ts`). A wrapper includes sources by relative path (a
  subworkflow reaches `../../../../modules/...`), so the whole `resources/wrappers` tree
  ships, content-addressed: uploaded once per version over SFTP, shared by every run, a
  changed tree lands beside the old one. `tests/` is left out except `tests/data` for a
  component whose default params point at it.
- **Run config.** `phi_remote.config` is passed with `-c`, so it outranks the wrapper's own
  config without touching any wrapper: `process.executor`, queue, `--account`, extra
  `sbatch` flags, `executor.queueSize`, `singularity.cacheDir`. The profile is the
  connection's runtime (`singularity` by default), which every wrapper already defines.
- **Inputs are cluster paths.** Local paths are meaningless there. Phi checks each
  `kind: input` exists on the cluster (a glob by its fixed directory) before launching and
  does not move data. `outdir` defaults to the run directory's `results/`; outputs stay on
  the cluster and are recorded with `location: 'remote'`.
- **Cluster check.** Before uploading, Phi runs a short script over SSH: it loads the login
  profile (so `module load` works), runs the connection's setup commands, and fails with the
  setting to change if Nextflow (or `sbatch` for Slurm) is missing. A missing container
  runtime only warns, since some sites provide it on compute nodes only. Failure reasons are
  written to the run's `nextflow.log`, which is what `wrapper_status` shows the agent.
- **A connection without HPC settings is refused**, so a run never falls back to executing
  every task on the login node.
- **Quitting and restarting.** On quit a remote run is *detached*, not cancelled; its record
  stays `running` with a `remote.json` snapshot (pid, run directory, log offset). On the next
  start `WrapperJobManager.adoptRemoteRuns` reattaches, replays the local log to rebuild
  progress, and carries on from the saved offset. A run that cannot be reattached is `lost`
  (outcome unknown), never `failed`; likewise when contact is lost for ~12 polls in a row
  (each preceded by a reconnect attempt).
- **Cancel** signals the head process (SIGTERM to the process group, or `scancel` for a Slurm
  head job); Nextflow then cancels its own Slurm jobs. SIGKILL after a grace period.
- **The head job's allocation** (`sbatch` controller). Defaults are 1 CPU, 4G and a 2-day time
  limit, using the connection's queue and account; the connection's head-job options
  (`controllerOptions`, e.g. `--time=7-00:00:00 --mem=8G`) come last and override them. The
  time limit is set on purpose: a partition default of an hour or two would silently kill a
  pipeline that is still running. A job Slurm ends for time or memory is reported with that
  state and the flag to change, not just "failed".
- **Who decides success under `sbatch`.** `launch.sh` records Nextflow's exit code in
  `exit_code` and exits with it, and that file outranks the scheduler's view (a Nextflow
  failure would otherwise be a `COMPLETED` job). Only when the job died without recording one
  (time limit, out of memory, cancelled) does the scheduler's state decide, read as `squeue`,
  then `scontrol show job`, then `sacct` (many clusters have no working `sacct`).

Verified with real Nextflow and Docker against a "remote" that is a local bash session
(modules, subworkflows with cross-directory includes, a two-run reattach), including the
`sbatch` controller against a stand-in Slurm (`tests/helpers/fakeSlurm.ts`) that really runs
the submitted script in the background and can really cancel it; unit tests cover the rest.
**Not verified against an external SSH host or a real Slurm cluster**: real `sbatch`/`scontrol`
output, nested submission, and Nextflow's own Slurm executor rest on the stand-in and on the
generated config having been accepted by `nextflow config`. The `squeue → scontrol → sacct`
status order was carried over from the older `SbatchRunner`, which was tried on a real cluster.

Since the original composition implementation, the remote project work added a
Phi-owned OpenSSH host profile and connection test, same-name remote file/command
tools, bounded remote result browsing and previews, and explicit result download.
These paths have automated coverage but still need external-host acceptance. Remote
input files are not uploaded automatically. PBS/LSF and cleanup of orphaned Slurm
task jobs after a hard-killed Nextflow head remain outside this implementation.
The validation record is in `docs/roadmap/remote-e01-validation.md`.
