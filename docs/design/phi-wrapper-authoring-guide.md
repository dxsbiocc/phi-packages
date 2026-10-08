# Phi Wrapper Authoring Guide (Phase 1)

Date: 2026-09-10

This is a practical guide for writing a `wrapper.yaml` and running it locally during
Phase 1. For the product rationale and phased scope, see
[phi-wrapper-product-prd.md](phi-wrapper-product-prd.md); for the full target manifest
shape and how each Phase 1 field maps to later phases, see
[phi-wrapper-technical-design.md](phi-wrapper-technical-design.md). This guide only
covers what Phase 1 actually validates and runs.

## Two Ways A Wrapper Exists In Phase 1

- **Bundled**: ships inside the app, installed automatically at startup. You can't add
  one of these from the UI — it's the app's own bundled wrapper package (see
  `resources/wrappers/`).
- **Custom**: a local folder you point Phi at. Open the Wrappers sidebar, click the
  **+** button (top-left, next to refresh), and choose the folder containing your
  `wrapper.yaml`. Custom wrappers show up in the catalog and can be inspected, but
  never become a default agent tool — the agent only gets tools for bundled wrappers
  (see the technical design's Trust And Registry section for why).

There is no CLI yet (`phi-wrapper install ...` is Phase 2) and no registry (Phase 3) —
Phase 1's only path onto disk is one of the two above.

## Minimal `wrapper.yaml`

```yaml
phiWrapperVersion: 1
id: acme/tools/my-wrapper # namespace/group/short-id — exactly three segments
shortId: my-wrapper
name: My Wrapper
version: 1.0.0
summary: One sentence describing what this runs.

runtime:
  minVersion: 1.0.0
  maxVersion: 1.x

resourceClass: light # light | standard | heavy | hpc — heavy/hpc need explicit
# acknowledgement before submit in Phase 1 (there's no remote to redirect to yet)

engine:
  type: nextflow # only nextflow is executable in Phase 1; other values parse but can't run
  entrypoint: main.nf
  profiles:
    - id: local # a "local" profile is required — Phase 1 has no other executor
      executor: local

inputs:
  - id: reads
    type: fastq_reads
    required: true

parameters:
  schema: # a real JSON Schema — validated with ajv before a plan is ever created
    type: object
    required: [reads]
    properties:
      reads:
        type: string

outputs: # at least one required
  - id: report
    label: Report
    type: html
    path: results/report.html
    primary: true

resources:
  defaults:
    cpus: 4
    memory: 8 GB
```

This is genuinely everything `manifest.ts` requires. Everything below is optional but
worth knowing about.

## Field-By-Field Notes

- **`id`** must match `^[a-z0-9-]+/[a-z0-9-]+/[a-z0-9-]+$` — lowercase, exactly three
  `/`-separated segments. This is the canonical identifier — it's what the sidebar,
  the plan card, and the run history show as the primary label, not `name`. `name` is
  free text you control; nothing about it is verified, so don't rely on it to tell two
  wrappers apart.
- **`version`** must be valid SemVer.
- **`runtime.minVersion`/`maxVersion`** are checked against Phi's own wrapper runtime
  version (`PHI_WRAPPER_RUNTIME_VERSION` in `manifest.ts`, currently `1.0.0`).
  `maxVersion` accepts a `N.x` wildcard form.
- **`engine.profiles`** must include one entry with `id: local`. Declaring `docker`,
  `slurm`, etc. profiles is fine — they're just not selectable yet (Phase 1's plan
  creation always resolves to `local`).
- **`steps`** (optional) declares the workflow DAG the chat plan card and Wrappers
  sidebar render as a diagram — omit it and you get a generic
  inputs → wrapper → outputs fallback instead of a real shape. Each step needs `id`,
  `label`, and `dependsOn` (an array of other step ids, `[]` for a starting step). Make
  `id` match the corresponding Nextflow **process name** — once a run is submitted, the
  executor's weblog listener attributes live per-step state (`running` /
  `completed` / `failed`) by matching Nextflow's process name against this id
  (case-insensitive; a tagged trace name like `"fastqc (sample1)"` matches on the part
  before the parenthesis). An id that doesn't match anything just stays `pending`
  forever — it doesn't error, but you won't get a live view.
- **`inputs[].samplesheet`** (optional): for a `fastq_reads` input with
  `layout: paired_end`, a glob value like `data/*_{R1,R2}.fastq.gz` gets resolved into
  a samplesheet (pairing files by `_R1`/`_R2` or `_1`/`_2` in the filename) and
  persisted alongside the plan. Non-glob values just resolve as a literal path.
- **`outputs[].path`** is resolved relative to the plan's output directory
  (`plan.cwd` + `plan.outputDir`) once a run actually executes — declare it as the
  path your pipeline publishes to under `params.outdir` (the local executor injects
  `outdir` into `params.json` automatically if your manifest doesn't already declare
  it as a parameter).
- **`verification`**, **`registryStatus`**, **`source`** are parsed and preserved but
  never read to decide trust — trust tier comes from _how the wrapper got onto disk_
  (bundled vs. custom), never from anything inside the manifest itself. Don't bother
  filling these in for a Phase 1 custom wrapper; they don't do anything yet.
- **`environment`**, **`summaries`**, **`permissions`**, **`tests`**, **`license`**,
  **`citations`**, **`dataPolicy`** are all parsed and carried through untouched but
  not yet acted on by anything in Phase 1 (no container provisioning, no permission
  enforcement, no test runner). They're safe to include for forward-compatibility with
  Phase 2/3, or safe to omit.

## Running It Locally

1. Add the wrapper as `custom` via the Wrappers sidebar's **+** button.
2. Ask the agent to prepare a plan (only works for `bundled` wrappers by default — see
   below) — or, for now, custom wrappers are inspect-only in chat since they don't
   register as an agent tool. To actually submit a custom wrapper's plan today, that
   path isn't wired into the UI yet; the practical way to exercise a custom wrapper
   during development is to promote it to the fixture directory temporarily and treat
   it as `bundled` for local testing, or drive `createWrapperRunPlan`/
   `submitWrapperRunPlan` directly from a test/script against `src/main/agent/wrappers`.
3. Local execution needs `nextflow` on `PATH` (and `docker` if your chosen profile
   declares `containerRuntime: docker`) — the executor's doctor check
   (`checkLocalDoctor` in `doctor.ts`) fails the run immediately with a clear reason if
   either is missing, rather than attempting a command that will fail confusingly.
4. Run output lands under `~/.phi/wrappers/runs/<runId>/` (params, logs,
   `outputs.json`, `summary.json`) plus your declared `outputDir` under the project
   directory you ran it against.

## Why The Agent Only Gets Tools For Bundled Wrappers

`wrapper_search` and `wrapper_inspect` see the whole catalog, bundled or custom, so
the agent can tell you a custom wrapper exists. But a per-wrapper `wrapper_<id>` tool
— the one that actually creates a run plan — is only generated for `bundled` entries.
This isn't a technical limitation; it's deliberate: a `custom` wrapper is something you
pointed Phi at yourself, with no verification at all, and Phase 1 doesn't yet have a
UI flow for "let me explicitly allow this specific custom wrapper as an agent tool for
this session." Building that flow is reasonable Phase 1/2 follow-up work, not
something to work around by weakening the default.
