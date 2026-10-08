# Skill contract

contractVersion: 1.1.0

This contract defines a Phi skill: its files, its frontmatter, how its scripts run
(`skill_run`), and how it declares **script tools**. It builds on the
[environment contract](environment.md) and the [execution contract](execution.md). See
[runtime foundation](../design/phi-runtime-foundation.md) §5.1, §6, and §10.

## 1. Files

```text
<skill-name>/
  SKILL.md                 # required: frontmatter + instructions
  references/              # optional: documents the agent reads on demand
  scripts/                 # optional: programs run through skill_run or script tools
  assets/                  # optional: templates, examples, data files
  schemas/                 # optional: JSON Schemas for script tool outputs
  environment.yml          # optional: the skill's own environment (see § 3.2)
  locks/<platform>.txt     # required when environment.yml is present
```

The directory name equals `name`. Nothing is executed from `references/` or `assets/`.

## 2. Frontmatter

`SKILL.md` starts with YAML frontmatter. Standard [Agent Skills](https://agentskills.io/specification)
fields keep their names and meaning; Phi adds nothing at the top level. **Everything Phi-specific
lives under one `phi` key**, so the skill stays a valid Agent Skill for other tools and Phi
fields never collide with future standard fields.

### 2.1 Standard fields

| Field                      | Required | Rule                                                      |
| -------------------------- | -------- | --------------------------------------------------------- |
| `name`                     | yes      | `^[a-z0-9][a-z0-9-]{0,63}$`; equals the directory name    |
| `description`              | yes      | 1–1024 characters; what the skill does and when to use it |
| `license`                  | no       | string                                                    |
| `compatibility`            | no       | string                                                    |
| `metadata`                 | no       | map of string to string                                   |
| `allowed-tools`            | no       | string                                                    |
| `disable-model-invocation` | no       | boolean                                                   |
| `hide`                     | no       | boolean                                                   |

Keys already read by the runtime (`globs`, `alwaysApply`) are accepted. Any other top-level
key is reported as a **warning** and ignored by Phi: skills from the wider Agent Skills
ecosystem often carry vendor keys (for example `required_environment_variables`). Problems
inside the `phi` block are **errors**.

### 2.2 The `phi` block

```yaml
phi:
  environment: phi:python@1 # where scripts and script tools run
  attachTo: [main] # who gets this skill's script tools
  toolPrefix: scanpy # standalone skills with script tools only
  scripts: # script tools (§ 5)
    - name: qc
      description: Compute QC metrics for an .h5ad file and write a summary table
      run: [python, ./scripts/qc.py]
      args: { … JSON Schema … }
      approval: read
      output: ./schemas/qc-result.json
```

| Field         | Required                                                    | Rule                                                                                                                                                                                                                         |
| ------------- | ----------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `environment` | when the skill has scripts (§ 3.1)                          | an environment reference: `phi:<name>@<major>`, `plugin:<name>`, `project:<name>`, or `./environment.yml`                                                                                                                    |
| `attachTo`    | no                                                          | list of agent names (PascalCase) and/or `main`; default `[main]`. Script tools are registered for these agents only                                                                                                          |
| `toolPrefix`  | when `scripts` is set and the skill is not part of a plugin | `^[a-z][a-z0-9]{1,11}$`, not an engine-reserved prefix (`skill`, `env`, `http`, `wrapper`, `agent`, `db`, `mcp`); inside a plugin the plugin's `toolPrefix` is used and this field is rejected                               |
| `scripts`     | no                                                          | list of script tool declarations (§ 5)                                                                                                                                                                                       |
| `deprecated`  | no                                                          | 1–300 characters (1.1.0): the skill still works but is being replaced; the message names the replacement. The UI marks it on the Skills page and in the catalog, where it sorts last. Enablement and behaviour are unchanged |

## 3. Environments

### 3.1 When a skill needs one

A skill **has scripts** when `scripts/` contains any file ending in `.py`, `.R`, `.r`, `.sh`,
`.js`, `.mjs`, or `.pl`, or when `phi.scripts` is set. Such a skill must declare
`phi.environment`. A skill without scripts must not declare one.

### 3.2 Skill-local environments

`environment: ./environment.yml` gives the skill its own environment. The file follows the
environment contract, and `locks/<platform>.txt` must exist for every platform Phi ships.
Its envId uses scope `skill` with the skill name as owner:
`skill-<skill name>-<environment name>-<hash12>` (environment contract 1.1.0).

### 3.3 Resolution order (`skill_run` and script tools)

1. the skill's own `phi.environment`;
2. otherwise the environment bound to the calling agent session;
3. otherwise `phi:python@1`, and the result carries the warning
   `skill <name> declares no environment`.

An environment that is not built yet is never replaced by the host: the call fails with
`environment <ref> is not ready` and the UI offers to build it (size and progress shown).

## 4. `skill_run`

A core tool that runs one program from a skill's `scripts/` directory.

```ts
skill_run({ skill: string, script: string, args?: string[], cwd?: string })
```

| Parameter | Rule                                                                                                                |
| --------- | ------------------------------------------------------------------------------------------------------------------- |
| `skill`   | an installed, enabled skill name                                                                                    |
| `script`  | a path relative to the skill's `scripts/` directory; it must resolve inside `scripts/` (no `..`, no symlink escape) |
| `args`    | passed verbatim after the script path                                                                               |
| `cwd`     | project-relative directory, default the project root; must stay inside the project                                  |

The command is `<interpreter> <absolute script path> <args…>`, where the interpreter comes from
the extension: `.py` → `python`, `.R`/`.r` → `Rscript`, `.sh` → `bash`, `.js`/`.mjs` → `node`,
`.pl` → `perl`, each resolved inside the environment by `runInEnvironment`. Anything else is
rejected.

**Approval:** `skill_run` executes code, so in `ask` mode every call needs approval like
`bash`. **Result:** `envId`, `resolvedCommand`, `exitCode`, `terminated`, `stdout` and `stderr`
(each truncated per the execution contract), `durationMs`, and any warning from § 3.3.

## 5. Script tools

A script tool is a typed tool the engine generates from a `phi.scripts` entry. The engine
contains no code for any particular skill.

### 5.1 Declaration

| Field            | Required | Rule                                                                                                                                                                                                                                                                       |
| ---------------- | -------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `name`           | yes      | `^[a-z][a-z0-9_]{0,31}$`, unique within the skill                                                                                                                                                                                                                          |
| `description`    | yes      | 1–1024 characters, written for the model                                                                                                                                                                                                                                   |
| `run`            | yes      | argv array. The first element is a command resolved inside the environment (`python`, `Rscript`, …). Elements starting with `./` are paths relative to the skill directory and must resolve inside it                                                                      |
| `args`           | yes      | JSON Schema with `type: object`, `properties`, optional `required`, `additionalProperties: false`. Property names match `^[a-z][a-z0-9_]*$`. Supported property types: `string`, `number`, `integer`, `boolean`, `array` of those, and `string` with a `format` from § 5.3 |
| `approval`       | yes      | `read` or `write`                                                                                                                                                                                                                                                          |
| `output`         | no       | `./<path>` to a JSON Schema file inside the skill; when present, stdout is validated against it                                                                                                                                                                            |
| `timeoutSeconds` | no       | 1–86400, default 600                                                                                                                                                                                                                                                       |

### 5.2 Registration and naming

The tool is registered as `<toolPrefix>_<name>` for every agent in `attachTo`. Its parameter
schema is `args`. Two tools with the same final name are an installation error.

### 5.3 Path formats

| Format         | Meaning                                                                                                                                |
| -------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| `input-path`   | an existing file or directory; resolved against the call's project directory; its real path must be inside the project                 |
| `project-path` | a file or directory to be written; resolved against the project directory; must be inside the project; its parent directory must exist |

The engine replaces these arguments with absolute paths before running the program. Paths outside
the project are rejected before anything runs.

### 5.4 Invocation

Arguments become flags in the order of `args.properties`:

| Value                   | Passed as                                    |
| ----------------------- | -------------------------------------------- |
| string, number, integer | `--<name> <value>`                           |
| boolean                 | `--<name>` when `true`, nothing when `false` |
| array                   | `--<name> <item>` repeated per item          |
| absent optional         | nothing                                      |

The program runs with `runInEnvironment` in the resolved environment, `cwd` = the project
directory, and `timeoutSeconds`.

### 5.5 Output and errors

- On success (exit code 0) the program writes **one JSON object** to stdout and nothing else;
  diagnostics go to stderr. The engine parses stdout, validates it against `output` when
  declared, and returns it to the model.
- On failure the program exits non-zero and may write `{"error": "<message>"}` to stdout. The
  engine reports that message, or the tail of stderr when there is none.
- An output that is not a single JSON object, or that fails `output` validation, is an error.
- An output may list files it produced under `artifacts: [<project path>, …]`; each needs a
  `<file>.phi-artifact.json` descriptor (artifact contract, step 4 of the implementation plan).
  Until that contract exists, the engine passes `artifacts` through unchanged.

### 5.6 Approval

`read` tools run without approval in `ask` mode only when every path argument is an
`input-path`; `write` tools, or any tool with a `project-path` argument, need approval in `ask`
mode.

## 6. Validation

`validateSkill(dir)` checks everything above and returns all problems, not just the first,
split into `errors` and `warnings`; a skill is valid when `errors` is empty. It runs
in `npm run lint` for every skill in the repository and when a skill is installed or added.

## 7. Versioning

`contractVersion` follows the content distribution design §4.4: minor versions are additive
only (new optional fields, new path formats, new interpreters); anything else needs an ADR and
a deprecation window.

## Changes

- **1.1.0** (2026-10-02): optional `phi.deprecated` message. Additive.
