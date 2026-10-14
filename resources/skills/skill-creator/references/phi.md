# Phi package rules

Use these rules when the skill is added to this Phi Packages repository. The normative text is `docs/contracts/skill.md` and `docs/contracts/environment.md`. The rest of this skill is the upstream skill-creator procedure.

## Layout

Standalone catalog skills live in `resources/skills/<name>/`. The directory name equals `name`.

```text
resources/skills/<name>/
  SKILL.md                 # required
  references/              # optional, read on demand
  scripts/                 # optional, run only through skill_run or script tools
  assets/                  # optional templates and data; never executed
  schemas/                 # optional JSON Schemas for script-tool output
```

A skill inside a plugin lives at `resources/plugins/<plugin>/skills/<name>/` and must be listed in that plugin's `components.skills`. `create-wrapper` stays in the Phi application. Do not add `resources/skills/create-wrapper/` here.

## Frontmatter

`SKILL.md` starts with YAML frontmatter. Standard Agent Skills fields stay at the top level. Everything Phi-specific goes under `phi:`.

| Field | Rule |
| --- | --- |
| `name` | Required. `^[a-z0-9][a-z0-9-]{0,63}$`, no leading, trailing, or doubled hyphen. Equals the directory name. |
| `description` | Required. 1–1024 characters. What it does and when to use it, including the words a user would actually say. |
| `license` | Optional string. |
| `compatibility` | Optional string. Only when the environment is unusual. |
| `metadata` | Optional map of string to string. `version` is the package version when it is semantic `MAJOR.MINOR.PATCH` with an optional prerelease; otherwise the catalog builder uses `1.0.0`. |
| `allowed-tools` | Optional space-separated string. |
| `phi` | Optional. Required only for the cases below. |

Do not put Phi keys beside `name`. Unknown top-level keys are warnings and are ignored.

A skill with no scripts must not declare `phi.environment`. A skill whose `scripts/` contains `.py`, `.R`, `.r`, `.sh`, `.js`, `.mjs`, or `.pl` must declare one. Prefer an official environment such as `phi:python@1` or `phi:r@1`. Use `./environment.yml` only with a lock for `darwin-arm64`, `darwin-x64`, and `linux-x64`.

Script tools (`phi.scripts`) also need `phi.toolPrefix` on a standalone skill. The prefix matches `^[a-z][a-z0-9]{1,11}$` and is not `skill`, `env`, `http`, `wrapper`, `agent`, `db`, or `mcp`. Inside a plugin, omit `toolPrefix` and use the plugin's prefix.

## Checklist

Reject the skill until every applicable line is true.

- Directory name equals `name`, and `name` matches the pattern above.
- `description` is one to 1024 characters and states both the job and the cues that should trigger it.
- `metadata`, when present, is a map of strings. No secrets, credentials, or machine-local paths.
- `references/` and `assets/` are not executed and are not treated as extra skills.
- `SKILL.md` links to supporting files with paths relative to the skill root.
- No script files and no `phi.scripts`: omit `phi` entirely, or omit `environment`.
- Any script file under `scripts/`: set `phi.environment` to `phi:<name>@<major>`, `plugin:<name>`, `project:<name>`, or `./environment.yml`.
- `./environment.yml` has `name`, `channels`, and `dependencies`, and `locks/darwin-arm64.txt`, `locks/darwin-x64.txt`, and `locks/linux-x64.txt` exist.
- `skill_run` paths stay inside `scripts/`. Do not tell the agent to execute `references/` or `assets/`.
- When `phi.scripts` is set, a standalone skill sets `phi.toolPrefix` and a plugin skill does not. Each entry has `name`, `description`, `run`, `args`, and `approval` (`read` or `write`). `args` is a JSON Schema object with `additionalProperties: false`. Successful stdout is one JSON object.
- `README.md` and `README.en.md` standalone skill counts match the number of directories in `resources/skills/`.
- Do not add a row to a local routing table. `find-skills` searches the public ecosystem with `npx skills find`.
- Do not add `phi-package.yaml` unless the package needs a title, dependency, or requirement the generated manifest cannot express.
- `node --test tests/skills.test.mjs tests/content-boundaries.test.mjs tests/content-icons.test.mjs` passes.
