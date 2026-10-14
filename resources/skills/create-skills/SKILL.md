---
name: create-skills
description: Author a new Phi skill that follows the skill contract. Use when the user wants to create, scaffold, or add a skill, write a SKILL.md, or package instructions, references, and scripts an agent can load. Covers frontmatter, directory layout, script environments, and the README skill count, including 新建技能 and 写技能.
allowed-tools: Read Write Edit Bash
metadata: {"version": "1.0.0", "skill-author": "Phi Packages"}
---

# Create skills

## Overview

Add one skill directory that Phi can validate and the catalog builder can package. The contract is `docs/contracts/skill.md` in the phi-packages repository. This skill is the authoring procedure; it does not relax that contract.

Run `find-skills` first. Extend an existing skill when the task is the same job. Create a skill when the procedure is reusable, specific, and not already described.

`create-wrapper` stays in the Phi application. Do not add `resources/skills/create-wrapper/` here.

## When to use this skill

Use it when the user asks to create, scaffold, or add a skill, or to turn a repeated workflow into instructions an agent can load. Use `find-skills` when they only want to locate one.

## Layout

Standalone catalog skills live here. The directory name is the skill name.

```text
resources/skills/<name>/
  SKILL.md                 # required
  references/              # optional, read on demand
  scripts/                 # optional, run only through skill_run or script tools
  assets/                  # optional templates and data; never executed
  schemas/                 # optional JSON Schemas for script-tool output
```

A skill inside a plugin lives at `resources/plugins/<plugin>/skills/<name>/` and must be listed in that plugin's `components.skills`.

Copy [assets/SKILL.template.md](assets/SKILL.template.md), replace the placeholders, and delete sections that do not apply. Keep `SKILL.md` under about 500 lines. Move long tables, field lists, and examples into `references/`.

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
| `phi` | Optional. Required only for the cases in the checklist. |

Do not put Phi keys beside `name`. Unknown top-level keys are warnings and are ignored.

A skill with no scripts must not declare `phi.environment`. A skill whose `scripts/` contains `.py`, `.R`, `.r`, `.sh`, `.js`, `.mjs`, or `.pl` must declare one. Prefer an official environment such as `phi:python@1` or `phi:r@1`. Use `./environment.yml` only with a lock for `darwin-arm64`, `darwin-x64`, and `linux-x64`.

Script tools (`phi.scripts`) also need `phi.toolPrefix` on a standalone skill. The prefix matches `^[a-z][a-z0-9]{1,11}$` and is not `skill`, `env`, `http`, `wrapper`, `agent`, `db`, or `mcp`. Inside a plugin, omit `toolPrefix` and use the plugin's prefix.

The full checklist is [references/checklist.md](references/checklist.md).

## After writing

1. Check the new directory against the checklist.
2. Set the standalone skill count in `README.md` and `README.en.md` to the number of directories in `resources/skills/`.
3. Add one row to `resources/skills/find-skills/references/routing.md` for each new catalog skill.
4. Do not add `phi-package.yaml` unless the package needs a title, dependency, or requirement the generated manifest cannot express. The builder derives `id`, `type: skill`, `version`, `title`, and `summary` from `SKILL.md`.
5. Do not commit secrets, `__pycache__/`, `.nextflow/`, `work/`, or `results/`.

## Report

Tell the user the skill path, the `name`, when it should trigger, whether it has scripts, and which README counts changed.
