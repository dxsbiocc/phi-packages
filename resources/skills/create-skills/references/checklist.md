# Skill checklist

Reject the skill until every applicable line is true. The normative text is `docs/contracts/skill.md` and `docs/contracts/environment.md`.

## Identity

- Directory name equals `name`.
- `name` matches `^[a-z0-9][a-z0-9-]{0,63}$` and does not start or end with `-` or contain `--`.
- `description` is one to 1024 characters and states both the job and the cues that should trigger it.
- `metadata`, when present, is a map of strings. A package version in `metadata.version` is `MAJOR.MINOR.PATCH` with an optional prerelease.
- No secrets, credentials, or machine-local paths.

## Files

- `SKILL.md` is the only required file.
- `references/` and `assets/` are not executed and are not treated as extra skills.
- `SKILL.md` links to supporting files with paths relative to the skill root, one level deep.
- The body teaches the procedure. It does not repeat the description's trigger list as a substitute for a precise description.

## Scripts and environments

- No script files and no `phi.scripts`: omit `phi` entirely, or omit `environment`.
- Any `.py`, `.R`, `.r`, `.sh`, `.js`, `.mjs`, or `.pl` file under `scripts/`: set `phi.environment`.
- Environment reference is one of `phi:<name>@<major>`, `plugin:<name>`, `project:<name>`, or `./environment.yml`.
- `./environment.yml` has `name`, `channels`, and `dependencies`, and `locks/darwin-arm64.txt`, `locks/darwin-x64.txt`, and `locks/linux-x64.txt` exist.
- `skill_run` paths stay inside `scripts/` (`skill`, `script`, optional `args` and `cwd`). Do not tell the agent to execute `references/` or `assets/`.

## Script tools

Only when `phi.scripts` is set:

- Standalone skill sets `phi.toolPrefix`. Plugin skill does not.
- Each entry has `name`, `description`, `run`, `args`, and `approval` (`read` or `write`).
- `name` matches `^[a-z][a-z0-9_]{0,31}$` and is unique in the skill.
- `args` is a JSON Schema object with `additionalProperties: false`. Property names match `^[a-z][a-z0-9_]*$`.
- `run` starts with an interpreter command (`python`, `Rscript`, `bash`, `node`, `perl`). Paths start with `./` and stay inside the skill.
- `output`, when set, is `./<path>` to a schema file inside the skill.
- Successful stdout is one JSON object. Diagnostics go to stderr.

## Repository

- New standalone skills are under `resources/skills/`.
- `README.md` and `README.en.md` standalone skill counts match the number of directories in `resources/skills/`.
- `resources/skills/find-skills/references/routing.md` has a row for the new skill.
- Plugin skills are listed in the plugin `phi-package.yaml` `components.skills`.
- `create-wrapper` was not copied into this repository.
- `node --test tests/skills.test.mjs tests/content-boundaries.test.mjs tests/content-icons.test.mjs` passes.
