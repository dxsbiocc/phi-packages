# Where skills live

Search these locations. Skip a path that does not exist.

## This repository

Standalone skills, one directory per skill:

```text
resources/skills/<name>/SKILL.md
```

Plugin skills:

```text
resources/plugins/<plugin>/skills/<name>/SKILL.md
```

The visualization plugin contributes `omics-visualization`. The Office plugin is not in this repository.

List frontmatter without reading every body:

```bash
for f in resources/skills/*/SKILL.md resources/plugins/*/skills/*/SKILL.md; do
  [ -f "$f" ] || continue
  awk '
    BEGIN { p = 0 }
    /^---$/ { p++; if (p == 2) exit }
    p == 1 && ($0 ~ /^(name|description):/) { print FILENAME ": " $0 }
  ' "$f"
done
```

A folded `description: >-` continues on indented lines until the next unindented key. Join those lines before judging a match.

## Installed on the machine

```text
~/.phi/packages/skill/<id>/<version>/SKILL.md
~/.phi/packages/plugin/<id>/<version>/skills/<name>/SKILL.md
```

One active version is installed per package. If several version directories exist, prefer the one the installer marked active; otherwise read the highest semantic version and say which path was used.

```bash
find ~/.phi/packages/skill ~/.phi/packages/plugin \
  -name SKILL.md -type f 2>/dev/null
```

## Catalog, not yet installed

Phi's official catalog is the `catalog-v1` release of this repository. Browsing it yields `index.json`. Skill rows have `"type": "skill"` and a `summary` taken from the skill description. An installed copy is what the agent can execute; a catalog row only says the skill can be installed.

`resources/skills/create-wrapper/` is intentionally absent. Wrapper authoring stays in the Phi application.

## What not to treat as a skill

- `references/`, `scripts/`, and `assets/` are files inside a skill. They are not separate skills.
- Connector packages under `resources/connectors/` are MCP servers, not skills. `literature-search` explains which of them to use for papers.
- Wrapper adapters under `resources/wrappers/` are Nextflow packages. Use the `nextflow` skill, or `create-wrapper` in the Phi app, for those.
