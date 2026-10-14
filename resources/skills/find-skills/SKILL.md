---
name: find-skills
description: Find and choose a Phi skill for a task. Use when the user asks which skill to use, what skills exist, or to search, list, or compare skills for single-cell analysis, plotting, chemistry, machine learning, survival, networks, sequencing files, Nextflow, document conversion, literature, writing, or a new workflow. Also use when the right skill is unclear before specialized work, including 查找技能, 用哪个技能, and skill discovery.
allowed-tools: Read Bash
metadata: {"version": "1.0.0", "skill-author": "Phi Packages"}
---

# Find skills

## Overview

Choose an installed skill from evidence in the skill directories, not from memory. A skill's `name` and `description` are the only fields loaded before it is opened, so match the user's task against those descriptions and then read the chosen `SKILL.md`.

The routing table in [references/routing.md](references/routing.md) is a hint for skills that ship with Phi Packages. A directory scan overrides it.

## When to use this skill

Use it when:

- The user asks what skill exists, which skill fits, or to list skills.
- A task might belong to a specialized skill (scanpy, nextflow, rdkit, and the rest) and the choice is not already obvious.
- Nothing in the current context names the skill to load.

Do not use it to author a skill. Hand that to `create-skills`. Do not use it to invent a workflow when a skill already covers the task.

## Search

Read [references/where-to-look.md](references/where-to-look.md) and search in this order:

1. Standalone skills in `resources/skills/<name>/SKILL.md` when this repository is the workspace.
2. Plugin skills in `resources/plugins/<plugin>/skills/<name>/SKILL.md`.
3. Installed packages under `~/.phi/packages/skill/` and `~/.phi/packages/plugin/`.
4. The signed catalog index, when the user has one, for skills that are published but not installed.

For each candidate, read only the YAML frontmatter first (`name`, `description`). Open the body of a skill after it is the best match, or when two descriptions still look tied.

Record the match as:

| Skill | Where it lives | Why it matches | Gaps |
| --- | --- | --- | --- |
| `scanpy` | `resources/skills/scanpy` | QC and clustering of scRNA-seq | Not for `.h5ad` format questions; use `anndata` |

## Choosing

- Prefer one skill. Name a second only when the task splits cleanly, such as `anndata` for the file format and `scanpy` for the analysis.
- Prefer an installed skill over a catalog entry that is not installed. Say when a skill must be installed before it can run.
- `create-wrapper` is an engine skill shipped with Phi, not a package in this repository. Recommend it for wrapper authoring, and do not look for it under `resources/skills/` here.
- If no description matches, say so. Offer `create-skills` when the user wants a new reusable procedure. Do not stretch a nearby library skill to cover an unrelated job.

## Report

Lead with the chosen skill name and the task it should take. Then give the path, the one-line reason, and any skill that should stay unused. Quote the description only when the choice needs justification.
