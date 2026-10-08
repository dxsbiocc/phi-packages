# Enablement contract

contractVersion: 1.1.0

Installing content and **using** it are separate. Only enabled items reach the main
agent, which keeps its prompt and tool routing small (the original goal of content
distribution). This contract defines the enablement state, how it is resolved, and
what "enabled" means for each kind of item. See
[content distribution design](../design/phi-content-distribution-design.md) §7 and
the [package contract](package.md).

## 1. Items

An item is identified by `<kind>:<id>`:

| Kind     | `id`               | Examples               |
| -------- | ------------------ | ---------------------- |
| `skill`  | the skill's `name` | `skill:scanpy`         |
| `plugin` | the plugin's `id`  | `plugin:visualization` |

Since 1.1.0 also `wrapper:<package id>` (wrapper contract § 5) and `mcp:<package id>`
(connector contract § 3).

## 2. State file

`~/.phi/state/enabled.json`, written atomically:

```json
{
  "version": 1,
  "global": { "skill:scanpy": true, "skill:pptx": false },
  "projects": {
    "<project real path>": { "skill:scanpy": false }
  }
}
```

A missing or invalid file is read as empty (and the problem is logged); unknown item
keys are kept as they are.

## 3. Resolution

For an item in a project, the first defined value wins:

1. the project's override (`projects[<real path of the project>]`);
2. the global value (`global`);
3. the **default** for the item's source:

| Source                                                                                                           | Default                           |
| ---------------------------------------------------------------------------------------------------------------- | --------------------------------- |
| core skills shipped with Phi that the engine depends on (`create-wrapper`, `nextflow`)                           | enabled; cannot be disabled       |
| other bundled skills (`resources/skills/`)                                                                       | disabled — added from the catalog |
| skills installed as packages, user-authored skills (`~/.phi/skills/`), project skills (`<project>/.phi/skills/`) | enabled                           |
| bundled plugins                                                                                                  | enabled                           |
| plugins installed as packages                                                                                    | enabled                           |
| wrapper packages (bundled or installed) and user-authored wrappers (1.1.0)                                       | enabled                           |
| connector packages the user added (1.1.0)                                                                        | enabled                           |

The core list lives in one place in the engine and is shown in the UI as "built-in".

## 4. Effect

- A **disabled skill** is not offered to the main agent (not in its skill list or
  prompt) and its script tools are not registered for the main agent.
- A **disabled plugin** behaves as plugin contract § 4 "disable": its agents, skills,
  and tools are hidden everywhere.
- A skill that a loaded agent lists in its own `skills:` stays available **to that
  agent** even when the skill is disabled for the main agent: an agent's toolbox is
  part of the agent. Its script tools follow their `attachTo` as usual.
- Changing enablement takes effect for new sessions; a running session keeps the
  set it started with.
- Enablement never deletes or rewrites content files.

## 5. Migration

On the first start with this contract (no `enabled.json` yet), every bundled skill
that the user's existing sessions have used is enabled globally, so nobody loses a
skill they rely on. "Used" means the skill name appears as a loaded/invoked skill in
the session history Phi keeps. Everything else starts at its default.

## 6. Versioning

`contractVersion` follows the content distribution design §4.4: minor versions are
additive only (new kinds, new sources with their defaults); anything else needs a
decision record and a deprecation window.

## Changes

- **1.1.0** (2026-10-02): `wrapper:` and `mcp:` items with their defaults. Additive.
