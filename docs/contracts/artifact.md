# Artifact contract

contractVersion: 1.0.0

An artifact is a result file a tool produced for the user (a figure, a table, a
report …) plus a small descriptor that says what it is and how it was made. Script
tools list the artifacts they wrote; the engine checks each descriptor and presents
the files in the conversation. The engine knows artifact kinds, never the domain
that produced them (decision 9 of the
[content distribution decision record](../decisions/content-distribution.md)). See
the [skill contract](skill.md) § 5.5.

## 1. Files

The descriptor sits beside the artifact and is named after it:

```text
figures/volcano.png
figures/volcano.png.phi-artifact.json
```

Both files are inside the project. The descriptor is UTF-8 JSON, at most 64 KiB, and
matches [artifact.schema.json](artifact.schema.json).

## 2. Descriptor

| Field             | Required | Rule                                                                                                                   |
| ----------------- | -------- | ---------------------------------------------------------------------------------------------------------------------- |
| `contractVersion` | yes      | `1.0.0`                                                                                                                |
| `kind`            | yes      | `figure`, `table`, `structure`, `molecule`, `network`, or `report`                                                     |
| `file`            | yes      | the artifact's file name; equals the descriptor's name without `.phi-artifact.json`                                    |
| `mediaType`       | yes      | the file's media type, for example `image/png`, `image/svg+xml`, `application/pdf`, `text/tab-separated-values`        |
| `title`           | yes      | 1–200 characters, shown to the user                                                                                    |
| `description`     | no       | at most 2000 characters                                                                                                |
| `provenance`      | yes      | § 2.1                                                                                                                  |
| `figure`          | no       | only for `kind: figure`: `{ widthPx?, heightPx?, format }`, `format` one of `png`, `pdf`, `svg`                        |
| `table`           | no       | only for `kind: table`: `{ rows?, columns? }` (non-negative integers; `columns` may instead be a list of column names) |

Unknown top-level keys are errors, except keys starting with `x-`, which are kept
and ignored (producer-specific data).

### 2.1 Provenance

| Field        | Required | Rule                                                                                       |
| ------------ | -------- | ------------------------------------------------------------------------------------------ |
| `createdAt`  | yes      | ISO 8601 timestamp with a time zone                                                        |
| `tool`       | yes      | the tool that made it, for example `viz_render`                                            |
| `skill`      | no       | the skill whose script made it                                                             |
| `inputs`     | no       | list of `{ path, sha256? }`: project-relative paths of the inputs that were read           |
| `script`     | no       | project-relative path of an edited script that produced it (for example a prepared plot.R) |
| `parameters` | no       | JSON object of the arguments that matter for reproducing it                                |

The engine adds `envId` from the run when it presents the artifact; producers do
not write it.

## 3. From a script tool to the conversation

1. A script tool's JSON output may contain `artifacts: [<project path>, …]` (skill
   contract § 5.5). Only the first 8 are presented; the rest are reported as a
   warning.
2. For each path the engine checks that the file exists inside the project and that
   its descriptor exists and is valid. A path that fails is dropped from
   presentation and reported to the model as a warning in the tool result, naming
   the problem; the tool call itself still succeeds.
3. The valid artifacts are presented in the conversation the same way as
   `present_files` (a card with the file, its title, and its kind), in the order
   listed. The tool result keeps `artifacts` so the model can refer to them.

## 4. Versioning

`contractVersion` follows the content distribution design §4.4: minor versions are
additive only (new kinds, new optional fields); anything else needs a decision record
and a deprecation window.
