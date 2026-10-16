# Registration checklist

Paths are relative to the installed **omics-visualization** skill root: the
sibling directory `../omics-visualization/` (plugin path
`skills/omics-visualization/`). Do not add `plot.R`, previews, or validators
under `create-template`.

## Files to add

```text
scripts/<family>/<slug>/
  plot.R
  example.tsv          # or nodes.tsv + links.tsv, etc.
  preview.png          # required by catalog
  preview.pdf          # recommended
```

## Catalog entry shape

In `references/catalog/<family>.yaml`:

```yaml
  - id: "family-slug"
    title: "Human title"
    aliases:
      - "中文名"
      - "english alias"
    source: "scripts/family/slug/plot.R"
    examples:
      - "scripts/family/slug/example.tsv"
    preview: "scripts/family/slug/preview.png"
    use_when: >
      When this geometry and claim are the right reading.
    avoid_when: >
      Point to near-miss template ids. Say what not to compute.
    input_shape: "one row per ...; required columns..."
    language: "r"
    status: "stable"
```

Also patch nearby templates' `avoid_when` so routers do not land on the wrong id.

## Contract entry shape

In `references/template_contracts.json` (keep JSON valid):

```json
{
  "id": "family-slug",
  "family": "family",
  "title": "Human title",
  "source": "scripts/family/slug/plot.R",
  "preview": "scripts/family/slug/preview.png",
  "preferred_shapes": ["paired_distribution"],
  "roles": {
    "required": ["pair_id", "category", "numeric_y"],
    "optional": ["target_entity"]
  },
  "intent_keywords": [
    "family-slug",
    "bezier paired boxplot",
    "贝塞尔连线箱线图"
  ],
  "use_when": "...",
  "avoid_when": "..."
}
```

Only use role names accepted by `scripts/validate_template_contracts.py`
(via `ROLE_HINTS` in `scripts/route_template.py`), both in the
omics-visualization skill. Common ones:

`source_entity`, `target_entity`, `category`, `numeric_y`, `numeric_x`,
`pair_id`, `weight`, `magnitude`, `signed_association`, `axis_entity`,
`focal_entity`, `significance`, `flow_value`

## Optional doc touch-ups

- `references/plots.yaml` — add Chinese/English family aliases users will say
- `SKILL.md` Easy mis-routes — one bullet when confusion is likely
  (the omics-visualization `SKILL.md`, not this skill)

## Commands

Run these with the working directory set to the installed omics-visualization
skill root (`../omics-visualization/` from this skill). The programs already
live in that skill; do not recreate them here.

```bash
# render
Rscript scripts/<family>/<slug>/plot.R \
  scripts/<family>/<slug>/example.tsv \
  scripts/<family>/<slug>/preview.png

# validate
python3 scripts/validate_template_contracts.py --json
python3 scripts/qa_single_plot.py scripts/<family>/<slug>/preview.png --json
```

Contracts must be `"ok": true` with an empty `errors` list.

## Extending an existing template

When fixing honesty (e.g. require `pair_id`):

1. Edit `plot.R` CONFIG + DATA PREPARATION
2. Update example tables
3. Refresh preview
4. Align catalog `input_shape` / `use_when` and contract `roles` / `use_when`
5. Re-validate contracts from the omics-visualization skill root

Do not mint a new id for a schema honesty fix.
