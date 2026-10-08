# Revise an existing figure

Use this path when the user asks to change a figure already made, including its
colors, labels, legend, size, or a small annotation. The project-local prepared
script is the source of truth. Keep the original template and data mapping.

1. Obtain the exact existing `plot.R`, input tables in render order, current
   output path, and requested change from the delegated task or prior report.
   Read the script and check those paths. If a path is missing, search only the
   known project output directory. If the source or inputs remain ambiguous,
   report the missing paths; do not start a new template selection.
2. For a color-only request, inspect the current CONFIG and color assignments,
   then read only the generated `references/palettes.yaml` or
   `references/palettes/colors.json`. Edit the smallest relevant lines in the
   project copy. Preserve data preparation, thresholds, group classification,
   labels, geometry, and export size unless the user also asked to change them.
3. Do not call `viz_route`, `viz_examples`, or `viz_prepare` for this local
   revision. Use `viz_render` on the edited script with the same input tables
   and the requested output destination. Inspect the rendered image and report
   the exact source and output paths, palette choice, and QA result.

If the user explicitly asks for a different chart type or template, return to
the new-figure route in `SKILL.md`. A missing editable script does not by itself
authorize recreating the figure from a new template.
