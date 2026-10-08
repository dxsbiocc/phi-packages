# Phi palettes

This directory is the single bundled source of named colors for Phi's main
agent and scientific visualization templates.

- `palettes.yaml` lists recommended palettes, their intended encodings, and
  default choices. The main agent's `palette_suggest` tool reads this index.
- `colors.json` contains the full named color catalog. Visualization's
  `palette_colors()` helper reads it when rendering a figure.

Keep palette IDs and color order stable when editing either file. The catalog
origin is recorded in `palettes.yaml`; do not copy these files into a user's
project for an ordinary plotting run.
