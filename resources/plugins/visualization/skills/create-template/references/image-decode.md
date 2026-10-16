# Image → template decode

Use when the user attaches a figure to reproduce as a catalog template in the
sibling omics-visualization skill. This is not a plot of the user's own table.

## Read the pixels

Open the image with the file Read tool. Note layout, legends, and encodings
yourself; do not trust an auto-generated caption alone.

## Fill this contract before coding

```text
Claim: <one sentence>
Family candidates: <graph|boxplot|heatmap|...>
Observation unit: <gene, sample-pair, edge, ...>

Marks:
  - <point|segment|ribbon|box|tile|...>

Encodings:
  - position: <cartesian|polar|dual-ring|...>
  - colour: <column + role>
  - size / linewidth / shape: <column or fixed>
  - facets / groups: <column>

Inputs:
  - table A: columns...
  - table B (optional): columns...

Supplied only (never compute in plot.R):
  - <weights, communities, p-values, layouts, ...>

Must NOT do in plot.R:
  - <DE, enrichment, clustering for membership, ...>
```

## Geometry cues → likely family

| Visual cue | Prefer |
|------------|--------|
| Nodes + edges, force/stress/circle | `graph` |
| Boxes ± pair links | `boxplot` (paired / bezier-point) |
| Matrix tiles / clustered heat | `heatmap` |
| Ranked bars / enrichment dots | `bar` |
| Genome rings / karyotype / lollipop | `ideogram` |
| Flows between stages | `sankey` / `graph-chord` |

If two families fit, pick the one whose **claim** matches; document why the
other is in `avoid_when`.

## Example data from a figure

- Recreate category labels and legend levels visible in the figure.
- Use plausible numeric ranges; do not fabricate significance stars unless the
  template only demos `ggpubr` labels on example data and documents that.
- For networks: every edge endpoint must exist in the node table.
- For pairs: every pair needs an explicit `pair_id` (or sibling equivalent),
  never “first n rows of group A vs group B”.

## Visual parity checklist

After the first preview:

- [ ] Same layout grammar (rings, facets, axes), not a vague lookalike
- [ ] Legend classes match supplied categories
- [ ] Size / width scales map to the claimed columns
- [ ] Labels readable at final size; no critical clipping
- [ ] Background/legend contrast acceptable (white panel if reference is white)
