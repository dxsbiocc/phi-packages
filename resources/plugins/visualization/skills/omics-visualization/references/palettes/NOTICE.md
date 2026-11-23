# Palette data notice

`colors.json` (and the index in `../palettes.yaml`) is a catalog of color
values collected from several sources. It is **not** covered by the repository's
MIT license (see the root `LICENSE`); each part stays under the terms of its
original source. The hex values are used as data for coloring figures.

| Category in `colors.json` | Source | Terms |
| --- | --- | --- |
| Collection as a whole (`Brand`, `Artwork`, `Qualitative`, `Concept`, `Diverging`, `Gradient`, `Quantitative`, `Aggregation`) | Loading Colors, <https://loading.io/color/feature/> (PlotDB Ltd.) | Its [license page](https://loading.io/license/) does not state which license the color palettes fall under. Not confirmed for redistribution. |
| ColorBrewer schemes (for example `Accent`, `Dark2`, `Paired`, `Set2`, `Blues`, `RdBu`, `Spectral`) | ColorBrewer, Cynthia Brewer, Mark Harrower and The Pennsylvania State University, <https://colorbrewer2.org> | Apache License 2.0 (as understood by the maintainers; confirm at the source). |
| CARTOColors schemes (for example `Safe`, `Bold`, `Vivid`, `Prism`, `Antique`, `Pastel`, `BluGrn`, `Burg`, `Emrld`) | CARTOColors, CARTO, <https://carto.com/carto-colors/> | CC BY 3.0, attribution required (as understood by the maintainers; confirm at the source). |
| `Brand` palettes | Colors used by the named companies and products | The names and marks belong to their owners. Their use here implies no affiliation or endorsement. |
| `Artwork` palettes | Colors sampled from the named films, animations and artworks | The works belong to their creators. Their use here implies no affiliation or endorsement. |

Palette values were also merged from OmicsAgent's frontend ColorPicker (same
project owner) and the `ldpalettepicker` library from PlotDB; check that
library's repository for its license.

If you redistribute `colors.json`, keep this notice with it. If you are a
rights holder and want an entry changed or removed, open an issue.

## Summary

`colors.json` collects palettes from several sources and is **not** covered by
the MIT license at the repository root; each part remains under the terms of its
original source. The loading.io [license page](https://loading.io/license/)
does not state which license tier the palette library falls under, so
redistribution is unconfirmed; ColorBrewer is used under its published Apache
2.0, and CARTOColors under its published CC BY 3.0 (attribution required), with
the original site as the authority. Brand and artwork palettes are only colour
sampling, with no affiliation or endorsement relationship to the companies or
works concerned.
