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

## 中文说明

`colors.json` 汇集了多个来源的配色，**不属于**仓库根目录 MIT 许可的范围，
各部分仍按原来源的条款使用。loading.io 的[许可页](https://loading.io/license/)
没有说明配色库属于哪一档许可，再分发未经确认；ColorBrewer 按其公布的 Apache 2.0、
CARTOColors 按其公布的 CC BY 3.0（需署名）使用，具体以原站为准。
品牌和作品配色只是取色，与相关公司或作品方无隶属或背书关系。
