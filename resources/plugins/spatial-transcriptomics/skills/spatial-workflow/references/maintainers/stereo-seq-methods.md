# Stereo-seq methods (maintainers)

Stereo-seq (BGI/华大) data formats and the methodological patterns shared by the
seven user-supplied application papers. The format section is verified against
`scripts/readers/stereoseq.py`; the analysis guidance is an original distillation
of established Stereo-seq practice and the papers' methods, not a copy of their
text. Where a specific per-paper parameter matters, confirm it against the
source PDF before citing it as a number — this note states the common pattern,
not invented exact values.

## The platform in one paragraph

Stereo-seq captures RNA on a DNA-nanoball (DNB) array with a centre-to-centre
spacing of ~500 nm (older chips ~715 nm). A single DNB (bin1) is sub-cellular, so
the **analysis unit is a choice**, not a given: aggregate N×N DNBs into a square
bin (binN), or use SAW's segmented cells (cellbin). bin20 ≈ 10 µm and bin50 ≈
25 µm at 500 nm spacing — comparable to a Visium spot (55 µm) but finer. This is
why the plugin makes `bin_size`/`cellbin` explicit decisions and flags an
auto-selected bin (`bin_auto`) for the agent to override.

## Formats the reader accepts

### Square-bin GEF (HDF5, `/geneExp/binN`)

- `gene` — compound dtype `(gene: S64, offset: uint, count: uint)`; gene *i* owns
  expression rows `[offset[i], offset[i]+count[i])`.
- `expression` — compound `(x: uint, y: uint, count: uint)`; x/y are DNB
  coordinates, `count` is the MID (UMI) count.
- A tissue GEF carries several bins (e.g. 1/5/10/20/50/100/150/200).
  `choose_bin` is deterministic: explicit `bin_size` > the only bin present >
  a ~10 µm preference (20, then 50/10/100) flagged `bin_auto`.
- The reader aggregates DNB coords to bin origins `(x // N) * N`, builds a
  sparse bin×gene matrix, and stores bin origins in `obsm["spatial"]`. Root
  attrs (`resolution`, minX/minY/maxX/maxY) are recorded in meta.

### Cellbin GEF (HDF5, `/cellBin`)

- Segmentation is produced **upstream by SAW cellCut**; the reader consumes it,
  never re-segments. `cellExp` links cells→genes→counts; `gene` maps IDs→names.
- Internal field names vary by SAW version, so the reader inspects the compound
  dtype and adapts (`cellID/cell/id`, `MIDCount/count/umi/exp`,
  `geneID/gene/geneIndex`), reporting the fields it used.
- Cell centroids come from `cell`/`cellBorder`/`cellCentroid`; without
  recognisable centroids the reader raises (spatial statistics need coordinates).

### GEM (SAW v8 TSV)

- `#` comment lines precede a tab-separated header.
- Square bin: `geneID, geneName, x, y, MIDCount, ExonCount`.
- Cellbin adds a `CellID` column; the reader groups by CellID and takes the first
  x/y as the centroid. `.gz` is handled transparently.

All paths return `(adata, meta)` with raw counts in `X` (caller adds the `counts`
layer) and coordinates in `obsm["spatial"]`; `finalize` asserts coordinates
exist and records `n_units`/`n_genes`/`coord_key`.

## Choosing the analysis unit (the biological decision)

| Question | Unit | Why |
| --- | --- | --- |
| Tissue domains / regional structure | larger square bin (bin50≈25µm, bin100) | enough counts per unit for stable clustering; a domain spans many cells |
| Cell-type composition / niches | small bin (bin20≈10µm) or cellbin | finer resolution separates adjacent cell types |
| Single-cell-level questions | cellbin (SAW segmented) | one cell per unit; deconvolution adds little here |
| Sub-cellular / localisation | bin1–bin5 | DNB-level; very shallow, QC and depth dominate |

Consequences the agent must state: smaller bins are shallower (lower
`min_genes`, more empty units at tissue edges); cellbin inherits SAW's
segmentation quality (merged/split cells are upstream artifacts, not biology);
the spatial adjacency graph should use `knn`/`radius` for irregular bins, not
`grid` (which assumes a regular lattice).

## The seven papers — methodological patterns

Each paper is a Stereo-seq application; together they exercise the plugin's
capability set. The points below are the shared methodological pattern each
represents, distilled at the level of common practice (confirm exact parameters
against the PDFs).

| # | Paper | System | Methodological pattern it anchors |
| --- | --- | --- | --- |
| 1 | Liu et al., Dev Cell 2022, doi:10.1016/j.devcel.2022.04.009 | zebrafish embryogenesis | large-field Stereo-seq across developmental stages; square-bin domain mapping + marker annotation of emerging tissues; cross-stage spatial comparison |
| 2 | Tan et al., Zool Res 2024, doi:10.24272/j.issn.2095-8137.2023.205 | dairy cow placenta | non-model organism (gene symbols/orthologs matter); bin-level QC + clustering + marker panels for placental cell types; deconvolution against a reference |
| 3 | Wu et al., Cell Discov 2024, doi:10.1038/s41421-024-00740-6 | mouse placentation | mouse `^mt-` QC regex; spatial domains + SVGs marking trophoblast invasion fronts; neighbourhood enrichment of maternal–foetal interfaces |
| 4 | Xu et al., Nat Genet 2024, doi:10.1038/s41588-024-01709-7 | mouse liver | zonation as a spatial gradient; SVGs (Moran's I) recovering periportal–pericentral axes; cellbin/small-bin resolution of hepatocyte zones |
| 5 | Kilfeather et al., Cell Rep 2024, doi:10.1016/j.celrep.2024.113784 | dopaminergic neurons | fine-resolution (cellbin/small bin) subtype mapping; marker-panel annotation of neuronal subtypes; spatial neighbourhood of niche cells |
| 6 | Lei et al., Cell Metab 2024, doi:10.1016/j.cmet.2024.01.003 | macaque hypothalamus | primate tissue, nuclei-rich; cell-type deconvolution + cell-cell communication (ligand-receptor) of metabolic signalling between neuronal/glial populations |
| 7 | Tsubosaka et al., Cell Rep 2023, 42(10) | stomach encyclopedia | reference atlas building; large-bin domain annotation + cross-region spatial differential; a curated marker panel per gastric cell type |

### What the plugin takes from them

- **Bin size is reported, never hidden.** Every paper's resolution choice maps to
  the `bin_size`/`cellbin` decision and the `bin_auto` flag.
- **Species-aware QC.** Mouse/placenta papers (2,3,4) need `^mt-`; the guarded mt
  gate and the absent-column warning exist for exactly this.
- **Non-model symbols.** Cow/macaque (2,6) remind us that panel genes and
  ligand-receptor symbols must match the object's identifier type; ligrec matches
  case-insensitively by orthologous symbol.
- **Spatial gradients as SVGs.** Liver zonation (4) is the canonical Moran's I
  use case; the SVG screen + spatial map confirmation is the workflow.
- **Interfaces as neighbourhood enrichment.** Placental maternal–foetal borders
  (3) and neuronal niches (5,6) are off-diagonal nhood z-scores.
- **Composition needs a reference.** Deconvolution (2,6) is NNLS against a
  single-cell reference from the single-cell plugin; proportions are estimates.
- **Communication is a screen.** Hypothalamic signalling (6) is ligrec with a
  curated resource; a hit is a hypothesis, confirmed with the spatial map.
- **Cross-region/condition comparison.** Stomach atlas (7) and developmental
  stages (1) map to `spatial-diff` with `group_column`, sample as the unit.

## Sources

Format layout: BGI SAW v8 output, verified against `scripts/readers/stereoseq.py`.
Methodological distillation: the seven papers above (user-supplied PDFs) plus
established Stereo-seq practice. No upstream text, figures or data are copied;
per-paper exact parameters should be confirmed against the source PDFs before
being cited as specific values. OSTA chapter mapping is in `osta-map.md`; the
per-step teaching notes are in `references/step-notes/`.