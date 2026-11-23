# Red-line cheatsheet (agent runtime reference)

Condenses the **hard red lines** and **common misreadings** of
[`analysis-playbook.md`](analysis-playbook.md) §1–§10 into one table for quick
reference at runtime. This table is a summary; the full Call/Read/Decide decision
procedure and the *why* remain governed by analysis-playbook.md; the per-step
teaching short notes are in [`step-notes/`](step-notes/).

## Global (applies at every stage)

| Hard red line | Common misreading |
| --- | --- |
| Expression of record = `layers["lognorm"]`; markers/panels/SVG/communication never use scaled values | Running tests/scoring on z-scores |
| Spatial statistics need `obsm["spatial"]` coordinates + an adjacency graph; without coordinates only non-spatial analysis is possible | Building the spatial graph on expression PCA/UMAP |
| Absence of `{"skipped": reason}` ≠ a clean negative result | Treating skipped as "not significant" |
| Model-based gaps (cell2location/tangram/stereopy/liana/scvi) → project environment, never install inside the read-only environment | `pip install` inside the plugin environment |
| Every threshold is a Visium/Stereo-seq starting point, not a default; read the distribution first, then decide, and state which value and why | A default that reached the output untouched = no decision was made |

## By stage

| Stage (tool) | Hard red line | Common misreading |
| --- | --- | --- |
| **§1 inspect** (`st_py_inspect`) | Platform/analysis unit named first; the Stereo-seq bin is chosen deliberately; no counts layer → no QC/seurat_v3/deconvolution; `sample`/`group_column` picks the donor/condition column | Accepting `bin_auto`; filtering silently at load time; treating a cluster column as the group column |
| **§2 QC** (`st_py_qc`) | `kept_fraction`<0.5 or zero units → stop and re-read; `pct_counts_mt` absent/median 0 → fix the regex/identifiers first; QC only on raw counts | Accepting default thresholds without reading the distribution; a never-fire mt gate taken as a pass; a pooled violin hiding low-quality edges |
| **§3 preprocess+graph** (`st_py_preprocess`) | The spatial graph is built on coordinates, not on expression PCA; never cluster in a 2D embedding; one cluster per sample = batch artefact; `neighbor_method` must be stated | counts-flavour HVG applied to log values; taking dims from the elbow only misses rare domains; `grid` used on irregular bins; a single resolution only |
| **§4 markers** (`st_py_markers`) | The corrected p is descriptive (double dipping), not inference; tests are not run on scaled values; `avg_log2FC` is not a portable effect size; marker = high detection AND low detection elsewhere | Marker p treated as significance; looking only at large FC without specificity; wrong identifiers make all markers vanish |
| **§5 annotate** (`st_py_annotate`) | Zero usable panel genes = identifier mismatch; report hit counts per label + unassigned; a spatial label = a region/dominant type, not a pure cell; coarse labels are not propagated into subtypes | A narrow-margin winner treated as identification; forced labelling hides unassigned; a generic panel used as tissue-specific |
| **§6 SVG/neighbourhoods** (`st_py_spatial`) | Autocorrelation/enrichment depend entirely on the graph + labels; the SVG p-value (normal approximation) is a screen, not mechanism; depth/batch genes score high too; **enrichment is adjacency, not interaction** | SVG significance treated as mechanism; a positive z-score treated as signal; trusting the ranking without the spatial map; trusting z with a weak null at small n |
| **§7 deconvolution** (`st_py_deconvolve`) | Proportions are NNLS **estimates, not cell counts**; an scRNA reference is required; no shared genes = identifier mismatch, fix it first; reference granularity/tissue match bounds the answer | Proportions treated as cell numbers; assuming missing types can be recovered; forcing deconvolution on a small bin/cellbin |
| **§8 communication** (`st_py_communication`) | **A screen: co-expression is not proof of signalling**; ligrec reads lognorm (`use_raw=False`); inherits the graph + labels; symbols are case-insensitive but species pairs need checking; the starter set is not exhaustive | A hit treated as proof of signalling; not checking that the ligand/receptor domains are adjacent; a human-symbol resource assumed validated in mouse |
| **§9 cross-sample differential** (`st_py_spatial`+`group_column`) | **The SAMPLE is the experimental unit, not the spot**; <2 replicates/group → no condition-level claim; a descriptive screen, not a differential significance test; differences may be driven by composition/depth; <2 conditions → skipped | Spots treated as replicates (pseudoreplication); spatial differences attributed directly to regulation; composition not checked |
| **§10 notebook** (`st_py_notebook`) | The renderer fills, you decide; every config value carries a source tag, `[default]` = no decision; `self_check=true` before delivery; script and notebook share thresholds | Defaults treated as decided; the two paths disagree on thresholds |

## One-sentence recollection

**counts/lognorm is the foundation, coordinates+graph are the spatial prerequisite; p-values are
mostly descriptive/screening, not inference; a label is a regional hypothesis, a proportion is an
estimate, communication is adjacency co-expression; the experimental unit is the sample; skipped ≠
negative, gaps go to a project environment.**