# metagenomics (metagenomic analysis) · 0.1.0

**Downstream analysis chains + knowledge layer** for two metagenomics routes:
read-based species/function profiling (Kraken2/Bracken, MetaPhlAn/HUMAnN
reports → matrix → rarefaction curves/diversity/differential/figures,
including SGB strain-level views and a functional differential chain) and the
assembly/MAG route (QUAST statistics → CheckM/CheckM2 quality tiers →
filtering → joined MAG table (CoverM reliability masking + relative
abundance) → prokka/prodigal/eggNOG annotation summaries → GTDB tree
distances).
Every heavy-computation step (QC, classification, assembly, binning, CheckM,
GTDB-Tk, dRep, CoverM, annotation) runs through existing wrappers; the plugin
scripts only do deterministic downstream statistics, in one private locked
Python environment.

Two metagenomics routes as downstream script chains plus a knowledge layer:
read-based profiling (reports → matrix → rarefaction/diversity/differential/
figures, with SGB strain-level views and functional differential chains) and
assembly/MAGs (QUAST stats → quality tiers → filtering → joined MAG table
with CoverM masking/normalization → annotation summaries → GTDB tree
distances).
All heavy steps run through existing wrappers; plugin scripts do
deterministic downstream statistics in one private locked Python environment.

## Components

| Component | Role |
|--|--|
| Agent `Metagenomics` | Route choice (read-based vs assembly/MAG), wrapper/plugin division of work, database disk and download confirmation, threshold decision rules |
| Skill `profiling-workflow` | read-based downstream: 7 `mg_*` tools (inspect/merge/rarefaction/diversity/differential/function/plot) |
| Skill `mag-workflow` | assembly/MAG downstream: 6 `mg_*` tools (assembly_stats/bin_qc/bin_filter/mag_table/mag_function/mag_tree) |
| Environment `metagenomics` | Python 3.13 + pandas/numpy/scipy/matplotlib-base, locks for three platforms |

Both skills share `plugin:metagenomics`. The spec and the three-platform
locks live in `environments/metagenomics/`. v0.1.0 has script chains only,
no notebook renderer.

## Routes (the agent's job)

- Existing Kraken2/Bracken/MetaPhlAn/HUMAnN reports → read-based downstream chain.
- Existing contigs/bins/CheckM tables, or the goal is genome-resolved (MAG) → assembly/MAG chain.
- Raw reads only → go upstream through the wrappers first (fastp → kraken2/bracken or
  metaphlan/humann3; or megahit/spades → binning → checkm2); the plugin starts from reports/tables.
- 16S/ITS amplicons → the ampliseq wrapper pipeline; neither chain applies.
- Viromes, pangenomes, strain tracking → the knowledge layer points the way (genomad/checkv, bowtie2 and other
  wrappers); not implemented in this plugin.

## Tools

read-based (prefix `mg`):

| Tool | Purpose |
|--|--|
| `mg_inspect` | Identify the input kind (Kraken2/Bracken/MetaPhlAn/MPA/HUMAnN/matrix/samplesheet; also recognizes QUAST, bin quality tables, prokka/prodigal/eggNOG outputs and Newick trees), reporting sample names and the follow-up steps available |
| `mg_merge` | Merge same-kind reports into a taxa × sample matrix (chosen rank; Kraken2/Bracken counts, MetaPhlAn relative abundance; rank T = SGB strain level, with an sgb_species mapping), returning `value_kind` |
| `mg_rarefaction` | Rarefaction curves for count matrices (seeded hypergeometric subsampling, mean observed richness ±sd, group coloring, flags for samples that have not plateaued); a sequencing-depth adequacy check |
| `mg_diversity` | α (observed/Chao1/Shannon/Simpson/InvSimpson), Bray-Curtis, PCoA, PERMANOVA (seeded permutations), with boxplots and ordination figures |
| `mg_differential` | Two-group Mann-Whitney U (log2FC + effect size r) / multi-group Kruskal-Wallis + BH FDR, volcano plot; the matrix can be species or functional groups (mg_function output, values=relative) |
| `mg_function` | HUMAnN-style functional tables: column normalization (relative/CPM), regroup via a mapping (UniRef90→KO etc.), stratified/unstratified split; normalized/regrouped tables feed mg_differential/mg_diversity/mg_plot |
| `mg_plot` | Composition stacked barplot (top N + Other, optionally group means), taxa heatmap (optional log10 and row/column hierarchical-clustering order) |

assembly/MAG (same prefix):

| Tool | Purpose |
|--|--|
| `mg_assembly_stats` | QUAST report.tsv → contigs/total length/N50/L50/GC% per assembly (using the >=0 bp rows) |
| `mg_bin_qc` | CheckM2/CheckM quality table → MIMAG tiers (HQ ≥90/≤5, MQ ≥50/≤10, thresholds can be passed in), tier table + completeness-contamination scatter plot |
| `mg_bin_filter` | Copy the bin FASTAs that pass the thresholds into a clean set + manifest (default MQ floor 50/10) |
| `mg_mag_table` | Left-join the quality table with the GTDB-Tk summary (split into d__…s__ columns), CoverM abundance and dRep Cdb clusters → mag-table.tsv + phylum-level barplot; masks CoverM rows with low covered_fraction (default 0.3) + optional per-sample relative abundance (genome_size weighted) |
| `mg_mag_function` | prokka directory / prodigal .faa(.gz) / eggNOG .emapper.annotations → per-bin gene counts, MIMAG structural evidence (5S/16S/23S + tRNA≥18, prediction level), bin×KO and bin×COG matrices + COG stacked barplot |
| `mg_mag_tree` | GTDB classify.tree (pure-Python Newick parsing) → patristic distance matrix between MAGs + each MAG's nearest GTDB reference and distance (`--summary` identifies user tips precisely) |

## Knowledge layer

Each skill has one `references/analysis-playbook.md` (the agent decision
handbook: call → read JSON → decide → red lines), each step has one
`references/step-notes/<step>.md` (purpose/principle/result interpretation/
pitfalls), `references/README.md` is an index by function, and
`maintainers/source-map.md` records provenance mapping, licenses and
citations. All original prose; no upstream text is copied.

## Division of work with wrappers

- Upstream and heavy computation: module wrappers for fastp, kneaddata, kraken2, bracken, metaphlan, humann3,
  megahit, spades, metabat2, maxbin2, concoct, dastool, checkm/checkm2,
  gtdbtk, quast, prodigal/prokka/eggnogmapper, abricate/rgi and others,
  plus the whole **nf-core/mag 5.5.0** pipeline wrapper (also expanded here from a 5-parameter starter).
  mag bundles prodigal/prokka gene prediction itself; eggNOG functional annotation is a
  separate wrapper step after mag.
- **dRep and CoverM are not part of mag 5.5.0** — they are separate wrapper steps after mag, and the
  playbook spells out the hand-off.
- Database size is a hard constraint: GTDB ~80+ GB, Kraken2 libraries tens of GB, MetaPhlAn/HUMAnN
  single-digit to low-teens GB; before any wrapper run that triggers a download, confirm disk space and user consent.

## Sources and licenses

- **EasyMetagenome v1.24 / EasyMicrobiome** (Liu lab, GPL-3.0; iMeta 2025,
  DOI 10.1002/imt2.70001): pipeline segmentation, QC parameters, dRep threshold conventions. Pipeline knowledge only.
- **Snakemake Workflow repository** (local): tool pairing and input/output shape cross-checks. No code copied.
- **sib-swiss/metagenomics-data-analysis-training** (GPL-3.0): decision-loop structure and
  interpretation caveats. **Ideas only; no text or code copied.**
- **cambiotraining/metagenomics** (CC BY 4.0): assembly-side practical details (CheckM's two-step
  usage, bowtie2 host removal, MIMAG quality discussion). Attribution in both source-map.md files.
- **nf-core/mag 5.5.0** (vendored wrapper): parameter names and defaults cross-checked against its
  nextflow_schema.json.

The statistics are original implementations: Shannon/Simpson/Chao1, Bray-Curtis, PCoA (Gower
centering), PERMANOVA (Anderson 2001), Mann-Whitney/Kruskal-Wallis (scipy), BH FDR,
rarefaction curves (sequential hypergeometric subsampling), heatmap clustering (average linkage, scipy).
The full citation list is in
[source-map.md](skills/profiling-workflow/references/maintainers/source-map.md).

## Boundaries

- No packages are installed in the managed environment; the environment has pandas/numpy/scipy/matplotlib only.
- No LEfSe LDA/cladograms (lefse + graphlan wrappers), Krona (krona
  wrapper), ANCOM-BC, or network plots — gaps are reported honestly, not pretended away.
- MIMAG HQ tiers here rest on completeness/contamination only; strict HQ also needs rRNA/tRNA evidence —
  `mg_mag_function` gives a prediction-level structural criterion from prokka annotations (mimag_structural_pass),
  not barrnap/tRNAscan-SE confirmation; both the tool caveat and the report must state the level.
- No tree rendering (the environment has no phylogenetics libraries); mg_mag_tree only gives the distance matrix and nearest references.
- No virome/pangenome analysis; 16S goes to ampliseq.