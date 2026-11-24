# metagenomics (metagenomic analysis) · 0.3.0

**Downstream analysis chains + knowledge layer** for three microbiome routes
sharing one private locked Python environment (`plugin:metagenomics`, Python
3.13 + pandas/numpy/scipy/matplotlib):

1. **Read-based profiling** (`profiling-workflow`, 9 `mg_*` tools): Kraken2 /
   Bracken / MetaPhlAn / HUMAnN reports → abundance matrix → rarefaction
   curves → alpha/beta diversity (Bray-Curtis, PCoA, PERMANOVA, **PERMDISP**,
   and **multi-covariate marginal PERMANOVA**) → differential abundance →
   **Spearman co-occurrence networks** → **leak-free cross-validated biomarker
   classification (mRMR + L2 logistic regression)** → functional normalization
   / regrouping → composition barplots and clustered heatmaps.
2. **Assembly / MAG** (`mag-workflow`, 7 `mg_*` tools): QUAST statistics →
   CheckM/CheckM2 MIMAG quality tiers → filtering → joined MAG table (GTDB-Tk
   + CoverM reliability masking/normalization + dRep clusters) →
   prokka/prodigal/eggNOG gene/KO/COG summaries → **AMR, virulence, stress,
   and antiSMASH BGC trait summaries (`mg_mag_traits`)** → GTDB Newick
   patristic tree distances.
3. **16S / ITS amplicon** (`amplicon-workflow`, 1 `mg_*` tool + shared
   `profiling-workflow` tools): ASV count table + rooted phylogenetic tree
   exported from `nf-core/ampliseq` or QIIME 2/DADA2 → exact **Weighted and
   Unweighted UniFrac** beta diversity, PCoA, PERMANOVA, and PERMDISP.

Every heavy-computation step (QC, host depletion, classification, denoising,
assembly, binning, CheckM2, GTDB-Tk, dRep, CoverM, Prokka/eggNOG,
AMRFinderPlus/ABRicate/RGI/antiSMASH) runs through existing wrappers; the
plugin scripts only perform deterministic downstream statistics and reporting.

## Components

| Component | Role |
| --- | --- |
| Agent `Metagenomics` | Route selection (read-based, assembly/MAG, 16S/ITS amplicon), wrapper/plugin division of work, database disk and download confirmation, threshold and red-line rules |
| Skill `profiling-workflow` | Read-based downstream: 9 `mg_*` tools (`inspect`, `merge`, `rarefaction`, `diversity`, `differential`, `function`, `plot`, `network`, `classify`) |
| Skill `mag-workflow` | Assembly/MAG downstream: 7 `mg_*` tools (`assembly_stats`, `bin_qc`, `bin_filter`, `mag_table`, `mag_function`, `mag_tree`, `mag_traits`) |
| Skill `amplicon-workflow` | 16S/ITS amplicon downstream: `mg_amplicon_unifrac` (Weighted & Unweighted UniFrac + PCoA/PERMANOVA/PERMDISP) + hand-off to `profiling-workflow` |
| Environment `metagenomics` | Python 3.13 + pandas/numpy/scipy/matplotlib-base, explicit locks for three platforms |

All three skills share `plugin:metagenomics`. The spec and three-platform locks
live in `environments/metagenomics/`. v0.2.0 ships deterministic script chains
with `contract-1.0.0` artifact descriptors (no notebook renderer yet).

## Routes (the agent's job)

- Existing Kraken2/Bracken/MetaPhlAn/HUMAnN reports → `profiling-workflow`.
- Existing contigs/bins/CheckM/GTDB/AMR/antiSMASH outputs, or genome-resolved
  (MAG) goals → `mag-workflow`.
- Existing ASV tables + rooted phylogenetic trees from `nf-core/ampliseq` or
  QIIME 2 → `amplicon-workflow` (and `profiling-workflow` for rarefaction,
  differential abundance, networks, and plots).
- Raw reads only → run upstream wrappers first (`fastp`/`kneaddata` →
  `kraken2`/`bracken` or `metaphlan`/`humann3`; `nf-core/mag`; or
  `nf-core/ampliseq`).
- Viromes and pangenomes → the knowledge layer points to wrappers
  (`genomad`/`checkv`, `virsorter2`/`dramv`, `roary`/`panaroo`); no dedicated
  plugin scripts yet.

## Tools

### Read-based (`profiling-workflow`, prefix `mg`)

| Tool | Purpose |
| --- | --- |
| `mg_inspect` | Identify report/table formats (Kraken2/Bracken/MetaPhlAn/MPA/HUMAnN/matrix/samplesheet, plus QUAST, CheckM, Prokka/Prodigal/eggNOG, and Newick trees) and list ready follow-up steps |
| `mg_merge` | Merge per-sample reports into a `taxa × sample` matrix at a chosen rank (Kraken2/Bracken counts; MetaPhlAn relative abundance; rank `T` = SGB strain level with `sgb_species` mapping) |
| `mg_rarefaction` | Seeded without-replacement rarefaction curves for integer count matrices (mean richness ± SD, group coloring, non-plateau flags) |
| `mg_diversity` | Alpha diversity (Observed/Chao1 for counts; Shannon/Simpson/InvSimpson), Bray-Curtis, PCoA, PERMANOVA, **PERMDISP** dispersion homogeneity test, and optional **`--covariates` multi-covariate marginal PERMANOVA** (categorical and numeric metadata) |
| `mg_differential` | Two-group Mann-Whitney U (`log2FC` + effect size `r`) or multi-group Kruskal-Wallis + BH FDR and volcano plot on taxa or functional tables |
| `mg_function` | HUMAnN-style functional tables: column normalization (`relative`/`cpm`), regrouping via mapping TSV (e.g., UniRef90 → KO), and stratified/unstratified splitting |
| `mg_plot` | Stacked composition barplot (top N + Other, optional group means) and taxa heatmap (optional `log10` + hierarchical clustering) |
| `mg_network` | Pairwise Spearman co-occurrence network with prevalence filtering, BH FDR control, hub node ranking (`degree`, `weighted_degree`), and network plot |
| `mg_classify` | Leak-free stratified k-fold cross-validated binary classification (in-fold `log1p` + z-score scaling, in-fold `mRMR`/`mannwhitney` selection, L2 logistic regression) with pooled OOF AUROC, Youden sensitivity/specificity, feature weights, and ROC plot |
| `mg_multivariable` | Compositionally-aware CLR + OLS multivariable differential abundance with covariate adjustment (`exposure + covariates`), Student-t p-values, BH FDR, and volcano plot |

### Assembly / MAG (`mag-workflow`, prefix `mg`)

| Tool | Purpose |
| --- | --- |
| `mg_assembly_stats` | QUAST `report.tsv` → contigs, total length, N50/L50, GC% per assembly (`>=0 bp` rows) |
| `mg_bin_qc` | CheckM2/CheckM quality table → MIMAG completeness/contamination tiers (HQ ≥90/≤5, MQ ≥50/≤10) + scatter plot |
| `mg_bin_filter` | Copy passing bin FASTAs into a clean set + `manifest.tsv` (default MQ floor 50/10) |
| `mg_mag_table` | Left-join quality + GTDB-Tk summary (`d__`…`s__`) + CoverM abundance (low `covered_fraction` masking + optional genome-size-weighted relative abundance) + dRep clusters → `mag-table.tsv` + phylum barplot |
| `mg_mag_function` | Prokka / Prodigal / eggNOG outputs → per-bin gene counts, prediction-level MIMAG structural evidence (`5S/16S/23S` + `tRNA≥18`), `bin × KO` and `bin × COG` matrices + COG barplot |
| `mg_mag_traits` | Downstream aggregator for `amrfinderplus`, `abricate` (AMR/VFDB), `rgi` (CARD), and `antismash` JSON → `mag-traits-summary.tsv`, `bin × AMR class` and `bin × BGC type` matrices + overview barplot |
| `mg_mag_metabolism` | MAG metabolic module completeness (10 built-in core KEGG modules or custom `--module_map`, plus `gapseq` `*-Pathways.tbl`) and CarveMe/gapseq SBML GEM complexity summaries (`metabolites`, `reactions`, `exchange_reactions`, `genes`) + completeness heatmap |
| `mg_mag_tree` | Pure-Python Newick parsing of GTDB `classify.tree` → pairwise patristic distance matrix and nearest GTDB reference per MAG |

### 16S / ITS Amplicon (`amplicon-workflow`, prefix `mg`)

| Tool | Purpose |
| --- | --- |
| `mg_amplicon_unifrac` | Exact Normalized Weighted UniFrac and Unweighted UniFrac distance matrices from an ASV count table and rooted Newick tree (strictly ID-aligned), plus PCoA, PERMANOVA, and PERMDISP |

## Division of work with wrappers

- Upstream and heavy computation stay on wrappers: `fastp`, `kneaddata`,
  `kraken2`, `bracken`, `metaphlan`, `humann3`, `megahit`, `spades`, `myloasm`,
  `metabat2`, `maxbin2`, `concoct`, `dastool`, `checkm`/`checkm2`, `gtdbtk`,
  `drep`, `coverm`, `quast`, `prodigal`/`prokka`/`eggnogmapper`,
  `amrfinderplus`/`abricate`/`rgi`/`antismash`, plus the full `nf-core/mag`
  and `nf-core/ampliseq` pipelines.
- Confirm disk space and user consent before any wrapper run that downloads
  large databases (GTDB ~80+ GB, Kraken2 tens of GB, MetaPhlAn/HUMAnN several GB).

## Boundaries

- Locked environment contains `pandas`, `numpy`, `scipy`, and `matplotlib-base`
  only; R-based multivariable differential abundance (`MaAsLin 3`, `ANCOM-BC`),
  LEfSe cladograms (`lefse` + `graphlan` wrappers), Krona (`krona` wrapper),
  and phylogenetic tree rendering remain outside this Python environment.
- Keep MetaPhlAn `--unknown_estimation` flags and database versions identical
  across all samples in a cohort, and match 16S taxonomy classifiers to the
  actual sequenced hypervariable region (e.g., V3–V4 vs V4).
- `mg_mag_traits` reports sequence-homology predictions from upstream wrappers,
  not phenotypic antimicrobial resistance or active secondary-metabolite
  expression.
