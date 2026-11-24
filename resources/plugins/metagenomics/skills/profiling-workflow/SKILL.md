---
name: profiling-workflow
description: >-
  Downstream analysis of read-based metagenomic profiles: inspect Kraken2 /
  Bracken / MetaPhlAn / MPA / HUMAnN reports, merge them into a taxa x sample
  matrix, check sequencing depth with rarefaction curves, compute alpha/beta
  diversity with PCoA and PERMANOVA, run non-parametric differential abundance
  (Wilcoxon / Kruskal-Wallis + BH FDR), normalize/regroup functional tables,
  and plot composition barplots and optionally clustered heatmaps. Use for
  宏基因组, metagenomics, 物种组成, 多样性, alpha, beta, PCoA, 差异丰度,
  rarefaction, 稀释曲线, kraken2, bracken, metaphlan, humann, 功能组成, or
  when per-sample taxonomic/functional reports need statistics and figures.
  Upstream classification (running kraken2/metaphlan/humann on reads) is done
  with wrappers, not here. The assembly/MAG counterpart is mag-workflow.
phi:
  environment: plugin:metagenomics
  attachTo:
    - Metagenomics
  scripts:
    - name: inspect
      description: >-
        Identify metagenomics profiling inputs (Kraken2 report, Bracken report,
        MetaPhlAn profile, MPA table, HUMAnN genefamilies/pathabundance,
        abundance matrix, samplesheet; also recognizes CheckM quality tables
        and QUAST reports as mag-route inputs). Reports per-file kind, sample
        names, ranks present, and which follow-up tools the input supports.
        Call this first on unfamiliar input.
      run: [python, ./scripts/project_inspect.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: A report/table file, or a directory of per-sample reports.
      approval: read
      output: ./schemas/inspect-result.json
      timeoutSeconds: 300
    - name: merge
      description: >-
        Merge per-sample taxonomic profiles into one taxa x sample TSV matrix
        at a chosen rank. Kraken2 uses clade-assigned reads, Bracken uses
        re-estimated reads (both counts); MetaPhlAn uses relative abundance
        (0-100). Rank T keeps MetaPhlAn/MPA SGB (strain-level) rows and
        returns an sgb_species mapping where the profile names the species.
        Sample names come from file stems unless metadata maps them.
        Pass the returned value_kind to mg_diversity and mg_differential.
      run: [python, ./scripts/merge_profiles.py]
      args:
        type: object
        additionalProperties: false
        required: [inputs, kind, output]
        properties:
          inputs:
            type: array
            items:
              type: string
              format: input-path
            description: Per-sample report files of the same kind.
          kind:
            type: string
            enum: [kraken2, bracken, metaphlan, mpa]
            description: Report format of all inputs.
          rank:
            type: string
            description: Taxonomic rank code (default S; D/K/P/C/O/F/G/S/T; T = SGB/strain, MetaPhlAn/MPA only).
          output:
            type: string
            format: project-path
            description: Output matrix TSV (rows taxa, columns samples).
          metadata:
            type: string
            format: input-path
            description: Optional CSV/TSV with sample names as first column; validates sample coverage.
          min_abundance:
            type: number
            description: Drop taxa whose max abundance across samples is below this.
          top:
            type: integer
            minimum: 1
            description: Keep only the N taxa with highest mean abundance.
      approval: write
      output: ./schemas/merge-result.json
      timeoutSeconds: 600
    - name: rarefaction
      description: >-
        Observed-richness rarefaction curves for a count matrix: seeded
        without-replacement subsampling to a depth grid, mean richness with
        spread band per sample, optional group coloring, and a flag for
        samples whose curves have not plateaued. Integer counts only
        (kraken2/bracken); rejects relative abundances and fractional
        estimates. Run between mg_merge and mg_diversity to check depth
        adequacy; never rarefy before differential testing.
      run: [python, ./scripts/rarefaction.py]
      args:
        type: object
        additionalProperties: false
        required: [input, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Taxa x sample integer-count matrix TSV (from mg_merge).
          output_dir:
            type: string
            format: project-path
            description: Directory for the curve table and figure.
          points:
            type: integer
            minimum: 2
            description: Number of depth grid points (default 10).
          iterations:
            type: integer
            minimum: 1
            description: Subsampling draws per depth (default 10).
          seed:
            type: integer
            description: RNG seed for subsampling (default 12345).
          metadata:
            type: string
            format: input-path
            description: Optional sample metadata; curves colored by group.
          group_col:
            type: string
            description: Metadata column defining groups (default group).
      approval: write
      output: ./schemas/rarefaction-result.json
      timeoutSeconds: 1800
    - name: diversity
      description: >-
        Alpha diversity (observed, Chao1, Shannon, Simpson, InvSimpson; the
        first two only for counts), Bray-Curtis beta diversity, PCoA
        ordination with variance explained, PERMANOVA and PERMDISP (homogeneity
        of multivariate dispersions) when a group column is given, plus optional
        multi-covariate marginal PERMANOVA (categorical or numeric metadata
        columns such as sex, age, BMI) via seeded permutations. Writes alpha
        table, distance matrix, PCoA coordinates, optional covariate-permanova
        table, boxplot and ordination figure.
      run: [python, ./scripts/diversity.py]
      args:
        type: object
        additionalProperties: false
        required: [input, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Taxa x sample matrix TSV (from mg_merge or equivalent).
          values:
            type: string
            enum: [counts, relative]
            description: What the matrix holds; relative skips observed/Chao1.
          metadata:
            type: string
            format: input-path
            description: Optional sample metadata with a group column and/or covariates.
          group_col:
            type: string
            description: Metadata column defining groups (default group).
          covariates:
            type: array
            items:
              type: string
            description: Optional metadata columns for marginal distance-based PERMANOVA (e.g. [group, sex, age, bmi]).
          output_dir:
            type: string
            format: project-path
            description: Directory for tables and figures.
          permutations:
            type: integer
            minimum: 99
            description: PERMANOVA/PERMDISP permutations (default 999).
          seed:
            type: integer
            description: RNG seed for permutations (default 12345).
      approval: write
      output: ./schemas/diversity-result.json
      timeoutSeconds: 1800
    - name: differential
      description: >-
        Differential abundance between groups on any feature x sample matrix
        (taxa from mg_merge, or functional groups from mg_function): two-sided
        Mann-Whitney U with tie correction for two groups (log2 fold change +
        effect size r), Kruskal-Wallis for more, Benjamini-Hochberg FDR across
        features. Writes a ranked table and a volcano plot. For LEfSe-style
        LDA/cladograms use the lefse and graphlan wrappers on the same matrix.
      run: [python, ./scripts/differential.py]
      args:
        type: object
        additionalProperties: false
        required: [input, metadata, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Feature x sample matrix TSV (taxa or functional groups).
          metadata:
            type: string
            format: input-path
            description: Sample metadata CSV/TSV with the group column.
          group_col:
            type: string
            description: Metadata column defining groups (default group).
          case:
            type: string
            description: Group label treated as case (with control, restricts to two groups).
          control:
            type: string
            description: Group label treated as control.
          values:
            type: string
            enum: [counts, relative]
            description: Matrix value kind; sets the fold-change pseudocount.
          output_dir:
            type: string
            format: project-path
            description: Directory for the table and volcano plot.
          fdr:
            type: number
            description: Significance threshold on BH FDR (default 0.05).
          fc_threshold:
            type: number
            description: Minimum absolute log2 fold change for two-group tests (default 1).
      approval: write
      output: ./schemas/differential-result.json
      timeoutSeconds: 1200
    - name: function
      description: >-
        Deterministic operations on HUMAnN-style functional tables
        (genefamilies, pathabundance, or any feature x sample table):
        normalize columns to relative abundance or CPM, regroup features via
        a two-column mapping TSV (e.g. UniRef90 -> KO), or keep only
        stratified / unstratified rows. Normalized/regrouped tables feed
        mg_differential (values=relative), mg_diversity and mg_plot. Plain
        column renormalization, not HUMAnN's stratification-aware renorm.
      run: [python, ./scripts/function_table.py]
      args:
        type: object
        additionalProperties: false
        required: [input, action, output]
        properties:
          input:
            type: string
            format: input-path
            description: Functional table TSV (features as rows, samples as columns).
          action:
            type: string
            enum: [normalize, regroup, stratified, unstratified]
            description: Operation to apply.
          output:
            type: string
            format: project-path
            description: Output table TSV.
          unit:
            type: string
            enum: [relative, cpm]
            description: Normalization unit (default relative).
          mapping:
            type: string
            format: input-path
            description: Two-column TSV feature -> group, required for regroup.
      approval: write
      output: ./schemas/function-result.json
      timeoutSeconds: 600
    - name: plot
      description: >-
        Composition figures from a taxa x sample matrix: stacked barplot of
        the top N taxa (plus Other), optionally a group-mean barplot when
        metadata is given, or a taxa x sample heatmap with optional log10
        transform and hierarchical clustering of rows/columns.
      run: [python, ./scripts/plot_composition.py]
      args:
        type: object
        additionalProperties: false
        required: [input, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Taxa x sample matrix TSV.
          output_dir:
            type: string
            format: project-path
            description: Directory for figures.
          kind:
            type: string
            enum: [barplot, heatmap]
            description: Figure type (default barplot).
          top:
            type: integer
            minimum: 1
            description: Number of taxa to show (default 20).
          metadata:
            type: string
            format: input-path
            description: Optional sample metadata for the group-mean barplot.
          group_col:
            type: string
            description: Metadata column defining groups (default group).
          trans:
            type: string
            enum: [none, log10]
            description: Heatmap transform (default none).
          pseudocount:
            type: number
            description: Pseudocount for log10 heatmap (default 1e-6).
          cluster:
            type: string
            enum: [none, rows, cols, both]
            description: >-
              Heatmap only: reorder by average-linkage hierarchical clustering
              (Euclidean on displayed values; default none). Pair with log10.
      approval: write
      output: ./schemas/plot-result.json
      timeoutSeconds: 600
    - name: network
      description: >-
        Spearman co-occurrence network analysis on a taxa or functional
        feature x sample matrix with prevalence filtering, Benjamini-Hochberg
        FDR control, hub node ranking (degree and weighted degree), and a
        co-occurrence network figure. Pass --features to restrict the network
        to differentially abundant biomarkers.
      run: [python, ./scripts/network.py]
      args:
        type: object
        additionalProperties: false
        required: [input, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Feature x sample abundance matrix TSV.
          output_dir:
            type: string
            format: project-path
            description: Directory for network-edges.tsv, network-nodes.tsv, and figure.
          min_prevalence:
            type: number
            description: Minimum fraction of samples with abundance > 0 (default 0.2).
          top:
            type: integer
            minimum: 2
            description: Maximum top-mean-abundance features to test when --features is omitted (default 100).
          min_abs_rho:
            type: number
            description: Minimum absolute Spearman correlation to retain an edge (default 0.4).
          fdr:
            type: number
            description: Maximum BH-adjusted q-value to retain an edge (default 0.05).
          features:
            type: array
            items:
              type: string
            description: Optional explicit subset of feature names (e.g. significant biomarkers).
      approval: write
      output: ./schemas/network-result.json
      timeoutSeconds: 1200
    - name: classify
      description: >-
        Leak-free stratified k-fold cross-validated binary classification
        (case vs control) on any taxonomic or functional feature x sample
        matrix: strictly in-fold log1p + z-score scaling, in-fold mRMR or
        Mann-Whitney feature selection, and L2-regularized logistic
        regression. Writes out-of-fold predictions, feature-importance table,
        and pooled ROC curve with AUROC and Youden sensitivity/specificity.
      run: [python, ./scripts/classify.py]
      args:
        type: object
        additionalProperties: false
        required: [input, metadata, case, control, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Feature x sample abundance matrix TSV (species, gene families, or KOs).
          metadata:
            type: string
            format: input-path
            description: Sample metadata CSV/TSV with the group column.
          group_col:
            type: string
            description: Metadata column defining groups (default group).
          case:
            type: string
            description: Group label treated as positive class (y=1).
          control:
            type: string
            description: Group label treated as negative class (y=0).
          output_dir:
            type: string
            format: project-path
            description: Directory for oof-predictions.tsv, feature-importance.tsv, and roc-curve.png.
          top_features:
            type: integer
            minimum: 1
            description: Number of features selected inside each training fold (default 15).
          selector:
            type: string
            enum: [mrmr, mannwhitney]
            description: In-fold feature selector (default mrmr).
          folds:
            type: integer
            minimum: 2
            description: Number of stratified CV folds (default 5).
          l2_alpha:
            type: number
            description: L2 regularization strength for logistic regression (default 1.0).
          seed:
            type: integer
            description: RNG seed for stratified fold splitting (default 12345).
      approval: write
      output: ./schemas/classify-result.json
      timeoutSeconds: 1200
    - name: multivariable
      description: >-
        Compositionally-aware multivariable differential abundance via Centered
        Log-Ratio (CLR) transformation + OLS linear regression with covariate
        adjustment (e.g., sex, age, BMI) and Benjamini-Hochberg FDR control.
        Writes multivariable-clr.tsv and multivariable-volcano.png.
      run: [python, ./scripts/multivariable.py]
      args:
        type: object
        additionalProperties: false
        required: [input, metadata, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Feature x sample abundance matrix TSV (taxa or functional groups).
          metadata:
            type: string
            format: input-path
            description: Sample metadata CSV/TSV containing exposure and covariate columns.
          exposure:
            type: string
            description: Primary variable of interest in metadata (categorical or numeric; default group).
          reference:
            type: string
            description: Optional reference level when exposure is categorical.
          covariates:
            type: array
            items:
              type: string
            description: Optional confounding metadata columns to adjust for (e.g. [sex, age, bmi]).
          min_prevalence:
            type: number
            description: Minimum fraction of samples with abundance > 0 (default 0.1).
          pseudocount:
            type: number
            description: Optional positive pseudocount for CLR transform (default auto = 0.5 * min non-zero).
          fdr:
            type: number
            description: BH FDR significance threshold (default 0.05).
          output_dir:
            type: string
            format: project-path
            description: Directory for multivariable-clr.tsv and volcano plot.
      approval: write
      output: ./schemas/multivariable-result.json
      timeoutSeconds: 1200
---

# Read-based metagenomic profiling — downstream chain

This skill owns the **deterministic downstream** of read-based profiling:
reports → matrix → diversity / differential / network / classification / functional tables / figures.
Running classifiers on raw reads (fastp, kneaddata, kraken2, bracken,
metaphlan, humann3, krakentools, krona) is the wrappers' job; the agent
orchestrates them with `wrapper_search` / `wrapper_inspect` /
`wrapper_plan_run` / `wrapper_run`.

## Before you start

1. Read [references/analysis-playbook.md](references/analysis-playbook.md) —
   the call → read → decide → red-line loop for every step.
2. Call `mg_inspect` on the input. Never guess the report format or the value
   kind (counts vs relative abundance); the inspect result tells you, and you
   pass it on explicitly.
3. Per-step principle and interpretation notes live in
   [references/step-notes/](references/step-notes/).

## Chain

```
reports ──mg_merge──▶ matrix ──┬─ mg_rarefaction  (depth adequacy, counts only)
                               ├─ mg_diversity   (alpha, Bray-Curtis, PCoA, PERMANOVA, PERMDISP, covariates)
                               ├─ mg_differential (Wilcoxon/KW + BH FDR, volcano)
                               ├─ mg_network      (Spearman co-occurrence network + hubs)
                               ├─ mg_classify     (in-fold mRMR + L2 logistic CV + ROC/AUROC)
                               ├─ mg_plot         (barplot / heatmap, optional clustering)
                               └─ (via lefse + graphlan wrappers for cladograms)
HUMAnN tables ──mg_function──▶ normalized / regrouped / unstratified
                               ──┬─ mg_differential (values=relative)
                                 ├─ mg_classify / mg_diversity / mg_plot
MetaPhlAn rank T ──mg_merge──▶ SGB matrix + sgb_species (strain-level view)
```

## Rules

- One matrix, one value kind. Mixing Kraken2 counts with MetaPhlAn relative
  abundances in one matrix is invalid; analyze them separately and compare
  conclusions, not numbers.
- **MetaPhlAn normalization consistency**: keep `--unknown_estimation` (and
  database version) identical across all samples in a cohort; mixing profiles
  normalized to total reads vs classified reads silently skews relative
  abundances and differential tests.
- **ID-aligned merging**: always align samples and features strictly by ID
  (as `mg_merge` and `mg_diversity` do), never by row order.
- Chao1 and observed richness need counts. On MetaPhlAn output pass
  `values: relative` and report that the two indices were skipped.
- Check `permdisp` whenever `permanova` is significant: if both are p < 0.05,
  report that group separation may partly reflect unequal dispersion.
- Rarefaction is a depth check, not a normalization: run it on counts before
  comparing richness, and never feed rarefied data to `mg_differential`.
- Thresholds (FDR, fold change, top N, PERMANOVA permutations) are the agent's
  decision, passed explicitly; defaults are starting points, not verdicts.
- Every JSON return is the audit record — cite its counts, do not recompute.
- Gaps to report honestly, not fake: LEfSe LDA/cladogram (lefse + graphlan
  wrappers), Krona (krona wrapper), ANCOM-BC / MaAsLin 3 (R-based multivariable
  abundance models not in the Python environment). See the playbook's red lines.

## Sources

Distilled from EasyMetagenome (iMeta 2025), the Snakemake Meta pipeline,
cambiotraining/metagenomics (CC BY 4.0) and sib-swiss
metagenomics-data-analysis-training (GPL-3.0; ideas only, no code copied).
Provenance and citations: [references/maintainers/source-map.md](references/maintainers/source-map.md).