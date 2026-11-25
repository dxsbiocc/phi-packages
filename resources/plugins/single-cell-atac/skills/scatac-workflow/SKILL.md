---
name: scatac-workflow
description: >-
  Run a single-cell ATAC-seq chromatin-accessibility workflow on 10x Cell Ranger
  ATAC fragments.tsv.gz, a BAM, a saved .h5ad or a peak/tile matrix: inspect the
  input and genome, QC (fragment count, duplication, mitochondrial fraction, TSS
  enrichment), scrublet doublet removal, a tile matrix, feature selection,
  spectral embedding (LSI-like), kNN/Leiden clustering, UMAP, per-cluster MACS3
  peak calling and a merged peak matrix, marker regions and optional pairwise
  differential accessibility, gene activity and cell-type annotation (marker
  panel or scRNA label transfer), motif enrichment, and sample integration
  (Harmony/Scanorama/MNN). Use for 单细胞ATAC, scATAC, ATAC-seq, 染色质可及性,
  表观基因组, snapatac2, fragments, TSSe, peak calling, 差异可及性, motif 富集,
  gene activity, .fragments.gz. Python route (snapatac2 2.10 + scanpy 1.12);
  chromVAR/ArchR/Signac (R) and episcanpy (numpy 2 incompatible) are reported
  gaps. Not for raw reads (cellranger-atac/chromap wrappers) or bulk ChIP/CUT&Tag
  (the chromatin-profiling sibling).
phi:
  environment: plugin:scatac
  attachTo:
    - SingleCellAtac
  scripts:
    - name: inspect
      description: >-
        Read a 10x ATAC fragments.tsv.gz, a BAM, a saved .h5ad or a peak/tile
        matrix and report its input kind, shape, genome/chrom_sizes, which QC
        metrics are present (n_fragment, frac_dup, frac_mito, tsse), the usable
        metadata columns, and which follow-up steps the object is ready for.
        Call this first on an unfamiliar input. With output, saves the loaded
        object as .h5ad so later steps can chain.
      run: [python, ./scripts/project_inspect.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: A 10x ATAC fragments.tsv.gz, a BAM, a saved .h5ad, or a peak/tile matrix (csv/tsv/h5).
          output:
            type: string
            format: project-path
            description: Where to save the loaded object as .h5ad, inside the project.
          genome:
            type: string
            description: Bundled chrom_sizes genome (hg38, mm10, hg19, mm39); default hg38. For another genome pass chrom_sizes.
          chrom_sizes:
            type: string
            format: input-path
            description: Two-column chrom<TAB>size file for a non-bundled genome; overrides genome.
          gene_anno:
            type: string
            format: input-path
            description: Gene annotation GTF/GFF3 to enable TSS enrichment; null skips it with a gap message.
          is_paired:
            type: boolean
            description: Fragments are paired-end (default true; 10x ATAC is paired).
          min_num_fragments:
            type: integer
            minimum: 0
            description: Minimum fragments per barcode at import (default 200); below this a barcode is dropped as empty.
      approval: read
      output: ./schemas/step-result.json
      timeoutSeconds: 600
    - name: qc
      description: >-
        Import fragments and compute per-cell QC (fragment count, duplication
        fraction, mitochondrial fraction, fragment-size distribution, and TSS
        enrichment when gene_anno is supplied), apply the given fragment-count
        and TSSe thresholds, and report how many cells each gate removed plus
        the metric distributions behind a sensible choice. Writes before/after
        QC panels. Doublet removal needs the feature matrix and runs in
        scatac_preprocess, so this module stays lightweight. Start from raw
        fragments.
      run: [python, ./scripts/qc_filter.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: A scATAC input (fragments.tsv.gz / BAM / .h5ad / peak matrix) with raw fragments.
          output:
            type: string
            format: project-path
            description: Where to save the filtered object as .h5ad.
          figures_dir:
            type: string
            format: project-path
            description: Directory for the before/after QC plots; omit to skip figures.
          genome:
            type: string
            description: Bundled chrom_sizes genome (default hg38).
          chrom_sizes:
            type: string
            format: input-path
            description: chrom<TAB>size file for a non-bundled genome.
          gene_anno:
            type: string
            format: input-path
            description: GTF/GFF3 to enable TSS enrichment; null skips it with a gap message.
          is_paired:
            type: boolean
            description: Fragments are paired-end (default true).
          min_num_fragments:
            type: integer
            minimum: 0
            description: Minimum fragments per barcode at import (default 200).
          min_counts:
            type: number
            description: Minimum fragments per cell (default 1000); read the n_fragment distribution before choosing.
          max_counts:
            type: number
            description: Maximum fragments per cell; null disables the upper gate (a high shoulder can be doublets).
          min_tsse:
            type: number
            description: Minimum TSS enrichment; null disables the gate (TSSe needs gene_anno). Read the tsse distribution first.
          max_tsse:
            type: number
            description: Maximum TSS enrichment; null disables the upper gate.
          quantiles:
            type: array
            items:
              type: number
            description: Quantiles to report per metric (default 0.01 and 0.99).
          sample_column:
            type: string
            description: obs column naming each cell's sample, for a per-sample breakdown.
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 1800
    - name: preprocess
      description: >-
        The QC chain plus a genome-wide tile matrix, feature selection, scrublet
        doublet removal, spectral embedding (the scATAC analogue of LSI/PCA), a
        kNN graph, Leiden clustering at one or more resolutions and UMAP. Tiles
        are used for clustering because they cover the genome uniformly; called
        peaks come later for differential accessibility. Reports cluster sizes,
        the selected feature count and the embeddings. Reads a scATAC input.
      run: [python, ./scripts/preprocess_cluster.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: A scATAC input (fragments.tsv.gz / BAM / .h5ad / peak matrix).
          output:
            type: string
            format: project-path
            description: Where to save the processed object as .h5ad.
          figures_dir:
            type: string
            format: project-path
            description: Directory for the QC, UMAP and cluster plots.
          genome:
            type: string
            description: Bundled chrom_sizes genome (default hg38).
          chrom_sizes:
            type: string
            format: input-path
            description: chrom<TAB>size file for a non-bundled genome.
          gene_anno:
            type: string
            format: input-path
            description: GTF/GFF3 to enable TSS enrichment; null skips it.
          blacklist:
            type: string
            format: input-path
            description: Blacklist BED; null uses the shipped ENCODE blacklist for hg38/mm10, else none.
          min_num_fragments:
            type: integer
            minimum: 0
            description: Minimum fragments per barcode at import (default 200).
          bin_size:
            type: integer
            description: Tile matrix bin size in bp (default 500); a biological decision about resolution.
          n_features:
            type: integer
            description: Tiles/peaks to select for the spectral embedding (default 500000).
          min_counts:
            type: number
            description: Minimum fragments per cell (default 1000).
          max_counts:
            type: number
            description: Maximum fragments per cell; null disables the upper gate.
          min_tsse:
            type: number
            description: Minimum TSS enrichment; null disables the gate.
          max_tsse:
            type: number
            description: Maximum TSS enrichment; null disables the upper gate.
          doublet_filter:
            type: boolean
            description: Run scrublet on the selected features and filter predicted doublets (default true).
          doublet_threshold:
            type: number
            description: Scrublet doublet probability threshold (default 0.5).
          expected_doublet_rate:
            type: number
            description: Expected doublet rate for scrublet (default 0.1); scales with loaded cell number.
          n_comps:
            type: integer
            description: Spectral embedding components (default 30).
          n_neighbors:
            type: integer
            description: Neighbors for the kNN graph on X_spectral (default 50).
          resolution:
            type: array
            items:
              type: number
            description: One or more Leiden resolutions, for example 1.0 (repeat the flag for several).
          run_umap:
            type: boolean
            description: Compute the UMAP embedding (default true).
          seed:
            type: integer
            description: Random seed (default 0, the snapatac2 default).
          n_jobs:
            type: integer
            description: Worker processes (default 1 for reproducibility; raise for large fragment files).
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 3600
    - name: peaks
      description: >-
        The clustered chain plus per-cluster MACS3 peak calling, a merged
        non-overlapping peak set (summits extended by peak_half_width) and a
        separate peak-count matrix. Peaks are the open-chromatin regions that
        differ between cell states; the peak matrix drives markers while tiles
        stay the embedding. MACS3 needs fragments, so a peak-matrix input
        without fragments reports a gap here. Reports peaks per cluster and the
        merged peak count, and writes a merged BED.
      run: [python, ./scripts/call_peaks.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: A scATAC input with fragments (fragments.tsv.gz / BAM / fragment-backed .h5ad).
          output:
            type: string
            format: project-path
            description: Where to save the object as .h5ad.
          tables_dir:
            type: string
            format: project-path
            description: Directory for the merged peak BED and the peak-count matrix.
          genome:
            type: string
            description: Bundled chrom_sizes genome (default hg38).
          chrom_sizes:
            type: string
            format: input-path
            description: chrom<TAB>size file for a non-bundled genome.
          blacklist:
            type: string
            format: input-path
            description: Blacklist BED; null uses the shipped ENCODE blacklist for hg38/mm10.
          group_by:
            type: string
            description: obs column defining the clusters for peak calling (default leiden).
          peak_qvalue:
            type: number
            description: MACS3 q-value cutoff (default 0.05).
          peak_half_width:
            type: integer
            description: Half-width in bp to extend each summit before merging (default 250).
          call_broad_peaks:
            type: boolean
            description: Call broad peaks for histone marks (default false; ATAC uses narrow peaks).
          bin_size:
            type: integer
            description: Tile bin size in bp for the clustering chain (default 500).
          n_features:
            type: integer
            description: Tiles to select for the embedding (default 500000).
          resolution:
            type: array
            items:
              type: number
            description: Leiden resolution(s) (default 1.0).
          seed:
            type: integer
            description: Random seed (default 0).
          n_jobs:
            type: integer
            description: Worker processes for MACS3 (default 1).
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 3600
    - name: markers
      description: >-
        The peak chain plus per-cluster marker regions (marker_regions tests each
        cluster against the rest) written as one table per cluster, and optional
        pairwise cluster-vs-cluster differential accessibility (diff_test, a
        logistic regression giving log fold change and significance). Pairwise
        diff is opt-in via run_diff because it is O(clusters^2); marker_regions
        is always cluster-vs-rest. Both run on the peak matrix when peaks were
        called, otherwise on tiles. The p-values are descriptive (clusters were
        defined from the same data); treat the ranking as evidence, not
        significance.
      run: [python, ./scripts/find_markers.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: A scATAC input with fragments.
          output:
            type: string
            format: project-path
            description: Where to save the object as .h5ad.
          tables_dir:
            type: string
            format: project-path
            description: Directory for per-cluster marker tables and any diff CSVs.
          group_by:
            type: string
            description: obs column to group by (default leiden).
          marker_pvalue:
            type: number
            description: p-value cutoff for marker_regions (default 0.01).
          run_diff:
            type: boolean
            description: Run pairwise cluster-vs-cluster diff_test (default false; O(clusters^2), slow).
          diff_min_log_fc:
            type: number
            description: Minimum log fold change for diff_test (default 0.25).
          diff_min_pct:
            type: number
            description: Minimum fraction of cells with accessibility for a region to enter diff_test (default 0.05).
          peak_qvalue:
            type: number
            description: MACS3 q-value cutoff (default 0.05).
          peak_half_width:
            type: integer
            description: Summit half-width before merging peaks (default 250).
          genome:
            type: string
            description: Bundled chrom_sizes genome (default hg38).
          blacklist:
            type: string
            format: input-path
            description: Blacklist BED; null uses the shipped ENCODE blacklist for hg38/mm10.
          resolution:
            type: array
            items:
              type: number
            description: Leiden resolution(s) (default 1.0).
          seed:
            type: integer
            description: Random seed (default 0).
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 3600
    - name: annotate
      description: >-
        The clustered chain plus a gene activity matrix (accessibility aggregated
        over gene bodies and promoters) and cell-type annotation, by scoring a
        marker panel on gene activity, transferring labels from a scRNA
        reference, or applying an explicit cluster->label mapping. Gene activity
        is a proxy for expression, not RNA, and needs a gene annotation; without
        gene_anno, reference_path or mapping_path the step reports the gap and
        continues (a script must not invent labels). Writes the label column
        into obs. The reference comes from the single-cell plugin.
      run: [python, ./scripts/annotate_clusters.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: A scATAC input with fragments.
          output:
            type: string
            format: project-path
            description: Where to save the annotated object as .h5ad.
          tables_dir:
            type: string
            format: project-path
            description: Directory for the gene activity matrix.
          figures_dir:
            type: string
            format: project-path
            description: Directory for the annotation dot plot.
          gene_anno:
            type: string
            format: input-path
            description: GTF/GFF3 for the gene activity matrix; null skips gene activity with a gap message.
          gene_upstream:
            type: integer
            description: Promoter upstream bp included in gene activity (default 2000).
          gene_body:
            type: boolean
            description: Include the gene body in gene activity (default true).
          panel_path:
            type: string
            format: input-path
            description: JSON object mapping cell type to an array of gene symbols; null uses a generic example panel.
          mapping_path:
            type: string
            format: input-path
            description: JSON object mapping every cluster to a label; overrides panel scoring and reference transfer.
          reference_path:
            type: string
            format: input-path
            description: scRNA reference .h5ad for label transfer; null skips transfer.
          reference_label_key:
            type: string
            description: obs column in the reference naming cell types (default cell_type).
          group_by:
            type: string
            description: Cluster column (default leiden).
          label_column:
            type: string
            description: obs column to write labels into (default cell_type).
          min_score:
            type: number
            description: Minimum panel score to accept a label (default 0).
          genome:
            type: string
            description: Bundled chrom_sizes genome (default hg38).
          resolution:
            type: array
            items:
              type: number
            description: Leiden resolution(s) (default 1.0).
          seed:
            type: integer
            description: Random seed (default 0).
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 3600
    - name: motif
      description: >-
        The marker chain plus motif enrichment on per-cluster marker regions
        (which TF motifs are over-represented in each cluster's open chromatin),
        with a binomial or hypergeometric test against a background. Needs a
        genome FASTA and a motif database; without genome_fasta the step reports
        the gap and continues. motif_db null uses snapatac2 cis_bp (a download).
        This is sequence-level prediction, not measured TF binding; chromVAR
        per-cell TF deviations are an R-route gap, not implemented here.
      run: [python, ./scripts/motif_enrichment.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: A scATAC input with fragments.
          output:
            type: string
            format: project-path
            description: Where to save the object as .h5ad.
          tables_dir:
            type: string
            format: project-path
            description: Directory for per-cluster motif enrichment CSVs.
          genome_fasta:
            type: string
            format: input-path
            description: Genome FASTA to extract region sequences; null skips motif with a gap message.
          motif_db:
            type: string
            format: input-path
            description: Motif database .meme; null uses snapatac2 cis_bp (a download).
          motif_method:
            type: string
            enum: [binomial, hypergeometric]
            description: Motif enrichment test (default binomial).
          motif_pvalue:
            type: number
            description: p-value cutoff for motif enrichment (default 1e-5).
          group_by:
            type: string
            description: obs column defining the clusters (default leiden).
          marker_pvalue:
            type: number
            description: p-value cutoff for the marker regions motif runs on (default 0.01).
          genome:
            type: string
            description: Bundled chrom_sizes genome (default hg38).
          resolution:
            type: array
            items:
              type: number
            description: Leiden resolution(s) (default 1.0).
          seed:
            type: integer
            description: Random seed (default 0).
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 3600
    - name: integrate
      description: >-
        The clustered chain plus multi-sample batch integration on the spectral
        embedding with Harmony, Scanorama or MNN. Needs a batch_column naming
        each cell's sample and at least two batches; otherwise the step reports
        the gap and continues. Integration is for joint clustering and
        visualization only — corrected values never enter peak calling or
        differential accessibility. The sample is the experimental unit.
      run: [python, ./scripts/integrate_samples.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: A scATAC input with fragments and a sample/batch obs column.
          output:
            type: string
            format: project-path
            description: Where to save the integrated object as .h5ad.
          batch_column:
            type: string
            description: obs column naming each cell's sample/batch (default sample).
          integrate_method:
            type: string
            enum: [harmony, scanorama, mnc]
            description: Batch integration method on X_spectral (default harmony).
          group_by:
            type: string
            description: obs column defining the clusters (default leiden).
          n_comps:
            type: integer
            description: Spectral embedding components (default 30).
          n_neighbors:
            type: integer
            description: Neighbors for the kNN graph (default 50).
          resolution:
            type: array
            items:
              type: number
            description: Leiden resolution(s) (default 1.0).
          genome:
            type: string
            description: Bundled chrom_sizes genome (default hg38).
          seed:
            type: integer
            description: Random seed (default 0).
      approval: write
      output: ./schemas/step-result.json
      timeoutSeconds: 3600
    - name: notebook
      description: >-
        Render a runnable scATAC analysis notebook into the project from the
        fragment library, then open it with Phi's notebook tools and run cells
        on demand. The template is rendered, not copied: placeholders are filled
        from the input data and the requested steps, and every tunable value
        lands in one config cell at the top that the analysis cells reference.
        Pick a preset with module (1 QC, 2 preprocess+cluster, 3 peaks+markers,
        4 gene activity+annotation, 5 motif+integration, 6 full chain) or name
        steps directly; dependencies are added automatically. inject adds an
        external Python script as a code cell; merge splices in another
        notebook, whole or a cell range. Starts from fragments.tsv.gz, a BAM, a
        .h5ad or a peak matrix: raw-read processing (cellranger-atac, chromap)
        is an upstream wrapper, not notebook cells. Use self_check to execute
        every cell in the output directory now and fail loudly if one breaks.
      run: [python, ./scripts/build_notebook.py]
      args:
        type: object
        additionalProperties: false
        required: [output, data]
        properties:
          module:
            type: integer
            enum: [1, 2, 3, 4, 5, 6]
            description: scATAC preset. Pass module or steps, not both.
          steps:
            type: array
            items:
              type: string
              enum: [setup, blacklist, load, qc-metrics, qc-plots, filter, features, embed, cluster, peaks, markers, gene-activity, annotate, motif, integrate, save]
            description: >-
              Steps to render instead of a preset; declared dependencies are
              added automatically. Pass module or steps, not both.
          output:
            type: string
            format: project-path
            description: Path of the .ipynb to write inside the project.
          data:
            type: string
            format: input-path
            description: A fragments.tsv.gz, BAM, .h5ad or peak/tile matrix the notebook starts from.
          language:
            type: string
            enum: [python]
            description: Notebook language and kernel (python only; the R/ArchR/Signac route is a future sibling).
          inject:
            type: array
            items:
              type: string
              format: input-path
            description: Python (.py) scripts to add verbatim as code cells, in order.
          merge:
            type: array
            items:
              type: string
            description: >-
              Notebooks to splice in after the rendered steps, in order. Whole
              notebook (path.ipynb) or a cell range (path.ipynb[0:6]). Must be python.
          config:
            type: string
            format: input-path
            description: JSON file of parameter overrides; keys mirror the parameter names in assets/steps.json.
          output_object:
            type: string
            description: Where the notebook saves its final object, relative to the notebook folder (default scatac.h5ad).
          sample_column:
            type: string
            description: obs column naming each cell's sample (default sample).
          genome:
            type: string
            description: Bundled chrom_sizes genome (hg38, mm10, hg19, mm39; default hg38).
          chrom_sizes:
            type: string
            format: input-path
            description: chrom<TAB>size file for a non-bundled genome.
          gene_anno:
            type: string
            format: input-path
            description: GTF/GFF3 for TSSe and gene activity; null skips them with a gap message.
          genome_fasta:
            type: string
            format: input-path
            description: Genome FASTA for motif enrichment; null skips motif with a gap message.
          blacklist:
            type: string
            format: input-path
            description: Blacklist BED; null uses the shipped ENCODE blacklist for hg38/mm10.
          motif_db:
            type: string
            format: input-path
            description: Motif database .meme; null uses snapatac2 cis_bp (a download).
          min_num_fragments:
            type: integer
            description: Minimum fragments per barcode at import (default 200).
          bin_size:
            type: integer
            description: Tile bin size in bp (default 500).
          n_features:
            type: integer
            description: Tiles/peaks to select for the embedding (default 500000).
          min_counts:
            type: number
            description: Minimum fragments per cell (default 1000).
          min_tsse:
            type: number
            description: Minimum TSS enrichment; null disables the gate.
          doublet_filter:
            type: boolean
            description: Run scrublet doublet removal (default true).
          n_comps:
            type: integer
            description: Spectral embedding components (default 30).
          n_neighbors:
            type: integer
            description: Neighbors for the kNN graph (default 50).
          resolution:
            type: number
            description: Clustering resolution (default 1.0).
          run_umap:
            type: boolean
            description: Compute the UMAP embedding (default true).
          group_by:
            type: string
            description: obs column defining the clusters (default leiden).
          peak_qvalue:
            type: number
            description: MACS3 q-value cutoff (default 0.05).
          marker_pvalue:
            type: number
            description: p-value cutoff for marker_regions (default 0.01).
          run_diff:
            type: boolean
            description: Run pairwise differential accessibility (default false).
          reference_path:
            type: string
            format: input-path
            description: scRNA reference .h5ad for label transfer; null skips it.
          panel_path:
            type: string
            format: input-path
            description: Marker panel JSON for gene-activity scoring; null uses a generic example panel.
          batch_column:
            type: string
            description: obs column naming each cell's sample/batch (default sample).
          integrate_method:
            type: string
            enum: [harmony, scanorama, mnc]
            description: Batch integration method (default harmony).
          seed:
            type: integer
            description: Random seed for reproducible clustering and embedding (default 0).
          kernel:
            type: string
            description: Kernel name written to the notebook (default plugin-single-cell-atac-scatac).
          self_check:
            type: boolean
            description: Execute every cell once in the output directory and fail if one breaks.
      approval: write
      output: ./schemas/notebook-result.json
      timeoutSeconds: 1800
    - name: env_probe
      description: >-
        Report which scATAC Python packages the current environment provides,
        list the packages each workflow module needs, and say which modules can
        run here and which optional sections cannot. Call this before promising
        motif enrichment, gene activity or integration. The plugin environment
        is read-only and prebuilt, so a missing package needs a project
        environment, not an in-environment install.
      run: [python, ./scripts/env_probe.py]
      args:
        type: object
        additionalProperties: false
        properties:
          module:
            type: integer
            description: Restrict the report to this module number (1-6).
          write_panel:
            type: string
            format: project-path
            description: Write a starter marker-panel JSON here for the annotate step.
      approval: read
      output: ./schemas/env-probe.json
      timeoutSeconds: 300
---

# single-cell ATAC workflow

Run a reproducible scATAC analysis chain on one object at a time. Every step
reads a scATAC input (fragments.tsv.gz, BAM, .h5ad or peak/tile matrix) and
writes one `.h5ad`, so steps chain, can be re-run in isolation, and record the
parameters they used. This is the Python route of the scATAC plugin (snapatac2
2.10 + scanpy 1.12); the R/ArchR/Signac route is a future sibling and does not
exist here. snapatac2 is the single engine; the analysis structure and teaching
text are an original distillation of the snapatac2 documentation (MIT), and the
code is original to this plugin and runs in the private `plugin:scatac`
environment.

## Before the first step

Call `scatac_inspect` on an unfamiliar input: it reports the input kind
(fragment / BAM / .h5ad / peak matrix), shape, genome/chrom_sizes, which QC
metrics are present, and which steps the object is ready for. Call
`scatac_env_probe` before promising motif enrichment, gene activity or
integration. The plugin environment ships `snapatac2`, `scanpy`, `anndata`,
`polars`, `pyarrow`, `macs3`, `leidenalg`, `python-igraph`, `umap-learn`,
`scikit-learn`, `scipy`, `numpy`, `pandas`, `h5py`, `plotly`, `matplotlib`,
`kaleido` and `ipykernel`. It does **not** ship `episcanpy` (breaks under
numpy 2), the R suites `chromVAR` / `ArchR` / `Signac`, or the pip-only
`scvi-tools` / `magic-impute`. Report a gap as a fact and name the next step (a
project environment via `env_request`); never install packages inside the
read-only environment.

## The chain

1. `scatac_inspect` — what input/genome is this? Which steps is it ready for?
2. `scatac_qc` — fragment count, duplication, mitochondrial fraction, TSSe, thresholds, distributions.
3. `scatac_preprocess` — tile matrix, feature selection, scrublet, spectral embedding, Leiden clusters, UMAP.
4. `scatac_peaks` — per-cluster MACS3 peaks, a merged peak set and a peak-count matrix.
5. `scatac_markers` — marker regions per cluster; optional pairwise differential accessibility (run_diff).
6. `scatac_annotate` — gene activity and labels from a panel, a scRNA reference or an explicit mapping.
7. `scatac_motif` — motif enrichment on per-cluster marker regions.
8. `scatac_integrate` — multi-sample batch integration on the spectral embedding.

Each tool re-runs its chain from the input you give it (snapatac2 objects are
file-backed and heavy to hand off mid-chain), so pass the same raw input and the
same thresholds to keep results consistent. Skip a step only when its result is
already in the object; `scatac_inspect` says which results are present. Peaks
need fragments; gene activity and TSSe need `gene_anno`; motif needs
`genome_fasta`; annotation needs a panel, reference or mapping; integration
needs a batch column with at least two samples. A missing resource is reported
as a gap and the chain continues.

## Rendering a notebook

`scatac_notebook` renders the fragment library (`assets/fragments/`) into one
finished notebook in the project directory. It is the same deterministic
renderer as the single-cell and spatial skills, driven by this skill's
`steps.json`: steps, parameter values, and extra sources arrive as arguments,
and placeholders are filled from them and from the data path. It does not choose
steps or thresholds and does not interpret the user's request — those decisions
are yours (see Rules).

- `module` renders a preset (1-6); `steps` renders the listed steps plus their
  declared dependencies (`assets/steps.json`).
- Every tunable value lands in the first code cell (the config cell) with a
  source tag; the analysis cells reference those variables.
- `inject` adds a script verbatim as a code cell; `merge` splices in another
  notebook whole or in part (`path.ipynb[0:6]`). Sources must be Python — one
  kernel per notebook.
- `self_check: true` executes every cell in the output directory before
  delivery, with the same relative paths the user will run with.

## Rules

- **The analysis unit is the cell (barcode), not a spot/bin.** scATAC has no
  spatial coordinates; its "space" is the genome — accessibility lives on
  peaks/tiles. One barcode holds one chromatin sample; a doublet inflates
  accessibility, so scrublet matters. Read `scatac_inspect`'s kind report and
  say which input you analysed.
- **The renderer fills; the agent decides.** Which steps to render, which
  thresholds to pass, which tile bin and peak width, whether to call peaks
  per-cluster, how to edit the marker panel, whether to run pairwise diff —
  these are judgements you make from the user's question and the data, then pass
  as arguments. A default that reached the notebook untouched is a decision you
  did not make.
- **One reading, both paths.** Pass the same thresholds read from `scatac_qc` to
  the script chain (as flags) and to `scatac_notebook` (as flags or a `config`
  JSON), so the two paths agree by construction, not by luck.
- **Thresholds are arguments, not guesses.** Read `distribution` from
  `scatac_qc` and pass a value the data supports; say which number you chose and
  why. A threshold that removes most cells is a warning (`kept_fraction` below
  half), not a result. TSSe is the ATAC-specific signal-to-noise metric: a cell
  with enough fragments but low TSSe can still be noise.
- **Raw accessibility comes first.** Peak calling (MACS3), differential
  accessibility and scrublet run on fragments/counts. Corrected, scaled or
  integrated values never enter peak calling or differential tests; integration
  is for joint clustering and visualization only.
- **Tiles for clustering, peaks for comparison.** The tile matrix bins the
  genome uniformly for the embedding; called peaks are the features of record
  for marker regions and differential accessibility. Tile bin size and peak
  half-width are biological decisions about resolution — say what you chose and
  why.
- **backed AnnData must be closed.** snapatac2 objects are backed on HDF5 and
  auto-sync; every backed object is closed before the process exits (the script
  chain saves and closes; the notebook ends with a close cell). snapatac2
  figures are plotly Figures written as interactive `.html`, not matplotlib
  Axes; gene-activity downstream uses scanpy/matplotlib.
- **Gene activity is a proxy, not measured RNA.** It is inferred from
  promoter/gene-body accessibility, sparse and biased. Use it for annotation and
  marker visualization, not differential-expression claims; annotation is
  inference, not observation.
- **Label transfer needs a scRNA reference.** The reference comes from the
  single-cell plugin (shared counts+obs or pseudobulk handoff). State its
  cell-type granularity and source; without a reference do not claim label
  transfer — only gene-activity marker inference, labelled as such.
- **Motif enrichment is prediction, not binding.** It needs a genome FASTA and a
  motif database; report which database (e.g. JASPAR/cis_bp) and background you
  used. Motif hits are sequence-level predictions requiring experimental
  validation. chromVAR per-cell TF deviations are an R-route gap, not shipped.
- **Cross-sample comparisons need replicates.** The sample is the experimental
  unit, not the cell. Fewer than two biological replicates per condition means no
  condition-level claim. Do not over-correct biology away during integration.

## Reading the results

The JSON each tool returns is the audit record: shapes before and after, the
filters that fired with their cell counts, the doublet rate, cluster sizes,
spectral variance, peaks per cluster and the merged peak count, marker/diff
rankings, motif enrichment, integration info, and every note about what could
not run. Quote those numbers in your answer instead of re-deriving them. A step
carrying a "skipped" note did not run — never read its absence as a clean
result.

For the per-stage decision procedure — which tool to call, which JSON fields and
figures to read, how to set thresholds and peak/tile strategy, and the red lines
— read [analysis-playbook.md](references/analysis-playbook.md).
[references/red-line-cheatsheet.md](references/red-line-cheatsheet.md) compresses
each red line into a quick self-check table.
[references/step-notes/](references/step-notes/) carries the per-step
learner-facing short version (the notebook's teaching cells), and
[references/maintainers/](references/maintainers/) the implementation notes.