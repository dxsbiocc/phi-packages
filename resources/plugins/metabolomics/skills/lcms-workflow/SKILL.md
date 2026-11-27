---
name: lcms-workflow
description: >-
  Downstream analysis of untargeted LC-MS feature tables: inspect xcms /
  MZmine / MS-DIAL / MetaboAnalyst exports, QC-sample RSD filtering,
  QC-RLSC LOESS injection-order drift correction (mb_drift_correct),
  imputation, normalization (TIC / median / PQN / internal standard) and log
  transform, PCA and PLS-DA with VIP and permutation Q2, univariate tests
  (Welch / Wilcoxon / Kruskal-Wallis + BH FDR) with volcano plots,
  deterministic m/z annotation (isotopes, adducts, neutral losses,
  pseudospectra), MS2 library matching (MSP / MGF cosine), pathway
  enrichment (ID hypergeometric or mummichog-style m/z permutation), and
  clustered heatmaps. Use for 代谢组学, metabolomics, LC-MS, 特征表, xcms,
  mzmine, ms-dial, metaboanalyst, QC RSD, PQN, PLS-DA, VIP, 加合物, adduct,
  MS2, 通路富集, mummichog, or when a feature x sample intensity table needs
  statistics and figures. Upstream peak picking (msconvert, xcms, MZmine,
  MS-DIAL) is not in this plugin; report it as a gap.
phi:
  environment: plugin:metabolomics
  attachTo:
    - Metabolomics
  scripts:
    - name: inspect
      description: >-
        Identify metabolomics inputs: feature tables (xcms, MZmine, MS-DIAL,
        MetaboAnalyst, or generic feature x sample), sample metadata, MSP/MGF
        spectral libraries, pathway mappings, compound mass databases, m/z or
        ID lists. Reports per-file kind, feature/sample counts, mz/rt columns
        and units, missing rate, QC/blank samples, and which follow-up tools
        the input supports. Call this first on unfamiliar input.
      run: [python, ./scripts/project_inspect.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: A table/library file, or a directory of such files.
      approval: read
      output: ./schemas/inspect-result.json
      timeoutSeconds: 300
    - name: qc
      description: >-
        QC-based feature filtering: keep features whose RSD across QC/pool
        samples is below a threshold (default 30%) and whose detection rate
        across biological samples is at least a minimum (default 0.8).
        Requires QC samples (names containing qc/pool); fails without them.
        Writes the filtered matrix, a per-feature QC report, and an RSD
        distribution figure. Run on raw intensities before mb_preprocess.
      run: [python, ./scripts/qc_filter.py]
      args:
        type: object
        additionalProperties: false
        required: [input, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Feature x sample intensity table (mz/rt columns optional).
          output_dir:
            type: string
            format: project-path
            description: Directory for the filtered matrix, report, and figure.
          rsd_threshold:
            type: number
            minimum: 0
            description: Maximum QC RSD in percent (default 30).
          detection_rate:
            type: number
            minimum: 0
            maximum: 1
            description: Minimum non-missing fraction across biological samples (default 0.8).
      approval: write
      output: ./schemas/qc-result.json
      timeoutSeconds: 600
    - name: drift_correct
      description: >-
        QC-RLSC (Quality Control-based Robust LOESS Signal Correction) for
        long LC-MS runs with injection-order drift: fits pure-NumPy tricube
        + bisquare IRLS LOESS per feature on pooled QC injections sorted by
        run order, rescues features whose QC RSD drops below cutoff, and writes
        the corrected matrix, per-feature RSD metrics, and PCA/RSD QC plot.
      run: [python, ./scripts/drift_correct.py]
      args:
        type: object
        additionalProperties: false
        required: [input, sample_sheet, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Feature x sample intensity table (mz/rt columns optional).
          sample_sheet:
            type: string
            format: input-path
            description: Metadata CSV/TSV with sample ID and injection_order columns.
          output_dir:
            type: string
            format: project-path
            description: Directory for drift-corrected matrix, metrics, and figure.
          qc_label:
            type: string
            description: Sample-type label identifying pooled QC injections (default QC).
          span:
            type: number
            minimum: 0.25
            maximum: 1.0
            description: LOESS smoothing span in [0.25, 1.0] (default 0.75).
          rsd_cutoff:
            type: number
            minimum: 0
            description: QC RSD cutoff in percent for precision rescue (default 30).
      approval: write
      output: ./schemas/drift-correct-result.json
      timeoutSeconds: 600
    - name: preprocess
      description: >-
        Deterministic preprocessing in fixed order impute -> normalize ->
        transform. Impute missing values (NaN and zero; none / half_min /
        min / zero), normalize samples (none / tic / median / pqn / is; pqn
        uses the median-of-ratios to a reference sample, is divides by an
        internal standard column), then transform (none / log2 / log10).
        Writes the processed matrix and a per-sample normalization report.
        Feed the result to mb_multivariate, mb_univariate, and mb_plot.
      run: [python, ./scripts/preprocess.py]
      args:
        type: object
        additionalProperties: false
        required: [input, output]
        properties:
          input:
            type: string
            format: input-path
            description: Feature x sample intensity table (from mb_qc or raw).
          output:
            type: string
            format: project-path
            description: Output processed matrix TSV.
          impute:
            type: string
            enum: [none, half_min, min, zero]
            description: Missing-value imputation (default none); NaN and zero both count as missing. half_min = half the feature's minimum observed value.
          normalize:
            type: string
            enum: [none, tic, median, pqn, is]
            description: Sample normalization (default none).
          is_col:
            type: string
            description: Internal-standard sample/feature column name; required when normalize=is.
          transform:
            type: string
            enum: [none, log2, log10]
            description: Log transform applied after normalization (default none).
      approval: write
      output: ./schemas/preprocess-result.json
      timeoutSeconds: 600
    - name: multivariate
      description: >-
        PCA (SVD) and PLS-DA (NIPALS, one-hot group Y) on a processed feature
        x sample matrix with scaling (none / uv / pareto). Reports PCA
        variance explained, PLS-DA R2Y/Q2 cum (LOO cross-validation, 7-fold
        when n>50), VIP scores, seeded label-permutation Q2 null, and score
        plots colored by group. OPLS-DA is not implemented; report it as a
        gap (ropls/mixOmics belong in a project environment).
      run: [python, ./scripts/multivariate.py]
      args:
        type: object
        additionalProperties: false
        required: [input, metadata, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Processed feature x sample matrix (from mb_preprocess).
          metadata:
            type: string
            format: input-path
            description: Sample metadata CSV/TSV with the group column.
          group_col:
            type: string
            description: Metadata column defining groups (default group).
          output_dir:
            type: string
            format: project-path
            description: Directory for scores, VIP, model stats, and figures.
          scaling:
            type: string
            enum: [none, uv, pareto]
            description: Feature scaling before fitting (default uv).
          n_components:
            type: integer
            minimum: 1
            description: Components for PCA display and the PLS-DA model (default 2).
          n_perm:
            type: integer
            minimum: 0
            description: Label permutations for the Q2 null (default 99; 0 skips).
          seed:
            type: integer
            description: RNG seed for permutations and folds (default 12345).
      approval: write
      output: ./schemas/multivariate-result.json
      timeoutSeconds: 1800
    - name: univariate
      description: >-
        Per-feature group comparison on a processed matrix: Welch t-test or
        two-sided Mann-Whitney U for two groups (case vs control),
        Kruskal-Wallis for more, Benjamini-Hochberg FDR across features,
        fold change (2^mean difference when logged), and a volcano plot.
        Pass logged=true when the matrix is log-transformed. Writes a ranked
        table and the volcano figure.
      run: [python, ./scripts/univariate.py]
      args:
        type: object
        additionalProperties: false
        required: [input, metadata, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Processed feature x sample matrix.
          metadata:
            type: string
            format: input-path
            description: Sample metadata CSV/TSV with the group column.
          group_col:
            type: string
            description: Metadata column defining groups (default group).
          case:
            type: string
            description: Group label treated as case (required for two-group tests).
          control:
            type: string
            description: Group label treated as control (required for two-group tests).
          test:
            type: string
            enum: [welch, wilcoxon, kruskal]
            description: Statistical test (default wilcoxon).
          logged:
            type: boolean
            description: True when input is log-transformed; FC = 2^mean difference (default false).
          output_dir:
            type: string
            format: project-path
            description: Directory for the table and volcano plot.
          fdr:
            type: number
            description: Significance threshold on BH FDR (default 0.05).
          fc_threshold:
            type: number
            description: Minimum absolute log2 fold change (default 1).
      approval: write
      output: ./schemas/univariate-result.json
      timeoutSeconds: 1200
    - name: annotate
      description: >-
        Deterministic m/z relationship annotation on a feature table with
        mz/rt columns: 13C isotopes (delta 1.003355), adducts ([M+H]+,
        [M+Na]+, [M+K]+, [M+NH4]+, [M-H]-, [M+Cl]-, [M+FA-H]-), and neutral
        losses (H2O, CO2) within a ppm tolerance and RT window; features
        sharing a neutral mass and RT neighborhood are grouped into
        pseudospectra. This is annotation (MSI level 3), not identification.
        Writes the relationship table and pseudospectrum groups.
      run: [python, ./scripts/annotate.py]
      args:
        type: object
        additionalProperties: false
        required: [input, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Feature table with mz and rt columns (vendor aliases recognized).
          output_dir:
            type: string
            format: project-path
            description: Directory for the annotation tables.
          mode:
            type: string
            enum: [pos, neg]
            description: Ionization mode; selects the adduct set (default pos).
          ppm:
            type: number
            minimum: 0
            description: Mass tolerance in ppm (default 10).
          rt_tol:
            type: number
            minimum: 0
            description: >-
              RT window in minutes (default 0.2). If the table's median RT
              exceeds 60 it is treated as seconds and the window converted,
              unless this value is passed explicitly.
      approval: write
      output: ./schemas/annotate-result.json
      timeoutSeconds: 600
    - name: ms2_match
      description: >-
        Match query MS2 spectra (MGF) against a spectral library (MSP or
        MGF; NIST/MoNA/RIKEN MSP dialects tolerated): precursor m/z within
        ppm tolerance, fragment matching within a Da tolerance, sqrt
        intensity transform, cosine similarity. Reports the best hits per
        query with matched-fragment counts. This is MSI level 2 (putative
        library match); level 1 needs an authentic standard (RT + MS2).
      run: [python, ./scripts/ms2_match.py]
      args:
        type: object
        additionalProperties: false
        required: [query, library, output_dir]
        properties:
          query:
            type: string
            format: input-path
            description: Query spectra MGF (TITLE = feature ID, PEPMASS = precursor m/z).
          library:
            type: string
            format: input-path
            description: Reference spectral library, MSP or MGF.
          output_dir:
            type: string
            format: project-path
            description: Directory for the match table.
          ppm:
            type: number
            minimum: 0
            description: Precursor m/z tolerance in ppm (default 10).
          fragment_tol:
            type: number
            minimum: 0
            description: Fragment m/z tolerance in Da (default 0.02).
          min_score:
            type: number
            minimum: 0
            maximum: 1
            description: Minimum cosine score to report a hit (default 0.7).
          top:
            type: integer
            minimum: 1
            description: Maximum hits kept per query (default 5).
      approval: write
      output: ./schemas/ms2-match-result.json
      timeoutSeconds: 1200
    - name: pathway
      description: >-
        Pathway enrichment without a built-in database. mode=ids: compound
        ID list against a user-supplied pathway -> compound mapping TSV,
        one-sided hypergeometric (Fisher) test with BH FDR, background = all
        compounds in the mapping. mode=masses: mummichog-style — significant
        m/z list matched to a user-supplied compound mass database via
        adduct hypotheses, then pathway overlap; seeded permutation of
        random m/z lists drawn from a background feature list gives
        empirical p-values with BH FDR. Writes the enrichment table and a
        top-pathway barplot.
      run: [python, ./scripts/pathway_enrich.py]
      args:
        type: object
        additionalProperties: false
        required: [mode, input, output_dir]
        properties:
          mode:
            type: string
            enum: [ids, masses]
            description: Enrichment strategy.
          input:
            type: string
            format: input-path
            description: Significant compound IDs (ids) or m/z values (masses), one per line or first TSV column.
          output_dir:
            type: string
            format: project-path
            description: Directory for the enrichment table and figure.
          mapping:
            type: string
            format: input-path
            description: Two-column TSV pathway -> compound (header optional); required for both modes.
          compound_db:
            type: string
            format: input-path
            description: TSV with compound_id and mass columns (name optional); required for mode=masses.
          background:
            type: string
            format: input-path
            description: All-feature m/z list for permutation draws; required for mode=masses.
          ion_mode:
            type: string
            enum: [pos, neg]
            description: Ionization mode for adduct hypotheses in mode=masses (default pos).
          ppm:
            type: number
            minimum: 0
            description: Mass tolerance in ppm for mode=masses (default 10).
          n_perm:
            type: integer
            minimum: 1
            description: Permutations for mode=masses empirical p (default 100).
          seed:
            type: integer
            description: RNG seed for permutations (default 12345).
      approval: write
      output: ./schemas/pathway-result.json
      timeoutSeconds: 1800
    - name: plot
      description: >-
        Clustered heatmap of a feature x sample matrix: top N features by
        variance (default 50), optional log10 transform with pseudocount,
        and average-linkage hierarchical clustering of rows/columns
        (Euclidean on displayed values). Pair clustering with log10 for
        intensity data.
      run: [python, ./scripts/plot_heatmap.py]
      args:
        type: object
        additionalProperties: false
        required: [input, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: Feature x sample matrix TSV.
          output_dir:
            type: string
            format: project-path
            description: Directory for the figure.
          top:
            type: integer
            minimum: 1
            description: Number of features to show, by variance (default 50).
          trans:
            type: string
            enum: [none, log10]
            description: Value transform (default none).
          pseudocount:
            type: number
            description: Pseudocount for log10 (default 1e-6).
          cluster:
            type: string
            enum: [none, rows, cols, both]
            description: Hierarchical clustering of rows/columns (default none).
      approval: write
      output: ./schemas/plot-result.json
      timeoutSeconds: 600
---

# Untargeted LC-MS metabolomics — downstream chain

This skill owns the **deterministic downstream** of untargeted LC-MS:
feature table → QC → preprocessing → statistics / annotation / figures.
Upstream work (msconvert conversion, xcms / MZmine / MS-DIAL peak picking
and alignment) is **not** in this plugin and has no wrapper here; report it
as a gap and ask the user for the exported feature table.

## Before you start

1. Read [references/analysis-playbook.md](references/analysis-playbook.md) —
   the call → read → decide → red-line loop for every step.
2. Call `mb_inspect` on the input. Never guess the export format, the RT
   unit (xcms seconds vs MZmine minutes), or which samples are QC; the
   inspect result tells you, and you pass decisions on explicitly.
3. Per-step principle and interpretation notes live in
   [references/step-notes/](references/step-notes/).

## Chain

```
feature table ──[mb_drift_correct]──▶ ──mb_qc──▶ filtered ──mb_preprocess──▶ processed matrix
  (QC-RLSC LOESS drift)      (QC RSD + detection)    (impute → normalize → transform)
processed ──┬─ mb_multivariate  (PCA, PLS-DA, VIP, permutation Q2)
            ├─ mb_univariate    (Welch/Wilcoxon/KW + BH FDR, volcano)
            └─ mb_plot          (clustered heatmap)
feature table (mz/rt) ──mb_annotate──▶ isotopes / adducts / pseudospectra
MS2 (MGF) + library (MSP/MGF) ──mb_ms2_match──▶ cosine hits (MSI level 2)
significant IDs + mapping ──mb_pathway(ids)──▶ hypergeometric + BH
significant m/z + db + background + mapping ──mb_pathway(masses)──▶ permutation + BH
```

## Rules

- QC filtering needs QC samples. Without them, do not fake a filter: report
  the gap and let the user decide (detection-rate-only filtering is a
  conscious choice, not a default).
- One matrix, one processing state. Statistics run on the processed matrix;
  annotation runs on the original mz/rt table. Never log-transform before
  QC RSD, and pass `logged: true` to `mb_univariate` after a log transform.
- Thresholds (QC RSD 30%, detection 0.8, FDR 0.05, |log2FC| 1, VIP > 1,
  10 ppm precursor / 20 ppm fragment starting points) are the agent's
  decision, passed explicitly; defaults are starting points, not verdicts.
- Identification levels follow Schymanski 2014: m/z relationships are
  level 3 (annotated), library MS2 matches are level 2 (putative); level 1
  requires an authentic standard measured on the same method. Say which
  level every claim rests on.
- No built-in compound/pathway/spectral databases. KEGG/HMDB mappings,
  compound mass lists, and MSP/MGF libraries come from the user; confirm
  licensing and provenance before use, and cite them in the report.
- Every JSON return is the audit record — cite its counts, do not recompute.
- Gaps to report honestly, not fake: upstream peak picking, OPLS-DA,
  CAMERA-style full pseudospectra with peak-shape correlation, SIRIUS /
  MetFrag / GNPS molecular networking, GC-MS retention indices, NMR, and
  targeted MRM quantification. See the playbook's red lines.

## Sources

Distilled from the RforMassSpectrometry *Metabolomics Book*
(CC-BY-4.0; rformassspectrometry.github.io/metaRbolomics-book) and standard
practice (xcms, CAMERA, mummichog, MetaboAnalyst). Provenance and citations:
[references/maintainers/source-map.md](references/maintainers/source-map.md).