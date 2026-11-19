# Differential expression

Use this Phi plugin together with the `local/modules/differential-expression-deseq2` wrapper. The plugin runs locally in `phi:python@1`. The wrapper supplies DESeq2 and its R dependencies.

1. Call `de_preflight` with the sample sheet, the user’s study design, both contrast levels and their column names. A paired comparison needs a sample-to-donor mapping. If the plugin returns `needs_input`, the main agent presents its options and waits. Once the design is valid, pass the returned wrapper options together with the approved counts, sample sheet and output directory to Wrapper.
2. Call `de_assess` with the published results TSV and its matching QC file. It verifies the checksum and independently counts rows with both an evaluable adjusted p-value and fold change. Hits satisfy `padj < padjThreshold` and `abs(log2FoldChange) >= log2fcThreshold`. A valid empty screen ends candidate-list ORA. An unreadable, inconsistent or entirely unevaluable result blocks interpretation.
3. For an approved enrichment request with hits, call `de_ora` with the same two files, the study’s GMT annotation and a project output path. It verifies the inputs again, tests eligible gene sets against the annotated evaluable universe, and writes a hypergeometric ORA table with Benjamini–Hochberg adjustment across all tested sets. The tool returns summary counts and the file reference.

Run the offline regression tests with `python3 -m unittest discover -s tests` from this directory.
