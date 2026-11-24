---
name: amplicon-workflow
description: >-
  Downstream phylogenetic and community analysis for 16S/ITS amplicon ASV
  tables and rooted phylogenetic trees exported from nf-core/ampliseq or
  QIIME 2/DADA2: exact Weighted and Unweighted UniFrac beta diversity,
  PCoA ordination, PERMANOVA and PERMDISP dispersion testing, plus hand-off
  guidance to profiling-workflow (mg_rarefaction, mg_differential, mg_network,
  mg_plot) and paired 16S + shotgun studies. Use for 16S, ITS, amplicon,
  ASV, UniFrac, DADA2, QIIME 2, or ampliseq downstream analysis. Upstream
  denoising, primer trimming, and taxonomic assignment run through the
  nf-core/ampliseq wrapper, not here.
phi:
  environment: plugin:metagenomics
  attachTo:
    - Metagenomics
  scripts:
    - name: amplicon_unifrac
      description: >-
        Compute exact Normalized Weighted UniFrac and Unweighted UniFrac
        distance matrices from an ASV x sample count table and a rooted Newick
        phylogenetic tree (e.g., nf-core/ampliseq or QIIME 2 export), matching
        ASVs strictly by identifier (never by row position). Runs PCoA,
        PERMANOVA, and PERMDISP on the chosen UniFrac metric and writes both
        distance matrices, PCoA coordinates, and ordination figure.
      run: [python, ./scripts/amplicon_unifrac.py]
      args:
        type: object
        additionalProperties: false
        required: [asv_table, tree, output_dir]
        properties:
          asv_table:
            type: string
            format: input-path
            description: ASV x sample count matrix TSV/CSV (first column = ASV ID).
          tree:
            type: string
            format: input-path
            description: Rooted phylogenetic tree in Newick format matching ASV IDs.
          metric:
            type: string
            enum: [weighted, unweighted]
            description: UniFrac metric used for PCoA, PERMANOVA, and PERMDISP (default weighted).
          metadata:
            type: string
            format: input-path
            description: Optional sample metadata CSV/TSV with a group column.
          group_col:
            type: string
            description: Metadata column defining groups (default group).
          output_dir:
            type: string
            format: project-path
            description: Directory for UniFrac matrices, PCoA coordinates, and figure.
          permutations:
            type: integer
            minimum: 99
            description: PERMANOVA/PERMDISP permutations (default 999).
          seed:
            type: integer
            description: RNG seed for permutations (default 12345).
      approval: write
      output: ./schemas/amplicon-unifrac-result.json
      timeoutSeconds: 1800
---

# 16S / ITS Amplicon — downstream chain

This skill bridges upstream amplicon processing (`nf-core/ampliseq` or QIIME 2
/ DADA2 exports) with phylogenetic beta diversity (**Weighted & Unweighted
UniFrac**) and the shared downstream tools in `profiling-workflow`
(`mg_rarefaction`, `mg_diversity`, `mg_differential`, `mg_network`,
`mg_classify`, `mg_plot`).

## Chain

```
raw 16S/ITS reads ──(nf-core/ampliseq wrapper)──▶ ASV table + taxonomy + rooted tree
ASV count table + rooted tree ──mg_amplicon_unifrac──▶ Weighted/Unweighted UniFrac + PCoA + PERMANOVA + PERMDISP
ASV / genus count matrix ──(profiling-workflow)──▶ mg_rarefaction / mg_differential / mg_network / mg_plot
```

## Red lines (amplicon-specific)

1. **Match the classifier to the sequenced hypervariable region**: verify the
   amplicon primer pair / region (e.g., V3–V4 vs V4 vs full-length 16S) before
   launching `nf-core/ampliseq` or QIIME 2 taxonomy assignment — applying a
   V4-only classifier to V3–V4 reads silently distorts genus-level abundances.
2. **Align tables strictly by ASV/taxon identifier**: never concatenate or
   merge per-sample ASV tables by row position; `mg_amplicon_unifrac` intersects
   ASV IDs with tree tips explicitly and reports `n_asvs_matched`.
3. **Rooted tree required**: UniFrac depends on root placement; pass a rooted
   phylogenetic tree (e.g., midpoint-rooted FastTree/IQ-TREE from `ampliseq` or
   QIIME 2 `rooted-tree.nwk`).
4. **Paired 16S + shotgun studies**: keep 16S genus/ASV matrices and shotgun
   species matrices in separate tables and compare ecological conclusions across
   modalities rather than merging counts across assays.
