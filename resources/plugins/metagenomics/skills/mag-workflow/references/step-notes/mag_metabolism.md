# `mg_mag_metabolism` — Metabolic pathway completeness and GEM summary

## Purpose

Evaluates metabolic potential across recovered MAGs from three complementary
upstream sources without adding new Python dependencies:
1. `--ko_matrix` (`mag-ko-matrix.tsv` from `mg_mag_function`) scored against
   10 curated core microbial KEGG modules (glycolysis, pentose phosphate, TCA,
   beta-oxidation, nitrogen fixation, denitrification, DNRA, sulfate reduction,
   Wood-Ljungdahl, cobalamin B12) or a custom `--module_map` TSV.
2. `--gapseq_pathways` (`*-Pathways.tbl` from `nf-core/gapseq`).
3. `--sbml` (genome-scale metabolic models `.xml` / `.xml.gz` from
   `nf-core/carveme` or `nf-core/gapseq`), parsed with Python stdlib XML to
   summarize metabolites, reactions, exchange reactions, and gene products.

## Interpretation and pitfalls

- Built-in KO module completeness measures fractional KO presence per module;
  for custom pathways (e.g. alkane degradation or SCFA synthesis), pass
  `--module_map` with `module_id`, `module_name`, and `ko_list`.
- SBML structural summaries do not execute Flux Balance Analysis (FBA) or
  verify growth-medium feasibility.
