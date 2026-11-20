# mb_pathway — enrichment without a built-in database

## Principle

Two strategies, both requiring user-supplied, licensed reference data:

- **ids**: classical over-representation. Input = putatively identified
  compound IDs; mapping = pathway→compound TSV. One-sided hypergeometric
  (Fisher) test per pathway: P(X ≥ k) with N = all mapped compounds,
  K = pathway size, n = input size. Background is the mapping's compound
  universe — it should reflect what the platform could detect.
- **masses**: mummichog-style (Li 2013) for mostly unidentified features.
  Each significant m/z maps to a *set* of candidate compounds via adduct
  hypotheses (pos: +H/+Na/+K/+NH4; neg: −H/+Cl/+FA−H) against a compound
  mass database, within ppm. Pathway overlap counts the union of candidate
  sets. Significance comes from permutation: draw random m/z lists of the
  same size from the **full detected background**, rematch, recount;
  empirical p = (1 + #null ≥ observed) / (1 + n_perm), then BH FDR.

## Why permutation over m/z, not compounds

m/z matching is ambiguous (isomers, adducts, isobars). Permuting at the
m/z level preserves that ambiguity structure in the null — permuting
compound IDs would pretend the matches were certain. This is the core
mummichog insight; the implementation here is original.

## Interpretation

- Read `n_matched_input`: many unmatched m/z means the database does not
  cover the chemistry (wrong adduct mode, missing lipids, wrong organism).
- Enrichment is a hypothesis generator. Keep MSI caveats attached: the
  input features are level 2/3 at best.
- Overlap counts matter as much as p: a pathway hit on 1 compound is weak
  regardless of significance.
- `masses` p-values are bounded below by 1/(1+n_perm); raise n_perm for
  finer resolution near the threshold.

## Pitfalls

- Background must be all detected features (masses) or the full mapping
  universe (ids) — never just the significant list.
- Mapping/database licensing: KEGG requires a license for bulk use; HMDB
  is CC-BY-NC; MetaboAnalyst libraries have their own terms. Confirm and
  cite; nothing is bundled.
- FELLA network diffusion and topology-aware methods are gaps.
- Same seed → identical permutations; cite it.