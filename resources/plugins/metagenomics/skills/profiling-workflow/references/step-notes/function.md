# Functional tables

**Purpose.** Turn HUMAnN-style feature tables into comparable, correctly
levelled functional profiles.

**Principle.** HUMAnN genefamilies are stratified: `UniRef90_X|g__Y.s__Z`
rows attribute function to species; unstratified rows are community totals.
Pathabundance rows are MetaCyc pathways. Regrouping sums features into
higher-level groups (KO, EC, CAZy family) via an explicit mapping file —
the mapping is data you supply, not something the tool infers.

**Interpretation.** After `normalize`, columns sum to 1 (relative) or 1e6
(CPM) — cross-sample comparison is now depth-independent. `UNMAPPED` /
`UNINTEGRATED` rows are real fractions of the community (reads without
annotation); a large UNMAPPED fraction limits functional conclusions.

**Pitfalls.** This normalize is plain column renormalization; HUMAnN's
renorm handles stratification differently — use the humann3 wrapper for
exact semantics. Stratified rows double-count if summed with unstratified
rows; pick one view per analysis. Regroup without a mapping fails loudly —
obtain the mapping (e.g. UniRef90→KO from the HUMAnN utility databases)
first.