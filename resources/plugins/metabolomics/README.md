# metabolomics (metabolomic analysis) · 0.2.0

**Downstream analysis chain + knowledge layer** for untargeted LC-MS feature
tables: vendor export inspection (xcms/MZmine/MS-DIAL/MetaboAnalyst) → QC
sample RSD/detection-rate filtering → imputation/normalization (TIC/median/
PQN/internal standard)/log transform → PCA and PLS-DA (VIP + seeded
permutation Q2) → univariate tests (Welch/Wilcoxon/Kruskal-Wallis + BH FDR +
volcano plots) → deterministic m/z annotation (isotopes/adducts/neutral
losses → pseudospectra) → MS2 library matching (MSP/MGF cosine) → pathway
enrichment (ID hypergeometric / mummichog-style m/z permutation) → clustered
heatmaps.
Upstream (msconvert conversion, xcms/MZmine/MS-DIAL peak picking and
alignment) is **not in this plugin**, and the repository has no corresponding
wrapper yet — the knowledge layer states the gap plainly and asks the user for
the exported table.
Plugin scripts do only deterministic downstream statistics and run in one
private locked Python environment.

Untargeted LC-MS feature-table downstream chain plus a knowledge layer:
vendor export inspection → QC RSD/detection filtering → imputation/
normalization (TIC/median/PQN/IS)/log transform → PCA & PLS-DA (VIP +
seeded permutation Q2) → univariate tests (Welch/Wilcoxon/KW + BH FDR,
volcano) → deterministic m/z annotation (isotopes/adducts/neutral losses →
pseudospectra) → MS2 library matching (MSP/MGF cosine) → pathway enrichment
(ID hypergeometric / mummichog-style m/z permutation) → clustered heatmaps.
Upstream peak picking (msconvert, xcms, MZmine, MS-DIAL) is out of scope and
has no wrapper here; the knowledge layer reports the gap honestly.
All scripts are deterministic downstream statistics in one private locked
Python environment.

## Components

| Component | Role |
|--|--|
| Agent `Metabolomics` | Route decision (feature-table downstream vs upstream gap), QC/threshold decisions, how to phrase MSI identification levels, database license confirmation |
| Skill `lcms-workflow` | 10 `mb_*` tools (inspect/qc/drift_correct/preprocess/multivariate/univariate/annotate/ms2_match/pathway/plot) |
| Environment `metabolomics` | Python 3.13 + pandas/numpy/scipy/matplotlib-base, three-platform locks (same stack as metagenomics) |

The spec and the three-platform locks are in `environments/metabolomics/`.
v0.1.0 has only the script chain, with no notebook renderer.

## Boundaries (the agent's responsibility)

- An existing feature table (+metadata/QC annotation) → this plugin's
  downstream chain.
- Only raw `.raw`/`.mzML` → upstream gap: ask the user for an
  xcms/MZmine/MS-DIAL export table; do not pretend to do peak picking.
- GC-MS feature tables work mechanically, but retention indices (RI) are not
  supported; NMR, targeted MRM and lipid species naming are out of scope.
- OPLS-DA, GSimp imputation, full CAMERA pseudospectra, molecular formula
  elucidation, SIRIUS/MetFrag/GNPS →
  the knowledge layer reports the gap instead of faking it.
- Pathway/compound/spectral libraries are **not built in**: KEGG/HMDB
  mappings, mass databases and MSP/MGF libraries come from the user; confirm
  the license and version before use.

## Knowledge layer

`skills/lcms-workflow/references/`: analysis-playbook (call → read the JSON →
decide → red lines), 9 step-notes, maintainers/source-map. The main source is
the RforMassSpectrometry *Metabolomics Book* (CC-BY-4.0; ideas and structure
borrowed, all prose original, no code copied); the two Zhihu tutorials
supplied by the user return 403 and are noted as inaccessible in the
source-map.

## Verification

Repository-level `tests/test_metabolomics_plugin.py` (structure + function +
counter-examples + determinism, synthetic fixtures with fixed seeds);
`node --test tests/*.test.mjs` covers the skill contracts inside the plugin.