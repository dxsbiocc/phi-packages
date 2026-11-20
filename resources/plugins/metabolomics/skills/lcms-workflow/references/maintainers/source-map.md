# Source map (maintainers)

What informed this skill, what was taken from each source, and licensing.
Nothing here is user-facing analysis guidance; that lives in the playbook
and step notes.

## Sources

### 1. RforMassSpectrometry Metabolomics Book (primary)

https://rformassspectrometry.github.io/metaRbolomics-book/ — the
RforMassSpectrometry community's metabolomics analysis book. Repository
license confirmed **CC-BY-4.0** (GitHub API license key `cc-by-4.0`).
Chapters consulted (2026-11-20):

- §2.1 Preprocessing: msconvert/ProteoWizard conversion and centroiding;
  xcms peak picking (matchedFilter, centWave, massifquant), m/z+RT density
  grouping, `featureValues`/`featureDefinitions` export shape; §2.1.7
  annotation: CAMERA correlation+peak-shape pseudospectra, CliqueMS,
  RAMClustR; Rdisop/enviPat formula assignment; masstrixR exact-mass
  lookup; the "exact-mass annotation ≠ identification" distinction.
- §2.2 MS2 matching: cosine/dot-product/X-Rank scoring (MSnbase,
  OrgMassSpecR, msPurity), MSP dialect variation, MGF (GNPS), MassBank.
- §2.5 Statistical analysis: profiling vs fingerprinting; ropls permutation
  validation and mixOmics VIP; QC-based normalization (crmn,
  MetNormalizer, NormalizeMets); left-censored imputation (GSimp);
  MetaboAnalystR's built-in mummichog.
- §2.7 Functional analysis: FELLA KEGG diffusion, mummichog (Li 2013),
  graphite/KEGGREST/RaMP mapping access.

- Used for: chain ordering (QC → preprocess → multivariate/univariate →
  annotation → MS2 → pathway), method choices and their caveats, MSI-level
  discipline, the gap list (OPLS-DA, CAMERA full, SIRIUS/MetFrag/GNPS,
  GSimp, drift correction).
- **Ideas and structure only; all prose here is original. No text or code
  copied.** CC-BY-4.0 attribution: RforMassSpectrometry community,
  *Metabolomics Book* (https://rformassspectrometry.github.io/metaRbolomics-book/).

### 2. Zhihu tutorials (requested, inaccessible)

The user asked to reference:
- https://zhuanlan.zhihu.com/p/1929478762258494316
- https://zhuanlan.zhihu.com/p/111735059

Both returned HTTP 403 on fetch (2026-11-20) and could not be located via
search. With the user's approval, the skill proceeded on the book +
standard practice alone. Nothing is attributed to these URLs; revisit if
they become accessible.

### 3. Classical method literature (idea-level citations)

- Smith et al. 2006, *BMC Bioinformatics* — xcms LC/MS preprocessing.
- Tautenhahn et al. 2008, *BMC Bioinformatics* — centWave peak detection.
- Kuhl et al. 2012, *Analytical Chemistry* — CAMERA adduct/isotope
  annotation (our mb_annotate is a simplified m/z+RT variant, original code).
- Schymanski et al. 2014, *Environ Sci Technol* — identification confidence
  levels (MSI levels 1–5 framing used throughout).
- Li et al. 2013, *PLoS Comput Biol* — mummichog m/z-level permutation
  enrichment (mb_pathway masses mode is an original implementation of the
  idea, not a port).
- Dieterle et al. 2006, *Analytical Chemistry* — PQN normalization.
- van den Berg et al. 2006, *BMC Genomics* — centering/scaling choices
  (uv/pareto).
- Westerhuis et al. 2008, *Metabolomics* — PLS-DA validation/permutation
  discipline; Chong & Jun 2005 — VIP definition.
- Dunn et al. 2011, *Nature Protocols*; Broadhurst et al. 2018,
  *Metabolomics* — QC sample design and RSD filtering conventions.
- Horai et al. 2010, *J Mass Spectrom* — MassBank; Wang et al. 2016,
  *Nature Methods* — GNPS/MGF ecosystem.
- Xia & Wishart 2016, *Current Protocols in Bioinformatics* —
  MetaboAnalyst workflow conventions.

No code copied from any of these; all scripts are original numpy/scipy/
pandas implementations (NIPALS PLS2, BH step-up, hypergeometric test,
greedy-fragment cosine, union-find pseudospectra).

## Environment

`plugin:metabolomics` mirrors the metagenomics plugin stack exactly
(Python 3.13 + pandas/numpy/scipy/matplotlib-base, conda-forge); locks were
solved with micromamba 2.1.1 dry-run against the tuna mirror and rewritten
to conda.anaconda.org. Only fontconfig differs from the metagenomics locks
(upstream 2.18.3 → 2.18.4 rebuild between solve dates).