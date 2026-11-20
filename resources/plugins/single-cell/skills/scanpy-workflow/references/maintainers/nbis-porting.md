# Porting the NBIS scanpy labs

Decisions behind the scanpy fragment library and script chain, in the same
spirit as the Seurat skill's scnotebooks-porting.md.

## Source and license

The step order and explanations follow the NBIS workshop-scRNAseq scanpy labs
(https://nbisweden.github.io/workshop-scRNAseq/, NBIS Sweden). The repository
LICENSE is CC-BY-4.0; the site footer's GPL-3 label contradicts it, so cite
the repo license. This repository is MIT and **no lab code is copied**: the
labs were downloaded as a step reference only, and every fragment and script
here is original code that re-implements the same analysis with the plugin's
conventions (config cell, JSON audit, no install cells). Attribution travels
with the library (`steps.json` `attribution`, notebook footer) as CC-BY
requires for the adapted explanatory material.

The former standalone `resources/skills/scanpy` skill (K-Dense, BSD-3-Clause)
was absorbed into this plugin: its workflow knowledge informed the script
chain, and its BSD-3 notice is compatible with MIT. Its host-install guidance
was deliberately dropped — plugin environments are private and read-only.

## Environment facts that shaped the port

- scanpy 1.13 requires python >= 3.13; the environment pins python 3.13 and
  scanpy 1.13 (current stable line).
- `bbknn`, `harmonypy`, `scrublet` live on **bioconda**, not conda-forge.
  bbknn has no python 3.13 build (py39–py312 only), so integration ships
  harmony + combat and bbknn is a reported gap, not a silent failure.
- `scikit-misc` ships because `seurat_v3` HVG selection needs its loess.
- Lock generation quirk: micromamba 2.9.0 reverse-aliases mirror URLs through
  `mirrored_channels` and then sorts channels alphabetically, which puts
  bioconda above conda-forge and excludes conda-forge's scanpy under strict
  priority. The locks were solved with mirror URLs as plain channels (order =
  priority, no `mirrored_channels`) and rewritten back to the official host.

## Analysis conventions (shared with the Seurat route)

- counts layer preserved at load; QC and pseudobulk start from it.
- HVG on counts, `lognorm` layer snapshot before scaling, PCA on the variable-gene mask
  instead of gene subsetting — the Python expression of "keep all genes,
  scale once, test on un-scaled values".
- scrublet runs pre-normalization, matching the labs.
- Panel scoring = panel mean − rest mean; `min_score` gate; `unassigned`
  fallback; identical gene panels across routes.
- Pseudobulk layout (`gene_id` first, `sample` sheet) matches the repository's
  DE wrappers, so both routes hand off identically.
- Thresholds are arguments; every filter reports removed cells; kept_fraction
  below half is a warning.

## What was deliberately not ported

- Colab install cells and Google Drive paths — the environment is prebuilt.
- The labs' exploratory dead ends (re-running PCA after subsetting, manual
  doublet gating by eye) — replaced by the deterministic chain.
- scVI / CellTypist / Scanorama / gseapy / spatial sections — packages or
  runtime model downloads the environment omits; `sc_py_env_probe` reports
  them as gaps with a project-environment next step.