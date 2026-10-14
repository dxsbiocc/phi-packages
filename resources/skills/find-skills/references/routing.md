# Routing hints

Confirm every row with a directory scan. Descriptions in `SKILL.md` win when they disagree with this table.

| Task looks like | Skill | Location |
| --- | --- | --- |
| Which skill exists or which one fits | `find-skills` | `resources/skills/find-skills` |
| Author or scaffold a skill | `create-skills` | `resources/skills/create-skills` |
| Papers, preprints, citations, trials | `literature-search` | `resources/skills/literature-search` |
| Plan inputs, commands, environment, outputs | `reproducible-analysis` | `resources/skills/reproducible-analysis` |
| Methods, results, or figure captions from outputs | `scientific-writing` | `resources/skills/scientific-writing` |
| `.h5ad` / AnnData objects and file layout | `anndata` | `resources/skills/anndata` |
| scRNA-seq QC, clustering, markers, DE | `scanpy` | `resources/skills/scanpy` |
| RNA velocity, latent time | `scvelo` | `resources/skills/scvelo` |
| SAM/BAM/CRAM, VCF/BCF, FASTA/FASTQ | `pysam` | `resources/skills/pysam` |
| Nextflow, nf-core, DSL2, nf-test | `nextflow` | `resources/skills/nextflow` |
| SMILES, fingerprints, substructure search | `rdkit` | `resources/skills/rdkit` |
| Graphs, centrality, communities | `networkx` | `resources/skills/networkx` |
| Classification, regression, clustering, pipelines | `scikit-learn` | `resources/skills/scikit-learn` |
| Censored survival, Cox, concordance | `scikit-survival` | `resources/skills/scikit-survival` |
| SHAP values and explanation plots | `shap` | `resources/skills/shap` |
| Statistical plots on data frames | `seaborn` | `resources/skills/seaborn` |
| Fine-grained Matplotlib figures | `matplotlib` | `resources/skills/matplotlib` |
| PDF, Office, HTML, or audio to Markdown | `markitdown` | `resources/skills/markitdown` |
| Publication omics figures from templates | `omics-visualization` | `resources/plugins/visualization/skills/omics-visualization` |
| Author a Nextflow wrapper | `create-wrapper` | Phi application only; not in this repository |

When a new standalone skill is added under `resources/skills/`, add one row here. The scan procedure in `where-to-look.md` stays the source of truth.
