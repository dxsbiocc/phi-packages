# Source map (maintainers) — mag-workflow

Shared provenance for both metagenomics skills lives in
[../../../profiling-workflow/references/maintainers/source-map.md](../../../profiling-workflow/references/maintainers/source-map.md)
(EasyMetagenome, Snakemake Workflow repo, sib-swiss GPL-3.0 ideas-only,
cambiotraining CC BY 4.0 attribution, tool citations). This file records
what is specific to the assembly/MAG route.

## Route-specific extraction

- **EasyMetagenome v1.24**: assembly–binning stage ordering
  (megahit/SPAdes → MetaBAT2/MaxBin2/CONCOCT → DAS Tool → CheckM → dRep →
  GTDB-Tk), dRep dereplication convention `-sa 0.95 -nc 0.30 -comp 50
  -con 10`, and the MQ-floor filtering habit. No code copied (GPL-3.0).
- **cambiotraining/metagenomics** (CC BY 4.0, attribution here): MaxBin2
  multi-reads usage, CheckM `taxonomy_wf` → `qa --tab_table -o 2` two-step
  pattern, quality-discussion framing, prokka/abricate as follow-up
  annotation/AMR steps. No verbatim text copied.
- **sib-swiss training** (GPL-3.0): MAG quality interpretation caveats
  (CheckM vs CheckM2 non-interchangeability, MIMAG rRNA/tRNA requirement,
  co-assembly grouping rules). Ideas only.
- **Vendored nf-core/mag 5.5.0** (`resources/wrappers/workflows/nf-core/mag/`):
  parameter names/defaults cross-checked against `nextflow_schema.json`
  ($defs groups: assembly_options, binning_options,
  bin_quality_check_options, taxonomic_profiling_options). Confirmed mag
  5.5.0 has **no dRep/CoverM** — the playbook places them as post-mag
  wrappers.

## Standards

- MIMAG: Bowers et al. 2017, Nat Biotechnol. doi:10.1038/nbt.3893
  (HQ ≥90/≤5 + rRNA/tRNA; MQ ≥50/≤10).
- ANI species boundary ~95%: standard practice reflected in dRep defaults
  and GTDB species clustering.

## Implementation notes

- Scripts are original; QUAST parsing selects >=0 bp rows to avoid cutoff
  ambiguity. CheckM2/CheckM detection is column-signature based.
- Bin identity convention: file stem without .fa/.fasta/.fna/.fsa —
  matches CheckM2 Name, CheckM Bin Id, GTDB user_genome, dRep genome.
- dRep cluster column preference: secondary_cluster (species-level at
  -sa 0.95), fallback cluster.