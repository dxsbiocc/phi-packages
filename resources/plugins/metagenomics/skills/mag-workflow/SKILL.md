---
name: mag-workflow
description: >-
  Downstream analysis of metagenome assembly and MAGs: summarize QUAST
  assembly statistics (N50/L50/total length), grade bin quality from
  CheckM/CheckM2 with MIMAG tiers, filter bins by completeness/contamination
  into a clean set, and join quality + GTDB-Tk taxonomy + CoverM abundance +
  dRep clusters into one MAG table with a phylum figure. Use for MAG,
  metagenome-assembled genome, 分箱, binning, assembly, 组装, N50, QUAST,
  CheckM, CheckM2, GTDB-Tk, dRep, CoverM, MIMAG, or when bin FASTAs and
  quality tables need grading, filtering and reporting. Running assemblers
  and binners (megahit/spades/metabat2/maxbin2/concoct/dastool/checkm2/
  gtdbtk) or the whole nf-core/mag pipeline is done with wrappers, not here.
  The read-based counterpart is profiling-workflow.
phi:
  environment: plugin:metagenomics
  attachTo:
    - Metagenomics
  scripts:
    - name: assembly_stats
      description: >-
        Parse a QUAST report.tsv into per-assembly statistics: contig count,
        total length, largest contig, N50/N75, L50/L75, GC%. Uses the >=0 bp
        (all-contig) rows when present. Call before deciding whether an
        assembly is worth binning.
      run: [python, ./scripts/assembly_stats.py]
      args:
        type: object
        additionalProperties: false
        required: [input]
        properties:
          input:
            type: string
            format: input-path
            description: QUAST report.tsv file.
      approval: read
      output: ./schemas/assembly-stats-result.json
      timeoutSeconds: 120
    - name: bin_qc
      description: >-
        Grade bins from a CheckM2 quality_report.tsv or CheckM qa --tab_table
        output into MIMAG-style tiers (HQ completeness>=90/contamination<=5,
        MQ >=50/<=10, else LQ; thresholds adjustable). Writes a per-bin
        quality table with tier column and a completeness-vs-contamination
        scatter. HQ here is completeness/contamination only — rRNA/tRNA
        evidence is not in these tables.
      run: [python, ./scripts/bin_qc.py]
      args:
        type: object
        additionalProperties: false
        required: [input, output_dir]
        properties:
          input:
            type: string
            format: input-path
            description: CheckM2 quality_report.tsv or CheckM qa tab_table TSV.
          output_dir:
            type: string
            format: project-path
            description: Directory for the quality table and scatter figure.
          hq_completeness:
            type: number
            description: HQ completeness threshold (default 90).
          hq_contamination:
            type: number
            description: HQ contamination threshold (default 5).
          mq_completeness:
            type: number
            description: MQ completeness threshold (default 50).
          mq_contamination:
            type: number
            description: MQ contamination threshold (default 10).
      approval: write
      output: ./schemas/bin-qc-result.json
      timeoutSeconds: 600
    - name: bin_filter
      description: >-
        Copy bins passing completeness/contamination thresholds from a bins
        directory into a clean output set with a manifest (default >=50%
        complete, <=10% contamination — the MIMAG medium-quality floor).
        Fails loudly if a quality-table bin has no FASTA. Feed the output
        directory to dRep/GTDB-Tk/CoverM wrappers or annotation.
      run: [python, ./scripts/bin_filter.py]
      args:
        type: object
        additionalProperties: false
        required: [bins_dir, quality, output_dir]
        properties:
          bins_dir:
            type: string
            format: input-path
            description: Directory of bin FASTAs (.fa/.fasta/.fna).
          quality:
            type: string
            format: input-path
            description: Quality table (mg_bin_qc output, CheckM2 or CheckM qa).
          output_dir:
            type: string
            format: project-path
            description: Directory receiving passing bins and manifest.tsv.
          min_completeness:
            type: number
            description: Minimum completeness percent (default 50).
          max_contamination:
            type: number
            description: Maximum contamination percent (default 10).
      approval: write
      output: ./schemas/bin-filter-result.json
      timeoutSeconds: 1800
    - name: mag_table
      description: >-
        Left-join a bin quality table with optional GTDB-Tk summary
        (user_genome/classification, split into d__/p__/.../s__ columns),
        CoverM genome abundance table, and dRep Cdb.csv clusters into one
        mag-table.tsv; draws a MAGs-per-phylum barplot when taxonomy is
        present. Missing optional inputs are reported in the joined list,
        never silently filled.
      run: [python, ./scripts/mag_table.py]
      args:
        type: object
        additionalProperties: false
        required: [quality, output_dir]
        properties:
          quality:
            type: string
            format: input-path
            description: Quality table (mg_bin_qc output, CheckM2 or CheckM qa).
          gtdb_summary:
            type: string
            format: input-path
            description: GTDB-Tk gtdbtk.*.summary.tsv.
          coverm:
            type: string
            format: input-path
            description: CoverM genome output TSV (Genome column + metrics).
          drep_clusters:
            type: string
            format: input-path
            description: dRep Cdb.csv (genome + secondary_cluster).
          output_dir:
            type: string
            format: project-path
            description: Directory for mag-table.tsv and the phylum figure.
      approval: write
      output: ./schemas/mag-table-result.json
      timeoutSeconds: 600
---

# Assembly / MAG route — downstream chain

This skill owns the **deterministic downstream** of the assembly route:
QUAST stats → bin quality grading → quality filtering → joined MAG table.
The heavy steps run through wrappers: either the whole **nf-core/mag**
pipeline (assembly + binning + DAS Tool + CheckM/BUSCO + GTDB-Tk in one run)
or individual module wrappers (megahit/spades → metabat2/maxbin2/concoct →
dastool → checkm2 → gtdbtk). **dRep and CoverM are not part of mag 5.5.0** —
run them as separate wrappers after mag, on the filtered bin set.

## Before you start

1. Read [references/analysis-playbook.md](references/analysis-playbook.md) —
   whole-pipeline vs stepwise choice, database sizes, MIMAG red lines.
2. Locate inputs with `mg_inspect` (profiling-workflow): it recognizes QUAST
   reports and bin quality tables too.
3. Per-tool notes: [references/step-notes/](references/step-notes/).

## Chain

```
QUAST report ──mg_assembly_stats──▶ decide: bin this assembly?
bin quality  ──mg_bin_qc──────────▶ tiers + scatter ──▶ choose thresholds
bin FASTAs + quality ──mg_bin_filter──▶ clean set + manifest
   clean set ──(dRep / GTDB-Tk / CoverM wrappers)──▶ summaries
quality + gtdb + coverm + drep ──mg_mag_table──▶ mag-table.tsv + phylum plot
```

## Rules

- MIMAG tiers from completeness/contamination alone are **not** full MIMAG
  HQ: the strict standard also needs 5S/16S/23S rRNA and ≥18 tRNAs. Say so
  in every report; the JSON caveat carries this.
- Thresholds are decisions: MQ floor (50/10) is a default, HQ (90/5) is the
  standard; EasyMetagenome-style dRep dereplication uses ANI 0.95 with
  `-comp 50 -con 10` — pass them explicitly to the dRep wrapper.
- GTDB-Tk needs its ~80+ GB reference database and CheckM2 its own database;
  confirm disk and download consent with the user before any wrapper run
  that fetches them.
- Every JSON return is the audit record — cite tier counts and pass/reject
  numbers from it, do not recompute.
- Gaps to report honestly: gene annotation (prokka/prodigal wrappers), AMR
  screening (abricate/rgi/hamronization wrappers), viral identification
  (genomad/checkv wrappers) are separate steps, not implemented here.

## Sources

Distilled from EasyMetagenome (iMeta 2025), cambiotraining/metagenomics
(CC BY 4.0), sib-swiss metagenomics-data-analysis-training (GPL-3.0; ideas
only) and the vendored nf-core/mag 5.5.0 wrapper schema. Provenance:
[references/maintainers/source-map.md](references/maintainers/source-map.md).