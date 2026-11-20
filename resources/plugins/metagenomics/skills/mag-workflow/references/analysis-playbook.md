# Assembly / MAG analysis playbook

The decision loop for the assembly route. Every step: **call → read the
JSON → decide → next step**, and stop at the red lines.

## 0. Route check

- Reads already classified (Kraken2/MetaPhlAn reports)? That is the
  read-based route — use profiling-workflow instead.
- 16S/ITS amplicons? Neither route; use the ampliseq wrapper workflow.
- Goal is strain tracking of one known genome? bowtie2 alignment to that
  reference (wrapper) beats de novo assembly.
- Goal is community genomes, novel taxa, or genome-resolved function?
  Continue here.

## 1. Upstream: whole pipeline vs stepwise

**Option A — nf-core/mag (wrapper `workflows/nf-core/mag`).** One run does
QC → assembly (MEGAHIT and/or SPAdes) → binning (MetaBAT2/MaxBin2/CONCOCT/
SemiBin/COMEBin/MetaBinner) → DAS Tool refinement → BUSCO/CheckM/CheckM2 →
GTDB-Tk → gene prediction. Prefer it for multi-sample studies: it handles
co-assembly groups, binning across samples (`binning_map_mode`), and
provenance. Inspect the wrapper first (`wrapper_inspect`); key decisions:

- `input`: samplesheet (sample,fastq_1,fastq_2 or long-read columns).
- Assemblers: `skip_spades`/`skip_megahit` — MEGAHIT is faster/leaner,
  metaSPAdes often yields longer contigs; running both doubles cost but
  gives DAS Tool more to work with. `coassemble_group` pools samples of the
  same group into one assembly (needed for low-depth samples; never pool
  unrelated environments — gut with skin).
- Host removal: `host_genome`/`host_fasta` for host-associated samples
  (human, plant); otherwise host reads become junk bins.
- Binning: `min_contig_size` (default 1500) — smaller contigs carry noisy
  coverage/tetranucleotide signal. Keep most binners on; DAS Tool
  (`refine_bins_dastool`, threshold 0.5) reconciles them.
- Quality: `run_checkm2` (fast, machine-learning) vs `run_checkm`
  (lineage-specific marker sets, the accepted gold standard); BUSCO runs
  by default with `busco_db_lineage` (auto or e.g. bacteria_odb10).
- Taxonomy: GTDB-Tk unless `skip_gtdbtk`; needs `gtdb_db` (~80+ GB —
  confirm disk and download consent first). `gtdbtk_min_completeness` 50 /
  `gtdbtk_max_contamination` 10 gate which bins get classified.

**Option B — stepwise module wrappers.** megahit or spades → quast →
metabat2/maxbin2/concoct → dastool → checkm2 or checkm → gtdbtk. Choose
when only part of the chain is missing (e.g. bins exist, quality does not),
or when mag's resource envelope does not fit. Same tools, same thresholds;
you own the orchestration and the audit trail.

Either way: assembly and binning are **wrapper runs**, never plugin
scripts. Read each wrapper's `wrapper.yaml` before planning.

## 2. Assembly quality — mg_assembly_stats

- **Call** on QUAST `report.tsv` (from quast wrapper or mag's quast output).
- **Read** N50/L50, total length, contig count per assembly.
- **Decide:** metagenome assemblies are fragmented by nature — N50 of a few
  kb is normal for complex communities; do not apply isolate-genome
  expectations. Total length far above the expected community genome sum
  suggests strain heterogeneity or contamination; far below suggests
  insufficient depth. Compare co-assemblies vs per-sample assemblies before
  choosing which to bin.
- **Red line:** do not claim "good assembly" from N50 alone; report the
  whole row and the biological context.

## 3. Bin quality — mg_bin_qc

- **Call** on CheckM2 `quality_report.tsv` or CheckM `qa --tab_table` TSV.
  CheckM classic usage is two steps: `taxonomy_wf` then `qa` (cambiotraining
  pattern); the wrapper exposes both.
- **Read** tier counts and the scatter. Bins cluster near the axes' origin
  when binning failed; a completeness–contamination trade-off tail suggests
  over-merged bins.
- **Decide thresholds with the user's goal:** genome-resolved ecology
  usually wants MQ+ (≥50/≤10); reference-grade claims want HQ (≥90/≤5).
  CheckM and CheckM2 values differ systematically (different marker models)
  — do not mix tables from both in one comparison.
- **Red line:** MIMAG HQ strictly requires 5S/16S/23S rRNA + ≥18 tRNAs;
  these tables cannot show that. Report "HQ by completeness/contamination"
  and, if full MIMAG is claimed, run rRNA/tRNA checks separately (barrnap/
  tRNAscan-SE — not in this environment; say so).

## 4. Filter — mg_bin_filter

- **Call** with bins_dir + quality + thresholds from step 3.
- **Read** n_passed / n_rejected. A pass rate near zero means the assembly
  or binning failed upstream — go back, do not lower thresholds to manufacture
  MAGs.
- **Decide:** the clean set feeds dRep (dereplicate), GTDB-Tk (if not already
  run), CoverM (abundance), annotation.
- **Red line:** thresholds must be stated in the report; filtering at 30/20
  to inflate counts is a decision to surface, not hide.

## 5. Post-mag wrappers (not part of mag 5.5.0)

- **dRep** on the clean set: species-level dereplication at ANI 0.95
  (EasyMetagenome convention `-sa 0.95 -nc 0.30 -comp 50 -con 10`); strain
  level 0.99. Output `Cdb.csv` feeds mg_mag_table.
- **CoverM genome** mode: map reads back to the (dereplicated) MAG set for
  mean coverage / covered fraction — relative abundance proxy. Covered
  fraction < ~0.3–0.5 of a genome means the abundance estimate is unreliable.
- **GTDB-Tk** if skipped in mag: classify_wf on the clean set; same database
  warning as step 1.

## 6. Report — mg_mag_table

- **Call** with quality + whatever of gtdb/coverm/drep exists.
- **Read** `joined` and `n_classified`: unclassified MAGs are novel-lineage
  candidates, not errors — but GTDB-Tk only classifies bins passing its
  completeness/contamination gate, so check the gate before calling
  something novel.
- **Deliver:** mag-table.tsv (one row per MAG: quality tier, taxonomy ranks,
  abundance, cluster), phylum barplot, and the JSON audit numbers.
- **Red lines:** do not equate dRep clusters with species without stating
  the ANI threshold; do not present CoverM mean coverage as absolute
  abundance without a spike-in or genome-length normalization discussion.

## Overall red lines

- No assembler/binner/CheckM/GTDB-Tk execution inside the plugin
  environment — wrappers only.
- No database downloads without disk check and user consent (GTDB ~80+ GB,
  CheckM2 several GB, BUSCO lineages hundreds of MB).
- No silent threshold changes; every filter/classification gate appears in
  the report.
- Annotation (prokka/prodigal), AMR (abricate/rgi/hamronization), viruses
  (genomad/checkv) are follow-up wrapper steps — report them as next steps,
  do not improvise them here.

Sources: [maintainers/source-map.md](maintainers/source-map.md).