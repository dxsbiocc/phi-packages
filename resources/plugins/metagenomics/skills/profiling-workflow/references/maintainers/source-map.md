# Source map (maintainers)

What informed this skill, what was taken from each source, and licensing.
Nothing here is user-facing analysis guidance; that lives in the playbook
and step notes.

## Sources

### 1. EasyMetagenome / EasyMicrobiome (local: ~/Downloads/Meta)

EasyMetagenome v1.24 (Liu Yong-Xin lab) — a four-stage shotgun pipeline:
QC (fastp q15/len50, optional 20M-read downsampling), read-based profiling
(kraken2+bracken / MetaPhlAn+HUMAnN), assembly–binning–MAG (megahit/SPAdes →
MetaBAT2/MaxBin2/CONCOCT → DAS Tool → CheckM → dRep `-sa 0.95 -nc 0.30
-comp 50 -con 10` → GTDB-Tk), and statistics/visualization (EasyStat,
EasyEVenn). Chinese course slides included practical parameter choices and
figure conventions.

- Used for: pipeline stage ordering, QC parameter defaults, dRep threshold
  conventions, the read-based vs assembly-based route split.
- Citation: Liu et al. EasyMetagenome: a user-friendly and flexible pipeline
  for shotgun metagenomics. *iMeta* 2025. doi:10.1002/imt2.70001.
- Related tools referenced by the course: CAMI2 benchmark
  (doi:10.1038/s41592-022-01431-4), Sangwan et al. recovering genomes from
  metagenomes (doi:10.1186/s40168-016-0154-5), EasyAmplicon
  (doi:10.1002/imt2.83), EasyAmplicon2 (doi:10.1002/advs.202512447), EVenn
  (doi:10.1002/imt2.184), dbCAN2 (doi:10.1093/nar/gkz935).
- No code copied. EasyMetagenome is GPL-3.0; we borrow workflow knowledge
  only.

### 2. Snakemake Workflow repo (local: ~/Documents/WorkSpace/VSCode/Workflow)

Snakemake NGS framework with Meta rules (kraken2, bracken, metaphlan,
humann, megahit, binning). Used for: confirming tool pairing
(kraken2→bracken, metaphlan→humann) and per-tool input/output file shapes.
No code copied; scripts here are original implementations.

### 3. sib-swiss/metagenomics-data-analysis-training

https://github.com/sib-swiss/metagenomics-data-analysis-training (sparse
checkout of Day1_am, Day1_pm, Day2, slides). Course covers read-based
profiling (Kraken2/Bracken, MetaPhlAn/HUMAnN), assembly + MAG recovery,
CheckM/GTDB-Tk quality and taxonomy, and community statistics.

- Used for: decision-loop structure, interpretation caveats (depth
  normalization, compositionality, PERMANOVA dispersion).
- License GPL-3.0: **ideas and workflow structure only; no text or code
  copied.**

### 4. cambiotraining/metagenomics

https://cambiotraining.github.io/metagenomics/ (repo cloned to
/tmp/metag-training/cambio during development). Cambridge Bioinformatics
Training metagenomics course: QC (FastQC/cutadapt/Trimmomatic), bowtie2
host removal, SPAdes/metaSPAdes assembly, MaxBin2 binning, CheckM/GTDB-Tk,
prokka annotation, abricate AMR screening.

- Used for: assembly-side practical details (flash merge, clumpify dedup
  rationale, `bowtie2 --un` host-read extraction, CheckM taxonomy_wf → qa
  two-step usage, MIMAG-style quality discussion).
- License CC BY 4.0: attribution required — this file and the playbook's
  source line provide it. No verbatim text copied.

### 5. Vendored nf-core/mag 5.5.0 wrapper

`resources/wrappers/workflows/nf-core/mag/` — used for the assembly/MAG
route's whole-pipeline option (mag-workflow skill) and for parameter
defaults cross-checked against its `nextflow_schema.json`.

## Statistical methods

Original implementations in `scripts/`: Shannon/Simpson/Chao1 (standard
ecological indices), Bray-Curtis (Bray & Curtis 1957), principal coordinates
analysis (Gower 1966), PERMANOVA (Anderson 2001,
doi:10.1111/j.1442-9993.2001.01070.pp.x), Mann-Whitney U / Kruskal-Wallis
(scipy), Benjamini-Hochberg FDR (Benjamini & Hochberg 1995).

## Tool citations (for report writing)

- Kraken2: Wood et al. 2019, Genome Biol. doi:10.1186/s13059-019-1891-0
- Bracken: Lu et al. 2017, PeerJ Comput Sci. doi:10.7717/peerj-cs.104
- MetaPhlAn 4: Blanco-Míguez et al. 2023, Nat Biotechnol.
  doi:10.1038/s41587-023-01688-w
- HUMAnN 3: Beghini et al. 2021, eLife. doi:10.7554/eLife.65088
- MEGAHIT: Li et al. 2015, Bioinformatics. doi:10.1093/bioinformatics/btv033
- metaSPAdes: Nurk et al. 2017, Genome Res. doi:10.1101/gr.213959.116
- MetaBAT2: Kang et al. 2019, PeerJ. doi:10.7717/peerj.7359
- MaxBin2: Wu et al. 2016, Bioinformatics. doi:10.1093/bioinformatics/btv638
- CONCOCT: Alneberg et al. 2014, Nat Methods. doi:10.1038/nmeth.3103
- DAS Tool: Sieber et al. 2018, Nat Microbiol. doi:10.1038/s41564-018-0171-1
- CheckM: Parks et al. 2015, Genome Res. doi:10.1101/gr.186072.114
- CheckM2: Chklovski et al. 2023, Nat Methods. doi:10.1038/s41592-023-01940-w
- dRep: Olm et al. 2017, ISME J. doi:10.1038/ismej.2017.126
- GTDB-Tk: Chaumeil et al. 2022, Bioinformatics. doi:10.1093/bioinformatics/btac672
- MIMAG standards: Bowers et al. 2017, Nat Biotechnol.
  doi:10.1038/nbt.3893