#!/usr/bin/env python3
"""Manually wire wrappers for the vendored nf-core subworkflows that the
structure-inferring generator (scripts/generate-wrappers.py) could not handle.

Each wrapper is a thin agent-facing adapter: wrapper/main.nf includes ONLY the
vendored subworkflow at ../main.nf (which itself references the already-migrated
module directories under resources/wrappers/modules via relative includes). No
process logic is re-implemented here. Optional channels become Channel.empty();
static reference/val inputs become Channel.value(...) tuples; defaults live in
params.json only. publishDir is global (matches the generated subworkflow
wrappers) so every nested module's outputs land in one results directory.

The six utils_* pipeline-infrastructure subworkflows are intentionally skipped:
they require a live workflow object / pipeline parameter context (nf-schema /
nf-validation plugins, version printing, igenomes reference resolution,
annotation-cache plumbing) and are not standalone analysis components.

Idempotent: skips any subworkflow that already has wrapper/wrapper.yaml.
"""
import json
import os

ROOT = os.path.join(os.path.dirname(__file__), "..", "resources", "wrappers", "subworkflows", "nf-core")
ROOT = os.path.abspath(ROOT)

SKIP = {
    "utils_nextflow_pipeline": "pipeline infrastructure: version/parameter dumping, needs pipeline context",
    "utils_nfcore_pipeline": "pipeline infrastructure: config validation, needs nextflow_cli_args",
    "utils_nfschema_plugin": "pipeline infrastructure: nf-schema plugin, requires a live workflow object",
    "utils_nfvalidation_plugin": "pipeline infrastructure: nf-validation plugin, requires a live workflow object",
    "utils_references": "pipeline infrastructure: igenomes/reference YAML resolution",
    "utils_annotation_cache": "pipeline infrastructure: annotation-cache plumbing tied to pipeline params/help text",
}

# param tuple: (name, kind, type, required, description, default, enum)
def P(name, kind, typ, req, desc, default=None, enum=None):
    return (name, kind, typ, req, desc, default, enum)

SPECS = []
def spec(rel, workflow, name, summary, params, body, flags=None, cpus=2, memory="4.GB", time="2.h"):
    SPECS.append(dict(rel=rel, workflow=workflow, name=name, summary=summary,
                      params=params, body=body, flags=flags or [], cpus=cpus, memory=memory, time=time))

READS = "    reads_ch = Channel\n        .fromFilePairs(params.reads, size: -1)\n        .map { id, reads -> [[id: id, single_end: reads.size() == 1], reads] }\n"

# ---------------------------------------------------------------- simple ones
spec("archive_extract", "ARCHIVE_EXTRACT", "Extract archive",
     "Extract a tar/zip archive and route extracted vs non-extracted files.",
     [P("archive", "input", "path", True, "Archive file to extract (tar.gz, zip, ...).")],
     "    archive_ch = Channel\n        .fromPath(params.archive, checkIfExists: true)\n        .map { f -> [[id: f.simpleName], f] }\n\n    ARCHIVE_EXTRACT(archive_ch)\n")

spec("faa_seqfu_seqkit", "FAA_SEQFU_SEQKIT", "Protein FASTA QC (seqfu/seqkit)",
     "Check and compute statistics on a protein FASTA with seqfu and seqkit.",
     [P("fasta", "input", "fasta", True, "Protein FASTA file."),
      P("skip_preprocessing", "option", "boolean", False, "Skip the preprocessing step.", False)],
     "    fasta_ch = Channel\n        .fromPath(params.fasta, checkIfExists: true)\n        .map { f -> [[id: f.simpleName], f] }\n\n    FAA_SEQFU_SEQKIT(fasta_ch, params.skip_preprocessing)\n")

spec("fasta_clean_faidx", "FASTA_CLEAN_FAIDX", "Clean + index FASTA",
     "Trim/renamed a reference FASTA and build fai, chromsizes, dict and sequence stats.",
     [P("fasta", "input", "fasta", True, "Reference FASTA to clean and index."),
      P("replace_dots", "option", "boolean", False, "Replace dots in headers with underscores.", False),
      P("get_chromsizes", "option", "boolean", False, "Emit a chromsizes file.", True),
      P("get_dict", "option", "boolean", False, "Emit a sequence dictionary.", True)],
     "    reference_ch = Channel\n        .fromPath(params.fasta, checkIfExists: true)\n        .map { f -> [[id: f.simpleName], f] }\n\n    FASTA_CLEAN_FAIDX(reference_ch, params.replace_dots, params.get_chromsizes, params.get_dict)\n")

spec("fasta_consensus_autocycler", "FASTA_CONSENSUS_AUTOCYCLER", "Consensus assembly (AutoCycler)",
     "Build a consensus assembly from grouped contigs with AutoCycler.",
     [P("contigs", "input", "path_glob", True, "Glob of contig FASTA files grouped into one consensus (e.g. 'contigs/*.fasta').")],
     "    contigs_ch = Channel\n        .fromPath(params.contigs, checkIfExists: true)\n        .collect()\n        .map { files -> [[id: 'consensus'], files] }\n\n    FASTA_CONSENSUS_AUTOCYCLER(contigs_ch)\n",
     flags=["all matched contigs are treated as one group; multi-group input needs curation"])

spec("fasta_index_dna", "FASTA_INDEX_DNA", "Build DNA aligner index",
     "Build a genome aligner index (bowtie2, bwamem, bwamem2, dragmap or snap).",
     [P("fasta", "input", "fasta", True, "Reference genome FASTA."),
      P("altliftover", "input", "path", False, "ALT liftover file (only used by snap)."),
      P("aligner", "option", "string", True, "Aligner to build the index for.", "bwamem2",
        ["bowtie2", "bwamem", "bwamem2", "dragmap", "snap"])],
     "    fasta_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true)])\n"
     "    alt_ch = params.altliftover ? Channel.value([[id: 'reference'], file(params.altliftover, checkIfExists: true)]) : Channel.empty()\n\n"
     "    FASTA_INDEX_DNA(fasta_ch, alt_ch, params.aligner)\n",
     cpus=4, memory="8.GB", time="4.h")

spec("fasta_index_bismark_bwameth", "FASTA_INDEX_BISMARK_BWAMETH", "Build methylation aligner index",
     "Build bismark and/or bwameth indexes for a reference FASTA.",
     [P("fasta", "input", "fasta", True, "Reference genome FASTA."),
      P("fai", "input", "path", False, "Pre-built FASTA index (.fai)."),
      P("bismark_index", "input", "path", False, "Pre-built bismark index directory (skips building)."),
      P("bwameth_index", "input", "path", False, "Pre-built bwameth index directory (skips building)."),
      P("aligner", "option", "string", True, "Aligner to index for.", "bismark", ["bismark", "bismark_hisat", "bwameth"]),
      P("use_mem2", "option", "boolean", False, "Use bwamem2 for the bwameth index.", False)],
     "    fasta_fai_ch = Channel.value([[id: 'reference'], [file(params.fasta, checkIfExists: true)], params.fai ? [file(params.fai, checkIfExists: true)] : []])\n"
     "    bismark_ch = params.bismark_index ? Channel.value([[id: 'reference'], [file(params.bismark_index, checkIfExists: true)]]) : Channel.empty()\n"
     "    bwameth_ch = params.bwameth_index ? Channel.value([[id: 'reference'], [file(params.bwameth_index, checkIfExists: true)]]) : Channel.empty()\n\n"
     "    FASTA_INDEX_BISMARK_BWAMETH(fasta_fai_ch, bismark_ch, bwameth_ch, params.aligner, params.use_mem2)\n",
     flags=["index directories are staged whole; verify companion semantics"], cpus=4, memory="8.GB", time="4.h")

spec("fasta_index_methylseq", "FASTA_INDEX_METHYLSEQ", "Build methylseq indexes",
     "Build the full set of methylseq aligner indexes (bismark/bwameth/bwamem) plus fai.",
     [P("fasta", "input", "fasta", True, "Reference genome FASTA."),
      P("fasta_index", "input", "path", False, "Pre-built FASTA index (.fai)."),
      P("bismark_index", "input", "path", False, "Pre-built bismark index directory."),
      P("bwameth_index", "input", "path", False, "Pre-built bwameth index directory."),
      P("bwamem_index", "input", "path", False, "Pre-built bwamem index directory."),
      P("aligner", "option", "string", True, "Aligner to index for.", "bismark", ["bismark", "bismark_hisat", "bwameth", "bwamem"]),
      P("collecthsmetrics", "option", "boolean", False, "Run picard CollectHsMetrics.", False),
      P("use_mem2", "option", "boolean", False, "Use bwamem2 for the bwameth index.", False)],
     "    fasta_ch = Channel.value([[id: 'reference'], [file(params.fasta, checkIfExists: true)]])\n"
     "    fai_ch = params.fasta_index ? Channel.value([[id: 'reference'], [file(params.fasta_index, checkIfExists: true)]]) : Channel.empty()\n"
     "    bismark_ch = params.bismark_index ? Channel.value([[id: 'reference'], [file(params.bismark_index, checkIfExists: true)]]) : Channel.empty()\n"
     "    bwameth_ch = params.bwameth_index ? Channel.value([[id: 'reference'], [file(params.bwameth_index, checkIfExists: true)]]) : Channel.empty()\n"
     "    bwamem_ch = params.bwamem_index ? Channel.value([[id: 'reference'], [file(params.bwamem_index, checkIfExists: true)]]) : Channel.empty()\n\n"
     "    FASTA_INDEX_METHYLSEQ(fasta_ch, fai_ch, bismark_ch, bwameth_ch, bwamem_ch, params.aligner, params.collecthsmetrics, params.use_mem2)\n",
     flags=["index directories are staged whole; verify companion semantics"], cpus=4, memory="8.GB", time="4.h")

spec("fasta_vclust_prefilter_align_cluster", "FASTA_VCLUST_PREFILTER_ALIGN_CLUSTER", "Cluster FASTA (vclust)",
     "Prefilter, align and cluster genome FASTAs with vclust at a chosen ANI/TAI threshold.",
     [P("fasta", "input", "path_glob", True, "Glob of genome FASTA files to cluster."),
      P("save_alignment", "option", "boolean", False, "Save the alignment output.", False),
      P("metric", "option", "string", False, "Clustering metric.", "ani", ["ani", "gani", "tani"]),
      P("tani", "option", "number", False, "Total ANI threshold.", 0.95),
      P("gani", "option", "number", False, "Global ANI threshold.", 0.95),
      P("ani", "option", "number", False, "ANI threshold.", 0.95)],
     "    fasta_ch = Channel\n        .fromPath(params.fasta, checkIfExists: true)\n        .collect()\n        .map { files -> [[id: 'vclust'], files] }\n\n"
     "    FASTA_VCLUST_PREFILTER_ALIGN_CLUSTER(fasta_ch, params.save_alignment, params.metric, params.tani, params.gani, params.ani)\n",
     flags=["all matched FASTAs treated as one group"], cpus=4, memory="8.GB", time="4.h")

spec("fastq_complexity_filter", "FASTQ_COMPLEXITY_FILTER", "Complexity filter FASTQ",
     "Filter low-complexity reads with prinseqplusplus, bbduk or fastp.",
     [P("reads", "input", "path_glob", True, "FASTQ glob (paired or single)."),
      P("complexity_filter_tool", "option", "string", True, "Tool to use.", "fastp", ["prinseqplusplus", "bbduk", "fastp"])],
     READS + "\n    FASTQ_COMPLEXITY_FILTER(reads_ch, params.complexity_filter_tool)\n")

spec("fastq_contam_seqtk_kraken", "FASTQ_CONTAM_SEQTK_KRAKEN", "Screen contamination (seqtk/kraken2)",
     "Subsample reads with seqtk and screen for contamination with kraken2.",
     [P("reads", "input", "path_glob", True, "FASTQ glob (paired or single)."),
      P("sample_size", "option", "string", True, "Number of reads to subsample.", "10000"),
      P("kraken2_db", "input", "path", True, "Path to the kraken2 database directory.")],
     READS + "    kraken2_db_ch = Channel.value(file(params.kraken2_db, checkIfExists: true))\n\n"
     "    FASTQ_CONTAM_SEQTK_KRAKEN(reads_ch, params.sample_size, kraken2_db_ch)\n",
     flags=["kraken2_db comment says string but module expects a staged DB; passed as path"])

spec("fastq_extract_kraken_krakentools", "FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS", "Extract taxon reads (kraken2)",
     "Classify reads with kraken2 and extract reads of given taxids with krakentools.",
     [P("reads", "input", "path_glob", True, "FASTQ glob (paired or single)."),
      P("db", "input", "path", True, "Path to the kraken2 database directory."),
      P("taxid", "option", "string", True, "Taxonomic ids to extract, separated by spaces.", "9606")],
     READS + "    db_ch = Channel.value(file(params.db, checkIfExists: true))\n\n"
     "    FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS(reads_ch, db_ch, params.taxid)\n")

spec("fastq_preprocess_seqkit", "FASTQ_PREPROCESS_SEQKIT", "Preprocess FASTQ (seqkit)",
     "Sanity-check, clean, replace and deduplicate FASTQ reads with seqkit.",
     [P("reads", "input", "path_glob", True, "FASTQ glob (paired or single)."),
      P("skip_seqkit_sana_pair", "option", "boolean", False, "Skip seqkit sana pair.", False),
      P("skip_seqkit_seq", "option", "boolean", False, "Skip seqkit seq.", False),
      P("skip_seqkit_replace", "option", "boolean", False, "Skip seqkit replace.", False),
      P("skip_seqkit_rmdup", "option", "boolean", False, "Skip seqkit rmdup.", False)],
     READS + "\n    FASTQ_PREPROCESS_SEQKIT(reads_ch, params.skip_seqkit_sana_pair, params.skip_seqkit_seq, params.skip_seqkit_replace, params.skip_seqkit_rmdup)\n")

spec("fastq_qc_stats", "FASTQ_QC_STATS", "FASTQ QC statistics",
     "Compute FastQC, seqfu, seqkit and seqtk statistics on FASTQ reads.",
     [P("reads", "input", "path_glob", True, "FASTQ glob (paired or single)."),
      P("skip_fastqc", "option", "boolean", False, "Skip FastQC.", False),
      P("skip_seqfu_check", "option", "boolean", False, "Skip seqfu check.", False),
      P("skip_seqfu_stats", "option", "boolean", False, "Skip seqfu stats.", False),
      P("skip_seqkit_stats", "option", "boolean", False, "Skip seqkit stats.", False),
      P("skip_seqtk_comp", "option", "boolean", False, "Skip seqtk comp.", False)],
     READS + "\n    FASTQ_QC_STATS(reads_ch, params.skip_fastqc, params.skip_seqfu_check, params.skip_seqfu_stats, params.skip_seqkit_stats, params.skip_seqtk_comp)\n")

spec("fastq_trim_fastp_fastqc", "FASTQ_TRIM_FASTP_FASTQC", "Trim FASTQ (fastp) + FastQC",
     "Trim adapters/quality with fastp and run FastQC on the trimmed reads.",
     [P("reads", "input", "path_glob", True, "FASTQ glob (paired or single)."),
      P("adapter_fasta", "input", "path", False, "Adapter FASTA for fastp."),
      P("save_trimmed_fail", "option", "boolean", False, "Save reads that fail trimming.", False),
      P("discard_trimmed_pass", "option", "boolean", False, "Discard reads that pass trimming.", False),
      P("save_merged", "option", "boolean", False, "Save merged reads (paired input).", False),
      P("skip_fastp", "option", "boolean", False, "Skip fastp trimming.", False),
      P("skip_fastqc", "option", "boolean", False, "Skip FastQC.", False)],
     "    reads_ch = Channel\n        .fromFilePairs(params.reads, size: -1)\n        .map { id, reads -> [[id: id, single_end: reads.size() == 1], reads, params.adapter_fasta ? file(params.adapter_fasta, checkIfExists: true) : []] }\n\n"
     "    FASTQ_TRIM_FASTP_FASTQC(reads_ch, params.save_trimmed_fail, params.discard_trimmed_pass, params.save_merged, params.skip_fastp, params.skip_fastqc)\n")

spec("mmseqs_fasta_cluster", "MMSEQS_FASTA_CLUSTER", "Cluster sequences (MMseqs2)",
     "Cluster protein/nucleotide sequences with MMseqs2 linclust or cluster.",
     [P("fasta", "input", "fasta", True, "Sequence FASTA to cluster."),
      P("clustering_tool", "option", "string", True, "MMseqs2 clustering workflow.", "linclust", ["linclust", "cluster"])],
     "    seq_ch = Channel\n        .fromPath(params.fasta, checkIfExists: true)\n        .map { f -> [[id: f.simpleName], f] }\n\n    MMSEQS_FASTA_CLUSTER(seq_ch, params.clustering_tool)\n",
     cpus=4, memory="8.GB", time="4.h")

spec("vcf_annotate_snpeff", "VCF_ANNOTATE_SNPEFF", "Annotate VCF (snpEff)",
     "Annotate a VCF with snpEff using a named database or local cache.",
     [P("vcf", "input", "path", True, "VCF file to annotate."),
      P("snpeff_db", "option", "string", True, "snpEff database version (e.g. GRCh38.99).", "GRCh38.99"),
      P("snpeff_cache", "input", "path", False, "Local snpEff cache directory (skips download).")],
     "    vcf_ch = Channel\n        .fromPath(params.vcf, checkIfExists: true)\n        .map { f -> [[id: f.simpleName], f] }\n"
     "    cache_ch = params.snpeff_cache ? Channel.value(file(params.snpeff_cache, checkIfExists: true)) : Channel.empty()\n\n"
     "    VCF_ANNOTATE_SNPEFF(vcf_ch, params.snpeff_db, cache_ch)\n")

spec("vcf_filter_bcftools_ensemblvep", "VCF_FILTER_BCFTOOLS_ENSEMBLVEP", "Filter VCF (bcftools/VEP)",
     "Filter a VCF with bcftools view and/or Ensembl VEP filter_vep.",
     [P("vcf", "input", "path", True, "VCF file to filter."),
      P("filter_vep_feature_file", "input", "path", False, "filter_vep feature/filter text file."),
      P("filter_with_bcftools", "option", "boolean", False, "Run bcftools view filtering.", True),
      P("filter_with_filter_vep", "option", "boolean", False, "Run filter_vep filtering.", False)],
     "    vcf_ch = Channel\n        .fromPath(params.vcf, checkIfExists: true)\n        .map { f -> [[id: f.simpleName], f] }\n"
     "    feature_ch = params.filter_vep_feature_file ? Channel.fromPath(params.filter_vep_feature_file, checkIfExists: true).map { f -> [[id: f.simpleName], f] } : Channel.empty()\n\n"
     "    VCF_FILTER_BCFTOOLS_ENSEMBLVEP(vcf_ch, feature_ch, params.filter_with_bcftools, params.filter_with_filter_vep)\n")

spec("vcf_gather_bcftools", "VCF_GATHER_BCFTOOLS", "Gather scattered VCFs (bcftools)",
     "Concatenate/gather scattered VCF shards per sample and index the result.",
     [P("vcfs", "input", "path_glob", True, "Glob of VCF shards to gather."),
      P("vcf_indexes", "input", "path_glob", False, "Glob of matching VCF index files (.tbi/.csi)."),
      P("common_meta", "option", "string", True, "Comma-separated meta fields used for grouping.", "id"),
      P("sort", "option", "boolean", False, "Sort the gathered output (needs a sort config).", False)],
     "    vcfs_ch = Channel\n        .fromPath(params.vcfs, checkIfExists: true)\n        .map { f -> [[id: f.simpleName], f, [], 1] }\n\n"
     "    VCF_GATHER_BCFTOOLS(vcfs_ch, params.common_meta.tokenize(','), params.sort)\n",
     flags=["single-group adapter; index/count companions and multi-sample grouping need curation"])

spec("homer_groseq", "HOMER_GROSEQ", "GRO-seq analysis (HOMER)",
     "Make tag directories, UCSC files, find peaks and convert to BED with HOMER.",
     [P("bam", "input", "bam", True, "Aligned BAM file."),
      P("fasta", "input", "fasta", True, "Reference genome FASTA (for HOMER tag directory)."),
      P("uniqmap", "input", "path", False, "Uniquely-mappable regions file or zip (optional).")],
     "    bam_ch = Channel\n        .fromPath(params.bam, checkIfExists: true)\n        .map { f -> [[id: f.simpleName], [f]] }\n"
     "    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))\n\n"
     "    HOMER_GROSEQ(bam_ch, fasta_ch, params.uniqmap ?: null)\n",
     flags=["uniqmap passed as a path string; verify zip handling"])

spec("gtf_hybridmerge_gffcompare", "GTF_HYBRIDMERGE_GFFCOMPARE", "Hybrid-merge GTF (gffcompare)",
     "Merge novel transcripts into a reference GTF keeping chosen gffcompare class codes.",
     [P("novel_gtf", "input", "path", True, "Novel transcript GTF."),
      P("reference_gtf", "input", "path", True, "Reference annotation GTF."),
      P("backbone_gtf", "input", "path", True, "Backbone GTF to merge into."),
      P("class_codes", "option", "string", True, "gffcompare class codes to keep (e.g. 'u' or 'u,j,i').", "u,j,i"),
      P("blacklist_bed", "input", "path", False, "Blacklist BED to intersect (skips if absent).")],
     "    novel_ch = Channel.fromPath(params.novel_gtf, checkIfExists: true).map { f -> [[id: 'hybrid'], f] }\n"
     "    ref_ch = Channel.value([[id: 'hybrid'], file(params.reference_gtf, checkIfExists: true)])\n"
     "    backbone_ch = Channel.value([[id: 'hybrid'], file(params.backbone_gtf, checkIfExists: true)])\n"
     "    blacklist_ch = params.blacklist_bed ? Channel.value([[id: 'hybrid'], file(params.blacklist_bed, checkIfExists: true)]) : Channel.empty()\n\n"
     "    GTF_HYBRIDMERGE_GFFCOMPARE(novel_ch, ref_ch, backbone_ch, params.class_codes, blacklist_ch)\n")

spec("fasta_hmmsearch_rank_fastas", "FASTA_HMMSEARCH_RANK_FASTAS", "Rank FASTAs by HMM hits",
     "Search HMM profiles against a FASTA, rank hits and extract matching sequences.",
     [P("hmms", "input", "path_glob", True, "Glob of HMM profile files."),
      P("fasta", "input", "fasta", True, "Single query FASTA file."),
      P("save_domtblout", "option", "boolean", False, "Also write the per-domain hit table.", False)],
     "    hmms_ch = Channel\n        .fromPath(params.hmms, checkIfExists: true)\n        .map { f -> [[id: f.simpleName], f] }\n"
     "    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))\n\n"
     "    FASTA_HMMSEARCH_RANK_FASTAS(hmms_ch, fasta_ch, params.save_domtblout)\n")

spec("fastq_align_dna", "FASTQ_ALIGN_DNA", "Align DNA reads",
     "Align FASTQ reads to a genome index with bowtie2/bwamem/bwamem2/dragmap/snap.",
     [P("reads", "input", "path_glob", True, "FASTQ glob (paired or single)."),
      P("index", "input", "path", True, "Aligner index directory."),
      P("fasta", "input", "fasta", True, "Reference genome FASTA."),
      P("aligner", "option", "string", True, "Aligner to use.", "bwamem2", ["bowtie2", "bwamem", "bwamem2", "dragmap", "snap"]),
      P("sort", "option", "boolean", False, "Sort the output BAM.", True)],
     READS + "    index_ch = Channel.value([[id: 'reference'], [file(params.index, checkIfExists: true)]])\n"
     "    fasta_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true)])\n\n"
     "    FASTQ_ALIGN_DNA(reads_ch, index_ch, fasta_ch, params.aligner, params.sort)\n",
     flags=["index staged as a single directory; aligner-specific index layout may need curation"], cpus=4, memory="8.GB", time="4.h")

# ------------------------------------------------------- align + dedup family
spec("fastq_align_bamcmp_bwa", "FASTQ_ALIGN_BAMCMP_BWA", "Align + compare contamination (bwa/bamcmp)",
     "Align reads to primary and contaminant indexes and filter with bamcmp.",
     [P("reads", "input", "path_glob", True, "FASTQ glob (paired or single)."),
      P("primary_index", "input", "path", True, "Primary bwa index directory."),
      P("contaminant_index", "input", "path", True, "Contaminant bwa index directory."),
      P("fasta", "input", "fasta", False, "Reference FASTA (optional)."),
      P("fai", "input", "path", False, "Reference FASTA index (optional).")],
     READS + "    primary_ch = Channel.value([[id: 'primary'], [file(params.primary_index, checkIfExists: true)]])\n"
     "    contaminant_ch = Channel.value([[id: 'contaminant'], [file(params.contaminant_index, checkIfExists: true)]])\n"
     "    fasta_fai_ch = params.fasta ? Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true)]) : Channel.empty()\n\n"
     "    FASTQ_ALIGN_BAMCMP_BWA(reads_ch, primary_ch, contaminant_ch, fasta_fai_ch)\n",
     cpus=4, memory="8.GB", time="4.h")

spec("fastq_align_chromap", "FASTQ_ALIGN_CHROMAP", "Align chromatin reads (chromap)",
     "Align FASTQ reads with chromap (chromatin accessibility / Hi-C).",
     [P("reads", "input", "path_glob", True, "FASTQ glob (paired or single)."),
      P("index", "input", "path", True, "chromap index directory."),
      P("fasta", "input", "fasta", True, "Reference genome FASTA."),
      P("fai", "input", "path", True, "Reference FASTA index (.fai)."),
      P("barcodes", "input", "path", False, "Barcode file (optional)."),
      P("whitelist", "input", "path", False, "Barcode whitelist (optional)."),
      P("chr_order", "input", "path", False, "Chromosome order file (optional)."),
      P("pairs_chr_order", "input", "path", False, "Pairs chromosome order file (optional)."),
      P("update_readgroups", "option", "boolean", False, "Update read groups post alignment.", False)],
     READS + "    index_ch = Channel.value([[id: 'reference'], [file(params.index, checkIfExists: true)]])\n"
     "    fasta_fai_ch = Channel.value([[id: 'reference'], [file(params.fasta, checkIfExists: true)], [file(params.fai, checkIfExists: true)]])\n"
     "    barcodes_ch = params.barcodes ? Channel.value(file(params.barcodes, checkIfExists: true)) : Channel.empty()\n"
     "    whitelist_ch = params.whitelist ? Channel.value(file(params.whitelist, checkIfExists: true)) : Channel.empty()\n"
     "    chr_order_ch = params.chr_order ? Channel.value(file(params.chr_order, checkIfExists: true)) : Channel.empty()\n"
     "    pairs_ch = params.pairs_chr_order ? Channel.value(file(params.pairs_chr_order, checkIfExists: true)) : Channel.empty()\n\n"
     "    FASTQ_ALIGN_CHROMAP(reads_ch, index_ch, fasta_fai_ch, barcodes_ch, whitelist_ch, chr_order_ch, pairs_ch, params.update_readgroups)\n",
     cpus=4, memory="8.GB", time="4.h")

spec("fastq_align_dedup_bismark", "FASTQ_ALIGN_DEDUP_BISMARK", "Align + dedup (bismark)",
     "Align bisulfite reads with bismark, deduplicate and optionally report cytosines.",
     [P("reads", "input", "path_glob", True, "FASTQ glob (paired or single)."),
      P("fasta", "input", "fasta", True, "Reference genome FASTA."),
      P("fai", "input", "path", True, "Reference FASTA index (.fai)."),
      P("bismark_index", "input", "path", True, "bismark index directory."),
      P("skip_deduplication", "option", "boolean", False, "Skip deduplication.", False),
      P("cytosine_report", "option", "boolean", False, "Run coverage2cytosine.", False)],
     READS + "    fasta_fai_ch = Channel.value([[id: 'reference'], [file(params.fasta, checkIfExists: true)], [file(params.fai, checkIfExists: true)]])\n"
     "    index_ch = Channel.value([[id: 'reference'], [file(params.bismark_index, checkIfExists: true)]])\n\n"
     "    FASTQ_ALIGN_DEDUP_BISMARK(reads_ch, fasta_fai_ch, index_ch, params.skip_deduplication, params.cytosine_report)\n",
     cpus=4, memory="8.GB", time="4.h")

spec("fastq_align_dedup_bwamem", "FASTQ_ALIGN_DEDUP_BWAMEM", "Align + dedup (bwa-mem)",
     "Align reads with bwa-mem (or parabricks GPU), sort, dedup and compute stats.",
     [P("reads", "input", "path_glob", True, "FASTQ glob (paired or single)."),
      P("fasta", "input", "fasta", True, "Reference genome FASTA."),
      P("fai", "input", "path", True, "Reference FASTA index (.fai)."),
      P("bwamem_index", "input", "path", True, "bwa-mem index directory."),
      P("skip_deduplication", "option", "boolean", False, "Skip deduplication.", False),
      P("use_gpu", "option", "boolean", False, "Use GPU parabricks fq2bam.", False),
      P("output_fmt", "option", "string", False, "parabricks output format.", "bam", ["bam", "cram"]),
      P("interval_file", "input", "path", False, "Interval file (optional)."),
      P("known_sites", "input", "path", False, "Known sites VCF (optional).")],
     READS + "    fasta_fai_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true)])\n"
     "    index_ch = Channel.value([[id: 'reference'], [file(params.bwamem_index, checkIfExists: true)]])\n"
     "    interval_ch = params.interval_file ? Channel.value([[id: 'reference'], [file(params.interval_file, checkIfExists: true)]]) : Channel.empty()\n"
     "    known_ch = params.known_sites ? Channel.value([[id: 'reference'], [file(params.known_sites, checkIfExists: true)]]) : Channel.empty()\n\n"
     "    FASTQ_ALIGN_DEDUP_BWAMEM(reads_ch, fasta_fai_ch, index_ch, params.skip_deduplication, params.use_gpu, params.output_fmt, interval_ch, known_ch)\n",
     flags=["GPU path requires parabricks container/license"], cpus=4, memory="8.GB", time="4.h")

spec("fastq_align_dedup_bwameth", "FASTQ_ALIGN_DEDUP_BWAMETH", "Align + dedup (bwameth)",
     "Align bisulfite reads with bwameth, sort, dedup and compute stats.",
     [P("reads", "input", "path_glob", True, "FASTQ glob (paired or single)."),
      P("fasta", "input", "fasta", True, "Reference genome FASTA."),
      P("fai", "input", "path", True, "Reference FASTA index (.fai)."),
      P("bwameth_index", "input", "path", True, "bwameth index directory."),
      P("skip_deduplication", "option", "boolean", False, "Skip deduplication.", False),
      P("use_gpu", "option", "boolean", False, "Use GPU alignment.", False)],
     READS + "    fasta_fai_ch = Channel.value([[id: 'reference'], [file(params.fasta, checkIfExists: true)], [file(params.fai, checkIfExists: true)]])\n"
     "    index_ch = Channel.value([[id: 'reference'], [file(params.bwameth_index, checkIfExists: true)]])\n\n"
     "    FASTQ_ALIGN_DEDUP_BWAMETH(reads_ch, fasta_fai_ch, index_ch, params.skip_deduplication, params.use_gpu)\n",
     cpus=4, memory="8.GB", time="4.h")

spec("fastq_decontaminate_deacon_hostile", "FASTQ_DECONTAMINATE_DEACON_HOSTILE", "Decontaminate FASTQ",
     "Remove host/contaminant reads with hostile or deacon.",
     [P("reads", "input", "path_glob", True, "FASTQ glob (paired or single)."),
      P("fasta", "input", "fasta", False, "Contaminant reference FASTA (hostile, optional)."),
      P("reference", "input", "path", False, "Pre-built reference directory (optional)."),
      P("reference_name", "option", "string", False, "Name for the reference.", "contaminant"),
      P("index_name", "option", "string", False, "deacon index name (optional)."),
      P("decontaminator", "option", "string", True, "Decontamination tool.", "hostile", ["hostile", "deacon"])],
     READS + "    fasta_ch = params.fasta ? Channel.value([[id: params.reference_name], [file(params.fasta, checkIfExists: true)]]) : Channel.empty()\n"
     "    reference_ch = params.reference ? Channel.value([params.reference_name, file(params.reference, checkIfExists: true)]) : Channel.empty()\n\n"
     "    FASTQ_DECONTAMINATE_DEACON_HOSTILE(reads_ch, fasta_ch, reference_ch, params.index_name ?: null, params.decontaminator)\n")

spec("fastq_removeadapters_merge", "FASTQ_REMOVEADAPTERS_MERGE", "Trim adapters / merge reads",
     "Trim adapters and optionally merge paired reads with a choice of tools.",
     [P("reads", "input", "path_glob", True, "FASTQ glob (paired or single)."),
      P("adapter_tool", "option", "string", True, "Adapter trimming tool.", "fastp",
        ["trimmomatic", "cutadapt", "trimgalore", "bbduk", "leehom", "fastp", "adapterremoval"]),
      P("custom_adapters_file", "input", "path", False, "Custom adapters FASTA (bbduk/fastp) or txt (adapterremoval)."),
      P("save_merged", "option", "boolean", False, "Return merged reads (fastp/adapterremoval).", False),
      P("fastp_discard_trimmed_pass", "option", "boolean", False, "fastp: discard trimmed-pass reads.", False),
      P("fastp_save_trimmed_fail", "option", "boolean", False, "fastp: save trimmed-fail reads.", False)],
     READS + "    adapters_ch = params.custom_adapters_file ? Channel.value(file(params.custom_adapters_file, checkIfExists: true)) : Channel.empty()\n\n"
     "    FASTQ_REMOVEADAPTERS_MERGE(reads_ch, params.adapter_tool, adapters_ch, params.save_merged, params.fastp_discard_trimmed_pass, params.fastp_save_trimmed_fail)\n")

spec("fasta_clean_fcs", "FASTA_CLEAN_FCS", "Clean assembly (FCS-adaptor)",
     "Screen and clean a genome assembly of adaptor contamination with NCBI FCS.",
     [P("fasta", "input", "fasta", True, "Genome assembly FASTA."),
      P("database", "input", "path", True, "FCS database directory or manifest."),
      P("ramdisk_path", "option", "string", False, "Path to a ramdisk for FCS (optional).")],
     "    fasta_ch = Channel\n        .fromPath(params.fasta, checkIfExists: true)\n        .map { f -> [[id: f.simpleName], [f]] }\n"
     "    db_ch = Channel.value(file(params.database, checkIfExists: true))\n\n"
     "    FASTA_CLEAN_FCS(fasta_ch, db_ch, params.ramdisk_path ?: null)\n",
     flags=["FCS may need a ramdisk and large memory"], cpus=4, memory="16.GB", time="4.h")

spec("fasta_gxf_busco_plot", "FASTA_GXF_BUSCO_PLOT", "BUSCO assessment + plots",
     "Run BUSCO on an assembly/annotation and plot completeness summaries.",
     [P("fasta", "input", "fasta", True, "Assembly or protein FASTA."),
      P("gxf", "input", "path", False, "Annotation GFF/GFF3/GTF (evaluated when mode=genome)."),
      P("mode", "option", "string", True, "BUSCO mode.", "genome", ["genome", "transcriptome", "proteins"]),
      P("lineages", "option", "string", True, "Comma-separated BUSCO lineage datasets.", "auto"),
      P("busco_lineages_path", "input", "path", False, "Local BUSCO lineages directory (optional)."),
      P("busco_config", "input", "path", False, "BUSCO config file (optional)."),
      P("busco_cleanup", "option", "boolean", False, "Remove BUSCO intermediate files.", True)],
     "    fasta_ch = Channel.fromPath(params.fasta, checkIfExists: true).map { f -> [[id: 'busco'], f] }\n"
     "    gxf_ch = params.gxf ? Channel.value([[id: 'busco'], file(params.gxf, checkIfExists: true)]) : Channel.empty()\n\n"
     "    FASTA_GXF_BUSCO_PLOT(fasta_ch, gxf_ch, params.mode, params.lineages.tokenize(','), params.busco_lineages_path ?: [], params.busco_config ?: [], params.busco_cleanup)\n",
     flags=["lineage download needs network unless busco_lineages_path provided"], cpus=4, memory="8.GB", time="8.h")

spec("vcf_annotate_ensemblvep", "VCF_ANNOTATE_ENSEMBLVEP", "Annotate VCF (Ensembl VEP)",
     "Annotate a VCF with Ensembl VEP using a cache or GTF.",
     [P("vcf", "input", "path", True, "VCF file to annotate."),
      P("custom_files", "input", "path_glob", False, "Optional custom annotation files glob."),
      P("fasta", "input", "fasta", False, "Reference FASTA (optional)."),
      P("genome", "option", "string", True, "Genome assembly (e.g. GRCh38).", "GRCh38"),
      P("species", "option", "string", True, "Species name (e.g. homo_sapiens).", "homo_sapiens"),
      P("cache_version", "option", "string", True, "VEP cache version.", "110"),
      P("cache", "input", "path", False, "Local VEP cache directory (optional)."),
      P("extra_files", "input", "path_glob", False, "Extra files glob (optional)."),
      P("gtf", "input", "path", False, "GTF file (mutually exclusive with cache)."),
      P("gtf_tbi", "input", "path", False, "GTF tabix index (with gtf).")],
     "    vcf_ch = Channel\n        .fromPath(params.vcf, checkIfExists: true)\n        .map { f -> [[id: f.simpleName], f, params.custom_files ? file(params.custom_files, checkIfExists: true) : []] }\n"
     "    fasta_ch = params.fasta ? Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true)]) : Channel.empty()\n"
     "    cache_ch = params.cache ? Channel.value([[id: 'cache'], file(params.cache, checkIfExists: true)]) : Channel.empty()\n"
     "    extra_ch = params.extra_files ? Channel.value(file(params.extra_files, checkIfExists: true)) : Channel.empty()\n"
     "    gtf_ch = params.gtf ? Channel.value([file(params.gtf, checkIfExists: true), file(params.gtf_tbi, checkIfExists: true)]) : Channel.empty()\n\n"
     "    VCF_ANNOTATE_ENSEMBLVEP(vcf_ch, fasta_ch, params.genome, params.species, params.cache_version, cache_ch, extra_ch, gtf_ch)\n",
     flags=["cache download needs network when no local cache provided"])

spec("vcf_annotate_ensemblvep_snpeff", "VCF_ANNOTATE_ENSEMBLVEP_SNPEFF", "Annotate VCF (VEP + snpEff)",
     "Annotate a VCF with Ensembl VEP and/or snpEff, scattering by sites per chunk.",
     [P("vcf", "input", "path", True, "VCF file to annotate."),
      P("vcf_tbi", "input", "path", True, "VCF tabix index."),
      P("custom_files", "input", "path_glob", False, "Optional custom annotation files glob."),
      P("fasta", "input", "fasta", False, "Reference FASTA (optional)."),
      P("vep_genome", "option", "string", True, "VEP genome assembly.", "GRCh38"),
      P("vep_species", "option", "string", True, "VEP species.", "homo_sapiens"),
      P("vep_cache_version", "option", "string", True, "VEP cache version.", "110"),
      P("vep_cache", "input", "path", False, "Local VEP cache directory (optional)."),
      P("vep_extra_files", "input", "path_glob", False, "VEP extra files glob (optional)."),
      P("snpeff_db", "option", "string", True, "snpEff database version.", "GRCh38.99"),
      P("snpeff_cache", "input", "path", False, "Local snpEff cache directory (optional)."),
      P("tools_to_use", "option", "string", True, "Comma-separated tools: ensemblvep,snpeff.", "ensemblvep,snpeff"),
      P("sites_per_chunk", "option", "integer", False, "Variants per scattered VCF chunk.", 10000),
      P("vep_gtf", "input", "path", False, "VEP GTF (mutually exclusive with vep_cache)."),
      P("vep_gtf_tbi", "input", "path", False, "VEP GTF tabix index.")],
     "    vcf_ch = Channel\n        .fromPath(params.vcf, checkIfExists: true)\n        .map { f -> [[id: f.simpleName], f, file(params.vcf_tbi, checkIfExists: true), params.custom_files ? file(params.custom_files, checkIfExists: true) : []] }\n"
     "    fasta_ch = params.fasta ? Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true)]) : Channel.empty()\n"
     "    vep_cache_ch = params.vep_cache ? Channel.value(file(params.vep_cache, checkIfExists: true)) : Channel.empty()\n"
     "    vep_extra_ch = params.vep_extra_files ? Channel.value(file(params.vep_extra_files, checkIfExists: true)) : Channel.empty()\n"
     "    snpeff_cache_ch = params.snpeff_cache ? Channel.value(file(params.snpeff_cache, checkIfExists: true)) : Channel.empty()\n"
     "    vep_gtf_ch = params.vep_gtf ? Channel.value([file(params.vep_gtf, checkIfExists: true), file(params.vep_gtf_tbi, checkIfExists: true)]) : Channel.empty()\n\n"
     "    VCF_ANNOTATE_ENSEMBLVEP_SNPEFF(vcf_ch, fasta_ch, params.vep_genome, params.vep_species, params.vep_cache_version, vep_cache_ch, vep_extra_ch, params.snpeff_db, snpeff_cache_ch, params.tools_to_use.tokenize(','), params.sites_per_chunk, vep_gtf_ch)\n",
     flags=["cache download needs network when no local cache provided"])

spec("vcf_extract_relate_somalier", "VCF_EXTRACT_RELATE_SOMALIER", "Relatedness (somalier)",
     "Extract somalier markers from VCFs and infer sample relatedness.",
     [P("vcfs", "input", "path_glob", True, "Glob of VCF files."),
      P("vcf_indexes", "input", "path_glob", True, "Glob of matching VCF indexes."),
      P("fasta", "input", "fasta", True, "Reference genome FASTA."),
      P("fai", "input", "path", True, "Reference FASTA index (.fai)."),
      P("somalier_sites", "input", "path", True, "somalier sites VCF."),
      P("ped", "input", "path", True, "Pedigree file."),
      P("sample_groups", "input", "path", False, "Sample groups txt (optional)."),
      P("common_id", "option", "string", False, "Common identifier for single-sample VCFs.")],
     "    vcfs_ch = Channel\n        .fromPath(params.vcfs, checkIfExists: true)\n        .map { f -> [[id: f.simpleName], f, [], 1] }\n"
     "    fasta_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true)])\n"
     "    fai_ch = Channel.value([[id: 'reference'], file(params.fai, checkIfExists: true)])\n"
     "    sites_ch = Channel.value(file(params.somalier_sites, checkIfExists: true))\n"
     "    peds_ch = Channel.value([[id: 'ped'], file(params.ped, checkIfExists: true)])\n"
     "    groups_ch = params.sample_groups ? Channel.value(file(params.sample_groups, checkIfExists: true)) : Channel.empty()\n\n"
     "    VCF_EXTRACT_RELATE_SOMALIER(vcfs_ch, fasta_ch, fai_ch, sites_ch, peds_ch, groups_ch, params.common_id ?: null)\n",
     flags=["VCF/index pairing and count val need curation for multi-sample input"])

spec("bam_variant_calling_mpileup_bcftools", "BAM_VARIANT_CALLING_MPILEUP_BCFTOOLS", "Call variants (mpileup/bcftools)",
     "Call variants from BAMs at given positions with bcftools mpileup, gather and optionally annotate.",
     [P("bam", "input", "bam", True, "Aligned BAM file."),
      P("bai", "input", "path", True, "BAM index."),
      P("posfile", "input", "path", True, "Positions file (CHR POS REF ALT)."),
      P("fasta", "input", "fasta", True, "Reference genome FASTA."),
      P("fai", "input", "path", True, "Reference FASTA index."),
      P("panel_id", "option", "string", False, "Panel id for the posfile meta.", "panel"),
      P("chr", "option", "string", False, "Chromosome for the posfile meta.", "chr22"),
      P("sample_merge_key", "option", "string", False, "Meta key for sample merging.", "id"),
      P("sample_merge_value", "option", "string", False, "Meta value for sample merging.", "all_samples"),
      P("region_gather_keys", "option", "string", False, "Comma-separated meta keys for region gather.", "panel_id,id"),
      P("sort_region_gather", "option", "boolean", False, "Sort region gather output.", False),
      P("annotate", "option", "boolean", False, "Annotate the gathered VCF.", False)],
     "    bam_ch = Channel.value([[id: params.sample_merge_value], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])\n"
     "    posfile_ch = Channel.value([[panel_id: params.panel_id, chr: params.chr], file(params.posfile, checkIfExists: true)])\n"
     "    fasta_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true)])\n\n"
     "    BAM_VARIANT_CALLING_MPILEUP_BCFTOOLS(bam_ch, posfile_ch, fasta_ch, [params.sample_merge_key], [params.sample_merge_value], params.region_gather_keys.tokenize(','), params.sort_region_gather, params.annotate)\n",
     flags=["single-sample/region adapter; scatter-gather meta keys need curation"])

spec("bam_variant_demix_boot_freyja", "BAM_VARIANT_DEMIX_BOOT_FREYJA", "Lineage demix (freyja)",
     "Demix SARS-CoV-2 lineages from a BAM with freyja (optionally bootstrapped).",
     [P("bam", "input", "bam", True, "Aligned BAM file."),
      P("fasta", "input", "fasta", True, "Reference genome FASTA."),
      P("skip_boot", "option", "boolean", False, "Skip bootstrapping.", True),
      P("repeats", "option", "integer", False, "Number of bootstrap repeats.", 100),
      P("db_name", "option", "string", False, "freyja database name.", "freyja"),
      P("barcodes", "input", "path", False, "freyja barcodes file (optional; else updated)."),
      P("lineages_meta", "input", "path", False, "Lineages metadata file (optional)."),
      P("lineages_topology", "input", "path", False, "Lineages topology file (optional).")],
     "    bam_ch = Channel\n        .fromPath(params.bam, checkIfExists: true)\n        .map { f -> [[id: f.simpleName], f] }\n"
     "    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))\n"
     "    barcodes_ch = params.barcodes ? Channel.value(file(params.barcodes, checkIfExists: true)) : Channel.empty()\n"
     "    meta_ch = params.lineages_meta ? Channel.value(file(params.lineages_meta, checkIfExists: true)) : Channel.empty()\n"
     "    topo_ch = params.lineages_topology ? Channel.value(file(params.lineages_topology, checkIfExists: true)) : Channel.empty()\n\n"
     "    BAM_VARIANT_DEMIX_BOOT_FREYJA(bam_ch, fasta_ch, params.skip_boot, params.repeats, params.db_name, barcodes_ch, meta_ch, topo_ch)\n",
     flags=["freyja update needs network when barcodes not provided"])

spec("fastq_create_umi_consensus_fgbio", "FASTQ_CREATE_UMI_CONSENSUS_FGBIO", "UMI consensus (fgbio)",
     "Build UMI consensus reads from FASTQ with fgbio (bwa-mem/bwa-mem2 alignment).",
     [P("reads", "input", "path_glob", True, "FASTQ glob (paired or single)."),
      P("fasta", "input", "fasta", True, "Reference genome FASTA."),
      P("fai", "input", "path", True, "Reference FASTA index (.fai)."),
      P("dict", "input", "path", True, "Reference sequence dictionary."),
      P("bwa_index", "input", "path", False, "Pre-built bwa index directory (optional)."),
      P("groupreadsbyumi_strategy", "option", "string", True, "fgbio grouping strategy.", "Adjacency", ["Identity", "Edit", "Adjacency", "Paired"]),
      P("aligner", "option", "string", True, "Aligner.", "bwa-mem", ["bwa-mem", "bwa-mem2"]),
      P("duplex", "option", "boolean", False, "Duplex UMI structure.", False),
      P("min_reads", "option", "string", True, "Min reads per consensus (int, or up to 3 space-separated for duplex).", "1"),
      P("min_baseq", "option", "integer", True, "Minimum base quality.", 10),
      P("max_base_error_rate", "option", "number", True, "Maximum base error rate.", 0.025)],
     READS + "    ref_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true), file(params.dict, checkIfExists: true)])\n"
     "    bwa_ch = params.bwa_index ? Channel.value([[id: 'reference'], file(params.bwa_index, checkIfExists: true)]) : Channel.empty()\n\n"
     "    FASTQ_CREATE_UMI_CONSENSUS_FGBIO(reads_ch, ref_ch, bwa_ch, params.groupreadsbyumi_strategy, params.aligner, params.duplex, params.min_reads, params.min_baseq, params.max_base_error_rate)\n",
     cpus=4, memory="8.GB", time="4.h")

spec("orftable_fasta_gtf_buildorfcatalogue", "ORFTABLE_FASTA_GTF_BUILDORFCATALOGUE", "Build ORF catalogue",
     "Normalise ribo-seq ORF tables and build a merged small-ORF catalogue.",
     [P("orf_table", "input", "path", True, "ORF table file."),
      P("caller", "option", "string", True, "ORF caller that produced the table.", "ribocode", ["ribocode", "ribotish", "ribotricer", "rpbp", "price"]),
      P("fasta", "input", "fasta", True, "Reference genome FASTA."),
      P("gtf", "input", "path", True, "Reference GTF."),
      P("collapse", "option", "boolean", False, "Cluster catalogue peptides by identity.", True)],
     "    orf_ch = Channel.value([[id: 'orf'], file(params.orf_table, checkIfExists: true), params.caller])\n"
     "    fasta_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true)])\n"
     "    gtf_ch = Channel.value([[id: 'reference'], file(params.gtf, checkIfExists: true)])\n\n"
     "    ORFTABLE_FASTA_GTF_BUILDORFCATALOGUE(orf_ch, fasta_ch, gtf_ch, params.collapse)\n",
     flags=["single-caller adapter; multi-caller routing needs curation"])

spec("bcl_demultiplex", "BCL_DEMULTIPLEX", "Demultiplex BCL",
     "Demultiplex Illumina BCL runs to FASTQ with bclconvert or bcl2fastq.",
     [P("samplesheet", "input", "path", True, "Illumina samplesheet CSV."),
      P("flowcell", "input", "path", True, "BCL run directory (flowcell)."),
      P("run_id", "option", "string", False, "Flowcell/run id.", "flowcell"),
      P("lane", "option", "integer", False, "Lane number.", 1),
      P("demultiplexer", "option", "string", True, "Demultiplexer to use.", "bclconvert", ["bclconvert", "bcl2fastq"])],
     "    flowcell_ch = Channel.value([[id: params.run_id, lane: params.lane], file(params.samplesheet, checkIfExists: true), file(params.flowcell, checkIfExists: true)])\n\n"
     "    BCL_DEMULTIPLEX(flowcell_ch, params.demultiplexer)\n",
     flags=["bcl2fastq/bclconvert containers are licensed Illumina tools"], cpus=4, memory="8.GB", time="4.h")

spec("fasta_newick_epang_gappa", "FASTA_NEWICK_EPANG_GAPPA", "Phylogenetic placement (EPA-ng/gappa)",
     "Place query sequences onto a reference phylogeny with EPA-ng and gappa.",
     [P("queryseqfile", "input", "fasta", True, "Query sequences FASTA."),
      P("refseqfile", "input", "fasta", True, "Reference sequences FASTA."),
      P("refphylogeny", "input", "path", True, "Reference phylogeny (newick)."),
      P("hmmfile", "input", "path", False, "HMM profile (optional; built if absent for hmmer)."),
      P("alignmethod", "option", "string", True, "Alignment method.", "hmmer", ["hmmer", "clustalo", "mafft"]),
      P("model", "option", "string", False, "Substitution model (optional).")],
     "    pp_ch = Channel.value([\n        meta: [id: 'epang'],\n        data: [\n            alignmethod: params.alignmethod,\n            queryseqfile: file(params.queryseqfile, checkIfExists: true),\n            refseqfile: file(params.refseqfile, checkIfExists: true),\n            refphylogeny: file(params.refphylogeny, checkIfExists: true),\n            hmmfile: params.hmmfile ? file(params.hmmfile, checkIfExists: true) : [],\n            model: params.model ?: null\n        ]\n    ])\n\n"
     "    FASTA_NEWICK_EPANG_GAPPA(pp_ch)\n",
     cpus=4, memory="8.GB", time="4.h")

spec("dia_proteomics_analysis", "DIA_PROTEOMICS_ANALYSIS", "DIA proteomics (DIA-NN)",
     "Quantify DIA proteomics runs with DIA-NN (library-free or spectral library).",
     [P("ms_file", "input", "path", True, "Mass-spec file (mzML)."),
      P("searchdb", "input", "fasta", True, "Search database FASTA."),
      P("expdesign", "input", "path", True, "Experimental design file."),
      P("enzyme", "option", "string", False, "Digestion enzyme.", "Trypsin"),
      P("fixed_mods", "option", "string", False, "Fixed modifications.", "Carbamidomethyl (C)"),
      P("variable_mods", "option", "string", False, "Variable modifications.", "Oxidation (M)"),
      P("precursor_tolerance", "option", "number", False, "Precursor tolerance.", 20),
      P("fragment_tolerance", "option", "number", False, "Fragment tolerance.", 20),
      P("precursor_tolerance_unit", "option", "string", False, "Precursor tolerance unit.", "ppm", ["ppm", "Da"]),
      P("fragment_tolerance_unit", "option", "string", False, "Fragment tolerance unit.", "ppm", ["ppm", "Da"]),
      P("random_preanalysis", "option", "boolean", False, "Enable random preanalysis.", False),
      P("random_preanalysis_n", "option", "integer", False, "Preanalysis sample count."),
      P("random_preanalysis_seed", "option", "integer", False, "Preanalysis seed."),
      P("scan_window", "option", "integer", False, "DIA-NN scan window.", 0),
      P("mass_acc_automatic", "option", "boolean", False, "Automatic mass accuracy.", True),
      P("scan_window_automatic", "option", "boolean", False, "Automatic scan window.", True),
      P("diann_debug", "option", "boolean", False, "Enable DIA-NN debug output.", False),
      P("pg_level", "option", "string", False, "Protein group level.", "genes"),
      P("speclib", "input", "path", False, "Pre-generated spectral library (optional)."),
      P("empirical_library", "input", "path", False, "Empirical library (optional)."),
      P("empirical_log", "input", "path", False, "Assembly log (optional).")],
     "    input_ch = Channel.value([[id: 'sample'], file(params.ms_file, checkIfExists: true), params.enzyme, params.fixed_mods, params.variable_mods, params.precursor_tolerance, params.fragment_tolerance, params.precursor_tolerance_unit, params.fragment_tolerance_unit])\n"
     "    searchdb_ch = Channel.value([[id: 'searchdb'], file(params.searchdb, checkIfExists: true)])\n"
     "    expdesign_ch = Channel.value([[id: 'experiment'], file(params.expdesign, checkIfExists: true)])\n"
     "    speclib_ch = params.speclib ? Channel.value(file(params.speclib, checkIfExists: true)) : Channel.empty()\n"
     "    empirical_ch = params.empirical_library ? Channel.value(file(params.empirical_library, checkIfExists: true)) : Channel.empty()\n"
     "    empirical_log_ch = params.empirical_log ? Channel.value(file(params.empirical_log, checkIfExists: true)) : Channel.empty()\n\n"
     "    DIA_PROTEOMICS_ANALYSIS(input_ch, searchdb_ch, expdesign_ch, [params.random_preanalysis, params.random_preanalysis_n, params.random_preanalysis_seed], params.scan_window, params.mass_acc_automatic, params.scan_window_automatic, params.diann_debug, params.pg_level, speclib_ch, empirical_ch, empirical_log_ch)\n",
     flags=["single-sample/experiment adapter; multi-file combinations need curation", "DIA-NN refuses conda profile"], cpus=8, memory="16.GB", time="8.h")

spec("cache_download_ensemblvep_snpeff", "CACHE_DOWNLOAD_ENSEMBLVEP_SNPEFF", "Download VEP + snpEff caches",
     "Download Ensembl VEP and snpEff annotation caches.",
     [P("vep_assembly", "option", "string", True, "VEP genome assembly.", "GRCh38"),
      P("vep_species", "option", "string", True, "VEP species.", "homo_sapiens"),
      P("vep_cache_version", "option", "string", True, "VEP cache version.", "110"),
      P("snpeff_db", "option", "string", True, "snpEff database version.", "GRCh38.99"),
      P("preflight_check", "option", "boolean", False, "Run VEP download preflight check.", False)],
     "    vep_info = Channel.value([[id: 'vep_cache'], params.vep_assembly, params.vep_species, params.vep_cache_version])\n"
     "    snpeff_info = Channel.value([[id: 'snpeff_cache'], params.snpeff_db])\n\n"
     "    CACHE_DOWNLOAD_ENSEMBLVEP_SNPEFF(vep_info, snpeff_info, params.preflight_check)\n",
     flags=["requires network access"], cpus=2, memory="4.GB", time="8.h")

# ---------------------------------------------------- imputation / phasing
def impute_meta_body():
    return "    def meta_panel = [id: params.panel_id, chr: params.chr]\n"

spec("bam_impute_quilt2", "BAM_IMPUTE_QUILT2", "Impute from BAM (QUILT2)",
     "Impute genotypes from low-coverage BAMs with QUILT2 against a reference panel.",
     [P("bam", "input", "bam", True, "Input BAM/CRAM."),
      P("bai", "input", "path", True, "BAM/CRAM index."),
      P("bampaths", "option", "string", False, "BAM path string for QUILT sample naming."),
      P("bamnames", "option", "string", False, "BAM sample names for QUILT."),
      P("panel_vcf", "input", "path", True, "Reference panel VCF."),
      P("panel_index", "input", "path", True, "Reference panel VCF index."),
      P("map", "input", "path", True, "Genetic map file."),
      P("panel_id", "option", "string", False, "Panel id for meta.", "panel"),
      P("chr", "option", "string", True, "Chromosome to impute.", "chr22"),
      P("start", "option", "integer", False, "Region start (optional)."),
      P("end", "option", "integer", False, "Region end (optional)."),
      P("fasta", "input", "fasta", False, "Reference FASTA (optional)."),
      P("fai", "input", "path", False, "Reference FASTA index (optional)."),
      P("n_gen", "option", "integer", False, "Generations since founding.", 100),
      P("buffer", "option", "integer", False, "Region buffer.", 100000),
      P("sample_id", "option", "string", False, "Sample id.", "sample")],
     impute_meta_body() +
     "    ch_input = Channel.value([[id: params.sample_id], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), params.bampaths ?: '', params.bamnames ?: ''])\n"
     "    ch_panel = Channel.value([meta_panel, file(params.panel_vcf, checkIfExists: true), file(params.panel_index, checkIfExists: true)])\n"
     "    ch_map = Channel.value([meta_panel, file(params.map, checkIfExists: true)])\n"
     "    ch_chunks = Channel.value([meta_panel, params.chr, params.start ?: [], params.end ?: []])\n"
     "    ch_fasta = params.fasta ? Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true)]) : Channel.empty()\n\n"
     "    BAM_IMPUTE_QUILT2(ch_input, ch_panel, ch_chunks, ch_map, ch_fasta, params.n_gen, params.buffer)\n",
     flags=["single-region adapter; multi-region scatter and bampaths/bamnames semantics need curation"], cpus=4, memory="8.GB", time="8.h")

spec("bam_impute_quilt", "BAM_IMPUTE_QUILT", "Impute from BAM (QUILT)",
     "Impute genotypes from low-coverage BAMs with QUILT (hap/legend/posfile panel).",
     [P("bam", "input", "bam", True, "Input BAM/CRAM."),
      P("bai", "input", "path", True, "BAM/CRAM index."),
      P("bampaths", "option", "string", False, "BAM path string for QUILT."),
      P("bamnames", "option", "string", False, "BAM sample names for QUILT."),
      P("hap", "input", "path", True, "Reference hap file."),
      P("legend", "input", "path", True, "Reference legend file."),
      P("posfile", "input", "path", True, "Reference positions file."),
      P("map", "input", "path", True, "Genetic map file."),
      P("panel_id", "option", "string", False, "Panel id for meta.", "panel"),
      P("chr", "option", "string", True, "Chromosome to impute.", "chr22"),
      P("start", "option", "integer", False, "Region start (optional)."),
      P("end", "option", "integer", False, "Region end (optional)."),
      P("fasta", "input", "fasta", False, "Reference FASTA (optional)."),
      P("fai", "input", "path", False, "Reference FASTA index (optional)."),
      P("n_gen", "option", "integer", False, "Generations since founding.", 100),
      P("buffer", "option", "integer", False, "Region buffer.", 100000),
      P("sample_id", "option", "string", False, "Sample id.", "sample")],
     impute_meta_body() +
     "    ch_input = Channel.value([[id: params.sample_id], [file(params.bam, checkIfExists: true)], [file(params.bai, checkIfExists: true)], params.bampaths ?: '', params.bamnames ?: ''])\n"
     "    ch_hlp = Channel.value([meta_panel, file(params.hap, checkIfExists: true), file(params.legend, checkIfExists: true), file(params.posfile, checkIfExists: true)])\n"
     "    ch_map = Channel.value([meta_panel, file(params.map, checkIfExists: true)])\n"
     "    ch_chunks = Channel.value([meta_panel, params.chr, params.start ?: [], params.end ?: []])\n"
     "    ch_fasta = params.fasta ? Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true)]) : Channel.empty()\n\n"
     "    BAM_IMPUTE_QUILT(ch_input, ch_hlp, ch_chunks, ch_map, ch_fasta, params.n_gen, params.buffer)\n",
     flags=["single-region adapter; multi-region scatter needs curation"], cpus=4, memory="8.GB", time="8.h")

spec("bam_impute_stitch", "BAM_IMPUTE_STITCH", "Impute from BAM (STITCH)",
     "Impute genotypes from low-coverage BAMs with STITCH.",
     [P("bam", "input", "bam", True, "Input BAM/CRAM."),
      P("bai", "input", "path", True, "BAM/CRAM index."),
      P("bampaths", "option", "string", False, "BAM path string for STITCH."),
      P("bamnames", "option", "string", False, "BAM sample names for STITCH."),
      P("posfile", "input", "path", True, "Reference positions file."),
      P("map", "input", "path", True, "Genetic map file."),
      P("panel_id", "option", "string", False, "Panel id for meta.", "panel"),
      P("chr", "option", "string", True, "Chromosome to impute.", "chr22"),
      P("start", "option", "integer", False, "Region start (optional)."),
      P("end", "option", "integer", False, "Region end (optional)."),
      P("fasta", "input", "fasta", False, "Reference FASTA (optional)."),
      P("fai", "input", "path", False, "Reference FASTA index (optional)."),
      P("k_val", "option", "integer", False, "STITCH k value.", 2),
      P("n_gen", "option", "integer", False, "Generations since founding.", 100),
      P("buffer", "option", "integer", False, "Region buffer.", 100000),
      P("seed", "option", "integer", False, "Random seed.", 1),
      P("sample_id", "option", "string", False, "Sample id.", "sample")],
     impute_meta_body() +
     "    ch_input = Channel.value([[id: params.sample_id], [file(params.bam, checkIfExists: true)], [file(params.bai, checkIfExists: true)], params.bampaths ?: '', params.bamnames ?: ''])\n"
     "    ch_posfile = Channel.value([meta_panel, file(params.posfile, checkIfExists: true)])\n"
     "    ch_map = Channel.value([meta_panel, file(params.map, checkIfExists: true)])\n"
     "    ch_chunks = Channel.value([meta_panel, params.chr, params.start ?: [], params.end ?: []])\n"
     "    ch_fasta = params.fasta ? Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true)]) : Channel.empty()\n\n"
     "    BAM_IMPUTE_STITCH(ch_input, ch_posfile, ch_chunks, ch_map, ch_fasta, params.k_val, params.n_gen, params.buffer, params.seed)\n",
     flags=["single-region adapter; multi-region scatter needs curation"], cpus=4, memory="8.GB", time="8.h")

spec("bam_vcf_impute_glimpse2", "BAM_VCF_IMPUTE_GLIMPSE2", "Impute VCF (GLIMPSE2)",
     "Impute a low-coverage VCF with GLIMPSE2 against a reference panel.",
     [P("vcf", "input", "path", True, "Input VCF."),
      P("vcf_index", "input", "path", True, "Input VCF index (.csi)."),
      P("infos", "option", "string", False, "Sample/info list (optional)."),
      P("ref_vcf", "input", "path", True, "Reference panel VCF."),
      P("ref_index", "input", "path", True, "Reference panel VCF index."),
      P("region", "option", "string", True, "Reference region label.", "chr22"),
      P("map", "input", "path", True, "Genetic map file."),
      P("panel_id", "option", "string", False, "Panel id for meta.", "panel"),
      P("chr", "option", "string", True, "Chromosome.", "chr22"),
      P("fasta", "input", "fasta", False, "Reference FASTA (optional)."),
      P("fai", "input", "path", False, "Reference FASTA index (optional)."),
      P("chunk", "option", "boolean", False, "Activate automatic chunking.", True),
      P("chunk_model", "input", "path", False, "Chunking model file (optional)."),
      P("splitreference", "option", "boolean", False, "Split reference step.", False),
      P("sample_id", "option", "string", False, "Sample id.", "sample")],
     impute_meta_body() +
     "    ch_input = Channel.value([[id: params.sample_id], file(params.vcf, checkIfExists: true), file(params.vcf_index, checkIfExists: true), params.infos ?: [], params.infos ?: []])\n"
     "    ch_ref = Channel.value([meta_panel, file(params.ref_vcf, checkIfExists: true), file(params.ref_index, checkIfExists: true), params.region])\n"
     "    ch_map = Channel.value([meta_panel, file(params.map, checkIfExists: true)])\n"
     "    ch_chunks = Channel.empty()\n"
     "    ch_fasta = params.fasta ? Channel.value([meta_panel, file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true)]) : Channel.empty()\n\n"
     "    BAM_VCF_IMPUTE_GLIMPSE2(ch_input, ch_ref, ch_chunks, ch_map, ch_fasta, params.chunk, params.chunk_model ?: null, params.splitreference)\n",
     flags=["chunk=true regenerates chunks; chunk_model needed for GLIMPSE2 chunking"], cpus=4, memory="8.GB", time="8.h")

spec("vcf_impute_glimpse", "VCF_IMPUTE_GLIMPSE", "Impute VCF (GLIMPSE1)",
     "Impute a low-coverage VCF with GLIMPSE against a reference panel.",
     [P("vcf", "input", "path", True, "Input VCF."),
      P("vcf_index", "input", "path", True, "Input VCF index (.csi)."),
      P("infos", "option", "string", False, "Sample/info list (optional)."),
      P("ref_vcf", "input", "path", True, "Reference panel VCF."),
      P("ref_index", "input", "path", True, "Reference panel VCF index."),
      P("region", "option", "string", True, "Reference region label.", "chr22"),
      P("map", "input", "path", True, "Genetic map file."),
      P("panel_id", "option", "string", False, "Panel id for meta.", "panel"),
      P("chr", "option", "string", True, "Chromosome.", "chr22"),
      P("chunk", "option", "boolean", False, "Activate automatic chunking.", True),
      P("sample_id", "option", "string", False, "Sample id.", "sample")],
     impute_meta_body() +
     "    ch_vcf = Channel.value([[id: params.sample_id], file(params.vcf, checkIfExists: true), file(params.vcf_index, checkIfExists: true), params.infos ?: []])\n"
     "    ch_ref = Channel.value([meta_panel, file(params.ref_vcf, checkIfExists: true), file(params.ref_index, checkIfExists: true), params.region])\n"
     "    ch_map = Channel.value([meta_panel, file(params.map, checkIfExists: true)])\n"
     "    ch_chunks = Channel.empty()\n\n"
     "    VCF_IMPUTE_GLIMPSE(ch_vcf, ch_ref, ch_chunks, ch_map, params.chunk)\n",
     flags=["chunk=true regenerates chunks from the reference panel"], cpus=4, memory="8.GB", time="8.h")

spec("vcf_impute_beagle5", "VCF_IMPUTE_BEAGLE5", "Impute VCF (BEAGLE5)",
     "Impute a VCF with BEAGLE5 against a reference panel.",
     [P("vcf", "input", "path", True, "Input VCF/BCF."),
      P("vcf_index", "input", "path", True, "Input VCF index."),
      P("panel_vcf", "input", "path", True, "Reference panel VCF."),
      P("panel_index", "input", "path", True, "Reference panel VCF index."),
      P("map", "input", "path", True, "Genetic map file."),
      P("region", "option", "string", True, "Region to impute (chr or chr:start-end).", "chr22"),
      P("panel_id", "option", "string", False, "Panel id for meta.", "panel"),
      P("chr", "option", "string", True, "Chromosome.", "chr22"),
      P("sample_id", "option", "string", False, "Sample id.", "sample")],
     impute_meta_body() +
     "    ch_input = Channel.value([[id: params.sample_id], file(params.vcf, checkIfExists: true), file(params.vcf_index, checkIfExists: true)])\n"
     "    ch_panel = Channel.value([meta_panel, file(params.panel_vcf, checkIfExists: true), file(params.panel_index, checkIfExists: true)])\n"
     "    ch_map = Channel.value([meta_panel, file(params.map, checkIfExists: true)])\n"
     "    ch_chunks = Channel.value([meta_panel, params.region])\n\n"
     "    VCF_IMPUTE_BEAGLE5(ch_input, ch_panel, ch_chunks, ch_map)\n",
     flags=["single-region adapter; multi-region scatter needs curation"], cpus=4, memory="8.GB", time="8.h")

spec("vcf_impute_minimac4", "VCF_IMPUTE_MINIMAC4", "Impute VCF (Minimac4)",
     "Impute a VCF with Minimac4 against a reference panel (VCF or MSAV).",
     [P("vcf", "input", "path", True, "Input VCF."),
      P("vcf_index", "input", "path", True, "Input VCF index."),
      P("panel_vcf", "input", "path", True, "Reference panel VCF/BCF or MSAV."),
      P("panel_index", "input", "path", True, "Reference panel index."),
      P("sites_vcf", "input", "path", True, "Sites-to-impute VCF."),
      P("sites_index", "input", "path", True, "Sites VCF index."),
      P("map", "input", "path", True, "Genetic map file."),
      P("region", "option", "string", True, "Region to impute.", "chr22"),
      P("panel_id", "option", "string", False, "Panel id for meta.", "panel"),
      P("chr", "option", "string", True, "Chromosome.", "chr22"),
      P("sample_id", "option", "string", False, "Sample id.", "sample")],
     impute_meta_body() +
     "    ch_input = Channel.value([[id: params.sample_id, chr: params.chr], file(params.vcf, checkIfExists: true), file(params.vcf_index, checkIfExists: true)])\n"
     "    ch_panel = Channel.value([meta_panel, file(params.panel_vcf, checkIfExists: true), file(params.panel_index, checkIfExists: true)])\n"
     "    ch_posfile = Channel.value([meta_panel, file(params.sites_vcf, checkIfExists: true), file(params.sites_index, checkIfExists: true)])\n"
     "    ch_map = Channel.value([meta_panel, file(params.map, checkIfExists: true)])\n"
     "    ch_chunks = Channel.value([meta_panel, params.region])\n\n"
     "    VCF_IMPUTE_MINIMAC4(ch_input, ch_panel, ch_posfile, ch_chunks, ch_map)\n",
     flags=["single-region adapter; multi-region scatter needs curation"], cpus=4, memory="8.GB", time="8.h")

spec("vcf_phase_shapeit5", "VCF_PHASE_SHAPEIT5", "Phase VCF (SHAPEIT5)",
     "Phase a VCF with SHAPEIT5 (common variants), optionally chunked.",
     [P("vcf", "input", "path", True, "Input VCF."),
      P("vcf_index", "input", "path", True, "Input VCF index."),
      P("ref_vcf", "input", "path", True, "Reference panel VCF."),
      P("ref_index", "input", "path", True, "Reference panel VCF index."),
      P("scaffold_vcf", "input", "path", True, "Scaffold VCF."),
      P("scaffold_index", "input", "path", True, "Scaffold VCF index."),
      P("map", "input", "path", True, "Genetic map file."),
      P("chunk_model", "input", "path", True, "Chunking model file (required when chunk=true)."),
      P("region", "option", "string", True, "Region label.", "chr22"),
      P("pedigree", "input", "path", False, "Pedigree file (optional)."),
      P("panel_id", "option", "string", False, "Sample/panel id for meta.", "sample"),
      P("chr", "option", "string", True, "Chromosome.", "chr22"),
      P("chunk", "option", "boolean", False, "Activate automatic chunking.", True)],
     "    def meta = [id: params.panel_id, chr: params.chr]\n"
     "    ch_vcf = Channel.value([meta, file(params.vcf, checkIfExists: true), file(params.vcf_index, checkIfExists: true), params.pedigree ? file(params.pedigree, checkIfExists: true) : [], params.region])\n"
     "    ch_ref = Channel.value([meta, file(params.ref_vcf, checkIfExists: true), file(params.ref_index, checkIfExists: true)])\n"
     "    ch_scaffold = Channel.value([meta, file(params.scaffold_vcf, checkIfExists: true), file(params.scaffold_index, checkIfExists: true)])\n"
     "    ch_map = Channel.value([meta, file(params.map, checkIfExists: true)])\n"
     "    ch_chunks = Channel.empty()\n"
     "    chunk_model_ch = Channel.value(file(params.chunk_model, checkIfExists: true))\n\n"
     "    VCF_PHASE_SHAPEIT5(ch_vcf, ch_chunks, ch_ref, ch_scaffold, ch_map, params.chunk, chunk_model_ch)\n",
     flags=["chunk=true regenerates chunks; ref/scaffold/map required for the join"], cpus=4, memory="8.GB", time="8.h")

# ---------------------------------------------------- somatic / telomere / PON
spec("bam_applybqsr", "BAM_APPLYBQSR", "Apply BQSR",
     "Apply GATK base quality score recalibration to a BAM/CRAM.",
     [P("bam", "input", "bam", True, "Input BAM/CRAM."),
      P("bai", "input", "path", True, "Input index."),
      P("recal", "input", "path", True, "BQSR recalibration table."),
      P("fasta", "input", "fasta", True, "Reference FASTA."),
      P("fai", "input", "path", True, "Reference FASTA index."),
      P("dict", "input", "path", True, "Reference dictionary."),
      P("intervals", "input", "path", False, "Interval file (optional)."),
      P("output_suffix", "option", "string", True, "Output format.", "bam", ["bam", "cram"]),
      P("sample_id", "option", "string", False, "Sample id.", "sample")],
     "    input_ch = Channel.value([[id: params.sample_id], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), file(params.recal, checkIfExists: true)])\n"
     "    references_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true), file(params.dict, checkIfExists: true)])\n"
     "    intervals_ch = params.intervals ? Channel.value([file(params.intervals, checkIfExists: true), 1]) : Channel.value([[], 0])\n\n"
     "    BAM_APPLYBQSR(input_ch, references_ch, intervals_ch, params.output_suffix)\n",
     flags=["single-interval adapter; scatter-gather needs curation"])

spec("bam_create_som_pon_gatk", "BAM_CREATE_SOM_PON_GATK", "Create somatic PON (GATK)",
     "Call normal samples with Mutect2 and build a GATK somatic panel of normals.",
     [P("bam", "input", "bam", True, "Normal-sample BAM/CRAM."),
      P("bai", "input", "path", True, "Input index."),
      P("fasta", "input", "fasta", True, "Reference FASTA."),
      P("fai", "input", "path", True, "Reference FASTA index."),
      P("dict", "input", "path", True, "Reference dictionary."),
      P("pon_norm", "option", "string", True, "Panel of normals name.", "pon"),
      P("intervals_gendb", "input", "path", False, "Interval file for GenomicsDBImport (optional)."),
      P("intervals", "input", "path", False, "Interval file for scattering (optional)."),
      P("sample_id", "option", "string", False, "Sample id.", "normal")],
     "    input_ch = Channel.value([[id: params.sample_id], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])\n"
     "    fasta_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true)])\n"
     "    fai_ch = Channel.value([[id: 'reference'], file(params.fai, checkIfExists: true), []])\n"
     "    dict_ch = Channel.value([[id: 'reference'], file(params.dict, checkIfExists: true)])\n"
     "    gendb_ch = params.intervals_gendb ? Channel.value(file(params.intervals_gendb, checkIfExists: true)) : Channel.empty()\n"
     "    intervals_num_ch = params.intervals ? Channel.value([file(params.intervals, checkIfExists: true), 1]) : Channel.value([[], 0])\n\n"
     "    BAM_CREATE_SOM_PON_GATK(input_ch, fasta_ch, fai_ch, dict_ch, params.pon_norm, gendb_ch, intervals_num_ch)\n",
     flags=["single-normal adapter; a real PON needs multiple normal samples", "GenomicsDBImport needs an interval file"], cpus=4, memory="8.GB", time="8.h")

spec("bam_tumor_only_somatic_variant_calling_gatk", "BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK", "Tumor-only somatic calling (GATK)",
     "Call somatic variants from a tumor-only BAM with GATK Mutect2 + contamination filtering.",
     [P("bam", "input", "bam", True, "Tumor BAM/CRAM."),
      P("bai", "input", "path", True, "Input index."),
      P("fasta", "input", "fasta", True, "Reference FASTA."),
      P("fai", "input", "path", True, "Reference FASTA index."),
      P("dict", "input", "path", True, "Reference dictionary."),
      P("germline_resource", "input", "path", True, "Germline resource VCF (e.g. gnomAD)."),
      P("germline_resource_tbi", "input", "path", True, "Germline resource index."),
      P("alleles", "input", "path", False, "Alleles VCF (optional)."),
      P("alleles_tbi", "input", "path", False, "Alleles VCF index (optional)."),
      P("panel_of_normals", "input", "path", False, "Panel of normals VCF (optional)."),
      P("panel_of_normals_tbi", "input", "path", False, "PON index (optional)."),
      P("interval_file", "input", "path", False, "Interval file (optional)."),
      P("sample_id", "option", "string", False, "Sample id.", "tumor")],
     "    input_ch = Channel.value([[id: params.sample_id], [file(params.bam, checkIfExists: true)], [file(params.bai, checkIfExists: true)], []])\n"
     "    fasta_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true)])\n"
     "    fai_ch = Channel.value([[id: 'reference'], file(params.fai, checkIfExists: true), []])\n"
     "    dict_ch = Channel.value([[id: 'reference'], file(params.dict, checkIfExists: true)])\n"
     "    alleles_ch = params.alleles ? Channel.value(file(params.alleles, checkIfExists: true)) : Channel.value([])\n"
     "    alleles_tbi_ch = params.alleles_tbi ? Channel.value(file(params.alleles_tbi, checkIfExists: true)) : Channel.value([])\n"
     "    germline_ch = Channel.value(file(params.germline_resource, checkIfExists: true))\n"
     "    germline_tbi_ch = Channel.value(file(params.germline_resource_tbi, checkIfExists: true))\n"
     "    pon_ch = params.panel_of_normals ? Channel.value(file(params.panel_of_normals, checkIfExists: true)) : Channel.value([])\n"
     "    pon_tbi_ch = params.panel_of_normals_tbi ? Channel.value(file(params.panel_of_normals_tbi, checkIfExists: true)) : Channel.value([])\n"
     "    interval_ch = params.interval_file ? Channel.value(file(params.interval_file, checkIfExists: true)) : Channel.value([])\n\n"
     "    BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK(input_ch, fasta_ch, fai_ch, dict_ch, alleles_ch, alleles_tbi_ch, germline_ch, germline_tbi_ch, pon_ch, pon_tbi_ch, interval_ch)\n",
     flags=["single-sample adapter; reference resource channels are plain file values"], cpus=4, memory="8.GB", time="8.h")

spec("bam_tumor_normal_somatic_variant_calling_strelka", "BAM_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING_STRELKA", "Tumor-normal somatic calling (Strelka)",
     "Call somatic SNVs/indels from matched tumor-normal CRAMs with Strelka.",
     [P("normal_cram", "input", "bam", True, "Normal CRAM/BAM."),
      P("normal_crai", "input", "path", True, "Normal index."),
      P("tumor_cram", "input", "bam", True, "Tumor CRAM/BAM."),
      P("tumor_crai", "input", "path", True, "Tumor index."),
      P("fasta", "input", "fasta", True, "Reference FASTA."),
      P("fai", "input", "path", True, "Reference FASTA index."),
      P("dict", "input", "path", False, "Reference dictionary (optional)."),
      P("manta_vcf", "input", "path", False, "Manta indel candidates VCF (optional)."),
      P("manta_tbi", "input", "path", False, "Manta VCF index (optional)."),
      P("intervals", "input", "path", False, "Interval BED.gz (optional)."),
      P("intervals_tbi", "input", "path", False, "Interval BED.gz index (optional)."),
      P("sample_id", "option", "string", False, "Sample id.", "tumor_normal")],
     "    cram_ch = Channel.value([[id: params.sample_id], file(params.normal_cram, checkIfExists: true), file(params.normal_crai, checkIfExists: true), file(params.tumor_cram, checkIfExists: true), file(params.tumor_crai, checkIfExists: true), params.manta_vcf ? file(params.manta_vcf, checkIfExists: true) : [], params.manta_tbi ? file(params.manta_tbi, checkIfExists: true) : []])\n"
     "    dict_ch = params.dict ? Channel.value([[id: 'reference'], file(params.dict, checkIfExists: true)]) : Channel.empty()\n"
     "    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))\n"
     "    fai_ch = Channel.value(file(params.fai, checkIfExists: true))\n"
     "    intervals_ch = params.intervals ? Channel.value([file(params.intervals, checkIfExists: true), file(params.intervals_tbi, checkIfExists: true), 1]) : Channel.value([[], [], 0])\n\n"
     "    BAM_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING_STRELKA(cram_ch, dict_ch, fasta_ch, fai_ch, intervals_ch)\n",
     flags=["single-pair adapter; interval scattering needs curation"], cpus=4, memory="8.GB", time="8.h")

spec("bam_telomere_estimation", "BAM_TELOMERE_ESTIMATION", "Telomere length/content estimation",
     "Estimate telomere length (telseq/telogator2) and optionally content (telomerehunter).",
     [P("bam", "input", "bam", True, "Aligned BAM."),
      P("bai", "input", "path", True, "BAM index."),
      P("data_type", "option", "string", True, "Read type.", "short", ["short", "long"]),
      P("bed", "input", "path", False, "Exome BED for telseq (optional)."),
      P("control_bam", "input", "bam", False, "Matched control BAM (optional)."),
      P("control_bai", "input", "path", False, "Control BAM index (optional)."),
      P("fasta", "input", "fasta", True, "Reference FASTA."),
      P("fai", "input", "path", True, "Reference FASTA index."),
      P("cytoband", "input", "path", False, "Cytoband file for telomerehunter (optional)."),
      P("run_telomerehunter", "option", "boolean", False, "Enable telomerehunter content profiling.", False),
      P("length_estimator", "option", "string", False, "Override length estimator (short/long)."),
      P("sample_id", "option", "string", False, "Sample id.", "sample")],
     "    reads_ch = Channel.value([[id: params.sample_id], params.data_type, file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), params.bed ? file(params.bed, checkIfExists: true) : []])\n"
     "    control_ch = params.control_bam ? Channel.value([[id: params.sample_id], file(params.control_bam, checkIfExists: true), file(params.control_bai, checkIfExists: true)]) : Channel.empty()\n"
     "    fasta_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true)])\n"
     "    cytoband_ch = Channel.value(params.cytoband ? [file(params.cytoband, checkIfExists: true)] : [[]])\n\n"
     "    BAM_TELOMERE_ESTIMATION(reads_ch, control_ch, fasta_ch, cytoband_ch, params.run_telomerehunter, params.length_estimator ?: null)\n",
     flags=["single-sample adapter"])

spec("snpclustering", "SNPCLUSTERING", "SNP clustering (PCA)",
     "Impute with BEAGLE5, run plink2 PCA and cluster samples (k-means/DBSCAN).",
     [P("vcf", "input", "path", True, "Input VCF."),
      P("vcf_index", "input", "path", True, "Input VCF index."),
      P("refpanel", "input", "path", False, "Reference panel (optional)."),
      P("genmap", "input", "path", False, "Genetic map (optional)."),
      P("region", "option", "string", False, "Region for BEAGLE (optional)."),
      P("n_pcs", "option", "integer", False, "Number of PCs.", 10),
      P("use_approx", "option", "boolean", False, "Use approximate PCA.", False),
      P("algorithm", "option", "string", False, "Clustering algorithm.", "kmeans", ["kmeans", "dbscan"]),
      P("n_clusters", "option", "integer", False, "Number of clusters (kmeans).", 3),
      P("dbscan_eps", "option", "number", False, "DBSCAN eps.", 0.5),
      P("dbscan_min_samples", "option", "integer", False, "DBSCAN min samples.", 5),
      P("sample_id", "option", "string", False, "Sample id.", "sample")],
     "    vcf_ch = Channel.value([[id: params.sample_id], file(params.vcf, checkIfExists: true), file(params.vcf_index, checkIfExists: true)])\n"
     "    refpanel_ch = params.refpanel ? Channel.value(file(params.refpanel, checkIfExists: true)) : Channel.empty()\n"
     "    genmap_ch = params.genmap ? Channel.value(file(params.genmap, checkIfExists: true)) : Channel.empty()\n\n"
     "    SNPCLUSTERING(vcf_ch, refpanel_ch, genmap_ch, params.region ?: null, params.n_pcs, params.use_approx, params.algorithm, params.n_clusters, params.dbscan_eps, params.dbscan_min_samples)\n",
     flags=["refpanel/genmap usage depends on BEAGLE configuration; verify semantics"])

# ---------------------------------------------------- mega composite
spec("fastq_shortreads_preprocess_qc", "FASTQ_SHORTREADS_PREPROCESS_QC", "Short-read preprocessing + QC",
     "Full short-read preprocessing: QC, UMI extract, adapter trim, complexity filter, dedup, decontamination.",
     [P("reads", "input", "path_glob", True, "FASTQ glob (paired or single)."),
      P("skip_fastqc", "option", "boolean", False, "Skip FastQC.", False),
      P("skip_seqfu_check", "option", "boolean", False, "Skip seqfu check.", False),
      P("skip_seqfu_stats", "option", "boolean", False, "Skip seqfu stats.", False),
      P("skip_seqkit_stats", "option", "boolean", False, "Skip seqkit stats.", False),
      P("skip_seqtk_comp", "option", "boolean", False, "Skip seqtk comp.", False),
      P("skip_seqkit_sana_pair", "option", "boolean", False, "Skip seqkit sana pair.", False),
      P("skip_seqkit_seq", "option", "boolean", False, "Skip seqkit seq.", False),
      P("skip_seqkit_replace", "option", "boolean", False, "Skip seqkit replace.", False),
      P("skip_seqkit_rmdup", "option", "boolean", False, "Skip seqkit rmdup.", False),
      P("skip_umitools_extract", "option", "boolean", False, "Skip UMI-tools extract.", True),
      P("umi_discard_read", "option", "integer", False, "UMI discard read (0/1/2).", 0),
      P("skip_adapterremoval", "option", "boolean", False, "Skip adapter removal.", False),
      P("adapter_tool", "option", "string", True, "Adapter trimming tool.", "fastp",
        ["trimmomatic", "cutadapt", "trimgalore", "bbduk", "leehom", "fastp", "adapterremoval"]),
      P("custom_adapters_file", "input", "path", False, "Custom adapters file (optional)."),
      P("save_merged", "option", "boolean", False, "Return merged reads.", False),
      P("fastp_discard_trimmed_pass", "option", "boolean", False, "fastp: discard trimmed-pass.", False),
      P("fastp_save_trimmed_fail", "option", "boolean", False, "fastp: save trimmed-fail.", False),
      P("skip_complexity_filtering", "option", "boolean", False, "Skip complexity filtering.", False),
      P("complexity_filter_tool", "option", "string", True, "Complexity filter tool.", "fastp", ["prinseqplusplus", "bbduk", "fastp"]),
      P("skip_deduplication", "option", "boolean", False, "Skip deduplication.", False),
      P("skip_decontamination", "option", "boolean", False, "Skip decontamination.", True),
      P("decontamination_fasta", "input", "fasta", False, "Decontamination reference FASTA (optional)."),
      P("decontamination_reference", "input", "path", False, "Decontamination reference dir (optional)."),
      P("decontamination_reference_name", "option", "string", False, "Decontamination reference name.", "contaminant"),
      P("decontamination_index_name", "option", "string", False, "deacon index name (optional)."),
      P("decontamination_tool", "option", "string", False, "Decontamination tool.", "hostile", ["hostile", "deacon"]),
      P("skip_final_concatenation", "option", "boolean", False, "Skip final concatenation.", False)],
     READS +
     "    adapters_ch = params.custom_adapters_file ? Channel.value(file(params.custom_adapters_file, checkIfExists: true)) : Channel.empty()\n"
     "    decontam_fasta_ch = params.decontamination_fasta ? Channel.value([[id: params.decontamination_reference_name], [file(params.decontamination_fasta, checkIfExists: true)]]) : Channel.empty()\n"
     "    decontam_ref_ch = params.decontamination_reference ? Channel.value([params.decontamination_reference_name, file(params.decontamination_reference, checkIfExists: true)]) : Channel.empty()\n\n"
     "    FASTQ_SHORTREADS_PREPROCESS_QC(reads_ch, params.skip_fastqc, params.skip_seqfu_check, params.skip_seqfu_stats, params.skip_seqkit_stats, params.skip_seqtk_comp, params.skip_seqkit_sana_pair, params.skip_seqkit_seq, params.skip_seqkit_replace, params.skip_seqkit_rmdup, params.skip_umitools_extract, params.umi_discard_read, params.skip_adapterremoval, params.adapter_tool, adapters_ch, params.save_merged, params.fastp_discard_trimmed_pass, params.fastp_save_trimmed_fail, params.skip_complexity_filtering, params.complexity_filter_tool, params.skip_deduplication, params.skip_decontamination, decontam_fasta_ch, decontam_ref_ch, params.decontamination_index_name ?: null, params.decontamination_tool, params.skip_final_concatenation)\n",
     flags=["27-input composite; verify each stage's companion semantics for multi-sample use"], cpus=4, memory="8.GB", time="4.h")


CONFIG_TMPL = """// Launcher config for the {rel} subworkflow wrapper adapter. Not part of the
// fixed wrapper/{{main.nf,params.json,wrapper.yaml}} contract — Nextflow picks
// this up automatically because it sits beside wrapper/main.nf.
//
// Global publishDir collects every nested module's outputs into one results
// directory (matches the generated subworkflow wrappers).

process {{
    cpus   = {cpus}
    memory = '{memory}'
    time   = '{time}'

    publishDir = [
        path: {{ "${{params.outdir}}/{proc_dir}" }},
        mode: 'copy'
    ]
}}

profiles {{
    docker {{
        docker.enabled = true
        docker.runOptions = '-u $(id -u):$(id -g)'
        conda.enabled = false
        singularity.enabled = false
    }}
    singularity {{
        singularity.enabled = true
        singularity.autoMounts = true
        docker.enabled = false
        conda.enabled = false
    }}
    conda {{
        conda.enabled = true
        docker.enabled = false
        singularity.enabled = false
    }}
}}
"""


def yaml_param(name, kind, typ, req, desc, enum):
    lines = [f"  {name}:", f"    kind: {kind}", f"    type: {typ}", f"    required: {'true' if req else 'false'}"]
    if enum:
        lines.append("    enum:")
        for e in enum:
            lines.append(f"      - {e}")
    lines.append(f"    description: {json.dumps(desc)}")
    return "\n".join(lines)


def emit(s):
    rel = s["rel"]
    wdir = os.path.join(ROOT, rel, "wrapper")
    existing = os.path.join(wdir, "wrapper.yaml")
    if os.path.exists(existing):
        with open(existing) as f:
            if "wire-manual-subworkflows.py" not in f.read():
                print("skip (foreign/curated)", rel)
                return
    proc_dir = s["workflow"].lower()
    comp_id = "nf-core/subworkflows/" + rel.replace("_", "-")
    params = s["params"] + [P("outdir", "output", "path", True, "Output directory for the subworkflow's results.", "results")]
    defaults = {}
    for (name, kind, typ, req, desc, default, enum) in params:
        defaults[name] = default
    yaml_lines = [
        "# Agent-facing wrapper adapter for the vendored " + rel + " subworkflow",
        "# at ../main.nf. Includes the subworkflow only — never another wrapper;",
        "# the subworkflow references the already-migrated module directories.",
        "# Follows docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.",
    ]
    if s["flags"]:
        yaml_lines.append("# REVIEW: " + "; ".join(s["flags"]))
    yaml_lines.append(f"id: {comp_id}")
    yaml_lines.append(f"name: {json.dumps(s['name'])}")
    yaml_lines.append(f"summary: {json.dumps(s['summary'])}")
    yaml_lines.append("")
    yaml_lines.append("params:")
    for p in params:
        yaml_lines.append(yaml_param(*p[:5], p[6]))
    yaml_lines.append("")
    yaml_lines.append("outputs:")
    yaml_lines.append("  results:")
    yaml_lines.append("    type: directory")
    yaml_lines.append("    path: '${outdir}/" + proc_dir + "'")
    yaml_lines.append("    primary: true")
    yaml_text = "\n".join(yaml_lines) + "\n"

    decls = []
    for (name, kind, typ, req, desc, default, enum) in params:
        decls.append(f"params.{name} = null" if default is None else f"params.{name} = {json.dumps(default)}")
    header = (
        "#!/usr/bin/env nextflow\n"
        "// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.\n"
        "// Includes the subworkflow only — never another wrapper; the subworkflow\n"
        "// itself references the already-migrated modules under resources/wrappers/modules.\n"
        "// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.\n"
        "nextflow.enable.dsl = 2\n\n"
        f"include {{ {s['workflow']} }} from '../main.nf'\n\n"
    )
    main_nf = header + "\n".join(decls) + "\n\nworkflow {\n" + s["body"] + "}\n"
    config = CONFIG_TMPL.format(rel=rel, proc_dir=proc_dir, cpus=s["cpus"], memory=s["memory"], time=s["time"])

    os.makedirs(wdir, exist_ok=True)
    with open(os.path.join(wdir, "wrapper.yaml"), "w") as f:
        f.write(yaml_text)
    with open(os.path.join(wdir, "main.nf"), "w") as f:
        f.write(main_nf)
    with open(os.path.join(wdir, "params.json"), "w") as f:
        f.write(json.dumps(defaults, indent=2) + "\n")
    with open(os.path.join(wdir, "nextflow.config"), "w") as f:
        f.write(config)
    print("wired", rel)


def main():
    for rel, reason in SKIP.items():
        print("skip (infra)", rel, "-", reason)
    for s in SPECS:
        emit(s)
    print(f"\n{len(SPECS)} specs, {len(SKIP)} infrastructure skips")


if __name__ == "__main__":
    main()