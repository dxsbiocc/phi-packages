#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTQ_SHORTREADS_PREPROCESS_QC } from '../main.nf'

params.reads = null
params.skip_fastqc = false
params.skip_seqfu_check = false
params.skip_seqfu_stats = false
params.skip_seqkit_stats = false
params.skip_seqtk_comp = false
params.skip_seqkit_sana_pair = false
params.skip_seqkit_seq = false
params.skip_seqkit_replace = false
params.skip_seqkit_rmdup = false
params.skip_umitools_extract = true
params.umi_discard_read = 0
params.skip_adapterremoval = false
params.adapter_tool = "fastp"
params.custom_adapters_file = null
params.save_merged = false
params.fastp_discard_trimmed_pass = false
params.fastp_save_trimmed_fail = false
params.skip_complexity_filtering = false
params.complexity_filter_tool = "fastp"
params.skip_deduplication = false
params.skip_decontamination = true
params.decontamination_fasta = null
params.decontamination_reference = null
params.decontamination_reference_name = "contaminant"
params.decontamination_index_name = null
params.decontamination_tool = "hostile"
params.skip_final_concatenation = false
params.outdir = "results"

workflow {
    reads_ch = Channel
        .fromFilePairs(params.reads, size: -1)
        .map { id, reads -> [[id: id, single_end: reads.size() == 1], reads] }
    adapters_ch = params.custom_adapters_file ? Channel.value(file(params.custom_adapters_file, checkIfExists: true)) : Channel.empty()
    decontam_fasta_ch = params.decontamination_fasta ? Channel.value([[id: params.decontamination_reference_name], [file(params.decontamination_fasta, checkIfExists: true)]]) : Channel.empty()
    decontam_ref_ch = params.decontamination_reference ? Channel.value([params.decontamination_reference_name, file(params.decontamination_reference, checkIfExists: true)]) : Channel.empty()

    FASTQ_SHORTREADS_PREPROCESS_QC(reads_ch, params.skip_fastqc, params.skip_seqfu_check, params.skip_seqfu_stats, params.skip_seqkit_stats, params.skip_seqtk_comp, params.skip_seqkit_sana_pair, params.skip_seqkit_seq, params.skip_seqkit_replace, params.skip_seqkit_rmdup, params.skip_umitools_extract, params.umi_discard_read, params.skip_adapterremoval, params.adapter_tool, adapters_ch, params.save_merged, params.fastp_discard_trimmed_pass, params.fastp_save_trimmed_fail, params.skip_complexity_filtering, params.complexity_filter_tool, params.skip_deduplication, params.skip_decontamination, decontam_fasta_ch, decontam_ref_ch, params.decontamination_index_name ?: null, params.decontamination_tool, params.skip_final_concatenation)
}
