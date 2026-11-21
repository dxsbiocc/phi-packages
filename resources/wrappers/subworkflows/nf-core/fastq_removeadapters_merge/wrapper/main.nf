#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTQ_REMOVEADAPTERS_MERGE } from '../main.nf'

params.reads = null
params.adapter_tool = "fastp"
params.custom_adapters_file = null
params.save_merged = false
params.fastp_discard_trimmed_pass = false
params.fastp_save_trimmed_fail = false
params.outdir = "results"

workflow {
    reads_ch = Channel
        .fromFilePairs(params.reads, size: -1)
        .map { id, reads -> [[id: id, single_end: reads.size() == 1], reads] }
    adapters_ch = params.custom_adapters_file ? Channel.value(file(params.custom_adapters_file, checkIfExists: true)) : Channel.empty()

    FASTQ_REMOVEADAPTERS_MERGE(reads_ch, params.adapter_tool, adapters_ch, params.save_merged, params.fastp_discard_trimmed_pass, params.fastp_save_trimmed_fail)
}
