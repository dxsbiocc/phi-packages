#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTQ_ALIGN_DEDUP_BISMARK } from '../main.nf'

params.reads = null
params.fasta = null
params.fai = null
params.bismark_index = null
params.skip_deduplication = false
params.cytosine_report = false
params.outdir = "results"

workflow {
    reads_ch = Channel
        .fromFilePairs(params.reads, size: -1)
        .map { id, reads -> [[id: id, single_end: reads.size() == 1], reads] }
    fasta_fai_ch = Channel.value([[id: 'reference'], [file(params.fasta, checkIfExists: true)], [file(params.fai, checkIfExists: true)]])
    index_ch = Channel.value([[id: 'reference'], [file(params.bismark_index, checkIfExists: true)]])

    FASTQ_ALIGN_DEDUP_BISMARK(reads_ch, fasta_fai_ch, index_ch, params.skip_deduplication, params.cytosine_report)
}
