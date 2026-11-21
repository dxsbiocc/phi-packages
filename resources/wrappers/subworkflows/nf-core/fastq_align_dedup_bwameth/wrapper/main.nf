#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTQ_ALIGN_DEDUP_BWAMETH } from '../main.nf'

params.reads = null
params.fasta = null
params.fai = null
params.bwameth_index = null
params.skip_deduplication = false
params.use_gpu = false
params.outdir = "results"

workflow {
    reads_ch = Channel
        .fromFilePairs(params.reads, size: -1)
        .map { id, reads -> [[id: id, single_end: reads.size() == 1], reads] }
    fasta_fai_ch = Channel.value([[id: 'reference'], [file(params.fasta, checkIfExists: true)], [file(params.fai, checkIfExists: true)]])
    index_ch = Channel.value([[id: 'reference'], [file(params.bwameth_index, checkIfExists: true)]])

    FASTQ_ALIGN_DEDUP_BWAMETH(reads_ch, fasta_fai_ch, index_ch, params.skip_deduplication, params.use_gpu)
}
