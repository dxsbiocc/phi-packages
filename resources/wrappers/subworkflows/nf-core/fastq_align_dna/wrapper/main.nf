#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { FASTQ_ALIGN_DNA } from '../main.nf'

params.reads = null
params.index = null
params.fasta = null
params.aligner = "bwamem2"
params.sort = true
params.outdir = "results"

workflow {
    reads_ch = Channel
        .fromFilePairs(params.reads, size: -1)
        .map { id, reads -> [[id: id, single_end: reads.size() == 1], reads] }
    index_ch = Channel.value([[id: 'reference'], [file(params.index, checkIfExists: true)]])
    fasta_ch = Channel.value([[id: 'reference'], file(params.fasta, checkIfExists: true)])

    FASTQ_ALIGN_DNA(reads_ch, index_ch, fasta_ch, params.aligner, params.sort)
}
