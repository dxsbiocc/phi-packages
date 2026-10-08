#!/usr/bin/env nextflow
// Thin agent-facing adapter over ../main.nf (BOWTIE_ALIGN). A pre-built Bowtie
// index (`index`) is used as-is; without one, ../../build/main.nf builds it
// from `fasta` first.
nextflow.enable.dsl = 2

include { BOWTIE_ALIGN } from '../main.nf'
include { BOWTIE_BUILD } from '../../build/main.nf'

params.reads      = null
params.index      = null
params.fasta      = null
params.outdir     = null

workflow {
    if (!params.index && !params.fasta) {
        error "bowtie-align needs either `index` (a pre-built Bowtie index directory) or `fasta` to build one."
    }

    def reads = file(params.reads, checkIfExists: true)
    reads_ch = Channel.value([[id: reads.simpleName, single_end: true], reads])

    if (params.index) {
        index_ch = Channel.value([[id: 'bowtie_index'], file(params.index, checkIfExists: true)])
    } else {
        BOWTIE_BUILD(Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)]))
        index_ch = BOWTIE_BUILD.out.index
    }

    BOWTIE_ALIGN(reads_ch, index_ch, true)
}
