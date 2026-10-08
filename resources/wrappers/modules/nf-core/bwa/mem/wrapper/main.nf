#!/usr/bin/env nextflow
// Thin agent-facing adapter over ../main.nf (BWA_MEM). A pre-built BWA index
// (`index`) is used as-is; without one, ../../index/main.nf builds it from
// `fasta` first.
nextflow.enable.dsl = 2

include { BWA_MEM } from '../main.nf'
include { BWA_INDEX } from '../../index/main.nf'

params.reads1     = null
params.reads2     = null
params.index      = null
params.fasta      = null
params.outdir     = null

workflow {
    if (!params.index && !params.fasta) {
        error "bwa-mem needs either `index` (a pre-built BWA index directory) or `fasta` to build one."
    }

    def read1 = file(params.reads1, checkIfExists: true)
    reads_ch = Channel.value([
        [id: read1.simpleName.replaceAll(/[._-]R?1$/, ''), single_end: false],
        [read1, file(params.reads2, checkIfExists: true)]
    ])
    fasta_ch = params.fasta
        ? Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
        : Channel.value([[:], []])

    if (params.index) {
        index_ch = Channel.value([[id: 'bwa_index'], file(params.index, checkIfExists: true)])
    } else {
        BWA_INDEX(fasta_ch)
        index_ch = BWA_INDEX.out.index
    }

    BWA_MEM(reads_ch, index_ch, fasta_ch, false)
}
