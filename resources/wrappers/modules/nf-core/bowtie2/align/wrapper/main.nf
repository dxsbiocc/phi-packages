#!/usr/bin/env nextflow
// Thin agent-facing adapter over ../main.nf (BOWTIE2_ALIGN). A pre-built
// Bowtie2 index (`index`) is used as-is; without one, ../../build/main.nf
// builds it from `fasta` first.
nextflow.enable.dsl = 2

include { BOWTIE2_BUILD } from '../../build/main.nf'
include { BOWTIE2_ALIGN } from '../main.nf'

params.reads_1 = null
params.reads_2 = null
params.index   = null
params.fasta   = null
params.outdir  = null

workflow {
    if (!params.index && !params.fasta) {
        error "bowtie2-align needs either `index` (a pre-built Bowtie2 index directory) or `fasta` to build one."
    }

    if (params.index) {
        ch_index = Channel.value([[id: 'bowtie2_index'], file(params.index, checkIfExists: true)])
    } else {
        ch_fasta = Channel
            .fromPath(params.fasta, checkIfExists: true)
            .map { fasta -> [[id: fasta.baseName], fasta] }
        BOWTIE2_BUILD(ch_fasta)
        ch_index = BOWTIE2_BUILD.out.index.first()
    }

    def read1 = file(params.reads_1, checkIfExists: true)
    reads_ch = channel.of([
        [id: read1.simpleName.replaceAll(/[._-]R?1$/, ''), single_end: false],
        [read1, file(params.reads_2, checkIfExists: true)]
    ])

    ch_empty_fasta = channel.value([[:], []])

    BOWTIE2_ALIGN(reads_ch, ch_index, ch_empty_fasta, false, true)
}
