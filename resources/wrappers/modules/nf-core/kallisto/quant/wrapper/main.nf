#!/usr/bin/env nextflow
// Thin agent-facing adapter over ../main.nf (KALLISTO_QUANT). A pre-built
// kallisto index file (`index`) is used as-is; without one, ../../index/main.nf
// builds it from `fasta` (transcript sequences) first.
nextflow.enable.dsl = 2

include { KALLISTO_INDEX } from '../../index/main.nf'
include { KALLISTO_QUANT } from '../main.nf'

params.reads_1 = null
params.reads_2 = null
params.index   = null
params.fasta   = null
params.outdir  = null

workflow {
    if (!params.index && !params.fasta) {
        error "kallisto-quant needs either `index` (a pre-built kallisto .idx file) or `fasta` to build one."
    }

    if (params.index) {
        ch_index = Channel.value([[id: 'kallisto_index'], file(params.index, checkIfExists: true)])
    } else {
        ch_fasta = Channel
            .fromPath(params.fasta, checkIfExists: true)
            .map { fasta -> [[id: fasta.baseName], fasta] }
        KALLISTO_INDEX(ch_fasta)
        ch_index = KALLISTO_INDEX.out.index.first()
    }

    def read1 = file(params.reads_1, checkIfExists: true)
    reads_ch = channel.of([
        [id: read1.simpleName.replaceAll(/[._-]R?1$/, ''), single_end: false],
        [read1, file(params.reads_2, checkIfExists: true)]
    ])

    KALLISTO_QUANT(reads_ch, ch_index, [], [], '', '')
}
