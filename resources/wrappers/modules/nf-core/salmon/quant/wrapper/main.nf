#!/usr/bin/env nextflow
// Thin agent-facing adapter over ../main.nf (SALMON_QUANT). A pre-built salmon
// index (`index`) is used as-is; without one, ../../index/main.nf builds it
// from `transcript_fasta` first. SALMON_QUANT only reads the transcript FASTA
// in alignment mode, which this wrapper never uses, so it is optional once an
// index is given.
nextflow.enable.dsl = 2

include { SALMON_INDEX } from '../../index/main.nf'
include { SALMON_QUANT } from '../main.nf'

params.reads_1          = null
params.reads_2          = null
params.index            = null
params.transcript_fasta = null
params.gtf              = null
params.outdir           = null

workflow {
    if (!params.index && !params.transcript_fasta) {
        error "salmon-quant needs either `index` (a pre-built salmon index directory) or `transcript_fasta` to build one."
    }

    ch_transcript_fasta = params.transcript_fasta ? file(params.transcript_fasta, checkIfExists: true) : []

    if (params.index) {
        ch_index = Channel.value(file(params.index, checkIfExists: true))
    } else {
        SALMON_INDEX([], ch_transcript_fasta)
        ch_index = SALMON_INDEX.out.index
    }

    def read1 = file(params.reads_1, checkIfExists: true)
    reads_ch = channel.of([
        [id: read1.simpleName.replaceAll(/[._-]R?1$/, ''), single_end: false],
        [read1, file(params.reads_2, checkIfExists: true)]
    ])

    SALMON_QUANT(
        reads_ch,
        ch_index,
        file(params.gtf, checkIfExists: true),
        ch_transcript_fasta,
        false,
        ''
    )
}
