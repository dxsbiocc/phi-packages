#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BWAMETH_ALIGN } from '../main.nf'
include { BWAMETH_INDEX } from '../../index/main.nf'

params.reads      = null
params.fasta      = null
params.outdir     = null

workflow {
    reads_ch = Channel.value([[id: 'sample', single_end: true], file(params.reads, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])

    BWAMETH_INDEX(fasta_ch, false)
    BWAMETH_ALIGN(reads_ch, fasta_ch, BWAMETH_INDEX.out.index)
}
