#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BISCUIT_PILEUP } from '../main.nf'
include { BISCUIT_INDEX } from '../../index/main.nf'

params.bam1       = null
params.bai1       = null
params.bam2       = null
params.bai2       = null
params.fasta      = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], [file(params.bam1, checkIfExists: true), file(params.bam2, checkIfExists: true)],
                          [file(params.bai1, checkIfExists: true), file(params.bai2, checkIfExists: true)], [], []])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
    BISCUIT_INDEX(fasta_ch)

    BISCUIT_PILEUP(in_ch, fasta_ch, BISCUIT_INDEX.out.index)
}
