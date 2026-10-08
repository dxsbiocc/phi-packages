#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PICARD_COLLECTALIGNMENTSUMMARYMETRICS } from '../main.nf'


params.bam        = null
params.fasta      = null
params.outdir     = null

workflow {
    bam_ch   = Channel.value([[id: 'test', single_end: false], file(params.bam, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'genome', single_end: false], [file(params.fasta, checkIfExists: true)]])

    PICARD_COLLECTALIGNMENTSUMMARYMETRICS(bam_ch, fasta_ch)
}
