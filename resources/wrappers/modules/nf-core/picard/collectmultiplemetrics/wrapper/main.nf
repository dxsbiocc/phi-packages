#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PICARD_COLLECTMULTIPLEMETRICS } from '../main.nf'


params.bam        = null
params.bai        = null
params.fasta      = null
params.outdir     = null

workflow {
    bam_ch   = Channel.value([[id: 'test', single_end: false], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
    fai_ch   = Channel.value([[id: 'genome'], []])

    PICARD_COLLECTMULTIPLEMETRICS(bam_ch, fasta_ch, fai_ch)
}
