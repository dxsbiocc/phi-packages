#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SAMTOOLS_SORMADUP } from '../main.nf'

params.bam    = null
params.fasta  = null
params.fai    = null
params.outdir = null

workflow {
    bam_ch   = Channel.value([[id: 'test'], file(params.bam, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'fasta'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true)])

    SAMTOOLS_SORMADUP(bam_ch, fasta_ch)
}
