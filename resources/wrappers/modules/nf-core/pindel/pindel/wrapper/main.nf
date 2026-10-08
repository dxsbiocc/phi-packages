#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PINDEL_PINDEL } from '../main.nf'


params.bam        = null
params.bai        = null
params.fasta      = null
params.fai        = null
params.bed        = null
params.outdir     = null

workflow {
    bam_ch = Channel.value([[id: 'sample', single_end: false], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])
    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))
    fai_ch   = Channel.value(file(params.fai, checkIfExists: true))
    bed_ch   = Channel.value(file(params.bed, checkIfExists: true))

    PINDEL_PINDEL(bam_ch, fasta_ch, fai_ch, bed_ch)
}
