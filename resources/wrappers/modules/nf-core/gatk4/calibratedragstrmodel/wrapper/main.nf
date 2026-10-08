#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_CALIBRATEDRAGSTRMODEL } from '../main.nf'


params.bam        = null
params.bai        = null
params.strtable   = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    bam_ch = Channel.value([[id: 'test', single_end: false], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])
    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))
    fai_ch   = Channel.value(file(params.fai, checkIfExists: true))
    dict_ch  = Channel.value(file(params.dict, checkIfExists: true))
    str_ch = Channel.value(file(params.strtable, checkIfExists: true))

    GATK4_CALIBRATEDRAGSTRMODEL(bam_ch, fasta_ch, fai_ch, dict_ch, str_ch)
}
