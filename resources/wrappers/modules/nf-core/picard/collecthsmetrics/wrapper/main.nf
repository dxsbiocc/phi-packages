#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PICARD_COLLECTHSMETRICS } from '../main.nf'


params.bam        = null
params.bai        = null
params.baits      = null
params.targets    = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    input_ch = Channel.value([[id: 'test', single_end: false], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), file(params.baits, checkIfExists: true), file(params.targets, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
    fai_ch   = Channel.value([[id: 'genome'], file(params.fai, checkIfExists: true)])
    dict_ch  = Channel.value([[id: 'genome'], file(params.dict, checkIfExists: true)])
    gzi_ch   = Channel.value([[:], []])

    PICARD_COLLECTHSMETRICS(input_ch, fasta_ch, fai_ch, dict_ch, gzi_ch)
}
