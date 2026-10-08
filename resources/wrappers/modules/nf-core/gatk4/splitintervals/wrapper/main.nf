#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_SPLITINTERVALS } from '../main.nf'


params.intervals  = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    intervals_ch = Channel.value([[id: 'test'], file(params.intervals, checkIfExists: true)])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
    fai_ch   = Channel.value([[id: 'genome'], file(params.fai, checkIfExists: true)])
    dict_ch  = Channel.value([[id: 'genome'], file(params.dict, checkIfExists: true)])

    GATK4_SPLITINTERVALS(intervals_ch, fasta_ch, fai_ch, dict_ch)
}
