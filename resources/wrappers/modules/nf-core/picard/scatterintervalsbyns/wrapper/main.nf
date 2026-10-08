#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PICARD_SCATTERINTERVALSBYNS } from '../main.nf'


params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    fasta_ch = Channel.value([[id: 'test', single_end: false], file(params.fasta, checkIfExists: true)])
    fai_ch   = Channel.value([[id: 'test', single_end: false], file(params.fai, checkIfExists: true)])
    dict_ch  = Channel.value([[id: 'test', single_end: false], file(params.dict, checkIfExists: true)])

    PICARD_SCATTERINTERVALSBYNS(fasta_ch, fai_ch, dict_ch)
}
