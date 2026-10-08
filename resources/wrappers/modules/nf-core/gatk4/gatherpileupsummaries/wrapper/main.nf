#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_GATHERPILEUPSUMMARIES } from '../main.nf'


params.pileup     = null
params.dict       = null
params.outdir     = null

workflow {
    pileup_ch = Channel.value([[id: 'gathered', single_end: false], [file(params.pileup, checkIfExists: true)]])
    dict_ch   = Channel.value(file(params.dict, checkIfExists: true))

    GATK4_GATHERPILEUPSUMMARIES(pileup_ch, dict_ch)
}
