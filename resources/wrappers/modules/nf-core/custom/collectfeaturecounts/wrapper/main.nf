#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CUSTOM_COLLECTFEATURECOUNTS } from '../main.nf'


params.counts1    = null
params.counts2    = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result'], [file(params.counts1, checkIfExists: true), file(params.counts2, checkIfExists: true)]])

    CUSTOM_COLLECTFEATURECOUNTS(in_ch)
}
