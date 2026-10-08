#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CUSTOM_RESOLVETAXONOMY } from '../main.nf'


params.sequences  = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result'], [], file(params.sequences, checkIfExists: true), true])

    CUSTOM_RESOLVETAXONOMY(in_ch)
}
