#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CUSTOM_CLUSTERVISUALIZATION } from '../main.nf'


params.features   = null
params.clusters   = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result'], file(params.features, checkIfExists: true), file(params.clusters, checkIfExists: true)])

    CUSTOM_CLUSTERVISUALIZATION(in_ch)
}
