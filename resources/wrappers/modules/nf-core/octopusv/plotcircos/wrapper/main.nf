#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { OCTOPUSV_PLOTCIRCOS } from '../main.nf'


params.svcf       = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.svcf, checkIfExists: true), 'png'])

    OCTOPUSV_PLOTCIRCOS(in_ch, [])
}
