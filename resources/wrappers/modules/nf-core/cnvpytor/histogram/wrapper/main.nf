#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CNVPYTOR_HISTOGRAM } from '../main.nf'


params.pytor      = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.pytor, checkIfExists: true)])

    CNVPYTOR_HISTOGRAM(in_ch, [])
}
