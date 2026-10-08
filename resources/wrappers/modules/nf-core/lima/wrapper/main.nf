#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { LIMA } from '../main.nf'


params.ccs        = null
params.primers    = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.ccs, checkIfExists: true)])

    LIMA(in_ch, [file(params.primers, checkIfExists: true)])
}
