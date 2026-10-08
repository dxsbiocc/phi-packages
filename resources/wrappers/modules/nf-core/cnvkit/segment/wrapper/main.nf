#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CNVKIT_SEGMENT } from '../main.nf'


params.cnr        = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.cnr, checkIfExists: true)])

    CNVKIT_SEGMENT(in_ch)
}
