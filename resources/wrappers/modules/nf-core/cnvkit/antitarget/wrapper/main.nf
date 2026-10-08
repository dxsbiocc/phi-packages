#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CNVKIT_ANTITARGET } from '../main.nf'


params.targets    = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.targets, checkIfExists: true)])

    CNVKIT_ANTITARGET(in_ch)
}
