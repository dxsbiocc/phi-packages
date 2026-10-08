#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CNVKIT_FIX } from '../main.nf'


params.target     = null
params.antitarget = null
params.reference  = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.target, checkIfExists: true), file(params.antitarget, checkIfExists: true), file(params.reference, checkIfExists: true)])

    CNVKIT_FIX(in_ch)
}
