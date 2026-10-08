#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_ANALYZECOVARIATES } from '../main.nf'


params.before     = null
params.after      = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'test'], file(params.before, checkIfExists: true), file(params.after, checkIfExists: true), []])

    GATK4_ANALYZECOVARIATES(in_ch)
}
