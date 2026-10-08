#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CUSTOM_ORFNORMALISE } from '../main.nf'


params.orfs       = null
params.gtf        = null
params.caller     = 'ribocode'
params.outdir     = null

workflow {
    orfs_ch = Channel.value([[id: 'result'], file(params.orfs, checkIfExists: true), params.caller])
    gtf_ch  = Channel.value([[id: 'reference'], file(params.gtf, checkIfExists: true)])

    CUSTOM_ORFNORMALISE(orfs_ch, gtf_ch)
}
