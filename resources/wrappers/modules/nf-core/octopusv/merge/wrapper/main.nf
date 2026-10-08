#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { OCTOPUSV_MERGE } from '../main.nf'


params.svcf1      = null
params.svcf2      = null
params.strategy   = 'intersect'
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], [file(params.svcf1, checkIfExists: true), file(params.svcf2, checkIfExists: true)], params.strategy])

    OCTOPUSV_MERGE(in_ch)
}
