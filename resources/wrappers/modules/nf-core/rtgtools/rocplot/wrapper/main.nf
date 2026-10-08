#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { RTGTOOLS_ROCPLOT } from '../main.nf'
include { UNTAR } from '../../../untar/main.nf'
include { RTGTOOLS_VCFEVAL } from '../../vcfeval/main.nf'

params.sdf        = null
params.query      = null
params.query_tbi  = null
params.truth      = null
params.truth_tbi  = null
params.truth_bed  = null
params.regions_bed = null
params.outdir     = null

workflow {
    UNTAR(Channel.value([[id: 'test'], file(params.sdf, checkIfExists: true)]))
    query_ch = Channel.value([[id: 'sample'], file(params.query, checkIfExists: true), file(params.query_tbi, checkIfExists: true),
                              file(params.truth, checkIfExists: true), file(params.truth_tbi, checkIfExists: true),
                              file(params.truth_bed, checkIfExists: true), file(params.regions_bed, checkIfExists: true)])

    RTGTOOLS_VCFEVAL(query_ch, UNTAR.out.untar)
    RTGTOOLS_ROCPLOT(RTGTOOLS_VCFEVAL.out.weighted_roc)
}
