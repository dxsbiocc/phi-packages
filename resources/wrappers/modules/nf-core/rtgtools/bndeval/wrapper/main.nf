#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { RTGTOOLS_BNDEVAL } from '../main.nf'
include { RTGTOOLS_SVDECOMPOSE } from '../../svdecompose/main.nf'

params.sv         = null
params.sv_tbi     = null
params.regions_bed = null
params.outdir     = null

workflow {
    truth_in = Channel.value([[id: 'sample'], file(params.sv, checkIfExists: true), file(params.sv_tbi, checkIfExists: true)])
    RTGTOOLS_SVDECOMPOSE(truth_in)
    query = Channel.value([[id: 'sample'], file(params.sv, checkIfExists: true), file(params.sv_tbi, checkIfExists: true)])
    regions = Channel.value([[id: 'sample'], file(params.regions_bed, checkIfExists: true)])

    RTGTOOLS_BNDEVAL(RTGTOOLS_SVDECOMPOSE.out.vcf.join(RTGTOOLS_SVDECOMPOSE.out.index).join(query).join(regions))
}
