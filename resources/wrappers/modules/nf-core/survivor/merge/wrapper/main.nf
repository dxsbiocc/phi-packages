#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SURVIVOR_MERGE } from '../main.nf'


params.vcf1       = null
params.vcf2       = null
params.max_distance = '0.2'
params.min_callers = 1
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample', single_end: false], [file(params.vcf1, checkIfExists: true), file(params.vcf2, checkIfExists: true)]])

    SURVIVOR_MERGE(in_ch, params.max_distance, params.min_callers, 0, 0, 0, 0)
}
