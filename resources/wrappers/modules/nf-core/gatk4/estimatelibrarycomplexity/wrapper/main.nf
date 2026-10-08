#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_ESTIMATELIBRARYCOMPLEXITY } from '../main.nf'


params.bam        = null
params.dict       = null
params.outdir     = null

workflow {
    bam_ch   = Channel.value([[id: 'test', single_end: false], file(params.bam, checkIfExists: true)])
    empty_ch = Channel.value([])
    dict_ch  = Channel.value(file(params.dict, checkIfExists: true))

    GATK4_ESTIMATELIBRARYCOMPLEXITY(bam_ch, empty_ch, empty_ch, dict_ch)
}
