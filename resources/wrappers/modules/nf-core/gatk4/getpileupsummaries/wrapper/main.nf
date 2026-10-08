#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_GETPILEUPSUMMARIES } from '../main.nf'


params.bam        = null
params.bai        = null
params.intervals  = null
params.variants   = null
params.variants_tbi = null
params.outdir     = null

workflow {
    bam_ch   = Channel.value([[id: 'test'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), file(params.intervals, checkIfExists: true)])
    ref_ch   = Channel.value([[], []])
    var_ch   = Channel.value(file(params.variants, checkIfExists: true))
    varidx_ch = Channel.value(file(params.variants_tbi, checkIfExists: true))

    GATK4_GETPILEUPSUMMARIES(bam_ch, ref_ch, ref_ch, ref_ch, var_ch, varidx_ch)
}
