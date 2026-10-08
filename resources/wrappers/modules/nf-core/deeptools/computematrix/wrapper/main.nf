#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// Composes the vendored deeptools/bamcoverage module
// (../../bamcoverage/main.nf) so the agent can compute a signal matrix
// from a raw BAM in one step, without a separately-built bigWig. See
// docs/design/phi-wrapper-agent-composition-design.md sections 1 and 3.
nextflow.enable.dsl = 2

include { DEEPTOOLS_BAMCOVERAGE }  from '../../bamcoverage/main.nf'
include { DEEPTOOLS_COMPUTEMATRIX } from '../main.nf'

params.bam    = null
params.bai    = null
params.bed    = null
params.outdir = null

workflow {
    bam_ch = Channel.value([[id: 'test', single_end: false], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])
    empty_ch = Channel.value([])
    blacklist_ch = Channel.value([[id: 'no_blacklist'], []])

    DEEPTOOLS_BAMCOVERAGE(bam_ch, empty_ch, empty_ch, blacklist_ch)

    bigwig_ch = DEEPTOOLS_BAMCOVERAGE.out.bigwig
    bed_ch    = Channel.value(file(params.bed, checkIfExists: true))

    DEEPTOOLS_COMPUTEMATRIX(bigwig_ch, bed_ch)
}
