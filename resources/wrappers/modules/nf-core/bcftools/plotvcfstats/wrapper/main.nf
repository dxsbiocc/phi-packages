#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BCFTOOLS_PLOTVCFSTATS } from '../main.nf'
include { BCFTOOLS_STATS } from '../../stats/main.nf'

params.vcf        = null
params.tbi        = null
params.outdir     = null

workflow {
    vcf_ch = Channel.value([[id: 'result', single_end: false], file(params.vcf, checkIfExists: true), file(params.tbi, checkIfExists: true)])
    empty_ch = Channel.value([[:], []])

    BCFTOOLS_STATS(vcf_ch, empty_ch, empty_ch, empty_ch, empty_ch, empty_ch)
    BCFTOOLS_PLOTVCFSTATS(BCFTOOLS_STATS.out.stats)
}
