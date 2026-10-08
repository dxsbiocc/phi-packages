#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BCFTOOLS_ROHVIZ } from '../main.nf'
include { BCFTOOLS_ROH } from '../../roh/main.nf'

params.vcf        = null
params.tbi        = null
params.outdir     = null

workflow {
    vcf_ch = Channel.value([[id: 'result', single_end: false], file(params.vcf, checkIfExists: true), file(params.tbi, checkIfExists: true)])
    vcf_only = Channel.value([[id: 'result'], file(params.vcf, checkIfExists: true)])

    BCFTOOLS_ROH(vcf_ch, [[], []], [], [], [], [])
    BCFTOOLS_ROHVIZ(BCFTOOLS_ROH.out.roh, vcf_only, [], [])
}
