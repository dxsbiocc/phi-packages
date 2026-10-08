#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BCFTOOLS_REHEADER } from '../main.nf'


params.vcf        = null
params.fai        = null
params.outdir     = null

workflow {
    vcf_ch = Channel.value([[id: 'result', single_end: false], file(params.vcf, checkIfExists: true), [], []])
    fai_ch = Channel.value([[id: 'genome'], file(params.fai, checkIfExists: true)])

    BCFTOOLS_REHEADER(vcf_ch, fai_ch)
}
