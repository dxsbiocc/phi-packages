#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_MERGEVCFS } from '../main.nf'

params.vcf1   = null
params.vcf2   = null
params.dict   = null
params.outdir = null

workflow {
    vcf_ch  = Channel.value([[id: 'merged'], [file(params.vcf1, checkIfExists: true), file(params.vcf2, checkIfExists: true)]])
    dict_ch = Channel.value([[id: 'genome'], file(params.dict, checkIfExists: true)])

    GATK4_MERGEVCFS(vcf_ch, dict_ch)
}
