#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BCFTOOLS_PLUGINSETGT } from '../main.nf'


params.vcf        = null
params.tbi        = null
params.target_gt  = 'a'
params.new_gt     = 'p'
params.outdir     = null

workflow {
    vcf_ch = Channel.value([[id: 'result', single_end: false], file(params.vcf, checkIfExists: true), file(params.tbi, checkIfExists: true)])

    BCFTOOLS_PLUGINSETGT(vcf_ch, params.target_gt, params.new_gt, [], [])
}
