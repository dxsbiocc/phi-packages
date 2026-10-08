#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PLINK_LD } from '../main.nf'


params.vcf        = null
params.snpfile    = null
params.outdir     = null

workflow {
    vcf_ch = Channel.value([[id: 'sample'], file(params.vcf, checkIfExists: true)])
    snp_ch = Channel.value([[id: 'sample'], file(params.snpfile, checkIfExists: true)])

    PLINK_LD([[id: 'none'], [], [], []], vcf_ch, [[id: 'none'], []], snp_ch)
}
