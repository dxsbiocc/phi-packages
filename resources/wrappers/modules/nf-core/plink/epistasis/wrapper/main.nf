#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PLINK_EPISTASIS } from '../main.nf'


params.vcf        = null
params.phe        = null
params.outdir     = null

workflow {
    vcf_ch = Channel.value([[id: 'sample'], file(params.vcf, checkIfExists: true)])
    phe_ch = Channel.value([[id: 'sample'], file(params.phe, checkIfExists: true)])

    PLINK_EPISTASIS([[id: 'none'], [], [], []], vcf_ch, [[id: 'none'], []], phe_ch)
}
