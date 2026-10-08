#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BCFTOOLS_CONCAT } from '../main.nf'

params.vcf1     = null
params.vcf2     = null
params.vcf1_csi = null
params.vcf2_csi = null
params.outdir   = null

workflow {
    input_ch = Channel.value([
        [id: 'concat'],
        [file(params.vcf1, checkIfExists: true), file(params.vcf2, checkIfExists: true)],
        [file(params.vcf1_csi, checkIfExists: true), file(params.vcf2_csi, checkIfExists: true)]
    ])

    BCFTOOLS_CONCAT(input_ch)
}
