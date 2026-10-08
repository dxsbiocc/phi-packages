#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BCFTOOLS_INDEX } from '../main.nf'

params.vcf    = null
params.outdir = null

workflow {
    vcf_ch = Channel.value([[id: 'test'], file(params.vcf, checkIfExists: true)])

    BCFTOOLS_INDEX(vcf_ch)
}
