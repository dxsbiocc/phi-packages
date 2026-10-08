#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BCFTOOLS_VIEW } from '../main.nf'

params.vcf     = null
params.tbi     = null
params.outdir  = null

workflow {
    // meta.id must not equal the input VCF's own basename — see the same
    // fix in bcftools/norm's wrapper for why.
    vcf_ch = Channel.value([[id: 'test_view'], file(params.vcf, checkIfExists: true), file(params.tbi, checkIfExists: true)])

    empty_ch = Channel.value([])

    BCFTOOLS_VIEW(vcf_ch, empty_ch, empty_ch, empty_ch)
}
