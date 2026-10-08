#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BCFTOOLS_STATS } from '../main.nf'

params.vcf    = null
params.outdir = null

workflow {
    vcf_ch = Channel
        .fromPath(params.vcf, checkIfExists: true)
        .map { vcf -> [[id: vcf.simpleName], vcf, []] }

    empty_ch = Channel.value([[:], []])

    BCFTOOLS_STATS(vcf_ch, empty_ch, empty_ch, empty_ch, empty_ch, empty_ch)
}
