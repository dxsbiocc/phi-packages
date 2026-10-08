#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { TRUVARI_CONSISTENCY } from '../main.nf'


params.vcf1       = null
params.vcf2       = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], [file(params.vcf1, checkIfExists: true), file(params.vcf2, checkIfExists: true)]])

    TRUVARI_CONSISTENCY(in_ch)
}
