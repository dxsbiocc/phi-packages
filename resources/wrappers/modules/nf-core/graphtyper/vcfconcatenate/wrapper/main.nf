#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GRAPHTYPER_VCFCONCATENATE } from '../main.nf'


params.vcf1       = null
params.vcf2       = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'test_all'], [file(params.vcf1, checkIfExists: true), file(params.vcf2, checkIfExists: true)]])

    GRAPHTYPER_VCFCONCATENATE(in_ch)
}
