#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BCFTOOLS_ANNOTATE } from '../main.nf'

params.vcf             = null
params.vcf_tbi         = null
params.annotations     = null
params.annotations_tbi = null
params.outdir          = null

workflow {
    input_ch = Channel.value([
        [id: 'test_annotated', single_end: false],
        file(params.vcf, checkIfExists: true),
        file(params.vcf_tbi, checkIfExists: true),
        file(params.annotations, checkIfExists: true),
        file(params.annotations_tbi, checkIfExists: true),
        [], [], []
    ])

    BCFTOOLS_ANNOTATE(input_ch)
}
