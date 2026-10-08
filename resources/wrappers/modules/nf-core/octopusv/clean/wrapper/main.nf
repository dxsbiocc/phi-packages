#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { OCTOPUSV_CLEAN } from '../main.nf'


params.vcf        = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.vcf, checkIfExists: true), []])

    OCTOPUSV_CLEAN(in_ch)
}
