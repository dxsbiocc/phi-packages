#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SOMALIER_RELATE } from '../main.nf'


params.normal     = null
params.tumour     = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'cohort', single_end: false], [file(params.normal, checkIfExists: true), file(params.tumour, checkIfExists: true)], []])

    SOMALIER_RELATE(in_ch, [])
}
