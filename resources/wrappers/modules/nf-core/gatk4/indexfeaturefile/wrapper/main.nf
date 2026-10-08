#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_INDEXFEATUREFILE } from '../main.nf'

params.feature_file = null
params.outdir       = null

workflow {
    feature_ch = Channel.value([[id: 'test'], file(params.feature_file, checkIfExists: true)])

    GATK4_INDEXFEATUREFILE(feature_ch)
}
