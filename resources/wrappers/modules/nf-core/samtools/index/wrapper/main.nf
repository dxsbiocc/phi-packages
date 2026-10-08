#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf — see
// docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SAMTOOLS_INDEX } from '../main.nf'

params.input  = null
params.outdir = null

workflow {
    Channel
        .fromPath(params.input, checkIfExists: true)
        .map { file -> [[id: file.baseName], file] }
        .set { input_ch }

    SAMTOOLS_INDEX(input_ch)
}
