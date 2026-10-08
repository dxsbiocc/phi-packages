#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf — see
// docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { HISAT2_EXTRACTSPLICESITES } from '../main.nf'

params.gtf    = null
params.outdir = null

workflow {
    Channel
        .fromPath(params.gtf, checkIfExists: true)
        .map { gtf -> [[id: gtf.baseName], gtf] }
        .set { gtf_ch }

    HISAT2_EXTRACTSPLICESITES(gtf_ch)
}
