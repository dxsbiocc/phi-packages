#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PICARD_POSITIONBASEDDOWNSAMPLESAM } from '../main.nf'


params.bam        = null
params.fraction   = '0.5'
params.outdir     = null

workflow {
    bam_ch = Channel.value([[id: 'test'], file(params.bam, checkIfExists: true), params.fraction])

    PICARD_POSITIONBASEDDOWNSAMPLESAM(bam_ch)
}
