#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BEDTOOLS_BAMTOBED } from '../main.nf'


params.bam        = null
params.outdir     = null

workflow {
    bam_ch = Channel.value([[id: 'result'], file(params.bam, checkIfExists: true)])

    BEDTOOLS_BAMTOBED(bam_ch)
}
