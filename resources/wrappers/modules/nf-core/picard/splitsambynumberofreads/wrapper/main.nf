#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PICARD_SPLITSAMBYNUMBEROFREADS } from '../main.nf'


params.bam        = null
params.reads_per_file = 50
params.outdir     = null

workflow {
    bam_ch = Channel.value([[id: 'test', single_end: false], file(params.bam, checkIfExists: true)])
    ref_ch = Channel.value([[], [], []])

    PICARD_SPLITSAMBYNUMBEROFREADS(bam_ch, ref_ch, params.reads_per_file, [], [])
}
