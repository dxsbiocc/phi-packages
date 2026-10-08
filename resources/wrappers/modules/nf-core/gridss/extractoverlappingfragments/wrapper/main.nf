#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GRIDSS_EXTRACTOVERLAPPINGFRAGMENTS } from '../main.nf'


params.bam        = null
params.bai        = null
params.bed        = null
params.outdir     = null

workflow {
    bam_ch = Channel.of([[id: 'sample'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])
    bed_ch = Channel.value([[id: 'targets'], file(params.bed, checkIfExists: true)])

    GRIDSS_EXTRACTOVERLAPPINGFRAGMENTS(bam_ch, bed_ch)
}
