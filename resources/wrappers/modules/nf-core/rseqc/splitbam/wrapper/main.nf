#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { RSEQC_SPLITBAM } from '../main.nf'

params.bam    = null
params.bai    = null
params.bed    = null
params.outdir = null

workflow {
    bam_ch = Channel.value([[id: 'test'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])
    bed_ch = Channel.value([[id: 'test_bed'], file(params.bed, checkIfExists: true)])

    RSEQC_SPLITBAM(bam_ch, bed_ch)
}
