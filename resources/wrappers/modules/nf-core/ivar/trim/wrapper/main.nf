#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { IVAR_TRIM } from '../main.nf'


params.bam        = null
params.bai        = null
params.bed        = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])

    IVAR_TRIM(in_ch, file(params.bed, checkIfExists: true))
}
