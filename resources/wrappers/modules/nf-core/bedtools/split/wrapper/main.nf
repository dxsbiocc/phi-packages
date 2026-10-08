#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BEDTOOLS_SPLIT } from '../main.nf'


params.bed        = null
params.count      = 2
params.outdir     = null

workflow {
    bed_ch = Channel.value([[id: 'result'], file(params.bed, checkIfExists: true), params.count])

    BEDTOOLS_SPLIT(bed_ch)
}
