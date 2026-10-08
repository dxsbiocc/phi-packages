#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BEDTOOLS_MERGE } from '../main.nf'

params.bed    = null
params.outdir = null

workflow {
    bed_ch = Channel.value([[id: 'merged'], file(params.bed, checkIfExists: true)])

    BEDTOOLS_MERGE(bed_ch)
}
