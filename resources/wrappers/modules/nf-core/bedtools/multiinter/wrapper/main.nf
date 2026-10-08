#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BEDTOOLS_MULTIINTER } from '../main.nf'


params.bed        = null
params.bed2       = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'result', single_end: false], [file(params.bed, checkIfExists: true), file(params.bed2, checkIfExists: true)]])

    BEDTOOLS_MULTIINTER(in_ch, [])
}
