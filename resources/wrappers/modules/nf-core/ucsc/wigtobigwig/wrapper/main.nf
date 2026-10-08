#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { UCSC_WIGTOBIGWIG } from '../main.nf'

params.wig    = null
params.sizes  = null
params.outdir = null

workflow {
    wig_ch   = Channel.value([[id: 'test', single_end: false], file(params.wig, checkIfExists: true)])
    sizes_ch = Channel.value(file(params.sizes, checkIfExists: true))

    UCSC_WIGTOBIGWIG(wig_ch, sizes_ch)
}
