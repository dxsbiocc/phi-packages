#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BEDTOOLS_UNIONBEDG } from '../main.nf'


params.bedgraph   = null
params.bed        = null
params.sizes      = null
params.outdir     = null

workflow {
    in_ch    = Channel.value([[id: 'result'], [file(params.bedgraph, checkIfExists: true), file(params.bed, checkIfExists: true)]])
    sizes_ch = Channel.value([[id: 'sizes'], file(params.sizes, checkIfExists: true)])

    BEDTOOLS_UNIONBEDG(in_ch, sizes_ch)
}
