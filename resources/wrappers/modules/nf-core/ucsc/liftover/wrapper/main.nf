#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { UCSC_LIFTOVER } from '../main.nf'

params.bed    = null
params.chain  = null
params.outdir = null

workflow {
    bed_ch   = Channel.value([[id: 'test'], file(params.bed, checkIfExists: true)])
    chain_ch = Channel.value(file(params.chain, checkIfExists: true))

    UCSC_LIFTOVER(bed_ch, chain_ch)
}
