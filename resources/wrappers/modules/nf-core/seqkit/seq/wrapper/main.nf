#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SEQKIT_SEQ } from '../main.nf'

params.fastx  = null
params.outdir = null

workflow {
    fastx_ch = Channel
        .fromPath(params.fastx, checkIfExists: true)
        .map { f -> [[id: "${f.simpleName}_seq", single_end: false], f] }

    SEQKIT_SEQ(fastx_ch)
}
