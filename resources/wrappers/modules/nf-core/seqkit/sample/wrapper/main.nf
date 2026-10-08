#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SEQKIT_SAMPLE } from '../main.nf'

params.fastx      = null
params.proportion = '0.5'
params.outdir     = null

workflow {
    fastx_ch = Channel
        .fromPath(params.fastx, checkIfExists: true)
        .map { f -> [[id: "${f.simpleName}_sampled"], f] }

    SEQKIT_SAMPLE(fastx_ch)
}
