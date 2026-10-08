#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SEQKIT_TRANSLATE } from '../main.nf'

params.fastx  = null
params.outdir = null

workflow {
    fastx_ch = Channel
        .fromPath(params.fastx, checkIfExists: true)
        .map { f -> [[id: "${f.simpleName}_translated"], f] }

    SEQKIT_TRANSLATE(fastx_ch)
}
