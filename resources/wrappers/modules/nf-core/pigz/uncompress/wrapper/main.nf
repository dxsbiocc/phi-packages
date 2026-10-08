#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PIGZ_UNCOMPRESS } from '../main.nf'

params.input  = null
params.outdir = null

workflow {
    input_ch = Channel
        .fromPath(params.input, checkIfExists: true)
        .map { f -> [[id: f.simpleName], f] }

    PIGZ_UNCOMPRESS(input_ch)
}
