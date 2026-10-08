#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SEQKIT_REPLACE } from '../main.nf'

params.fastx   = null
params.pattern = 'A'
params.replacement = 'N'
params.out_ext = ''
params.outdir  = null

workflow {
    // meta.id must not equal the input's own basename: the module defaults
    // the output name/extension to match the input, and an identical name
    // makes seqkit refuse to run ("input and output files cannot be the same").
    fastx_ch = Channel
        .fromPath(params.fastx, checkIfExists: true)
        .map { fastx -> [[id: "${fastx.simpleName}_replaced"], fastx] }

    SEQKIT_REPLACE(fastx_ch, params.out_ext)
}
