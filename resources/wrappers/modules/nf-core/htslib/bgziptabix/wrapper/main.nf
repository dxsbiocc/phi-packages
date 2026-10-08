#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { HTSLIB_BGZIPTABIX } from '../main.nf'

params.input      = null
params.action     = 'compress'
params.make_index = true
params.out_ext    = 'vcf'
params.outdir     = null

workflow {
    input_ch = Channel
        .fromPath(params.input, checkIfExists: true)
        .map { f -> [[id: f.simpleName], f, [], []] }

    HTSLIB_BGZIPTABIX(input_ch, params.action, params.make_index, params.out_ext)
}
