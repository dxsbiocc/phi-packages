#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { FQ_GENERATE } from '../main.nf'

params.outdir = null

workflow {
    meta_ch = Channel.value([id: 'synthetic'])

    FQ_GENERATE(meta_ch)
}
