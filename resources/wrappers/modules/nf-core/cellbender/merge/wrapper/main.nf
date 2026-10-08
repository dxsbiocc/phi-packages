#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// cellbender_h5 defaults to a pre-computed CellBender output from
// nf-core's own test-datasets rather than composing with
// cellbender/removebackground: that module's smoke run takes long enough
// (a training loop, even at --epochs 5) that chaining it here would make
// every merge smoke test pay that cost too, for a merge step that only
// reads the finished .h5 file.
nextflow.enable.dsl = 2

include { CELLBENDER_MERGE } from '../main.nf'

params.filtered           = null
params.unfiltered         = null
params.cellbender_h5      = null
params.output_layer_name  = ''
params.outdir             = null

workflow {
    merge_ch = Channel.value([
        [id: 'test'],
        file(params.filtered, checkIfExists: true),
        file(params.unfiltered, checkIfExists: true),
        file(params.cellbender_h5, checkIfExists: true)
    ])

    CELLBENDER_MERGE(merge_ch, params.output_layer_name)
}
