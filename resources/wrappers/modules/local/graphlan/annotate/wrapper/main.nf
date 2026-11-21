#!/usr/bin/env nextflow
// Agent-facing wrapper over the local module at ../main.nf, ported from the
// Snakemake wrapper in the Workflow project (wrappers/graphlan/annotate, MIT, dengxsh).
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GRAPHLAN_ANNOTATE } from '../main.nf'

params.tree = null
params.annotation = null
params.args = null
params.outdir = null

workflow {
    GRAPHLAN_ANNOTATE(Channel.value([
        [id: 'graphlan'],
        file(params.tree, checkIfExists: true),
        file(params.annotation, checkIfExists: true)
    ]))
}
