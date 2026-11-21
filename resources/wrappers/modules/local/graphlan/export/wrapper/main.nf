#!/usr/bin/env nextflow
// Agent-facing wrapper over the local module at ../main.nf, ported from the
// Snakemake wrapper in the Workflow project (wrappers/graphlan/export, MIT, dengxsh).
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GRAPHLAN_EXPORT } from '../main.nf'

params.taxonomy = null
params.args = null
params.outdir = null

workflow {
    GRAPHLAN_EXPORT(Channel.value([[id: 'graphlan'], file(params.taxonomy, checkIfExists: true)]))
}
