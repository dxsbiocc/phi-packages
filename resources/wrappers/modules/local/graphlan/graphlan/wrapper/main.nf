#!/usr/bin/env nextflow
// Agent-facing wrapper over the local module at ../main.nf, ported from the
// Snakemake wrapper in the Workflow project (wrappers/graphlan/graphlan, MIT, dengxsh).
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GRAPHLAN_GRAPHLAN } from '../main.nf'

params.xml = null
params.format = null
params.args = null
params.outdir = null

workflow {
    GRAPHLAN_GRAPHLAN(Channel.value([[id: 'graphlan'], file(params.xml, checkIfExists: true)]), params.format)
}
