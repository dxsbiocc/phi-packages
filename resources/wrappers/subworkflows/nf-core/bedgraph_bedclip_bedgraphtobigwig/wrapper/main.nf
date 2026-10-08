#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1. Includes
// the subworkflow only — never another wrapper (section 6).
nextflow.enable.dsl = 2

include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG } from '../main.nf'

params.bedgraph = null
params.sizes    = null
params.outdir   = null

workflow {
    bedgraph_ch = Channel
        .fromPath(params.bedgraph, checkIfExists: true)
        .map { bedgraph -> [[id: bedgraph.simpleName], bedgraph] }

    sizes_ch = Channel.fromPath(params.sizes, checkIfExists: true)

    BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG(bedgraph_ch, sizes_ch)
}
