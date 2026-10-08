#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { UCSC_BEDCLIP } from '../main.nf'

params.bedgraph = null
params.sizes    = null
params.outdir   = null

workflow {
    bedgraph_ch = Channel
        .fromPath(params.bedgraph, checkIfExists: true)
        .map { bedgraph -> [[id: bedgraph.simpleName], bedgraph] }

    sizes_ch = Channel.fromPath(params.sizes, checkIfExists: true)

    UCSC_BEDCLIP(bedgraph_ch, sizes_ch)
}
