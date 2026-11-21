#!/usr/bin/env nextflow
// Agent-facing wrapper over the local module at ../main.nf, ported from the
// Snakemake wrapper in the Workflow project (wrappers/lefse/res, MIT, dengxsh).
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { LEFSE_PLOT_RES } from '../main.nf'

params.res = null
params.args = null
params.outdir = null

workflow {
    LEFSE_PLOT_RES(Channel.value([[id: 'lefse'], file(params.res, checkIfExists: true)]))
}
