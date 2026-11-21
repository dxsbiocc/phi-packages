#!/usr/bin/env nextflow
// Agent-facing wrapper over the local module at ../main.nf, ported from the
// Snakemake wrapper in the Workflow project (wrappers/lefse/run, MIT, dengxsh).
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { LEFSE_RUN } from '../main.nf'

params.formatted = null
params.args = null
params.outdir = null

workflow {
    LEFSE_RUN(Channel.value([[id: 'lefse'], file(params.formatted, checkIfExists: true)]))
}
