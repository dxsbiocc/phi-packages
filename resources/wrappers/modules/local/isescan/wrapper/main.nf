#!/usr/bin/env nextflow
// Agent-facing wrapper over the local module at ../main.nf, ported from the
// Snakemake wrapper in the Workflow project (wrappers/isescan, MIT, dengxsh).
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { ISESCAN } from '../main.nf'

params.contigs = null
params.args = null
params.outdir = null

workflow {
    ISESCAN(Channel.value([[id: 'isescan'], file(params.contigs, checkIfExists: true)]))
}
