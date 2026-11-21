#!/usr/bin/env nextflow
// Agent-facing wrapper over the local module at ../main.nf, ported from the
// Snakemake wrapper in the Workflow project (wrappers/pysradb, MIT, dengxsh).
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PYSRADB } from '../main.nf'

params.command = null
params.query = null
params.args = null
params.outdir = null

// Requires network access to NCBI/ENA.
workflow {
    PYSRADB(Channel.value([[id: 'pysradb'], params.command, params.query]))
}
