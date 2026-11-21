#!/usr/bin/env nextflow
// Agent-facing wrapper over the local module at ../main.nf, ported from the
// Snakemake wrapper in the Workflow project (wrappers/lefse/format, MIT, dengxsh).
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { LEFSE_FORMAT } from '../main.nf'

params.table = null
params.args = null
params.outdir = null

// The calling pipeline prepares the class/subclass/subject header rows
// (the Snakemake wrapper's group_line logic); this module only runs the formatter.
workflow {
    LEFSE_FORMAT(Channel.value([[id: 'lefse'], file(params.table, checkIfExists: true)]))
}
