#!/usr/bin/env nextflow
// Agent-facing wrapper over the local module at ../main.nf, ported from the
// Snakemake wrapper in the Workflow project (wrappers/virsorter2, MIT, dengxsh).
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { VIRSORTER2 } from '../main.nf'

params.contigs = null
params.db_dir = null
params.mode = null
params.include_groups = null
params.min_length = null
params.args = null
params.outdir = null

workflow {
    VIRSORTER2(
        Channel.value([[id: 'virsorter2'], file(params.contigs, checkIfExists: true)]),
        params.db_dir ? file(params.db_dir, checkIfExists: true) : [],
        params.mode,
        params.include_groups ?: '',
        params.min_length
    )
}
