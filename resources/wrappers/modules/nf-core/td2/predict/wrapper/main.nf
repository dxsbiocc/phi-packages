#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { TD2_PREDICT } from '../main.nf'

params.fasta = null
params.orfs_dir = null
params.args = null
params.outdir = null

workflow {
    TD2_PREDICT(Channel.value([
        [id: 'td2'],
        file(params.fasta, checkIfExists: true),
        file(params.orfs_dir, checkIfExists: true)
    ]))
}
