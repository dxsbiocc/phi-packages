#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { INTERPROSCAN } from '../main.nf'

params.fasta = null
params.interproscan_database = null
params.args = null
params.outdir = null

workflow {
    INTERPROSCAN(
        Channel.value([[id: 'interproscan'], file(params.fasta, checkIfExists: true)]),
        file(params.interproscan_database, checkIfExists: true)
    )
}
