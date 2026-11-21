#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { KAIJU_MKFMI } from '../main.nf'

params.fasta = null
params.nodes_dmp = null
params.names_dmp = null
params.keep_intermediate = null
params.args = null
params.outdir = null

workflow {
    KAIJU_MKFMI(
        Channel.value([[id: 'kaiju', single_end: true], file(params.fasta, checkIfExists: true)]),
        file(params.nodes_dmp, checkIfExists: true),
        file(params.names_dmp, checkIfExists: true),
        params.keep_intermediate
    )
}
