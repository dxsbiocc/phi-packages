#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { FOLDDISCO_QUERY } from '../main.nf'

params.pdb = null
params.query = null
params.index = null
params.structures = null
params.args = null
params.outdir = null

workflow {
    structures_ch = Channel.fromPath("${params.structures}/*").collect()

    FOLDDISCO_QUERY(
        Channel.value([[id: 'query'], file(params.pdb, checkIfExists: true), params.query]),
        structures_ch.map { s -> [[id: 'index'], file(params.index, checkIfExists: true), s] }
    )
}
