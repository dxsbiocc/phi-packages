#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { SOURCEPREDICT } from '../main.nf'

params.kraken_parse = null
params.sources = null
params.labels = null
params.taxa_sqlite = null
params.taxa_sqlite_traverse_pkl = null
params.args = null
params.outdir = null

workflow {
    SOURCEPREDICT(
        Channel.value([[id: 'sourcepredict'], file(params.kraken_parse, checkIfExists: true)]),
        file(params.sources, checkIfExists: true),
        file(params.labels, checkIfExists: true),
        file(params.taxa_sqlite, checkIfExists: true),
        file(params.taxa_sqlite_traverse_pkl, checkIfExists: true)
    )
}
