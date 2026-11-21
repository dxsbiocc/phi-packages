#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { MMSEQS_CREATETSV } from '../main.nf'

params.db_result = null
params.db_query = null
params.db_target = null
params.args = null
params.outdir = null

workflow {
    MMSEQS_CREATETSV(
        Channel.value([[id: 'createtsv'], file(params.db_result, checkIfExists: true)]),
        Channel.value([[id: 'query'], params.db_query ? file(params.db_query, checkIfExists: true) : []]),
        Channel.value([[id: 'target'], file(params.db_target, checkIfExists: true)])
    )
}
