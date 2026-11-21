#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { MERQURY_HAPMERS } from '../main.nf'

params.child_meryl = null
params.maternal_meryl = null
params.paternal_meryl = null
params.args = null
params.outdir = null

workflow {
    MERQURY_HAPMERS(
        Channel.value([[id: 'hapmers'], file(params.child_meryl, checkIfExists: true)]),
        file(params.maternal_meryl, checkIfExists: true),
        file(params.paternal_meryl, checkIfExists: true)
    )
}
