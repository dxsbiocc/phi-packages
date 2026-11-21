#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { KRONA_KTIMPORTTAXONOMY } from '../main.nf'

params.report = null
params.taxonomy = null
params.args = null
params.outdir = null

workflow {
    KRONA_KTIMPORTTAXONOMY(
        Channel.value([[id: 'krona'], file(params.report, checkIfExists: true)]),
        file(params.taxonomy, checkIfExists: true)
    )
}
