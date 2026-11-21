#!/usr/bin/env nextflow
// Agent-facing wrapper over the local module at ../main.nf, ported from the
// Snakemake wrapper in the Workflow project (wrappers/greenhill, MIT, dengxsh).
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GREENHILL } from '../main.nf'

params.hic = null
params.aware = null
params.bubble = null
params.cph = null
params.long = null
params.hic_type = null
params.args = null
params.outdir = null

workflow {
    if (params.aware && params.cph) {
        error "greenhill: 'aware' and 'cph' cannot be used at the same time."
    }
    def files = { p -> p ? file(p, checkIfExists: true).listFiles().findAll { it.isFile() }.toList() : [] }

    GREENHILL(
        Channel.value([[id: 'greenhill'], file(params.hic, checkIfExists: true)]),
        params.aware ? file(params.aware, checkIfExists: true) : [],
        files(params.bubble),
        params.cph ? file(params.cph, checkIfExists: true) : [],
        files(params.long),
        params.hic_type
    )
}
