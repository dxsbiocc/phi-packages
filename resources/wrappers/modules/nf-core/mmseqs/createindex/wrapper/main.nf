#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { MMSEQS_CREATEINDEX } from '../main.nf'
include { UNTAR } from '../../../untar/main.nf'

params.db_archive = null
params.outdir     = null

workflow {
    UNTAR(Channel.value([[id: 'test'], file(params.db_archive, checkIfExists: true)]))
    MMSEQS_CREATEINDEX(UNTAR.out.untar)
}
