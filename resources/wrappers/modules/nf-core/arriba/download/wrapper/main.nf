#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { ARRIBA_DOWNLOAD } from '../main.nf'


params.genome     = 'GRCh38'
params.outdir     = null

workflow {
    ARRIBA_DOWNLOAD(params.genome)
}
