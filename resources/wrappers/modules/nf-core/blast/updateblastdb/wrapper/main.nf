#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BLAST_UPDATEBLASTDB } from '../main.nf'


params.name       = 'mito'
params.outdir     = null

workflow {
    BLAST_UPDATEBLASTDB(Channel.value([[id: 'result'], params.name]))
}
