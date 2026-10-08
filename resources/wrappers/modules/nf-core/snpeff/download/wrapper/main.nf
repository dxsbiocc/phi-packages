#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SNPEFF_DOWNLOAD } from '../main.nf'

params.db     = 'NC_045512.2'
params.outdir = null

workflow {
    db_ch = Channel.value([[id: params.db], params.db])

    SNPEFF_DOWNLOAD(db_ch)
}
