#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { ZIP } from '../main.nf'

params.file   = null
params.outdir = null

workflow {
    files_ch = Channel.value([[id: 'zipped'], [file(params.file, checkIfExists: true)]])

    ZIP(files_ch)
}
