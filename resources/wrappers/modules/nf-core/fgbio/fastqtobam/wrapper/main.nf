#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { FGBIO_FASTQTOBAM } from '../main.nf'


params.reads1     = null
params.reads2     = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], [file(params.reads1, checkIfExists: true), file(params.reads2, checkIfExists: true)]])

    FGBIO_FASTQTOBAM(in_ch)
}
