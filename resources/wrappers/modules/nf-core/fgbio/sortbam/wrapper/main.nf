#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { FGBIO_SORTBAM } from '../main.nf'


params.bam        = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'sample'], file(params.bam, checkIfExists: true)])

    FGBIO_SORTBAM(in_ch)
}
