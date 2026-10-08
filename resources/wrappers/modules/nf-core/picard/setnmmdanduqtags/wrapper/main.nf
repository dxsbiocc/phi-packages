#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PICARD_SETNMMDANDUQTAGS } from '../main.nf'


params.bam        = null
params.fasta      = null
params.outdir     = null

workflow {
    bam_ch = Channel.value([[id: 'test', single_end: false], file(params.bam, checkIfExists: true)])
    ref_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])

    PICARD_SETNMMDANDUQTAGS(bam_ch, ref_ch)
}
