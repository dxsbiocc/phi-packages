#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CNVKIT_BATCH } from '../main.nf'


params.tumor      = null
params.tumor_bai  = null
params.normal     = null
params.normal_bai = null
params.fasta      = null
params.baits      = null
params.outdir     = null

workflow {
    bam_ch = Channel.value([[id: 'sample'], file(params.tumor, checkIfExists: true), file(params.tumor_bai, checkIfExists: true),
                            file(params.normal, checkIfExists: true), file(params.normal_bai, checkIfExists: true)])
    fasta_ch = Channel.value([[:], file(params.fasta, checkIfExists: true), []])
    targets_ch = Channel.value([[:], file(params.baits, checkIfExists: true)])
    ref_ch = Channel.value([[:], []])

    CNVKIT_BATCH(bam_ch, fasta_ch, targets_ch, ref_ch, false)
}
