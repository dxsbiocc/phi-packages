#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CUSTOM_ORFCOLLAPSE } from '../main.nf'


params.bed12      = null
params.catalogue  = null
params.orf_to_gene = null
params.aa_fasta   = null
params.cluster    = null
params.outdir     = null

workflow {
    in_ch = Channel.value([[id: 'cohort'], file(params.bed12, checkIfExists: true), file(params.catalogue, checkIfExists: true),
                           file(params.orf_to_gene, checkIfExists: true), file(params.aa_fasta, checkIfExists: true), file(params.cluster, checkIfExists: true)])

    CUSTOM_ORFCOLLAPSE(in_ch)
}
