#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Composes star/genomegenerate directly (vendored module, not a wrapper) so
// callers only supply a genome FASTA/GTF instead of a pre-built STAR index
// — same reasoning as the star/align and align_star wrappers.
nextflow.enable.dsl = 2

include { STAR_GENOMEGENERATE } from '../../genomegenerate/main.nf'
include { STAR_STARSOLO } from '../main.nf'

params.fasta       = null
params.gtf         = null
params.cb_read     = null
params.cdna_read   = null
params.umi_len     = 12
params.solo_type   = 'CB_UMI_Simple'
params.outdir      = null

workflow {
    fasta_ch = Channel.value([[id: 'genome'], [file(params.fasta, checkIfExists: true)]])
    gtf_ch   = Channel.value([[id: 'genome'], [file(params.gtf, checkIfExists: true)]])

    STAR_GENOMEGENERATE(fasta_ch, gtf_ch)

    reads_ch = Channel.of([
        [id: 'test', umi_len: params.umi_len],
        params.solo_type,
        [file(params.cb_read, checkIfExists: true), file(params.cdna_read, checkIfExists: true)]
    ])

    whitelist_ch = Channel.value(file('NO_FILE'))

    STAR_STARSOLO(reads_ch, whitelist_ch, STAR_GENOMEGENERATE.out.index)
}
