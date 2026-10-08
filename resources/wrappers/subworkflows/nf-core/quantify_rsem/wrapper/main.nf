#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Composes rsem/preparereference (vendored module, not a wrapper) so
// callers only supply a genome FASTA/GTF instead of a pre-built RSEM index
// — same reasoning as the rsem/calculateexpression module wrapper.
nextflow.enable.dsl = 2

include { RSEM_PREPAREREFERENCE } from '../../../../modules/nf-core/rsem/preparereference/main.nf'
include { QUANTIFY_RSEM } from '../main.nf'

params.reads_1     = null
params.reads_2     = null
params.fasta       = null
params.gtf         = null
params.samplesheet = null
params.id          = 'gene_id'
params.extra       = 'gene_name'
params.skip_merge  = false
params.outdir      = null

workflow {
    RSEM_PREPAREREFERENCE(
        file(params.fasta, checkIfExists: true),
        file(params.gtf, checkIfExists: true)
    )

    reads_ch = Channel.of([
        [id: 'test', single_end: false],
        [file(params.reads_1, checkIfExists: true), file(params.reads_2, checkIfExists: true)]
    ])

    samplesheet_ch = Channel.value([[id: 'samplesheet'], file(params.samplesheet, checkIfExists: true)])

    gtf_ch = Channel.of(file(params.gtf, checkIfExists: true))

    QUANTIFY_RSEM(
        samplesheet_ch,
        reads_ch,
        RSEM_PREPAREREFERENCE.out.index,
        gtf_ch,
        params.id,
        params.extra,
        params.skip_merge
    )
}
