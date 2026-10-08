#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Composes star/genomegenerate directly (vendored module, not a wrapper) so
// callers only supply a genome FASTA/GTF instead of a pre-built STAR index
// — same reasoning as the star/align and align_star wrappers. Unlike the
// module's own nf-test (which pins STAR_GENOMEGENERATE to an old 2.6.1d
// container to reproduce a genuinely legacy-format index), this wrapper
// builds the index with the current STAR version, so the upgrade rewrite
// is a documented no-op here — it still proves the channel wiring and
// process run end-to-end; only the "was the string actually rewritten"
// behavior isn't exercised.
nextflow.enable.dsl = 2

include { STAR_GENOMEGENERATE } from '../../../nf-core/star/genomegenerate/main.nf'
include { STAR_GENOMEPARAMS_UPGRADE } from '../main.nf'

params.fasta  = null
params.gtf    = null
params.outdir = null

workflow {
    fasta_ch = Channel.value([[id: 'genome'], [file(params.fasta, checkIfExists: true)]])
    gtf_ch   = Channel.value([[id: 'genome'], [file(params.gtf, checkIfExists: true)]])

    STAR_GENOMEGENERATE(fasta_ch, gtf_ch)

    index_ch = STAR_GENOMEGENERATE.out.index.map { meta, index -> [[id: 'star_index'], index] }

    STAR_GENOMEPARAMS_UPGRADE(index_ch)
}
