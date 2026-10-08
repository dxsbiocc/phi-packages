#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
//
// QUANTIFY_PSEUDO_ALIGNMENT expects a ready index, so this first builds one
// with SALMON_INDEX or KALLISTO_INDEX — the setup chain nf-core's own
// subworkflow test uses (see ../tests/main.nf.test). Those modules are
// included directly (section 6: never include another wrapper). The
// subworkflow branches on `pseudo_aligner` itself; only the matching index
// is built here.
//
// The fragment-length `val` inputs must not be null (Nextflow rejects that), so
// unset ones are passed as '' — Kallisto only reads them for single-end data.
//
// Two explicit read-file params instead of one glob: Nextflow rejects glob
// patterns for https:// sources. `reads_2` is optional (single-end when unset).
nextflow.enable.dsl = 2

include { SALMON_INDEX             } from '../../../../modules/nf-core/salmon/index/main.nf'
include { KALLISTO_INDEX           } from '../../../../modules/nf-core/kallisto/index/main.nf'
include { QUANTIFY_PSEUDO_ALIGNMENT } from '../main.nf'

params.reads_1             = null
params.reads_2             = null
params.index               = null
params.transcript_fasta    = null
params.gtf                 = null
params.samplesheet         = null
params.outdir              = null
params.pseudo_aligner      = 'salmon'
params.sample_id           = 'sample'
params.gtf_id_attribute    = 'gene_id'
params.gtf_extra_attribute = 'gene_name'
params.lib_type            = 'A'
params.kallisto_fraglen    = null
params.kallisto_fraglen_sd = null

workflow {
    ch_transcript_fasta = channel.value(file(params.transcript_fasta, checkIfExists: true))
    ch_gtf              = channel.value(file(params.gtf, checkIfExists: true))

    if (params.index) {
        def index = file(params.index, checkIfExists: true)
        // KALLISTO_QUANT's index input is tuple(meta, index) — unlike Salmon's bare path.
        ch_index = params.pseudo_aligner == 'salmon'
            ? channel.value(index)
            : channel.value([[id: 'kallisto_index'], index])
    } else if (params.pseudo_aligner == 'salmon') {
        SALMON_INDEX([], ch_transcript_fasta)
        ch_index = SALMON_INDEX.out.index
    } else {
        KALLISTO_INDEX(ch_transcript_fasta.map { fasta -> [[id: fasta.baseName], fasta] })
        // KALLISTO_QUANT's index input is tuple(meta, index) — unlike Salmon's bare path.
        ch_index = KALLISTO_INDEX.out.index
    }

    def single_end = !params.reads_2
    def reads      = [file(params.reads_1, checkIfExists: true)]
    if (!single_end) reads << file(params.reads_2, checkIfExists: true)
    ch_reads = channel.of([[id: params.sample_id, single_end: single_end], reads])

    ch_samplesheet = params.samplesheet
        ? channel.value([[id: 'samplesheet'], file(params.samplesheet, checkIfExists: true)])
        : channel.value([[:], []])

    QUANTIFY_PSEUDO_ALIGNMENT(
        ch_samplesheet,
        ch_reads,
        ch_index,
        ch_transcript_fasta,
        ch_gtf,
        params.gtf_id_attribute,
        params.gtf_extra_attribute,
        params.pseudo_aligner,
        false,
        params.lib_type,
        params.kallisto_fraglen ?: '',
        params.kallisto_fraglen_sd ?: '',
        false
    )
}
