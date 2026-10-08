#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored DESeq2 differential-expression
// module at ../main.nf (nf-core/modules `deseq2/differential`, as used by
// nf-core/differentialabundance). See
// docs/design/phi-wrapper-agent-composition-design.md section 1.
//
// ../main.nf's `input:` block takes four tuples: contrast spec, the
// samplesheet+counts matrix, and two optional file inputs (control genes for
// spike-in size-factor estimation, transcript lengths for length-scaled
// offsets) that this adapter does not expose — pass `[[], []]` for each,
// exactly as the module's own tests do to skip them.
nextflow.enable.dsl = 2

include { DESEQ2_DIFFERENTIAL } from '../main.nf'
include { CUSTOM_FILTERDIFFERENTIALTABLE } from '../../../../nf-core/custom/filterdifferentialtable/main.nf'

params.counts            = null
params.samplesheet       = null
params.contrast_variable = null
params.reference_level   = null
params.target_level      = null
params.formula           = null
params.comparison        = null
params.sample_id_col     = 'sample'
params.gene_id_col       = 'gene_id'
params.prefix            = 'de'
// A gene is differentially expressed when |log2FC| >= log2fc_threshold and padj < padj_threshold.
params.padj_threshold    = 0.05
params.log2fc_threshold  = 1
params.outdir            = null

workflow {
    contrast_ch = Channel.of(
        tuple(
            [id: params.prefix, variable: params.contrast_variable],
            params.contrast_variable,
            params.reference_level,
            params.target_level,
            params.formula,
            params.comparison
        )
    )

    matrix_ch = Channel.of(
        tuple(
            [id: params.prefix],
            file(params.samplesheet, checkIfExists: true),
            file(params.counts, checkIfExists: true)
        )
    )

    // Optional control-genes / transcript-lengths inputs, deliberately not
    // exposed as params (no test data to verify the spike-in / length-offset
    // paths) — this is the module's own convention for "skip me".
    spikes_ch  = Channel.of([[], []])
    lengths_ch = Channel.of([[], []])

    DESEQ2_DIFFERENTIAL(contrast_ch, matrix_ch, spikes_ch, lengths_ch)

    // The filter module takes a linear fold change and applies it to |log2FoldChange|.
    CUSTOM_FILTERDIFFERENTIALTABLE(
        DESEQ2_DIFFERENTIAL.out.results,
        Channel.value(['log2FoldChange', Math.pow(2, params.log2fc_threshold as double), '>=']),
        Channel.value(['padj', params.padj_threshold, '<'])
    )
}
