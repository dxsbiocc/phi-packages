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
// 'independent' keeps the existing unpaired bulk behaviour; 'paired' adds an
// additive donor block (~ donor + condition) and requires donor_col.
params.study_design      = 'independent'
params.donor_col         = null
params.outdir            = null

// The paired-matrix contract, evaluated before the DE process is dispatched so
// a missing/incomplete mapping becomes a user question even without preflight.
def pairsComplete(rows, sampleCol, groupCol, donorCol, reference, target) {
    if (!rows || ![sampleCol, groupCol, donorCol].every { rows[0].containsKey(it) }) {
        return false
    }
    def samples = rows.collect { it[sampleCol] }
    if (samples.any { !it?.trim() } || samples.toSet().size() != samples.size()) {
        return false
    }
    def paired = rows.findAll { it[groupCol] in [reference, target] }
    def donors = paired.groupBy { it[donorCol]?.trim() }
    return donors.size() >= 2 && !donors.any { donor, records ->
        !donor || donor.toLowerCase() in ['na', 'n/a', 'nan', 'null'] ||
            records.size() != 2 || records.collect { it[groupCol] }.toSet() != [reference, target].toSet()
    }
}

// Static, domain-owned decision; no sample values or log content. CREATE_NEW
// guarantees an existing question is never overwritten.
def publishPairingQuestion(outdir, prefix) {
    def payload = groovy.json.JsonOutput.toJson([
        schemaVersion: 1,
        question: [
            header: 'Paired design',
            question: 'The paired DESeq2 run stopped before inference because the sample-to-donor mapping is missing or incomplete. How should the analysis proceed?',
            options: [
                [label: 'Correct the mapping',
                 description: 'Provide a samplesheet with one sample per group for each of at least two donors, then launch a fresh run into a new output directory.'],
                [label: 'Descriptive report only',
                 description: 'Skip inferential statistics and report descriptive summaries of the samples without pairing.']
            ]
        ]
    ])
    def folder = file(outdir).resolve('deseq2_differential')
    java.nio.file.Files.createDirectories(folder)
    java.nio.file.Files.writeString(folder.resolve(prefix + '.deseq2.input.json'), payload, java.nio.file.StandardOpenOption.CREATE_NEW)

}

workflow {
    if (!(params.prefix ==~ /^[A-Za-z0-9_-][A-Za-z0-9_.-]{0,127}$/)) {
        error "Use a short filename prefix with letters, numbers, dots, underscores or hyphens."
    }
    if (params.study_design == 'paired' && !params.donor_col) {
        publishPairingQuestion(params.outdir, params.prefix)
        error "Paired DESeq2 requires the sample-to-donor mapping before inference."
    }
    if (!(params.study_design in ['independent', 'paired'])) {
        error "Confirm whether the study is independent or paired before running DESeq2."
    }
    if (params.study_design == 'paired' && (params.formula || params.comparison ||
        !params.donor_col || !(params.donor_col ==~ /^[A-Za-z.][A-Za-z0-9._]*$/) ||
        params.donor_col in [params.sample_id_col, params.contrast_variable] ||
        params.reference_level == params.target_level)) {
        error "Paired DESeq2 needs a donor column distinct from sample and condition, two different groups and the simple contrast. Ask for the missing design information."
    }
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

    sheet = file(params.samplesheet, checkIfExists: true)
    matrix = tuple([id: params.prefix], sheet, file(params.counts, checkIfExists: true))
    if (params.study_design == 'paired') {
        sep = sheet.name.toLowerCase().endsWith('.csv') ? ',' : '\t'
        matrix_ch = Channel.of(sheet).splitCsv(header: true, sep: sep, quote: '"').toList().map { rows ->
            if (!pairsComplete(rows, params.sample_id_col, params.contrast_variable, params.donor_col, params.reference_level, params.target_level)) {
                publishPairingQuestion(params.outdir, params.prefix)
                throw new IllegalArgumentException("Paired DESeq2 requires a complete sample-to-donor mapping before inference.")
            }
            matrix
        }
    } else {
        matrix_ch = Channel.of(matrix)
    }

    // Optional control-genes / transcript-lengths inputs, deliberately not
    // exposed as params (no test data to verify the spike-in / length-offset
    // paths) — this is the module's own convention for "skip me".
    spikes_ch  = Channel.of([[], []])
    lengths_ch = Channel.of([[], []])

    // Filtered gene tables are written by DESEQ2_DIFFERENTIAL itself
    // (|log2FoldChange| >= log2fc_threshold and padj < padj_threshold).
    DESEQ2_DIFFERENTIAL(contrast_ch, matrix_ch, spikes_ch, lengths_ch)
}