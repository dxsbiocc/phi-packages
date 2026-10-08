#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
//
// FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS has ~30 inputs. The ones not exposed
// in wrapper.yaml are fixed here: no BBSplit, no rRNA removal, no UMI, the
// Salmon index is built inside the subworkflow (make_salmon_index = true),
// and the strandedness thresholds keep nf-core's defaults. Genome files are
// only needed for `strandedness: auto`; without them empty channels are passed
// so the Salmon branch never runs.
//
// Two explicit read-file params instead of one glob: Nextflow rejects glob
// patterns for https:// sources. `reads_2` is optional (single-end when unset).
nextflow.enable.dsl = 2

include { FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS } from '../main.nf'

params.reads_1           = null
params.reads_2           = null
params.fasta             = null
params.transcript_fasta  = null
params.gtf               = null
params.outdir            = null
params.strandedness      = 'auto'
params.trimmer           = 'fastp'
params.skip_trimming     = false
params.skip_fastqc       = false
params.skip_linting      = false
params.min_trimmed_reads = 10000
params.sample_id         = 'sample'

// An optional genome file: a value channel when given, otherwise an empty channel.
def optionalFile(path) {
    return path ? channel.value(file(path, checkIfExists: true)) : channel.empty()
}

def percent(value) {
    return value == null ? '' : String.format('%.2f', value as double)
}

workflow {
    if (params.strandedness == 'auto' && !(params.fasta && params.transcript_fasta && params.gtf)) {
        error "strandedness 'auto' needs fasta, transcript_fasta and gtf (or set strandedness to forward, reverse or unstranded)"
    }

    def single_end = !params.reads_2
    def reads      = [file(params.reads_1, checkIfExists: true)]
    if (!single_end) reads << file(params.reads_2, checkIfExists: true)
    ch_reads = channel.of([[id: params.sample_id, single_end: single_end, strandedness: params.strandedness], reads])

    FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS(
        ch_reads,
        optionalFile(params.fasta),            // ch_fasta
        optionalFile(params.transcript_fasta), // ch_transcript_fasta
        optionalFile(params.gtf),              // ch_gtf
        [],                                    // ch_salmon_index (built inside)
        [],                                    // ch_sortmerna_index
        [],                                    // ch_bowtie2_index
        [],                                    // ch_bbsplit_index
        [],                                    // ch_rrna_fastas
        true,                                  // skip_bbsplit
        params.skip_fastqc,
        params.skip_trimming,
        true,                                  // skip_umi_extract
        params.skip_linting,
        true,                                  // make_salmon_index
        false,                                 // make_sortmerna_index
        false,                                 // make_bowtie2_index
        params.trimmer,
        params.min_trimmed_reads,
        false,                                 // save_trimmed (fastp: also keep reads failing filters)
        false,                                 // fastp_merge
        false,                                 // remove_ribo_rna
        'sortmerna',                           // ribo_removal_tool (unused)
        false,                                 // with_umi
        0,                                     // umi_discard_read
        false,                                 // save_merged_fastq
        0.8,                                   // stranded_threshold
        0.1                                    // unstranded_threshold
    )

    // The inferred strandedness lives only in the meta map — write it out.
    FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.reads
        .map { meta, _fastq ->
            def a = meta.salmon_strand_analysis
            [meta.id, meta.strandedness, percent(a?.forwardFragments), percent(a?.reverseFragments), percent(a?.unstrandedFragments)].join('\t')
        }
        .collectFile(
            name: 'strandedness.tsv',
            storeDir: "${params.outdir}/fastq_qc_trim_filter_setstrandedness",
            seed: 'sample\tstrandedness\tforward_pct\treverse_pct\tunstranded_pct',
            newLine: true
        )
}
