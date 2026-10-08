#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
//
// FASTQ_REMOVE_RRNA builds its own SortMeRNA / Bowtie2 index when told to, so
// this only includes the subworkflow (never another wrapper, section 6) and
// turns on the matching `make_*_index` flag. The rRNA FASTA is needed by
// sortmerna/bowtie2 only, so its absence is rejected up front for those.
//
// Two explicit read-file params instead of one glob: Nextflow rejects glob
// patterns for https:// sources. `reads_2` is optional (single-end when unset).
nextflow.enable.dsl = 2

include { FASTQ_REMOVE_RRNA } from '../main.nf'

params.reads_1           = null
params.reads_2           = null
params.rrna_fasta        = null
params.outdir            = null
params.ribo_removal_tool = 'sortmerna'
params.sample_id         = 'sample'

workflow {
    def tool = params.ribo_removal_tool
    if (tool != 'ribodetector' && !params.rrna_fasta) {
        error "ribo_removal_tool '${tool}' needs rrna_fasta (only ribodetector works without one)"
    }

    def single_end = !params.reads_2
    def reads      = [file(params.reads_1, checkIfExists: true)]
    if (!single_end) reads << file(params.reads_2, checkIfExists: true)
    ch_reads = channel.of([[id: params.sample_id, single_end: single_end], reads])

    ch_rrna_fastas = tool == 'ribodetector'
        ? channel.empty()
        : channel.of(file(params.rrna_fasta, checkIfExists: true))

    FASTQ_REMOVE_RRNA(
        ch_reads,
        ch_rrna_fastas,
        [],                        // ch_sortmerna_index (built inside)
        [],                        // ch_bowtie2_index (built inside)
        tool,
        tool == 'sortmerna',       // make_sortmerna_index
        tool == 'bowtie2'          // make_bowtie2_index
    )
}
