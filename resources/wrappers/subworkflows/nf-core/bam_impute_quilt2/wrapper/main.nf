#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { BAM_IMPUTE_QUILT2 } from '../main.nf'

params.bam = null
params.bai = null
params.bampaths = null
params.bamnames = null
params.panel_vcf = null
params.panel_index = null
params.map = null
params.panel_id = "panel"
params.chr = "chr22"
params.start = null
params.end = null
params.fasta = null
params.fai = null
params.n_gen = 100
params.buffer = 100000
params.sample_id = "sample"
params.outdir = "results"

workflow {
    def meta_panel = [id: params.panel_id, chr: params.chr]
    ch_input = Channel.value([[id: params.sample_id], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true), params.bampaths ?: '', params.bamnames ?: ''])
    ch_panel = Channel.value([meta_panel, file(params.panel_vcf, checkIfExists: true), file(params.panel_index, checkIfExists: true)])
    ch_map = Channel.value([meta_panel, file(params.map, checkIfExists: true)])
    ch_chunks = Channel.value([meta_panel, params.chr, params.start ?: [], params.end ?: []])
    ch_fasta = Channel.value(params.fasta ? [[id: 'reference'], file(params.fasta, checkIfExists: true), params.fai ? file(params.fai, checkIfExists: true) : []] : [[id: 'reference'], [], []])

    BAM_IMPUTE_QUILT2(ch_input, ch_panel, ch_chunks, ch_map, ch_fasta, params.n_gen, params.buffer)
}
