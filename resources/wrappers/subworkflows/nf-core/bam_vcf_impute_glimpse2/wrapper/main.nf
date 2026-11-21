#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { BAM_VCF_IMPUTE_GLIMPSE2 } from '../main.nf'

params.vcf = null
params.vcf_index = null
params.infos = null
params.ref_vcf = null
params.ref_index = null
params.region = "chr22"
params.map = null
params.panel_id = "panel"
params.chr = "chr22"
params.fasta = null
params.fai = null
params.chunk = true
params.chunk_model = null
params.splitreference = false
params.sample_id = "sample"
params.outdir = "results"

workflow {
    def meta_panel = [id: params.panel_id, chr: params.chr]
    ch_input = Channel.value([[id: params.sample_id], file(params.vcf, checkIfExists: true), file(params.vcf_index, checkIfExists: true), [], params.infos ?: []])
    ch_ref = Channel.value([meta_panel, file(params.ref_vcf, checkIfExists: true), file(params.ref_index, checkIfExists: true), params.region])
    ch_map = Channel.value([meta_panel, file(params.map, checkIfExists: true)])
    ch_chunks = Channel.empty()
    ch_fasta = Channel.value(params.fasta ? [[id: 'reference'], file(params.fasta, checkIfExists: true), params.fai ? file(params.fai, checkIfExists: true) : []] : [[], [], []])

    BAM_VCF_IMPUTE_GLIMPSE2(ch_input, ch_ref, ch_chunks, ch_map, ch_fasta, params.chunk, params.chunk_model ?: null, params.splitreference)
}
