#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { VCF_IMPUTE_GLIMPSE } from '../main.nf'

params.vcf = null
params.vcf_index = null
params.infos = null
params.ref_vcf = null
params.ref_index = null
params.region = "chr22"
params.map = null
params.panel_id = "panel"
params.chr = "chr22"
params.chunk = true
params.sample_id = "sample"
params.outdir = "results"

workflow {
    def meta_panel = [id: params.panel_id, chr: params.chr]
    ch_vcf = Channel.value([[id: params.sample_id], file(params.vcf, checkIfExists: true), file(params.vcf_index, checkIfExists: true), params.infos ?: []])
    ch_ref = Channel.value([meta_panel, file(params.ref_vcf, checkIfExists: true), file(params.ref_index, checkIfExists: true), params.region])
    ch_map = Channel.value([meta_panel, file(params.map, checkIfExists: true)])
    ch_chunks = Channel.empty()

    VCF_IMPUTE_GLIMPSE(ch_vcf, ch_ref, ch_chunks, ch_map, params.chunk)
}
