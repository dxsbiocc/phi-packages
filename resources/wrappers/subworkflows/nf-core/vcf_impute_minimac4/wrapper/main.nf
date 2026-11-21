#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { VCF_IMPUTE_MINIMAC4 } from '../main.nf'

params.vcf = null
params.vcf_index = null
params.panel_vcf = null
params.panel_index = null
params.sites_vcf = null
params.sites_index = null
params.map = null
params.region = "chr22"
params.panel_id = "panel"
params.chr = "chr22"
params.sample_id = "sample"
params.outdir = "results"

workflow {
    def meta_panel = [id: params.panel_id, chr: params.chr]
    ch_input = Channel.value([[id: params.sample_id, chr: params.chr], file(params.vcf, checkIfExists: true), file(params.vcf_index, checkIfExists: true)])
    ch_panel = Channel.value([meta_panel, file(params.panel_vcf, checkIfExists: true), file(params.panel_index, checkIfExists: true)])
    ch_posfile = Channel.value([meta_panel, file(params.sites_vcf, checkIfExists: true), file(params.sites_index, checkIfExists: true)])
    ch_map = Channel.value([meta_panel, file(params.map, checkIfExists: true)])
    ch_chunks = Channel.value([meta_panel, params.region])

    VCF_IMPUTE_MINIMAC4(ch_input, ch_panel, ch_posfile, ch_chunks, ch_map)
}
