#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { VCF_PHASE_SHAPEIT5 } from '../main.nf'

params.vcf = null
params.vcf_index = null
params.ref_vcf = null
params.ref_index = null
params.scaffold_vcf = null
params.scaffold_index = null
params.map = null
params.chunk_model = null
params.region = "chr22"
params.pedigree = null
params.panel_id = "sample"
params.chr = "chr22"
params.chunk = true
params.outdir = "results"

workflow {
    def meta = [id: params.panel_id, chr: params.chr]
    ch_vcf = Channel.value([meta, file(params.vcf, checkIfExists: true), file(params.vcf_index, checkIfExists: true), params.pedigree ? file(params.pedigree, checkIfExists: true) : [], params.region])
    ch_ref = Channel.value([meta, file(params.ref_vcf, checkIfExists: true), file(params.ref_index, checkIfExists: true)])
    ch_scaffold = Channel.value([meta, file(params.scaffold_vcf, checkIfExists: true), file(params.scaffold_index, checkIfExists: true)])
    ch_map = Channel.value([meta, file(params.map, checkIfExists: true)])
    ch_chunks = Channel.empty()
    chunk_model_ch = Channel.value(file(params.chunk_model, checkIfExists: true))

    VCF_PHASE_SHAPEIT5(ch_vcf, ch_chunks, ch_ref, ch_scaffold, ch_map, params.chunk, chunk_model_ch)
}
