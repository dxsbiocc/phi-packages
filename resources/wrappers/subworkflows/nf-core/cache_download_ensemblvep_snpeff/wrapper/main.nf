#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { CACHE_DOWNLOAD_ENSEMBLVEP_SNPEFF } from '../main.nf'

params.vep_assembly = "GRCh38"
params.vep_species = "homo_sapiens"
params.vep_cache_version = "110"
params.snpeff_db = "GRCh38.99"
params.preflight_check = false
params.outdir = "results"

workflow {
    vep_info = Channel.value([[id: 'vep_cache'], params.vep_assembly, params.vep_species, params.vep_cache_version])
    snpeff_info = Channel.value([[id: 'snpeff_cache'], params.snpeff_db])

    CACHE_DOWNLOAD_ENSEMBLVEP_SNPEFF(vep_info, snpeff_info, params.preflight_check)
}
