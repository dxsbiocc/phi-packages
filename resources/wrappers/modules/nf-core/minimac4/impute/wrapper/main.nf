#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { MINIMAC4_IMPUTE } from '../main.nf'
include { MINIMAC4_COMPRESSREF } from '../../compressref/main.nf'

params.target     = null
params.target_csi = null
params.ref        = null
params.ref_csi    = null
params.sites      = null
params.sites_csi  = null
params.map        = null
params.outdir     = null

workflow {
    MINIMAC4_COMPRESSREF(Channel.value([[id: 'input'], file(params.ref, checkIfExists: true), file(params.ref_csi, checkIfExists: true)]))
    target_ch = Channel.value([[id: 'NA12878', chr: 'chr22'], file(params.target, checkIfExists: true), file(params.target_csi, checkIfExists: true),
                               file(params.sites, checkIfExists: true), file(params.sites_csi, checkIfExists: true), file(params.map, checkIfExists: true)])
    in_ch = target_ch.combine(MINIMAC4_COMPRESSREF.out.msav)
        .map { meta, target_vcf, target_index, sites_vcf, sites_index, map, metaRef, ref -> [meta, target_vcf, target_index, ref, sites_vcf, sites_index, map, 'chr22'] }

    MINIMAC4_IMPUTE(in_ch)
}
