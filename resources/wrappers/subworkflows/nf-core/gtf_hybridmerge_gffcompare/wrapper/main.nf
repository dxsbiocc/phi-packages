#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { GTF_HYBRIDMERGE_GFFCOMPARE } from '../main.nf'

params.novel_gtf = null
params.reference_gtf = null
params.backbone_gtf = null
params.class_codes = "u,j,i"
params.blacklist_bed = null
params.outdir = "results"

workflow {
    novel_ch = Channel.fromPath(params.novel_gtf, checkIfExists: true).map { f -> [[id: 'hybrid'], f] }
    ref_ch = Channel.value([[id: 'hybrid'], file(params.reference_gtf, checkIfExists: true)])
    backbone_ch = Channel.value([[id: 'hybrid'], file(params.backbone_gtf, checkIfExists: true)])
    blacklist_ch = params.blacklist_bed ? Channel.value([[id: 'hybrid'], file(params.blacklist_bed, checkIfExists: true)]) : Channel.empty()

    GTF_HYBRIDMERGE_GFFCOMPARE(novel_ch, ref_ch, backbone_ch, params.class_codes, blacklist_ch)
}
