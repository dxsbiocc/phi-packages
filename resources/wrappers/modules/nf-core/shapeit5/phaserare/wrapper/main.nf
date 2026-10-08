#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SHAPEIT5_PHASERARE } from '../main.nf'
include { SHAPEIT5_PHASECOMMON } from '../../phasecommon/main.nf'
include { BCFTOOLS_INDEX } from '../../../bcftools/index/main.nf'

params.vcf        = null
params.csi        = null
params.region     = 'chr22'
params.outdir     = null

workflow {
    scaffold_in = Channel.value([[id: 'scaffold', single_end: false], file(params.vcf, checkIfExists: true), file(params.csi, checkIfExists: true), [], params.region, [], [], [], [], []])
    SHAPEIT5_PHASECOMMON(scaffold_in)
    BCFTOOLS_INDEX(SHAPEIT5_PHASECOMMON.out.phased_variant)

    rare_in = Channel.value([[id: 'input', single_end: false], file(params.vcf, checkIfExists: true), file(params.csi, checkIfExists: true), [], params.region])
        .join(SHAPEIT5_PHASECOMMON.out.phased_variant.map { meta, vcf -> [[id: 'input'], vcf] })
        .join(BCFTOOLS_INDEX.out.index.map { meta, csi -> [[id: 'input'], csi] })
        .map { meta, vcf, csi, samples, region, scaffold, scaffold_csi -> [meta, vcf, csi, samples, region, scaffold, scaffold_csi, params.region, []] }

    SHAPEIT5_PHASERARE(rare_in, [])
}
