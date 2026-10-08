#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GLIMPSE_SAMPLE } from '../main.nf'
include { GLIMPSE_PHASE } from '../../phase/main.nf'
include { BCFTOOLS_INDEX } from '../../../bcftools/index/main.nf'

params.vcf        = null
params.csi        = null
params.ref        = null
params.ref_csi    = null
params.map        = null
params.outdir     = null

workflow {
    ch_sample = Channel.of('NA12878 2').collectFile(name: 'sampleinfos.txt')
    regions = Channel.of(['chr21:16570000-16610000', 'chr21:16570000-16610000'])
    input_vcf = Channel.value([[id: 'input'], file(params.vcf, checkIfExists: true), file(params.csi, checkIfExists: true)])
    ref_panel = Channel.value([file(params.ref, checkIfExists: true), file(params.ref_csi, checkIfExists: true)])
    ch_map = Channel.value([file(params.map, checkIfExists: true)])
    phase_in = input_vcf.combine(ch_sample).combine(regions)
        .map { meta, vcf, index, sample, regionI, regionO -> [[id: meta.id, region: regionI], vcf, index, sample, regionI, regionO] }
        .combine(ref_panel).combine(ch_map)
    GLIMPSE_PHASE(phase_in)
    BCFTOOLS_INDEX(GLIMPSE_PHASE.out.phased_variants)

    GLIMPSE_SAMPLE(GLIMPSE_PHASE.out.phased_variants.join(BCFTOOLS_INDEX.out.index))
}
