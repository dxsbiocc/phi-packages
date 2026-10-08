#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1. Includes
// the subworkflow only — never another wrapper (section 6).
nextflow.enable.dsl = 2

include { SAMTOOLS_INDEX } from '../../../../modules/nf-core/samtools/index/main.nf'
include { BAM_RSEQC } from '../main.nf'

params.bam            = null
params.bed            = null
params.rseqc_modules  = 'bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication,tin'
params.outdir         = null

workflow {
    bam_ch = Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: bam.simpleName], bam] }

    SAMTOOLS_INDEX(bam_ch)

    bam_bai_ch = bam_ch.join(SAMTOOLS_INDEX.out.index)
        .map { meta, bam, bai -> [meta, [bam, bai]] }

    bed_ch = Channel.value(file(params.bed, checkIfExists: true))

    modules_list = params.rseqc_modules.split(',') as List

    BAM_RSEQC(bam_bai_ch, bed_ch, modules_list)
}
