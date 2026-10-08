#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1. Includes
// the subworkflow only — never another wrapper (section 6).
nextflow.enable.dsl = 2

include { SAMTOOLS_INDEX } from '../../../../modules/nf-core/samtools/index/main.nf'
include { BAM_QC_RNASEQ } from '../main.nf'

params.bam             = null
params.gtf             = null
params.bed             = null
params.fasta           = null
params.fai              = null
params.biotypes_header = null
params.single_end      = false
params.strandedness    = 'reverse'
params.biotype         = 'gene_biotype'
params.tools           = 'preseq,biotype_qc,qualimap,dupradar,rseqc_bam_stat,rseqc_inner_distance,rseqc_infer_experiment,rseqc_junction_annotation,rseqc_junction_saturation,rseqc_read_distribution,rseqc_read_duplication,rseqc_tin'
params.outdir          = null

workflow {
    bam_ch = Channel
        .fromPath(params.bam, checkIfExists: true)
        .map { bam -> [[id: bam.simpleName, single_end: params.single_end, strandedness: params.strandedness], bam] }

    SAMTOOLS_INDEX(bam_ch)

    bam_bai_ch = bam_ch.join(SAMTOOLS_INDEX.out.index)

    gtf_ch = Channel.value([[id: 'annotation'], file(params.gtf, checkIfExists: true)])

    bed_ch = Channel.value(file(params.bed, checkIfExists: true))

    fasta_fai_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true)])

    biotypes_header_ch = Channel.value([[:], file(params.biotypes_header, checkIfExists: true)])

    tools_list = params.tools.split(',') as List

    BAM_QC_RNASEQ(bam_bai_ch, gtf_ch, bed_ch, fasta_fai_ch, biotypes_header_ch, tools_list, params.biotype)
}
