#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CTATSPLICING_STARTOCANCERINTRONS } from '../main.nf'
include { STARFUSION_BUILD } from '../../../starfusion/build/main.nf'
include { CTATSPLICING_PREPGENOMELIB } from '../../prepgenomelib/main.nf'
include { STAR_ALIGN } from '../../../star/align/main.nf'
include { SAMTOOLS_INDEX } from '../../../samtools/index/main.nf'

params.fasta      = null
params.gtf        = null
params.fusion_lib = null
params.pfam       = null
params.dfam_hmm   = null
params.dfam_h3f   = null
params.dfam_h3i   = null
params.dfam_h3m   = null
params.dfam_h3p   = null
params.cancer_introns = null
params.reads1     = null
params.reads2     = null
params.outdir     = null

workflow {
    STARFUSION_BUILD(Channel.value([[id: 'minigenome_fasta'], file(params.fasta, checkIfExists: true)]),
                     Channel.value([[id: 'minigenome_gtf'], file(params.gtf, checkIfExists: true)]),
                     file(params.fusion_lib, checkIfExists: true), 'homo_sapiens', file(params.pfam, checkIfExists: true),
                     [file(params.dfam_hmm, checkIfExists: true), file(params.dfam_h3f, checkIfExists: true), file(params.dfam_h3i, checkIfExists: true),
                      file(params.dfam_h3m, checkIfExists: true), file(params.dfam_h3p, checkIfExists: true)],
                     'https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/AnnotFilterRule.pm')
    CTATSPLICING_PREPGENOMELIB(STARFUSION_BUILD.out.reference, [file(params.cancer_introns, checkIfExists: true)])
    reads_ch = Channel.value([[id: 'test', single_end: true], [file(params.reads1, checkIfExists: true), file(params.reads2, checkIfExists: true)]])
    idx_ch = CTATSPLICING_PREPGENOMELIB.out.reference.map { meta, reference -> [meta, "${reference}/ref_genome.fa.star.idx"] }
    gtf_ch = CTATSPLICING_PREPGENOMELIB.out.reference.map { meta, reference -> [meta, "${reference}/ref_annot.gtf"] }
    STAR_ALIGN(reads_ch, idx_ch, gtf_ch, false)
    SAMTOOLS_INDEX(STAR_ALIGN.out.bam_sorted_aligned)

    in_ch = STAR_ALIGN.out.spl_junc_tab.join(STAR_ALIGN.out.junction, by: [0]).join(STAR_ALIGN.out.bam_sorted_aligned, by: [0]).join(SAMTOOLS_INDEX.out.index, by: [0])

    CTATSPLICING_STARTOCANCERINTRONS(in_ch, CTATSPLICING_PREPGENOMELIB.out.reference)
}
