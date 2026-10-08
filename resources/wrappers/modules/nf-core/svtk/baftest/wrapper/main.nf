#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { SVTK_BAFTEST } from '../main.nf'
include { MANTA_GERMLINE } from '../../../manta/germline/main.nf'
include { SVTK_VCF2BED } from '../../vcf2bed/main.nf'

params.cram       = null
params.crai       = null
params.fasta      = null
params.fai        = null
params.baf        = null
params.baf_tbi    = null
params.outdir     = null

workflow {
    cram_ch = Channel.value([[id: 'sample'], file(params.cram, checkIfExists: true), file(params.crai, checkIfExists: true), [], []])
    fasta_ch = Channel.value([[id: 'genome'], file(params.fasta, checkIfExists: true)])
    fai_ch = Channel.value([[id: 'genome'], file(params.fai, checkIfExists: true)])
    config_ch = Channel.of('[manta]', 'enableRemoteReadRetrievalForInsertionsInGermlineCallingModes = 0').collectFile(name: 'manta_options.ini', newLine: true)

    MANTA_GERMLINE(cram_ch, fasta_ch, fai_ch, config_ch)
    SVTK_VCF2BED(MANTA_GERMLINE.out.diploid_sv_vcf.combine(MANTA_GERMLINE.out.diploid_sv_vcf_tbi, by: 0).map { meta, vcf, tbi -> [[id: 'sample'], vcf, tbi] })
    baf_files = Channel.value([[id: 'sample'], file(params.baf, checkIfExists: true), file(params.baf_tbi, checkIfExists: true)])
    batch = Channel.of('sample\tgroup\tbatch', 'sample\tsample\tsample').collectFile(name: 'batch.txt', newLine: true, sort: true)

    SVTK_BAFTEST(SVTK_VCF2BED.out.bed.combine(baf_files, by: 0).combine(batch))
}
