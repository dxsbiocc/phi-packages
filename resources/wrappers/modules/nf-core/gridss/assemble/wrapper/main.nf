#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GRIDSS_ASSEMBLE } from '../main.nf'
include { BWA_INDEX } from '../../../bwa/index/main.nf'
include { GRIDSS_PREPROCESS } from '../../preprocess/main.nf'

params.bam        = null
params.bai        = null
params.fasta      = null
params.fai        = null
params.config     = null
params.outdir     = null

workflow {
    bam_ch = Channel.of([[id: 'sample'], file(params.bam, checkIfExists: true), file(params.bai, checkIfExists: true)])
    fasta_ch = Channel.of([[id: 'fasta'], file(params.fasta, checkIfExists: true), file(params.fai, checkIfExists: true)])
    BWA_INDEX(Channel.value([[id: 'fasta'], file(params.fasta, checkIfExists: true)]))
    ref_ch = fasta_ch.join(BWA_INDEX.out.index)
    config_ch = Channel.value([[id: 'gridss_config'], file(params.config, checkIfExists: true)])
    GRIDSS_PREPROCESS(bam_ch, ref_ch)
    GRIDSS_ASSEMBLE(bam_ch.join(GRIDSS_PREPROCESS.out.preprocess_dir), ref_ch, config_ch)
}
