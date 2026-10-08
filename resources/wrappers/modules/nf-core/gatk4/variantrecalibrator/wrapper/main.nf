#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { GATK4_VARIANTRECALIBRATOR } from '../main.nf'


params.vcf        = null
params.tbi        = null
params.hapmap     = null
params.hapmap_tbi = null
params.omni       = null
params.omni_tbi   = null
params.g1000      = null
params.g1000_tbi  = null
params.dbsnp      = null
params.dbsnp_tbi  = null
params.fasta      = null
params.fai        = null
params.dict       = null
params.outdir     = null

workflow {
    vcf_ch  = Channel.value([[id: 'test'], file(params.vcf, checkIfExists: true), file(params.tbi, checkIfExists: true)])
    res_ch  = Channel.value([file(params.hapmap, checkIfExists: true), file(params.omni, checkIfExists: true), file(params.g1000, checkIfExists: true), file(params.dbsnp, checkIfExists: true)])
    tbi_ch  = Channel.value([file(params.hapmap_tbi, checkIfExists: true), file(params.omni_tbi, checkIfExists: true), file(params.g1000_tbi, checkIfExists: true), file(params.dbsnp_tbi, checkIfExists: true)])
    labels  = [
        '--resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz',
        '--resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz',
        '--resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.hg38.vcf.gz',
        '--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.hg38.vcf.gz'
    ]
    fasta_ch = Channel.value(file(params.fasta, checkIfExists: true))
    fai_ch   = Channel.value(file(params.fai, checkIfExists: true))
    dict_ch  = Channel.value(file(params.dict, checkIfExists: true))

    GATK4_VARIANTRECALIBRATOR(vcf_ch, res_ch, tbi_ch, labels, fasta_ch, fai_ch, dict_ch)
}
