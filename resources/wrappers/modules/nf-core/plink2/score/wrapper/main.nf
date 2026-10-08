#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PLINK2_SCORE } from '../main.nf'
include { PLINK2_VCF } from '../../vcf/main.nf'

params.vcf        = null
params.scorefile  = null
params.outdir     = null

workflow {
    vcf_ch = Channel.value([[id: 'test', single_end: false], [file(params.vcf, checkIfExists: true)]])
    PLINK2_VCF(vcf_ch)
    in_ch = PLINK2_VCF.out.pgen.join(PLINK2_VCF.out.psam).join(PLINK2_VCF.out.pvar)
    PLINK2_SCORE(in_ch, file(params.scorefile, checkIfExists: true))

}
