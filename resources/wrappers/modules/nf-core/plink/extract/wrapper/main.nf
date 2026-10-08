#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { PLINK_EXTRACT } from '../main.nf'
include { PLINK_VCF } from '../../vcf/main.nf'
include { GAWK } from '../../../gawk/main.nf'

params.vcf        = null
params.outdir     = null

workflow {
    vcf_ch = Channel.value([[id: 'sample', single_end: false], file(params.vcf, checkIfExists: true)])
    PLINK_VCF(vcf_ch)
    bbf_ch = PLINK_VCF.out.bed.join(PLINK_VCF.out.bim).join(PLINK_VCF.out.fam)
    GAWK(PLINK_VCF.out.bim.map { meta, bim -> [meta, [bim]] }, [], false)

    PLINK_EXTRACT(bbf_ch.join(GAWK.out.output))
}
