#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { BCFTOOLS_PLUGINFILLTAGS } from '../main.nf'


params.vcf        = null
params.csi        = null
params.outdir     = null

workflow {
    vcf_ch = Channel.value([[id: 'result', region: 'chr22:16570065-16609999'], file(params.vcf, checkIfExists: true), file(params.csi, checkIfExists: true)])

    BCFTOOLS_PLUGINFILLTAGS(vcf_ch, [], [], [])
}
