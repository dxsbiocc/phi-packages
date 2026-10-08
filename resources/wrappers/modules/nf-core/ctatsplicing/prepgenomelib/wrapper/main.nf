#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
nextflow.enable.dsl = 2

include { CTATSPLICING_PREPGENOMELIB } from '../main.nf'
include { STARFUSION_BUILD } from '../../../starfusion/build/main.nf'

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
params.outdir     = null

workflow {
    STARFUSION_BUILD(Channel.value([[id: 'minigenome_fasta'], file(params.fasta, checkIfExists: true)]),
                     Channel.value([[id: 'minigenome_gtf'], file(params.gtf, checkIfExists: true)]),
                     file(params.fusion_lib, checkIfExists: true), 'homo_sapiens', file(params.pfam, checkIfExists: true),
                     [file(params.dfam_hmm, checkIfExists: true), file(params.dfam_h3f, checkIfExists: true), file(params.dfam_h3i, checkIfExists: true),
                      file(params.dfam_h3m, checkIfExists: true), file(params.dfam_h3p, checkIfExists: true)],
                     'https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/AnnotFilterRule.pm')

    CTATSPLICING_PREPGENOMELIB(STARFUSION_BUILD.out.reference, [file(params.cancer_introns, checkIfExists: true)])
}
