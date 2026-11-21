#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { CELLRANGER_MULTI } from '../main.nf'

params.gex_fastqs = null
params.gex_reference = null
params.vdj_fastqs = null
params.ab_fastqs = null
params.beam_fastqs = null
params.cmo_fastqs = null
params.crispr_fastqs = null
params.gex_frna_probeset = null
params.gex_targetpanel = null
params.vdj_reference = null
params.vdj_primer_index = null
params.fb_reference = null
params.beam_antigen_panel = null
params.beam_control_panel = null
params.cmo_reference = null
params.cmo_barcodes = null
params.cmo_barcode_assignment = null
params.frna_sampleinfo = null
params.ocm_barcodes = null
params.skip_renaming = null
params.args = null
params.outdir = null

// REVIEW: only the GEX-only and full-library combinations are exercised upstream;
// other partial combinations are not runtime-verified.
workflow {
    def fq  = { p -> p ? file(p, checkIfExists: true).listFiles().findAll { it.isFile() }.toList() : [] }
    def ref = { p -> p ? file(p, checkIfExists: true) : [] }

    CELLRANGER_MULTI(
        [id: 'multi'],
        Channel.value([[id: 'gex'],    fq(params.gex_fastqs),    '']),
        Channel.value([[id: 'vdj'],    fq(params.vdj_fastqs),    '']),
        Channel.value([[id: 'ab'],     fq(params.ab_fastqs),     '']),
        Channel.value([[id: 'beam'],   fq(params.beam_fastqs),   '']),
        Channel.value([[id: 'cmo'],    fq(params.cmo_fastqs),    '']),
        Channel.value([[id: 'crispr'], fq(params.crispr_fastqs), '']),
        Channel.value([[id: 'gex_reference'], ref(params.gex_reference)]),
        ref(params.gex_frna_probeset),
        ref(params.gex_targetpanel),
        ref(params.vdj_reference),
        ref(params.vdj_primer_index),
        ref(params.fb_reference),
        ref(params.beam_antigen_panel),
        ref(params.beam_control_panel),
        ref(params.cmo_reference),
        ref(params.cmo_barcodes),
        ref(params.cmo_barcode_assignment),
        ref(params.frna_sampleinfo),
        ref(params.ocm_barcodes),
        params.skip_renaming
    )
}
