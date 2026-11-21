#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { MULTIQCSAV } from '../main.nf'

params.run_dir = null
params.extra_multiqc_files = null
params.multiqc_config = null
params.multiqc_logo = null
params.replace_names = null
params.sample_names = null
params.args = null
params.outdir = null

workflow {
    def run = file(params.run_dir, checkIfExists: true)
    def xml = []
    def interop = []
    run.eachFileRecurse { f ->
        if (f.fileName.toString() == 'RunInfo.xml') xml << f
        if (f.parent?.fileName?.toString() == 'InterOp' && f.fileName.toString().endsWith('.bin')) interop << f
    }
    def extra = params.extra_multiqc_files
        ? file(params.extra_multiqc_files, checkIfExists: true).listFiles().findAll { it.isFile() }.toList()
        : []
    def opt = { p -> p ? file(p, checkIfExists: true) : [] }

    MULTIQCSAV(Channel.value([
        [id: 'sav'], xml, interop, extra,
        opt(params.multiqc_config), opt(params.multiqc_logo),
        opt(params.replace_names), opt(params.sample_names)
    ]))
}
