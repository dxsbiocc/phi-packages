#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored subworkflow at ../main.nf.
// Includes the subworkflow only — never another wrapper; the subworkflow
// itself references the already-migrated modules under resources/wrappers/modules.
// See docs/contracts/wrapper.md. Manually wired by scripts/wire-manual-subworkflows.py.
nextflow.enable.dsl = 2

include { DIA_PROTEOMICS_ANALYSIS } from '../main.nf'

params.ms_file = null
params.searchdb = null
params.expdesign = null
params.enzyme = "Trypsin"
params.fixed_mods = "Carbamidomethyl (C)"
params.variable_mods = "Oxidation (M)"
params.precursor_tolerance = 20
params.fragment_tolerance = 20
params.precursor_tolerance_unit = "ppm"
params.fragment_tolerance_unit = "ppm"
params.random_preanalysis = false
params.random_preanalysis_n = null
params.random_preanalysis_seed = null
params.scan_window = 0
params.mass_acc_automatic = true
params.scan_window_automatic = true
params.diann_debug = false
params.pg_level = "genes"
params.speclib = null
params.empirical_library = null
params.empirical_log = null
params.outdir = "results"

workflow {
    input_ch = Channel.value([[id: 'sample'], file(params.ms_file, checkIfExists: true), params.enzyme, params.fixed_mods, params.variable_mods, params.precursor_tolerance, params.fragment_tolerance, params.precursor_tolerance_unit, params.fragment_tolerance_unit])
    searchdb_ch = Channel.value([[id: 'searchdb'], file(params.searchdb, checkIfExists: true)])
    expdesign_ch = Channel.value([[id: 'experiment'], file(params.expdesign, checkIfExists: true)])
    speclib_ch = params.speclib ? Channel.value(file(params.speclib, checkIfExists: true)) : Channel.empty()
    empirical_ch = params.empirical_library ? Channel.value(file(params.empirical_library, checkIfExists: true)) : Channel.empty()
    empirical_log_ch = params.empirical_log ? Channel.value(file(params.empirical_log, checkIfExists: true)) : Channel.empty()

    DIA_PROTEOMICS_ANALYSIS(input_ch, searchdb_ch, expdesign_ch, [params.random_preanalysis, params.random_preanalysis_n, params.random_preanalysis_seed], params.scan_window, params.mass_acc_automatic, params.scan_window_automatic, params.diann_debug, params.pg_level, speclib_ch, empirical_ch, empirical_log_ch)
}
