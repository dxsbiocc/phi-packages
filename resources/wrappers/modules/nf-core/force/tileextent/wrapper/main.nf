#!/usr/bin/env nextflow
// Thin agent-facing adapter over the vendored module at ../main.nf.
// See docs/design/phi-wrapper-agent-composition-design.md section 1.
// Manually wired: the module's staging pattern is too exotic for the generator.
nextflow.enable.dsl = 2

include { FORCE_TILEEXTENT } from '../main.nf'

params.aoi = null
params.datacube_definition = null
params.shapefile_dbf = null
params.shapefile_prj = null
params.shapefile_shx = null
params.outdir = null

workflow {
    def shp = { p -> p ? file(p, checkIfExists: true) : [] }

    FORCE_TILEEXTENT(
        file(params.aoi, checkIfExists: true),
        file(params.datacube_definition, checkIfExists: true),
        shp(params.shapefile_dbf), shp(params.shapefile_prj), shp(params.shapefile_shx)
    )
}
