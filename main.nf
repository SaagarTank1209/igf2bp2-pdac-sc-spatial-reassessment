include { ENV_SMOKE_TEST; GPU_SMOKE_TEST } from './modules/local/env_smoke_test.nf'
include { LOAD_QC } from './modules/local/load_qc.nf'
include { ANNOTATE_SC } from './modules/local/annotate_sc.nf'
include { STATE_PROGRAM } from './modules/local/state_program.nf'
include { SPATIAL_TESTS } from './modules/local/spatial_tests.nf'
include { POSTHOC_DIAGNOSTICS } from './modules/local/posthoc_diagnostics.nf'
include { FIGURES } from './modules/local/figures.nf'

params.step    = 'smoke'
params.datadir = "/scratch/${System.getenv('USER')}/pdac_spatial"

workflow {
    if (params.step == 'smoke') {
        ENV_SMOKE_TEST('main')
        if (params.smoke_gpu) {
            GPU_SMOKE_TEST()
        }
    }
    else if (params.step == 'load_qc') {
        LOAD_QC("${params.datadir}/counts", "${params.datadir}/raw")
    }
    else if (params.step == 'annotate') {
        ANNOTATE_SC(file("${params.outdir}/01_load_qc"))
    }
    else if (params.step == 'state') {
        STATE_PROGRAM(file("${params.outdir}/02_annotation"))
    }
    else if (params.step == 'spatial') {
        SPATIAL_TESTS(file("${params.outdir}/01_load_qc"))
    }
    else if (params.step == 'posthoc') {
        POSTHOC_DIAGNOSTICS(file("${params.outdir}/01_load_qc"))
    }
    else if (params.step == 'figures') {
        FIGURES(file("${params.outdir}/01_load_qc"), file("${params.outdir}/02_annotation"),
                file("${params.outdir}/03_state"), file("${params.outdir}/04_spatial"),
                file("${params.outdir}/05_posthoc"))
    }
    else if (params.step == 'all') {
        LOAD_QC("${params.datadir}/counts", "${params.datadir}/raw")
        ANNOTATE_SC(LOAD_QC.out)
        STATE_PROGRAM(ANNOTATE_SC.out)
        SPATIAL_TESTS(LOAD_QC.out)
        POSTHOC_DIAGNOSTICS(LOAD_QC.out)
        FIGURES(LOAD_QC.out, ANNOTATE_SC.out, STATE_PROGRAM.out, SPATIAL_TESTS.out, POSTHOC_DIAGNOSTICS.out)
    }
    else {
        error "Unknown --step ${params.step}"
    }
}