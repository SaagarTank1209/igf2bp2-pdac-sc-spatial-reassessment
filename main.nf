include { ENV_SMOKE_TEST; GPU_SMOKE_TEST } from './modules/local/env_smoke_test.nf'
include { LOAD_QC } from './modules/local/load_qc.nf'
include { ANNOTATE_SC } from './modules/local/annotate_sc.nf'
include { STATE_PROGRAM } from './modules/local/state_program.nf'

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
        ANNOTATE_SC("${params.outdir}/01_load_qc/h5ad")
    }
    else if (params.step == 'state') {
        STATE_PROGRAM("${params.outdir}/02_annotation")
    }
    else {
        error "Unknown --step ${params.step}"
    }
}