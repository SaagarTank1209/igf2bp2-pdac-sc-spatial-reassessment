include { ENV_SMOKE_TEST; GPU_SMOKE_TEST } from './modules/local/env_smoke_test.nf'
include { LOAD_QC } from './modules/local/load_qc.nf'

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
    else {
        error "Unknown --step ${params.step}"
    }
}