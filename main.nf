include { ENV_SMOKE_TEST; GPU_SMOKE_TEST } from './modules/local/env_smoke_test.nf'

workflow {
    ENV_SMOKE_TEST('main')
    if (params.smoke_gpu) {
        GPU_SMOKE_TEST()
    }
}