process ENV_SMOKE_TEST {
    tag "${env_name}"
    publishDir "${params.outdir}/smoke_test", mode: 'copy'

    input:
    val env_name

    output:
    path "smoke_${env_name}.txt"

    script:
    """
    python -c "import sys, importlib.metadata as m; print('python', sys.version.split()[0]); [print(p, m.version(p)) for p in ('numpy','scanpy','anndata')]" > smoke_${env_name}.txt
    hostname >> smoke_${env_name}.txt
    """
}

process GPU_SMOKE_TEST {
    label 'gpu'
    publishDir "${params.outdir}/smoke_test", mode: 'copy'

    output:
    path "smoke_gpu.txt"

    script:
    """
    python -c "import torch, cell2location, scvi; print('imports OK'); print('torch', torch.__version__); ok=torch.cuda.is_available(); print('gpu visible', ok); print(torch.cuda.get_device_name(0) if ok else 'none')" > smoke_gpu.txt
    hostname >> smoke_gpu.txt
    """
}