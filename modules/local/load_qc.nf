process LOAD_QC {
    cpus 2
    memory '16 GB'
    time '1h'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val counts_dir
    val raw_dir

    output:
    path "01_load_qc"

    script:
    """
    python ${projectDir}/notebooks/01_load_qc.py \\
        --counts_dir ${counts_dir} \\
        --raw_dir ${raw_dir} \\
        --outdir 01_load_qc
    """
}