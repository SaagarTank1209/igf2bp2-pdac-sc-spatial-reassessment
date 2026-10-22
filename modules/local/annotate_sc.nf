process ANNOTATE_SC {
    cpus 4
    memory '16 GB'
    time '1h'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val h5ad_dir

    output:
    path "02_annotation"

    script:
    """
    python ${projectDir}/notebooks/02_sc_annotation.py \\
        --h5ad_dir ${h5ad_dir} \\
        --outdir 02_annotation
    """
}