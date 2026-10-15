process QUANT_QC {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' ?
        'docker://ghcr.io/andreyhgl/deseq2:main' :
        'ghcr.io/andreyhgl/deseq2:main'}"

    input:
    tuple val(meta), path(metadata), path(tx2gene), path(quant, stageAs: 'quant/*')

    output:
    tuple val(meta), path("${prefix}_pca.pdf"),     emit: pca
    tuple val(meta), path("${prefix}_samples.tsv"), emit: samples
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    quant_qc.R \\
        --metadata ${metadata} \\
        --tx2gene ${tx2gene} \\
        --quant_dir quant \\
        --filename ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-deseq2: \$(Rscript -e 'cat(as.character(packageVersion("DESeq2")))')
        r-ggplot2: \$(Rscript -e 'cat(as.character(packageVersion("ggplot2")))')
        r-patchwork: \$(Rscript -e 'cat(as.character(packageVersion("patchwork")))')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_pca.pdf ${prefix}_samples.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: stub
        bioconductor-deseq2: stub
        r-ggplot2: stub
        r-patchwork: stub
    END_VERSIONS
    """
}