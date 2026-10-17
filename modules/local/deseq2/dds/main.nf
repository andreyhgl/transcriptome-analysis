process DESEQ2_DDS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container "${workflow.containerEngine == 'singularity' ?
        'docker://ghcr.io/andreyhgl/deseq2-env:main' :
        'ghcr.io/andreyhgl/deseq2-env:main'}"

    input:
    tuple val(meta), path(metadata), path(tx2gene), path(quant, stageAs: 'quant/*')

    output:
    tuple val(meta), path("${prefix}.Rds"), emit: rds
    tuple val(meta), path("${prefix}.log"), emit: log
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -o pipefail

    deseq2_dds.R \\
        --filename ${prefix} \\
        --metadata ${metadata} \\
        --tx2gene ${tx2gene} \\
        --quant_dir quant \\
        ${args} \\
        2>&1 | tee ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-deseq2: \$(Rscript -e 'cat(as.character(packageVersion("DESeq2")))')
        bioconductor-tximport: \$(Rscript -e 'cat(as.character(packageVersion("tximport")))')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.Rds ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: stub
        bioconductor-deseq2: stub
        bioconductor-tximport: stub
    END_VERSIONS
    """
}