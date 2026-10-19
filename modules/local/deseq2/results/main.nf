process DESEQ2_RESULTS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container "${workflow.containerEngine == 'singularity' ?
        'docker://ghcr.io/andreyhgl/deseq2-env:main' :
        'ghcr.io/andreyhgl/deseq2-env:main'}"

    input:
    tuple val(meta), path(rds)

    output:
    tuple val(meta), path("${prefix}_*_vs_*.tsv"),  emit: results
    tuple val(meta), path("${prefix}_summary.tsv"), emit: summary
    tuple val(meta), path("${prefix}_results.log"), emit: log
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -o pipefail

    deseq2_results.R \\
        --dds ${rds} \\
        ${args} \\
        2>&1 | tee ${prefix}_results.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-deseq2: \$(Rscript -e 'cat(as.character(packageVersion("DESeq2")))')
        bioconductor-apeglm: \$(Rscript -e 'cat(as.character(packageVersion("apeglm")))')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_groupA_vs_control.tsv
    touch ${prefix}_summary.tsv ${prefix}_results.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: stub
        bioconductor-deseq2: stub
        bioconductor-apeglm: stub
    END_VERSIONS
    """
}
