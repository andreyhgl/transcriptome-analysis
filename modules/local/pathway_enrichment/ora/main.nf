process PATHWAY_ENRICHMENT_ORA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/../environment.yml"
    container "${workflow.containerEngine == 'singularity' ?
        'docker://ghcr.io/andreyhgl/enrichment-env:main' :
        'ghcr.io/andreyhgl/enrichment-env:main'}"

    input:
    tuple val(meta), path(results)
    path annotation

    output:
    tuple val(meta), path("${prefix}.tsv"),  emit: tsv
    tuple val(meta), path("${prefix}.xlsx"), emit: xlsx
    tuple val(meta), path("${prefix}.Rds"),  emit: rds
    tuple val(meta), path("${prefix}.log"),  emit: log
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def anno = annotation ? "--annotation ${annotation}" : ''
    prefix = task.ext.prefix ?: "${meta.id}_ora"
    """
    set -o pipefail

    pathway_enrichment.R \\
        --results ${results} \\
        --contrast ${meta.contrast} \\
        --prefix ${prefix} \\
        ${anno} \\
        ${args} \\
        2>&1 | tee ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-clusterprofiler: \$(Rscript -e 'cat(as.character(packageVersion("clusterProfiler")))')
        bioconductor-reactomepa: \$(Rscript -e 'cat(as.character(packageVersion("ReactomePA")))')
        r-openxlsx: \$(Rscript -e 'cat(as.character(packageVersion("openxlsx")))')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_ora"
    """
    touch ${prefix}.tsv ${prefix}.xlsx ${prefix}.Rds ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: stub
        bioconductor-clusterprofiler: stub
        bioconductor-reactomepa: stub
        r-openxlsx: stub
    END_VERSIONS
    """
}