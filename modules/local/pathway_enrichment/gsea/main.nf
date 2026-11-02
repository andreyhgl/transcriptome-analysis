process PATHWAY_ENRICHMENT_GSEA {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container "${workflow.containerEngine == 'singularity' ?
        'docker://ghcr.io/andreyhgl/enrichment-env:main' :
        'ghcr.io/andreyhgl/enrichment-env:main'}"

    input:
    tuple val(meta), path(results)
    path annotation // optional: [] when no annotation table is given

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
    prefix = task.ext.prefix ?: "${meta.id}_gsea"
    """
    pathway_enrichment_gsea.R \\
        --results ${results} \\
        --contrast ${meta.contrast} \\
        --prefix ${prefix} \\
        --cores ${task.cpus} \\
        ${anno} \\
        ${args} \\
        2>&1 | tee ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-clusterprofiler: \$(Rscript -e 'cat(as.character(packageVersion("clusterProfiler")))')
        bioconductor-reactomepa: \$(Rscript -e 'cat(as.character(packageVersion("ReactomePA")))')
        bioconductor-fgsea: \$(Rscript -e 'cat(as.character(packageVersion("fgsea")))')
        r-openxlsx: \$(Rscript -e 'cat(as.character(packageVersion("openxlsx")))')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_gsea"
    """
    touch ${prefix}.tsv ${prefix}.xlsx ${prefix}.Rds ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: stub
        bioconductor-clusterprofiler: stub
        bioconductor-reactomepa: stub
        bioconductor-fgsea: stub
        r-openxlsx: stub
    END_VERSIONS
    """
}