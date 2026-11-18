process GO_THEMES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' ?
        'docker://ghcr.io/andreyhgl/enrichment-env:c237546' :
        'ghcr.io/andreyhgl/enrichment-env:c237546'}"

    input:
    // every contrast's pathway_enrichment_ora.R Rds, collected
    tuple val(meta), path(rds, stageAs: 'enrichment/*')

    output:
    tuple val(meta), path("${prefix}.Rds"), emit: rds
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    tuple val(meta), path("${prefix}.log"), emit: log
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_go_themes"
    """
    pathway_enrichment_go_themes.R \\
        --enrichment_dir enrichment \\
        --prefix ${prefix} \\
        ${args} \\
        2>&1 | tee ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-gosemsim: \$(Rscript -e 'cat(as.character(packageVersion("GOSemSim")))')
        bioconductor-go.db: \$(Rscript -e 'cat(as.character(packageVersion("GO.db")))')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_go_themes"
    """
    touch ${prefix}.Rds ${prefix}.tsv ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: stub
        bioconductor-gosemsim: stub
        bioconductor-go.db: stub
    END_VERSIONS
    """
}