process QUARTO_REPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' ?
        'docker://ghcr.io/andreyhgl/deseq2-env:main' :
        'ghcr.io/andreyhgl/deseq2-env:main'}"

    input:
    tuple val(meta), path(rds), path(results, stageAs: 'deseq2/*'), path(summary)
    path annotation
    path enrichment, stageAs: 'enrichment/*'
    path quarto_report_file
    path quarto_meta_file, stageAs: 'quarto_meta.yml'
    path quarto_css_file

    output:
    tuple val(meta), path("${prefix}.html"), emit: html
    tuple val(meta), path("*.{pdf,txt}"),    emit: figures, optional: true
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def anno = annotation ? "-P annotation:${annotation}" : ''
    prefix = task.ext.prefix ?: "${meta.id}_report"
    """
    # TMPDIR is used by Quarto.makeTempDirSync(), set manually to avoid
    # HPC issues
    export TMPDIR="./cache/"

    # Set environment variables needed for Quarto rendering
    export XDG_DATA_HOME="./cache/"
    export XDG_CACHE_HOME="./cache/"
    export XDG_RUNTIME_DIR="./cache/"

    # Everything is staged next to the .qmd, so the report finds its inputs
    # and pandoc its CSS and figures relative to the task directory
    quarto render ${quarto_report_file} \\
        --metadata-file ${quarto_meta_file} \\
        --css ${quarto_css_file} \\
        --output ${prefix}.html \\
        -P data_dir:"\$PWD" \\
        -P results_dir:deseq2 \\
        -P dds:${rds} \\
        -P summary:${summary} \\
        -P enrichment_dir:enrichment \\
        ${anno} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-deseq2: \$(Rscript -e 'cat(as.character(packageVersion("DESeq2")))')
        r-ggvenn: \$(Rscript -e 'cat(as.character(packageVersion("ggvenn")))')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_report"
    """
    touch ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: stub
        r-base: stub
        bioconductor-deseq2: stub
        r-ggvenn: stub
    END_VERSIONS
    """
}
