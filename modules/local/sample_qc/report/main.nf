process SAMPLE_QC_REPORT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' ?
        'docker://ghcr.io/andreyhgl/deseq2-env:f070b01' :
        'ghcr.io/andreyhgl/deseq2-env:f070b01'}"

    input:
    tuple val(meta), path(rds)
    path quarto_report_file
    path quarto_meta_file, stageAs: 'quarto_meta.yml'
    path quarto_css_file

    output:
    tuple val(meta), path("${prefix}.html"), emit: html
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_report"
    """
    # Quarto writes caches and temp files here; keep them inside the task dir
    # to avoid HPC issues (TMPDIR is used by Quarto.makeTempDirSync())
    export XDG_DATA_HOME="./cache/"
    export XDG_CACHE_HOME="./cache/"
    export XDG_RUNTIME_DIR="./cache/"
    export TMPDIR="./cache/"

    quarto render ${quarto_report_file} \\
        --metadata-file ${quarto_meta_file} \\
        --css ${quarto_css_file} \\
        -P qc:${rds} \\
        ${args} \\
        --output ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        r-ggplot2: \$(Rscript -e 'cat(as.character(packageVersion("ggplot2")))')
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
        r-ggplot2: stub
    END_VERSIONS
    """
}
