process SAMPLE_QC_REPORT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // TODO: set a container with quarto, ggplot2, patchwork, ggrepel, scales,
    // knitr, rmarkdown and reactable (the QUARTO_REPORT image has them all)
    // container "docker://ghcr.io/<org>/quarto-report:<tag>"

    input:
    tuple val(meta), path(rds)
    path quarto_report_file
    // Not _quarto.yml: Quarto would also load that as the project config
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
    # The staged files are symlinks to the pipeline's assets/. Quarto resolves
    # symlinks and would write its output next to the originals, so replace
    # them with real copies
    for f in ${quarto_report_file} ${quarto_meta_file} ${quarto_css_file}; do
        cp --remove-destination "\$(readlink -f "\$f")" "\$f"
    done

    # Quarto (deno) makes a temp folder at start-up. The host's TMPDIR (e.g. a
    # job scratch folder) is passed to the container but not mounted in it,
    # so use a folder in the task directory, which always is
    export TMPDIR="\$PWD/.tmp"
    mkdir -p "\$TMPDIR"

    # Quarto and deno need a writable cache; HOME may not exist in containers
    export XDG_CACHE_HOME="\$PWD/.cache"
    export XDG_DATA_HOME="\$PWD/.local/share"
    [ -w "\${HOME:-/nonexistent}" ] || export HOME="\$PWD"

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
