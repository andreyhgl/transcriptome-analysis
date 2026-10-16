process DESEQ2_GENE_PLOTS {
    tag "${meta.id}"
    label 'process_parallel'

    conda "${moduleDir}/../environment.yml"
    container "${workflow.containerEngine == 'singularity' ?
        'docker://ghcr.io/andreyhgl/deseq2:main' :
        'ghcr.io/andreyhgl/deseq2:main'}"

    input:
    tuple val(meta), path(rds), path(results, stageAs: 'results/*')
    path annotation // optional: [] when no annotation table is given

    output:
    tuple val(meta), path("${prefix}_gene_plots.pdf"), emit: plots
    tuple val(meta), path("${prefix}_gene_plots.tsv"), emit: genes
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def anno = annotation ? "--annotation ${annotation}" : ''
    def batch_size = task.ext.batch_size ?: 500
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # 1. Render one PDF per gene, in parallel, into ${prefix}_gene_plots/
    deseq2_gene_plots.R \\
        --dds ${rds} \\
        --results_dir results \\
        --filename ${prefix} \\
        --cores ${task.cpus} \\
        ${anno} \\
        ${args}

    # 2. Merge them into one PDF in rank order (sort -V: 2 before 10). qpdf
    #    merges in batches so the argument list never gets too long.
    find ${prefix}_gene_plots -name '*.pdf' | sort -V > pages.txt

    mkdir batches
    split -l ${batch_size} pages.txt batches/list_

    batch_pdfs=()
    for list in batches/list_*; do
        mapfile -t pages < "\$list"
        echo "Merging \$(basename "\$list"): \${#pages[@]} pages"

        qpdf --empty --pages "\${pages[@]}" -- "\${list}.pdf"
        batch_pdfs+=("\${list}.pdf")
    done

    qpdf --empty --pages "\${batch_pdfs[@]}" -- ${prefix}_gene_plots.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-deseq2: \$(Rscript -e 'cat(as.character(packageVersion("DESeq2")))')
        r-ggplot2: \$(Rscript -e 'cat(as.character(packageVersion("ggplot2")))')
        r-ggsignif: \$(Rscript -e 'cat(as.character(packageVersion("ggsignif")))')
        r-ggrepel: \$(Rscript -e 'cat(as.character(packageVersion("ggrepel")))')
        qpdf: \$(qpdf --version | head -n 1 | sed 's/^qpdf version //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_gene_plots.pdf ${prefix}_gene_plots.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: stub
        bioconductor-deseq2: stub
        r-ggplot2: stub
        r-ggsignif: stub
        r-ggrepel: stub
        qpdf: stub
    END_VERSIONS
    """
}