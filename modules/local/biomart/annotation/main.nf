process BIOMART_ANNOTATION {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/bioconductor-biomart:2.68.0--r45hdfd78af_0' :
        'quay.io/biocontainers/bioconductor-biomart:2.68.0--r45hdfd78af_0'}"

    input:
    tuple val(meta), val(organism), val(version)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    tuple val(meta), path("${prefix}.log"), emit: log
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def release = version ? "--version ${version}" : ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    biomart_annotation.R \\
        --organism '${organism}' \\
        ${release} \\
        --filename ${prefix} \\
        ${args} \\
        2>&1 | tee ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-biomart: \$(Rscript -e 'cat(as.character(packageVersion("biomaRt")))')
        ensembl: ${version ?: 'current'}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'gene_id\\tgene_name\\tchr\\tstart\\tend\\tgene_type\\tdescription\\n' > ${prefix}.tsv
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: stub
        bioconductor-biomart: stub
        ensembl: ${version ?: 'current'}
    END_VERSIONS
    """
}