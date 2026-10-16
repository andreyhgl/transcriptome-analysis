//
// Gene annotation table for the gene plots, fetched from Ensembl BioMart and
// pinned to one Ensembl release (BIOMART_ANNOTATION). Columns: gene_id,
// gene_name, chr, start, end, gene_type, description.
//
// Needs internet access to ensembl.org from the node that runs the task.
//

include { BIOMART_ANNOTATION } from '../../../modules/local/biomart/annotation'

workflow GENE_ANNOTATION {

    take:
    organism // value: species, e.g. 'mus_musculus' or 'Mus musculus'
    version  // value: Ensembl release, e.g. 113; null for the current release

    main:
    ch_versions = channel.empty()

    if (!version) {
        log.warn "GENE_ANNOTATION: no Ensembl version given, using the " +
            "current release. Gene IDs may not match the salmon index."
    }
    else if (!(version.toString() ==~ /\d+/)) {
        error "GENE_ANNOTATION: version must be an Ensembl release number " +
            "such as 113, got '${version}'"
    }

    def org_id = organism.toString().toLowerCase().replaceAll(/[^a-z0-9]+/, '_')
    def meta = [id: "ensembl_${org_id}_v${version ?: 'current'}"]

    BIOMART_ANNOTATION(channel.of([meta, organism, version]))
    ch_versions = ch_versions.mix(BIOMART_ANNOTATION.out.versions)

    // A value channel, so every process that reads it gets the same table
    ch_annotation = BIOMART_ANNOTATION.out.tsv
        .map { _meta, tsv -> tsv }
        .first()

    emit:
    annotation = ch_annotation // channel: path(tsv)
    versions   = ch_versions   // channel: path(versions.yml)
}
