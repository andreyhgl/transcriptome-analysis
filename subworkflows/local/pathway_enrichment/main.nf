//
// Pathway enrichment: over-representation of the significant genes of each
// DESeq2 contrast in GO, KEGG and Reactome, run as one task per contrast
//

include { PATHWAY_ENRICHMENT_ORA } from '../../../modules/local/pathway_enrichment/ora'

workflow PATHWAY_ENRICHMENT {

    take:
    ch_results    // channel: [ meta, [ tsv ] ], one table per contrast
    ch_annotation // channel: path(annotation) or [], as a value channel

    main:
    ch_versions = channel.empty()

    // One item per contrast. The contrast is the file name without the
    // results prefix: DDS_summer_2022_vs_control.tsv => summer_2022_vs_control
    ch_contrasts = ch_results
        .flatMap { meta, tsvs ->
            [tsvs].flatten().collect { tsv ->
                def contrast = tsv.baseName.replaceFirst("^${meta.id}_", '')
                [meta + [id: meta.id + '_' + contrast, contrast: contrast], tsv]
            }
        }

    PATHWAY_ENRICHMENT_ORA(ch_contrasts, ch_annotation)
    ch_versions = ch_versions.mix(PATHWAY_ENRICHMENT_ORA.out.versions.first())

    emit:
    tsv      = PATHWAY_ENRICHMENT_ORA.out.tsv  // channel: [ meta, tsv ]
    xlsx     = PATHWAY_ENRICHMENT_ORA.out.xlsx // channel: [ meta, xlsx ]
    rds      = PATHWAY_ENRICHMENT_ORA.out.rds  // channel: [ meta, Rds ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
