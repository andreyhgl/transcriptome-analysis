/*
 * Pathway enrichment of each DESeq2 contrast in GO, KEGG and Reactome, run
 * as one task per contrast and method:
 *   - ORA: over-representation of the significant genes (Up, Down, Both)
 *   - GSEA: all genes ranked by one score, no significance cutoff
 * Either method can be switched off with params.enrichment_methods.
 * The enriched GO terms of all ORA tasks are then grouped into themes
 * (GO_THEMES), once, so a theme is the same in every contrast.
 */

include { PATHWAY_ENRICHMENT_ORA  } from '../../../modules/local/pathway_enrichment/ora'
include { PATHWAY_ENRICHMENT_GSEA } from '../../../modules/local/pathway_enrichment/gsea'
include { GO_THEMES               } from '../../../modules/local/pathway_enrichment/go_themes'

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

    PATHWAY_ENRICHMENT_ORA ( ch_contrasts, ch_annotation )
    ch_versions = ch_versions.mix(PATHWAY_ENRICHMENT_ORA.out.versions.first())

    PATHWAY_ENRICHMENT_GSEA ( ch_contrasts, ch_annotation )
    ch_versions = ch_versions.mix(PATHWAY_ENRICHMENT_GSEA.out.versions.first())

    // All ORA results together, with the meta of the results input:
    // [ [id: 'DDS'], [ Rds, ... ] ]. Without ORA results nothing is emitted,
    // so GO_THEMES doesn't run
    ch_themes_input = ch_results
        .map { meta, _tsvs -> meta }
        .combine(
            PATHWAY_ENRICHMENT_ORA.out.rds
                .map { _meta, rds -> rds }
                .collect()
                .map { files -> [files] }
        )

    GO_THEMES ( ch_themes_input )
    ch_versions = ch_versions.mix ( GO_THEMES.out.versions )

    emit:
    ora_tsv   = PATHWAY_ENRICHMENT_ORA.out.tsv   // channel: [ meta, tsv ]
    ora_xlsx  = PATHWAY_ENRICHMENT_ORA.out.xlsx  // channel: [ meta, xlsx ]
    ora_rds   = PATHWAY_ENRICHMENT_ORA.out.rds   // channel: [ meta, Rds ]
    gsea_tsv  = PATHWAY_ENRICHMENT_GSEA.out.tsv  // channel: [ meta, tsv ]
    gsea_xlsx = PATHWAY_ENRICHMENT_GSEA.out.xlsx // channel: [ meta, xlsx ]
    gsea_rds  = PATHWAY_ENRICHMENT_GSEA.out.rds  // channel: [ meta, Rds ]
    go_themes = GO_THEMES.out.rds                // channel: [ meta, Rds ]
    versions  = ch_versions                      // channel: [ versions.yml ]
}