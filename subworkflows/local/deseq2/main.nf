/*
 * Differential expression with DESeq2: 
 *  > QUANT_QC:       quality-check the salmon quantifications with a PCA
 *  > DESEQ2_DDS:     import quant-files with tximport, filter genes and fit the
 *                     model
 *  > DESEQ2_RESULTS: extract every group-vs-reference contrast with 
 *                     apeglm-shrunken fold changes
 *
 * Subworkflow emits csv with genes and their statistical values for subsequent
 * analysis
 */

include { QUANT_QC          } from '../../../modules/local/quant_qc'
include { DESEQ2_DDS        } from '../../../modules/local/deseq2/dds'
include { DESEQ2_RESULTS    } from '../../../modules/local/deseq2/results'
include { DESEQ2_GENE_PLOTS } from '../../../modules/local/deseq2/gene_plots'

workflow DESEQ2 {

    take:
    ch_metadata // channel: path(metadata.csv)
    ch_tx2gene  // channel: path(tx2gene.tsv)
    ch_quant    // channel: path(sample dir with quant.sf), one item per sample

    main:
    ch_versions = channel.empty()

    // All samples go into one dataset: [ meta, metadata, tx2gene, [ dirs ] ]
    ch_input = ch_metadata
        .combine(ch_tx2gene)
        .combine(ch_quant.collect().map { dirs -> [dirs] })
        .map { metadata, tx2gene, dirs ->
            [[id: 'dds'], metadata, tx2gene, dirs]
        }

    QUANT_QC(ch_input)
    ch_versions = ch_versions.mix(QUANT_QC.out.versions.first())
    
    DESEQ2_DDS(ch_input)
    ch_versions = ch_versions.mix(DESEQ2_DDS.out.versions.first())

    DESEQ2_RESULTS(DESEQ2_DDS.out.rds)
    ch_versions = ch_versions.mix(DESEQ2_RESULTS.out.versions.first())

    ch_plot_input = DESEQ2_DDS.out.rds.join(DESEQ2_RESULTS.out.results)

    DESEQ2_GENE_PLOTS(ch_plot_input, ch_annotation)
    ch_versions = ch_versions.mix(DESEQ2_GENE_PLOTS.out.versions.first())

    emit:
    qc_pca     = QUANT_QC.out.pca             // channel: [ meta, pdf ]
    qc_samples = QUANT_QC.out.samples         // channel: [ meta, tsv ]
    rds        = DESEQ2_DDS.out.rds           // channel: [ meta, Rds ]
    results    = DESEQ2_RESULTS.out.results   // channel: [ meta, [ tsv ] ]
    summary    = DESEQ2_RESULTS.out.summary   // channel: [ meta, tsv ]
    gene_plots = DESEQ2_GENE_PLOTS.out.plots  // channel: [ meta, dir ]
    versions   = ch_versions                  // channel: [ versions.yml ]
}
