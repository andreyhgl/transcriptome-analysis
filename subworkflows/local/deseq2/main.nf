/*
 * Differential expression with DESeq2: import salmon quantifications with
 * tximport, filter genes and fit the model (DESEQ2_DDS), then extract every
 * group-vs-reference contrast with apeglm-shrunken fold changes
 * (DESEQ2_RESULTS)
*/

include { DESEQ2_DDS     } from '../../../modules/local/deseq2/dds'
//include { DESEQ2_RESULTS } from '../../../modules/local/deseq2/results'

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

    DESEQ2_DDS(ch_input)
    ch_versions = ch_versions.mix(DESEQ2_DDS.out.versions.first())

    //DESEQ2_RESULTS(DESEQ2_DDS.out.rds)
    //ch_versions = ch_versions.mix(DESEQ2_RESULTS.out.versions.first())

    emit:
    rds      = DESEQ2_DDS.out.rds          // channel: [ meta, Rds ]
    //results  = DESEQ2_RESULTS.out.results  // channel: [ meta, [ tsv ] ]
    //summary  = DESEQ2_RESULTS.out.summary  // channel: [ meta, tsv ]
    versions = ch_versions                 // channel: [ versions.yml ]
}
