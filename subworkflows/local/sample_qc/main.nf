/*
 * Sample quality control, independent of the experimental design: per-sample
 * metrics, PCA, clustering and metadata checks on the salmon quantifications
 * (SAMPLE_QC_METRICS), summarised in a log and a Quarto report
 * (SAMPLE_QC_REPORT)
*/

include { SAMPLE_QC_METRICS } from '../../../modules/local/sample_qc/metrics'
include { SAMPLE_QC_REPORT  } from '../../../modules/local/sample_qc/report'

workflow SAMPLE_QC {

    take:
    ch_samples         // channel: path(sample dir with quant.sf), one item per sample
    ch_metadata        // channel: path(metadata.csv)
    ch_tx2gene         // channel: path(tx2gene.tsv)
    quarto_report_file // path: sample_qc.qmd
    quarto_meta_file   // path: Quarto metadata (format, theme, execute options)
    quarto_css_file    // path: stylesheet

    main:
    ch_versions = channel.empty()

    // All samples go into one check: [ meta, metadata, tx2gene, [ dirs ] ]
    ch_input = ch_metadata
        .combine(ch_tx2gene)
        .combine(ch_samples.collect().map { dirs -> [dirs] })
        .map { metadata, tx2gene, dirs ->
            [[id: 'sample_qc'], metadata, tx2gene, dirs]
        }

    SAMPLE_QC_METRICS ( ch_input )
    ch_versions = ch_versions.mix ( SAMPLE_QC_METRICS.out.versions.first() )

    SAMPLE_QC_REPORT (
        SAMPLE_QC_METRICS.out.rds,
        quarto_report_file,
        quarto_meta_file,
        quarto_css_file
    )
    ch_versions = ch_versions.mix ( SAMPLE_QC_REPORT.out.versions.first() )

    emit:
    qc_log   = SAMPLE_QC_METRICS.out.log     // channel: [ meta, log ]
    report   = SAMPLE_QC_REPORT.out.html     // channel: [ meta, html ]
    samples  = SAMPLE_QC_METRICS.out.samples // channel: [ meta, tsv ]
    rds      = SAMPLE_QC_METRICS.out.rds     // channel: [ meta, Rds ]
    versions = ch_versions                   // channel: [ versions.yml ]
}