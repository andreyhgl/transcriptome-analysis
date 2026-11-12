#!/usr/bin/env nextflow

// ~~ Import processes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
include { SAMPLE_QC          } from './subworkflows/local/sample_qc'
include { GENE_ANNOTATION    } from './subworkflows/local/gene_annotation'
include { DESEQ2             } from './subworkflows/local/deseq2'
include { PATHWAY_ENRICHMENT } from './subworkflows/local/pathway_enrichment'
include { QUARTO_REPORT      } from './modules/local/quarto_report'
include { paramsSummary; completionSummary } from './subworkflows/local/utils_pipeline'

workflow {

    main:
    // pipeline parameters 
    log.info(paramsSummary(params))

    ch_metadata = channel.fromPath(params.metadata, checkIfExists: true)
    ch_tx2gene  = channel.fromPath(params.tx2gene, checkIfExists: true)
    ch_quarto_meta = file(params.quarto_meta_file, checkIfExists: true)
    ch_sample_qc_meta = file(params.sample_qc_meta_file, checkIfExists: true)
    ch_quarto_css  = file(params.quarto_css_file, checkIfExists: true)

    // One directory per sample; the folder name is the sample id
    ch_quant = channel
        .fromPath("${params.quant_dir}/*/quant.sf", checkIfExists: true)
        .map { quant -> quant.parent }

    if (!(params.species in ['human', 'mouse'])) {
        error "--species must be 'human' or 'mouse', not '${params.species}'"
    }

    // Gene annotation for the gene plots:
    if (params.gene_annotation) {
        // table included
        ch_annotation = channel
            .fromPath(params.gene_annotation, checkIfExists: true)
            .collect()
    }
    else {
        // fetched from Ensembl BioMart, w/ parameters:
        // params.species
        // params.ensembl_version
        GENE_ANNOTATION(params.species, params.ensembl_version)
        ch_annotation = GENE_ANNOTATION.out.annotation
    }

    if (!params.skip_sample_qc) {
        SAMPLE_QC (
            ch_quant,
            ch_metadata,
            ch_tx2gene,
            file(params.sample_qc_report_file, checkIfExists: true),
            ch_quarto_meta,
            file(params.sample_qc_dashboard_file, checkIfExists: true),
            file(params.sample_qc_dashboard_meta, checkIfExists: true),
            ch_quarto_css
        )
    }

    DESEQ2 ( ch_metadata, ch_tx2gene, ch_quant, ch_annotation )

    PATHWAY_ENRICHMENT ( DESEQ2.out.results, ch_annotation )

    // Pathway enrichment emits ORA and GSEA tables for every contrast
    ch_enrichment = PATHWAY_ENRICHMENT.out.ora_rds
        .mix(PATHWAY_ENRICHMENT.out.gsea_rds)
        .map { meta, rds -> rds }
        .toList()

    QUARTO_REPORT (
        DESEQ2.out.rds
            .join(DESEQ2.out.results)
            .join(DESEQ2.out.summary),
        ch_annotation,
        ch_enrichment,
        file(params.quarto_report_file, checkIfExists: true),
        ch_quarto_meta,
        ch_quarto_css
    )

    onComplete:
    log.info(completionSummary(params, workflow))
}