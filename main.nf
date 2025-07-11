#!/usr/bin/env Nextflow

/*
 * Github: https://github.com/andreyhgl/transcriptome-analysis
 */

// ~~ Import processes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
//include { PCA_PLOTS           } from './modules/pca_plots.nf'
include { ENSEMBL         } from './modules/ensembl/main'
include { EDGER_DGELIST   } from './modules/edgeR/DGEList/main'
include { EDGER_QC_PLOTS  } from './modules/edgeR/QC_plots/main'
include { DIFF_EXPRESSION } from './modules/diff_expression/main'
//include { LONGTABLE       } from './modules/longtable.nf'
//include { DMG_TABLE           } from './modules/diffmeth_tables.nf'
//include { GENE_ONTOLOGY       } from './modules/gene_ontology.nf'
//include { SUPPLEMENTARY_EXCEL } from './modules/supplementary.nf'
//include { SUPPLEMENTARY_PLOTS } from './modules/supplementary.nf'
//include { WRAPPER             } from './modules/wrapper.nf'
//include { REPORT              } from './modules/report.nf'

// ~~ Channels ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *

log.info \
  """
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Outdir                    : ${params.outdir}
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  """
  .stripIndent(true)

// ~~ Workflow ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *

workflow {

  ENSEMBL( params.reference_genome )

  ch_ensembl_dataset  = ENSEMBL.out.ENSEMBL_DATASET

  ch_quantfiles       = Channel.fromPath(params.quant_path)
  ch_metadata         = Channel.fromPath(params.metadata)
  ch_tx2gene          = Channel.fromPath(params.tx2gene)

  EDGER_DGELIST (
    ch_metadata,
    ch_quantfiles,
    ch_tx2gene
  )

  ch_DGEList          = EDGER_DGELIST.out.DGEList

  EDGER_QC_PLOTS (
    ch_DGEList
  )

  DIFF_EXPRESSION (
    ch_DGEList,
    ch_ensembl_dataset
  )

/*
  LONGTABLE (
    ch_DGEList,
    ch_ensembl_dataset
  )







  DMR_TABLE (
    METHYLKIT.out.collect(),
    ch_ensembl_dataset,
    ch_cpgislands_GRCm39,
    ch_refseq_UCSC_GRCm39,
    genomic_features
  )

  DMG_TABLE (
    DMR_TABLE.out.DMR_TABLES.collect(),
    ch_ensembl_dataset,
    DMG_table_output,
    genomic_features
  )

  GENE_ONTOLOGY (
    DMR_TABLE.out.DMR_TABLES.collect(),
    ch_ensembl_dataset,
    genomic_features
  )

  SUPPLEMENTARY_EXCEL (
    DMR_TABLE.out.DMR_TABLES.collect(),
    genomic_features
  )

  SUPPLEMENTARY_PLOTS (
    ch_metadata,
    DMR_TABLE.out.DMR_TABLES.collect(),
    BETAVALUES.out.BETAVALUES.collect(),
    genomic_features
  )

  WRAPPER ( // concatinate all tables w/ genomic features
    DMR_TABLE.out.DMR_TABLES.collect(),
    DMG_TABLE.out.DMG_TABLES.collect(),
    BETAVALUES.out.BETAVALUES.collect(),
    GENE_ONTOLOGY.out.GO_TABLES.collect()
  )

  REPORT (
    WRAPPER.out.collect(),
    ch_metadata,
    ch_ensembl_dataset,
    generations,
    genomic_features
  )

*/  
}

workflow.onComplete {
  def msg = """\
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Pipeline execution summary
    --------------------------
    Completed at     : ${workflow.complete}
    Duration         : ${workflow.duration}
    Success          : ${workflow.success}
    workDir          : ${workflow.workDir}
    exit status      : ${workflow.exitStatus}
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """
    .stripIndent()

  sendMail (
    to: "${params.email}",
    subject: 'Transcriptome analysis',
    body: msg
  )
}