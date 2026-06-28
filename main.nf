#!/usr/bin/env Nextflow

/*
 * Github: https://github.com/andreyhgl/transcriptome-analysis
 */

// ~~ Import processes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *

include { METADATA            } from './modules/metadata/main'
include { ENSEMBL             } from './modules/ensembl/main'
include { DESEQ2_DESEQ        } from './modules/DESeq2/DESeq/main'
include { DESEQ2_RESULTS      } from './modules/DESeq2/results/main'
include { DESEQ2_QC_PLOTS     } from './modules/DESeq2/QC_plots/main'
include { DESEQ2_DGE_PLOTS    } from './modules/DESeq2/DGE_plots/main'
include { PATHWAY_ENRICHMENT  } from './modules/pathway_enrichment/main'
include { QUARTO_NOTEBOOK     } from './modules/Quarto/notebook/main'
//include { QUARTO_REPORT       } from './modules/Quarto/report/main'

log.info \
  """
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   > Diff. analysis R-packages  : ${params.diff_analysis_package}
   > Outdir                     : ${params.outdir}
  
  ~~ Ensembl arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   > Genome                     : ${params.ensembl_genome}
   > Version                    : ${params.ensembl_version}
   > Species                    : ${params.species}
  
  ~~ Experiment arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   > Samples info               : ${params.sample_info}
   > Seq. files metadata?       : ${params.seqfiles_metadata}
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  """
  .stripIndent(true)

// ~~ Workflow ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *

workflow {

  /*
    Work-around for passing an 'empty' argument to a process
    input: path() cannot handle null, false or '', but and empty list is OK, []
  */

  sample_info = params.sample_info ? file ( params.sample_info, checkIfExists: true) : []

  seqfiles_metadata = params.seqfiles_metadata ? file ( params.seqfiles_metadata, checkIfExists: true) : []
    
  ENSEMBL (
    params.species,
    params.ensembl_version
  )

  ch_ensembl_table          = ENSEMBL.out.ensembl_table

  METADATA (
    sample_info,
    seqfiles_metadata,
    samples_to_remove
  )

  ch_metadata               = METADATA.out.metadata

  DESEQ2_DESEQ (
    ch_metadata,
    params.quant_files,
    params.tx2gene
  )

  ch_DDS                    = DESEQ2_DESEQ.out.DDS


  DESEQ2_RESULTS (
    ch_ensembl_table,
    ch_DDS
  )

  ch_genexp_table           = DESEQ2_RESULTS.out.genexp_table
//  ch_count                  = DESEQ2_RESULTS.out.unique_genes_count

  DESEQ2_QC_PLOTS (
    ch_DDS
  )

  DESEQ2_DGE_PLOTS (
    ch_DDS,
    ch_genexp_table
  )

  PATHWAY_ENRICHMENT (
    params.species,
    ch_ensembl_table,
    ch_DDS,
    ch_genexp_table
  )

  ch_pathway_enrichment     = PATHWAY_ENRICHMENT.out.pathway_enrichment_table

  QUARTO_NOTEBOOK (
    params.quarto_report_file,
    params.quarto_meta,
    params.quarto_css,
    ch_ensembl_table,
    ch_DDS,
    ch_genexp_table,
    ch_pathway_enrichment
  )

  //ch_count.splitText().view( num -> " \n> Number of unique genes: $num\n")

}

workflow.onComplete {

  if (workflow.success) {

    def summary = """
    ==============================================
     Pipeline Complete
    ==============================================
     Completed at             : ${workflow.complete}
     Duration                 : ${workflow.duration}
     Success                  : ${workflow.success}
     Exit status              : ${workflow.exitStatus}
     Work dir                 : ${workflow.workDir}
     Run name                 : ${workflow.runName}
     Outdir                   : ${params.outdir}
    ----------------------------------------------
     Parameters
    ----------------------------------------------
     Coverage files           : ${params.quant_files}
     Seq. files metadata      : ${params.seqfiles_metadata}
     Ensembl Genome           : ${params.ensembl_genome}
     Ensembl Version          : ${params.ensembl_version}
     Species                  : ${params.species}
    ----------------------------------------------
     Results
    ----------------------------------------------
    """.stripIndent()

    /*
     * Place holder, adjust ongoing project 
    

    // Read each log file and append
    def log_files = [
      "logs/table1.log",
      "logs/table2.tsv"
    ]

    log_files.each { rel_path ->
      def f = new File("${params.outdir}/${rel_path}")
      if (f.exists()) {
        summary += " [${f.name}]\n"
        summary += f.text + "\n"
      } else {
        summary += " [${rel_path}] not found\n\n"
      }
    }
    */


    summary += "==============================================\n"
    summary += " End of Report\n"
    summary += "==============================================\n"

    def timestamp = workflow.complete.format("yyyy-MM-dd_HH-mm-ss")
    def log_file = new File("${params.outdir}/logs/pipeline_summary_${timestamp}.log")
    log_file.text = summary
    log.info summary

  } else {
    log.error "Pipeline failed, exit: ${workflow.exitStatus}"
  }
}