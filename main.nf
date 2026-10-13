#!/usr/bin/env nextflow

// ~~ Import processes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
include { DESEQ2 } from './subworkflows/local/deseq2'

// ~~ verbose pipeline settings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
params.monochrome_logs = false

def c = [
  reset : '\u001B[0m',
  bold  : '\u001B[1m',
  dim   : '\u001B[2m',
  green : '\u001B[32m',
  yellow: '\u001B[33m',
  blue  : '\u001B[34m',
  purple: '\u001B[35m',
  cyan  : '\u001B[36m'
]

def col = { colour, text -> "${c[colour]}${text}${c.reset}" }
def arrow = col('dim', '>')
def sub = col('dim', '-')

// blank out all codes when colour is unwanted
if (params.monochrome_logs) c = c.collectEntries { k, v -> [k, ''] }

log.info \
  """
  ${col('purple', '---- Pipeline parameters ---------------------------------')}
  ${arrow} Outdir                 : ${col('green', params.outdir)}
  ${arrow} quant-files dir        : ${col('green', params.quant_files_dir)}
  ${arrow} tx2gene file           : ${col('green', params.tx2gene_file)}
  ${arrow} metadata               : ${col('green', params.metadata)}
  
  ${col('purple', '---- Ensembl Parameters ----------------------------------')}
  ${arrow} Genome                 : ${col('cyan', params.ensembl_genome)}
  ${arrow} Version                : ${col('cyan', params.ensembl_version)}
  ${arrow} Species                : ${col('cyan', params.species)}

  ${col('purple', '---- Experimental parameters -----------------------------')}
  ${arrow} ${col('bold', 'DESeq2')}
    ${sub} design                 : ${col('yellow', params.deseq2_design)}
    ${sub} min. count             : ${col('yellow', params.deseq2_min_count)}
    ${sub} max. zero fraction     : ${col('yellow', params.deseq2_max_zero_frac)}
  
  ${col('purple', '----------------------------------------------------------')}

  """
  .stripIndent(true)

workflow {

    main:
    ch_metadata = channel.fromPath(params.metadata, checkIfExists: true)
    ch_tx2gene  = channel.fromPath(params.tx2gene, checkIfExists: true)

    // One directory per sample; the folder name is the sample id
    ch_quant = channel
        .fromPath("${params.quant_dir}/*/quant.sf", checkIfExists: true)
        .map { quant -> quant.parent }

    DESEQ2(ch_metadata, ch_tx2gene, ch_quant)
}
