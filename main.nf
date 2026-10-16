#!/usr/bin/env nextflow

// ~~ Import processes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
include { GENE_ANNOTATION } from './subworkflows/local/gene_annotation'
include { DESEQ2          } from './subworkflows/local/deseq2'

// ~~ printoout pipeline parameters functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *

params.monochrome_logs = false

// ANSI colour codes, or empty strings when colour is unwanted
def ansiCodes(monochrome) {
    def codes = [
        reset : '\u001B[0m',
        bold  : '\u001B[1m',
        dim   : '\u001B[2m',
        green : '\u001B[32m',
        yellow: '\u001B[33m',
        blue  : '\u001B[34m',
        purple: '\u001B[35m',
        cyan  : '\u001B[36m'
    ]
    // toBoolean(): `--monochrome_logs false` can arrive as the string 'false'
    def off = monochrome.toString().toBoolean()
    return off ? codes.collectEntries { k, _v -> [k, ''] } : codes
}

def col(codes, colour, text) {
    return "${codes[colour]}${text}${codes.reset}"
}

def paramsSummary(monochrome) {
    def c = ansiCodes(monochrome)
    def arrow = col(c, 'dim', '>')
    def sub = col(c, 'dim', '-')

    return """
    ${col(c, 'purple', '---- Pipeline parameters ---------------------------------')}
    ${arrow} Outdir                  : ${col(c, 'green', params.outdir)}
    ${arrow} quant-files dir         : ${col(c, 'green', params.quant_dir)}
    ${arrow} tx2gene file            : ${col(c, 'green', params.tx2gene)}
    ${arrow} metadata                : ${col(c, 'green', params.metadata)}

    ${col(c, 'purple', '---- Ensembl Parameters ----------------------------------')}
    ${arrow} Version                 : ${col(c, 'cyan', params.ensembl_version)}
    ${arrow} Species                 : ${col(c, 'cyan', params.species)}

    ${col(c, 'purple', '---- Experimental parameters -----------------------------')}
    ${arrow} ${col(c, 'bold', 'DESeq2')}
      ${sub} design                : ${col(c, 'yellow', params.deseq2_design)}
      ${sub} min. count            : ${col(c, 'yellow', params.deseq2_min_count)}
      ${sub} max. zero fraction    : ${col(c, 'yellow', params.deseq2_max_zero_frac)}
      ${sub} alpha                 : ${col(c, 'yellow', params.deseq2_alpha)}
      ${sub} LFC threshold         : ${col(c, 'yellow', params.deseq2_lfc_threshold)}
      ${sub} QC fill               : ${col(c, 'yellow', params.qc_fill)}
      ${sub} QC PCA n_top          : ${col(c, 'yellow', params.qc_ntop)}
      ${sub} Gene plots top_n      : ${col(c, 'yellow', params.gene_plots_top_n)}
      ${sub} Gene annotation       : ${col(c, 'yellow', params.gene_annotation)}
      ${sub} Gene plots labels     : ${col(c, 'yellow', params.gene_plots_labels)}

    ${col(c, 'purple', '----------------------------------------------------------')}
    """.stripIndent(true)
}

workflow {

    main:
    log.info(paramsSummary(params.monochrome_logs))

    ch_metadata = channel.fromPath(params.metadata, checkIfExists: true)
    ch_tx2gene  = channel.fromPath(params.tx2gene, checkIfExists: true)

    // One directory per sample; the folder name is the sample id
    ch_quant = channel
        .fromPath("${params.quant_dir}/*/quant.sf", checkIfExists: true)
        .map { quant -> quant.parent }


    // Gene annotation for the gene plots:
    if (params.gene_annotation) {
        // a provided table
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


    DESEQ2(ch_metadata, ch_tx2gene, ch_quant, ch_annotation)
}
