#!/usr/bin/env nextflow

// ~~ Import processes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *
include { DESEQ2 } from './subworkflows/local/deseq2'

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
    ${arrow} Outdir                 : ${col(c, 'green', params.outdir)}
    ${arrow} quant-files dir        : ${col(c, 'green', params.quant_files_dir)}
    ${arrow} tx2gene file           : ${col(c, 'green', params.tx2gene_file)}
    ${arrow} metadata               : ${col(c, 'green', params.metadata)}

    ${col(c, 'purple', '---- Ensembl Parameters ----------------------------------')}
    ${arrow} Genome                 : ${col(c, 'cyan', params.ensembl_genome)}
    ${arrow} Version                : ${col(c, 'cyan', params.ensembl_version)}
    ${arrow} Species                : ${col(c, 'cyan', params.species)}

    ${col(c, 'purple', '---- Experimental parameters -----------------------------')}
    ${arrow} ${col(c, 'bold', 'DESeq2')}
      ${sub} design                 : ${col(c, 'yellow', params.deseq2_design)}
      ${sub} min. count             : ${col(c, 'yellow', params.deseq2_min_count)}
      ${sub} max. zero fraction     : ${col(c, 'yellow', params.deseq2_max_zero_frac)}

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

    DESEQ2(ch_metadata, ch_tx2gene, ch_quant)
}
