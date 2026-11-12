// ~~ Print pipeline parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ *

// ANSI colour codes, or empty strings when colour is unwanted
def ansiCodes(monochrome) {
    def codes = [
        reset : '\u001B[0m',
        bold  : '\u001B[1m',
        dim   : '\u001B[2m',
        red   : '\u001B[31m',
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

// One section per block of the log: a title, a colour for the values and
// the rows as label: value (printed in this order)
def paramsSections(p) {
    return [
        [
            title : 'Input',
            colour: 'green',
            rows  : [
                'Outdir'         : p.outdir,
                'Quant-files dir': p.quant_dir,
                'tx2gene file'   : params.tx2gene,
                'Metadata'       : params.metadata
            ]
        ],
        [
            title : 'Quarto',
            colour: 'green',
            rows  : [
                'DE report'       : params.quarto_report_file,
                'Sample QC report': params.sample_qc_report_file,
                'Metadata'        : params.quarto_meta_file,
                'CSS'             : params.quarto_css_file
            ]
        ],
        [
            title : 'Ensembl',
            colour: 'cyan',
            rows  : [
                'Version': params.ensembl_version,
                'Species': params.species
            ]
        ],
        [
            title : 'Sample QC',
            colour: 'yellow',
            rows  : [
                'Skip': params.skip_sample_qc
            ]
        ],
        [
            title : 'DESeq2',
            colour: 'yellow',
            rows  : [
                'Design'            : params.deseq2_design,
                'Min. count'        : params.deseq2_min_count,
                'Max. zero fraction': params.deseq2_max_zero_frac,
                'Alpha'             : params.deseq2_alpha,
                'LFC threshold'     : params.deseq2_lfc_threshold,
                'QC fill'           : params.qc_fill,
                'QC PCA n_top'      : params.qc_ntop,
                'Gene plots top_n'  : params.gene_plots_top_n,
                'Gene plots labels' : params.gene_plots_labels,
                'Gene annotation'   : params.gene_annotation ?: 'Ensembl BioMart'
            ]
        ]
    ]
}

def paramsSummary(p) {
    def c = ansiCodes(p.monochrome_logs)
    def blocks = paramsSections(p).collect { s ->
        renderBlock(c, s.title, s.colour, s.rows)
    }
    return '\n' + (blocks + col(c, 'purple', '-' * 58)).join('\n\n') + '\n'
}

// A purple header line, then the rows as `> label : value`
def renderBlock(c, title, colour, rows) {
    def width = 58
    def label_width = 22
    def arrow = col(c, 'dim', '>')

    def header = col(c, 'purple', "---- ${title} ".padRight(width, '-'))
    def lines = rows.collect { label, value ->
        "${arrow} ${label.padRight(label_width)}: ${col(c, colour, value)}"
    }
    return ([header] + lines).join('\n')
}

// Status, duration, task counts and where to find the results; on failure
// also the error and the command to resume
def completionSummary(p, wf) {
    def c = ansiCodes(p.monochrome_logs)
    def stats = wf.stats

    def status = wf.success
        ? col(c, 'green', 'Completed successfully')
        : col(c, 'red', 'Failed')

    def rows = [
        'Status'  : status,
        'Run name': wf.runName,
        'Duration': wf.duration,
        'Tasks'   : "${stats.succeededCount} succeeded, " +
                    "${stats.cachedCount} cached, " +
                    "${stats.failedCount} failed",
        'Results' : p.outdir,
        'Work dir': wf.workDir
    ]
    if (stats.ignoredCount > 0) {
        rows['Ignored tasks'] = stats.ignoredCount
    }

    def blocks = [renderBlock(c, 'Pipeline finished', 'cyan', rows)]

    if (!wf.success) {
        def error = wf.errorMessage?.readLines()?.find { line -> line.trim() } ?: '-'
        def resume = wf.commandLine.contains('-resume')
            ? wf.commandLine
            : "${wf.commandLine} -resume"

        blocks << renderBlock(c, 'Error', 'red', [
            'Exit status': wf.exitStatus,
            'Message'    : error,
            'Resume with': resume
        ])
    }

    return '\n' + (blocks + col(c, 'purple', '-' * 58)).join('\n\n') + '\n'
}