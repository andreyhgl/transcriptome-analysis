# Analysis pipeline

## Processes

+ Produce quant files
+ [Check sample quality (`SAMPLE_QC`)](#Sample-quality-control)
+ Calculate differential gene expression
+ Add annotation to genes
+ Pathway enrichment:
  + ORA
  + GSEA

# Computational pipeline

> [!IMPORTANT]
> The pipeline loads with default parameters and looks inside folders `data/` and `results/` for metadata files and quantification files, respectively. Quantification files are generated with [nf-core/rnaseq](https://nf-co.re/rnaseq/latest) pipeline.

```sh
project/
├── assets/                  #│ Quarto files contained within assets/
│   ├── styles.css           #│ shared by all reports
│   ├── quarto_report/
│   │   ├── _quarto.yml
│   │   ├── transcriptome_report.qmd
│   │   └── sample_qc.qmd
│   └── quarto_dashboard/
│       ├── _quarto.yml
│       └── sample_qc_dashboard.qmd
│
├── bin/                     #│ Mostly R-scripts
│
├── conf/                    #│ pipeline configs
│   ├── base.config
│   └── modules.config
│
├── modules/local/           #│
├── subworkflows/local/      #│
│   ├── deseq2/              #│ Processes and workflows
│   ├── gene_annotation/     #│
│   ├── pathway_enrichment/  #│
│   └── sample_qc/           #│
│
├── main.nf
├── nextflow.config
│
├── data/                    #│ contains metadata for the project
├── results/                 #│ contains quantfiles, output from nf-core/rnaseq 
│
└── README.md
```

## Metadata

The pipeline requires a metadata file with columns: `sample`. Pipeline defaults to path: `data/metadata.csv`, but can be manually set: `--metadata 'path/to/metadata.csv'`. The DESeq2 formula design is decided by a column name and is passed with `--deseq2_design 'column-name'`.

```sh
nextflow run andreyhgl/transcriptome-analysis \
  --metadata 'path/to/metadata.csv' \
  --deseq2_design 'column-name'
```

```
# metadata.csv

sample,group,season,year,concentration
1,control,control,NA,0
2,winter,winter,2022,200
3,spring,spring,2022,200
4,summer,summer,2022,200
5,autumn,autumn,2022,200
6,winter,winter,2023,200
```

## Quantification files

Defaults to `--quant_dir results/star_salmon`

```sh
nextflow run andreyhgl/transcriptome-analysis \
  --metadata 'path/to/metadata.csv' \
  --quant_dir 'path/to/quant_files' \
  --deseq2_design 'column-name' 
```

## Sample quality control

Sample quality control checks samples prior to the differential expression analysis. The check is blind to experimental design (DESeq2 design: `~1`). The check runs in parallel to the main analysis and can be turned off with `--skip_sample_qc`. Every sample with both metadata and a quant-file is checked.

+ **Per-sample metrics:** library size, mapping rate (from salmon's `aux_info/meta_info.json`), detected genes, share of reads in the 10 most abundant genes, size factor, and median correlation with the other samples.
+ **PCA** on the vst counts of the most variable genes. Each principal component is regressed on every metadata column and technical metric, to show what drives it.
+ **Clustering** on Spearman correlations of the expressed genes.
+ **Confounded metadata:** pairs of categorical columns compared with [Cramér's V](https://en.wikipedia.org/wiki/Cram%C3%A9r%27s_V), and nested columns.
+ **Flags:** a sample is flagged when a metric is more than `--sample_qc_mad` (default: 3) median absolute deviations, and at least a minimum difference, from the median, on the bad side (low library size, mapping rate, detected genes or correlation; high top-10 share). The minimum differences are set in `bin/sample_qc.R` and stop very tight metrics from giving trivial flags, such as a mapping rate of 99.1 % against 99.5 %

A flag says that a sample differs from the others, not why, and no sample is removed automatically. To leave a sample out, remove it from the metadata and rerun with `-resume`.

The results are written to `<outdir>/sample_qc/`:

| File | Content |
|---|---|
| `sample_qc.log` | Summary printed by the script, with the flagged samples |
| `sample_qc_samples.tsv` | One row per sample: metadata, metrics, PCs and flags |
| `sample_qc.Rds` | All results, read by the report and the dashboard |
| `*.html` | Quarto report (`sample_qc.qmd`) and dashboard (`sample_qc_dashboard.qmd`) |

| Parameter | Default | |
|---|---|---|
| `--skip_sample_qc` | `false` | Skip the sample QC |
| `--sample_qc_sample_col` | `sample` | Metadata column with the sample ids |
| `--sample_qc_columns` | `''` | Metadata columns to test, comma separated; `''` = every column that separates the samples |
| `--sample_qc_group` | `''` | Column for the groups and colours of the report; `''` = the first |
| `--sample_qc_ntop` | `500` | Most variable genes in the PCA |
| `--sample_qc_min_count` | `10` | Reads a gene needs in half of the samples to be used for the correlations |
| `--sample_qc_mad` | `3` | Median absolute deviations from the median before a sample is flagged |

To preview the report or the dashboard locally:

```sh
quarto preview assets/quarto_dashboard/sample_qc_dashboard.qmd \
  -P qc:<outdir>/sample_qc/sample_qc.Rds
```

## Differential gene expression

## Genomic annotation

Categorises genes according to their genomic positions and overlapping regulatory regions; transposable elements, enhancers, imprinted genes.

## Pathway analysis

Over-representation analysis. Comparison of genesets from gene ontology, Reactome, KEGG. (clusterProfiler, ReactomePA).

## Report

Interactive quarto results report, 

# Run

To start the pipeline:

```sh
nextflow run andreyhgl/transcriptome-analysis \
  -profile apptainer \
  -params-file params.yml \
  -resume \
  -bg
```

All parameters comes with defaults, set inside `nextflow.config`

```groovy
params {
    // Input
    metadata  = 'data/metadata.csv'
    tx2gene   = 'results/star_salmon/salmon.merged.tx2gene.tsv'
    quant_dir = 'results/star_salmon'

    ensembl_version = 115
    species         = 'human'

    // Sample QC
    skip_sample_qc       = false
    sample_qc_sample_col = 'sample'
    sample_qc_columns    = ''
    sample_qc_group      = ''
    sample_qc_ntop       = 500
    sample_qc_min_count  = 10
    sample_qc_mad        = 3

    // DESeq2
    deseq2_design        = 'group'
    deseq2_min_count     = 10
    deseq2_max_zero_frac = 0.5
    deseq2_alpha         = 0.05
    deseq2_lfc_threshold = 0
```

## Environment

The pipeline uses docker containers:

+ sample quality control
+ differential gene expression
+ gene annotations
+ pathway analysis
+ quarto report