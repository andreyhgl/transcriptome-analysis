# Analysis pipeline

+ Produce quant files
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
├── assets/quarto_report/    #│ Quarto files contained within assets/
│   ├── _quarto.yml
│   ├── styles.css
│   └── report.qmd
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
│   └── pathway_enrichment/  #│
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

    // DESeq2
    deseq2_design        = 'group'
    deseq2_min_count     = 10
    deseq2_max_zero_frac = 0.5
    deseq2_alpha         = 0.05
    deseq2_lfc_threshold = 0
```

## Environment

The pipeline uses docker containers:

+ differential gene expression
+ gene annotations
+ pathway analysis
+ quarto report
