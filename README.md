Analysis pipeline, from quant-files => calculation of differential gene expression => gene annotation analysis => pathway analysis => report.

# Computational pipeline

The pipeline comes with some default parameters and will look inside folders `data/`, for metadata files and `results/` for quantification files, assuming the alignment and quantification has be performed with [nf-core/rnaseq](https://nf-co.re/rnaseq/latest) pipeline.

```css
project/
├── README.md
├── bin/
├── conf/
│   ├── base.config
│   └── modules.config
├── modules/
├── subworkflows/
│
├── data/ 						# metadata for the project
├── results/ 					# nf-core/rnaseq output
│
├── main.nf
└── nextflow.config
```

## Metadata

The pipeline requires a metadata file with columns: `sample`. Pipeline defaults to `data/metadata.csv`, but can be manually set: `--metadata 'path/to/metadata.csv'`.

The DESeq2 formula design is decided by a column name and is passed with `--deseq2_design 'column-name'`.

```sh
nextflow run andreyhgl/transcriptome-analysis \
	--metadata 'path/to/metadata.csv' \
	--deseq2_design 'column-name'
```

## Quantification files

Defaults to `--quant_dir results/star_salmon`

```sh
nextflow run andreyhgl/transcriptome-analysis \
	--metadata 'path/to/metadata.csv' \
	--quant_dir 'path/to/quant_files' \
	--deseq2_design 'column-name' 
```

## Diferrential gene expression

## Genomic annotation

Categories genes according to their genomic positions and overlapping regulatory regions; transposable elements, enhancers, imprinted genes.

## Pathway analysis

Over-representation analaysis. Comparison of genesets from gene ontology, Reactome, KEGG. (clusterProfiler, ReatcomePA).

## Report

Interactive quarto results report, 

# Run

Analysis pipeline is initiated via `bin/run_analysis.sh`. All parameters are set via CLI upon pipeline start. The parameters come with defaults.

To start pipeline:

```sh
bash bin/run_analysis.sh
```

Inside the script the nextflow pipeline is initiated:

```sh
nextflow run andreyhgl/transcriptome_analysis \
	-profile apptainer \
	-params params.yml \
	-resume \
	-bg
```

Parameters are set (or default loaded) into `params.yml` and included in the pipeline:

```yml
# params.yml
output: analysis
quant_files_dir:
ensembl_genome: human
ensembl_version: 115
```

## Environment

The pipeline uses docker containers:

+ differential gene expression
+ gene annotations
+ pathway analysis
+ quarto report
