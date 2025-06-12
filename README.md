[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![R](https://img.shields.io/badge/-script-276DC3.svg?style=flat&logo=R)](https://cran.r-project.org)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

# README

This repository holds a nextflow pipeline for analysing gene expression studies. The experimental design is build on multiple generations (F0, F1, F2, etc) and doses (0, 10, 100, 1000, etc). The pipeline expect quantification files (`quant.sf`) as input generated with Salmon and outputs tables of (1) **differentially expressed genes** and (2) **gene ontology** analysis results. These tables come in `Rds` file format containing all generations and doses in a single file, respectively. The pipeline also outputs the tables as `excel` files to be included as supplementary tables in a scientic report.

> [!NOTE]
> This pipeline is per default setup for the [mouse genome (GRCm39)](https://www.ensembl.org/Mus_musculus/)

<details>
  <summary>Quantification files</summary>

>To generate methylation coverage files from sequencing files refer to [nf-core/rnaseq pipeline](https://nf-co.re/rnaseq/latest)

</details>

<details>
  <summary>Differential expression</summary>

>Differential expression was analysed with the R-package [`edgeR`](https://bioconductor.org/packages/release/bioc/html/edgeR.html), which utilizes negative binomial distributions and generalized linear model as statistical method. `FDR < 0.05` was used for multiple testing correction (Benjamini-Hochberg qvalue). Default settings were used for most of the functions expect; `estimateDisp(robust = TRUE)` and `glmQLFit(robust = TRUE)`. Summaries findings in a long table w/ a significant gene as a unique row, add results from the differentail gene expression analysis as columns.

</details> 

<details>
  <summary>Gene ontology analysis</summary>

>To investigate if any biological functions, processes or pathways are enriched (over-represented) the _Over Representation Analysis (ORA)_ [Boyle et al., 2004](https://doi.org/10.1093/bioinformatics/bth456) method is used. ORA uses hypergeometric distribution and compares the differentially methylated genes with all genes in the dataset. The _p_-values are adjusted to _q_-values for multiple corretion (significance threshold `qvalue < 0.2`).

>Enrichment is analysed in three databases; (1) Gene Ontology (**GO**), (2) Kyoto Encyclopedia of Genes and Genomes (**KEGG**), and **Reactome** pathways. GO and KEGG enrichment are tested with the R-package [`clusterProfiler`](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html), [Yu et al., 2012](https://doi.org/10.1089/omi.2011.0118), [Wu et al., 2021](https://doi.org/10.1016/j.xinn.2021.100141). The reactome pathways are tested with the R-package [`ReactomePA`](https://bioconductor.org/packages/release/bioc/html/ReactomePA.html), [Yu et al., 2016](https://doi.org/10.1039/C5MB00663E). 

</details>

## The pipeline

The nextflow pipeline produce the following:

+ PCA plots (for sanity check of experiment)
+ DGEList object
+ Differentially expressed genes table
+ Gene ontology analysis
+ Supplementary files (plots, excel-tables)
+ Concatinated tables (for easy import for results report)

## Preparation

### Metadata

Setup the `data/metadata.csv` to look like the following:

```csv
gen,id,treatment,...
F0,F0_1,0,...
F0,F0_2,0,...
F1,F1_3,10,...
F2,F2_4,100,...
```

Each row represents a sample in the column order: generation, sample id and treatment/dose.

## Parameters

The pipeline accepts three parameters:

Experimental design:

+ generations (F0, F1, F2, etc)
+ doses (0, 10, 100, 1000, etc)
+ genomic features (CpG-sites, Promoters, CpG-islands)

## Reproducibility

For reproducibility this pipeline uses two singularity containers, which can be downloaded from the [Cloud Library](https://cloud.sylabs.io/library). The `RNAseq` container holds most of the R-packages used in the analysis, while `gene-ontology` container holds gene ontology related R-packages

```sh
# apptainer (instead of singulartiy) also works

IMAGE1='library://andreyhgl/singularity-r/rnaseq:latest'
IMAGE2='library://andreyhgl/singularity-r/gene-ontology:latest'

singularity pull ${IMAGE1}
singularity pull ${IMAGE2}
```

<details>
  <summary>Run interactively</summary>

To run scripts manually with the containers use the `exec` flag or run the script interactively with `shell`.

```sh
# execute script
singularity exec ${IMAGE} <scriptfile>

# run script interactively
singularity shell ${IMAGE}
$ Rscript <scriptfile>
```

</details>

## Run the pipeline

```sh
#!/bin/bash -l

export NXF_HOME=".nextflow/"

nextflow pull andreyhgl/transcriptome-analysis

nextflow run andreyhgl/transcriptome-analysis \
  --quant_path 'path-to-quant-files' \
  --metadata 'path-to-metadata.csv' \
  --tx2gene 'path-to-tx2gene.tsv' \
  --generation "F0,F1,F2" \
  --treatment "10,100" \
  -profile local \
  -resume
```