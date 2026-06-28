[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![R](https://img.shields.io/badge/-script-276DC3.svg?style=flat&logo=R)](https://cran.r-project.org)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

# README

This repository holds a nextflow pipeline for analysing **gene expression** studies. The the pipeline allows for experimental design with multiple doses (0, 10, 100, 1000, etc). The pipeline expect quantification files (`quant.sf`) with a transcript-to-gene index file (`tx2gene.tsv`), both generated with Salmon, and sample metadata file as input. The pipeline outputs tables of (1) **differentially expressed genes** and (2) **gene ontology** analysis results which are combined to a report generated with Quarto. The pipeline also outputs the tables as `excel` files to be included as supplementary tables in a scientic report.

> [!IMPORTANT]
> This pipeline is optimised for SLURM on a high-performance computing (HPC) cluster ([UPPMAX](https://docs.uppmax.uu.se/cluster_guides/pelle/)).

<details><summary>Quantification files</summary><br>

  To generate methylation coverage files from sequencing files refer to [nf-core/rnaseq pipeline](https://nf-co.re/rnaseq/latest)
  </details>

<details><summary>Differential expression</summary><br>

  Differential expression was analysed with the R-package [`DESeq2`](https://bioconductor.org/packages//release/bioc/html/DESeq2.html).
  </details>

<details><summary>Gene ontology analysis</summary><br>

  To investigate if any biological functions, processes or pathways are enriched (over-represented) the _Over Representation Analysis (ORA)_ [Boyle et al., 2004](https://doi.org/10.1093/bioinformatics/bth456) method is used. ORA uses hypergeometric distribution and compares the differentially methylated genes with all genes in the dataset. The _p_-values are adjusted to _q_-values for multiple corretion (significance threshold `qvalue < 0.2`).

  Enrichment is analysed in three databases; (1) Gene Ontology (**GO**), (2) Kyoto Encyclopedia of Genes and Genomes (**KEGG**), and **Reactome** pathways. GO and KEGG enrichment are tested with the R-package [`clusterProfiler`](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html), [Yu et al., 2012](https://doi.org/10.1089/omi.2011.0118), [Wu et al., 2021](https://doi.org/10.1016/j.xinn.2021.100141). The reactome pathways are tested with the R-package [`ReactomePA`](https://bioconductor.org/packages/release/bioc/html/ReactomePA.html), [Yu et al., 2016](https://doi.org/10.1039/C5MB00663E).
  </details>

## Reproducibility

### Run the pipeline

To start the pipeline, wrap the code below into a shell script. Make sure to change the paths to the correct files/folders.

```sh
NXF_HOME=".nextflow/"

nextflow pull andreyhgl/transcriptome-analysis

PROJECT='slurm-projID'
QUANT_FILES='path/to/quant_files'
TX2GENE_FILE='path/to/tx2gene.tsv'
SAMPLE_INFO_FILE='path/to/sample_info'
SPECIES='mouse/human'
ENSEMBL_VERSION=115

nextflow run main.nf \
  -profile uppmax \
  --project "$PROJECT" \
  --diff_analysis_package 'DESeq2' \
  --quant_files "$QUANT_FILES" \
  --tx2gene "$TX2GENE_FILE" \
  --sample_info "$SAMPLE_INFO_FILE" \
  --species "$SPECIES" \
  --ensembl_version "$ENSEMBL_VERSION"
```

<details><summary>Singularity containers</summary><br>

  For reproducibility this pipeline uses two singularity containers, which can be downloaded from the [Cloud Library](https://cloud.sylabs.io/library). The `RNAseq` container holds most of the R-packages used in the analysis, while `gene-ontology` container holds gene ontology related R-packages

  ```sh
  # apptainer/singulartiy also works

  IMAGE1='library://andreyhgl/singularity-r/rnaseq'
  IMAGE2='library://andreyhgl/singularity-r/gene-ontology'

  apptainer pull "$IMAGE1"
  apptainer pull "$IMAGE2"
  ```

  To run scripts manually with the containers use the `exec` flag or run the script interactively with `shell`.

  ```sh
  # execute script
  apptainer exec "$IMAGE" <scriptfile>

  # run script interactively
  apptainer shell "$IMAGE"
  $ Rscript <scriptfile>
  ```
  </details>

### The pipeline

The nextflow pipeline produce the following:

+ Ensembl database table containing gene annotations
+ Quality control plots: PCA, distance plots
+ Differentially expressed genes table
+ Gene ontology analysis
+ Supplementary files (plots, excel-tables)
+ Concatinated tables (for easy import for results report)
+ A final report that summaries the results

## Preparation

### Metadata

Setup the `metadata.csv` with each row representing a sample in the column order: sample id and treatment/dose:

```
id,dose
sample1,0
sample2,0
sample3,10
sample4,100
sample5,10
sample6,100

etc....
```