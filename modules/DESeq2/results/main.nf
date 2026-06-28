process DESEQ2_RESULTS {
  tag "Calculating DE genes"
  label 'process_medium'

  conda "${moduleDir}/environment.yml"
  container 'library://andreyhgl/singularity-r/rnaseq'

  //container "${ workflow.containerEngine == 'singularity' ?
  //  'docker://ghcr.io/karlssonlaboratory/methylkit-env:6b7f121' :
  //  'ghcr.io/karlssonlaboratory/methylkit-env:6b7f121' }"

  input:
  path ensembl
  path DDS

  output:
  path 'genexp_table.csv.gz', emit: genexp_table
  path 'significant_genes_S1.xlsx'

  script:
  """
  DESeq2_results.R
  """
}