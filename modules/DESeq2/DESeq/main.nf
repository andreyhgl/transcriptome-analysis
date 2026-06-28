process DESEQ2_DESEQ {
  tag "Building DESeqDataSet object"
  label 'process_medium'

  conda "${moduleDir}/environment.yml"
  container 'library://andreyhgl/singularity-r/rnaseq'

  //container "${ workflow.containerEngine == 'singularity' ?
  //  'docker://ghcr.io/karlssonlaboratory/methylkit-env:6b7f121' :
  //  'ghcr.io/karlssonlaboratory/methylkit-env:6b7f121' }"

  input:
  path metadata
  path quant_files
  path tx2gene

  output:
  path 'DDS.Rds', emit: DDS

  script:
  """
  DESeq2_DESeq.R
  """
}