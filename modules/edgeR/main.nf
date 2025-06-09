process EDGER {
  tag 
  //label

  conda "${moduleDir}/environment.yml"
  container 'library://andreyhgl/singularity-r/rnaseq:latest'

  time 1.h
  memory 15.GB
  cpus 2

  input:
  path metadata
  path quant_files
  path tx2gene

  output:
  path '*.Rds', emit: DGEList

  """
  edgeR.R
  """
}