process DIFF_EXPRESSION {
  tag 
  //label

  conda "${moduleDir}/environment.yml"
  container 'library://andreyhgl/singularity-r/rnaseq:latest'

  time 1.h
  memory 8.GB
  cpus 1

  input:
  path DGEList

  output:
  path 'DGE_table.csv', emit: DGE_table

  """
  diff_expression.R
  """
}