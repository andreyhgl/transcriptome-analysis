process ENSEMBL {
  tag "Building ensembl table: ${species} v. ${ensembl_version}"
  label 'process_single'
  
  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' ?
    'docker://ghcr.io/karlssonlaboratory/methylkit-env:6b7f121' :
    'ghcr.io/karlssonlaboratory/methylkit-env:6b7f121' }"

  input:
  val species
  val ensembl_version

  output:
  path 'ensembl_table.csv.gz', emit: ensembl_table

  script:
  """
  # set environment variables for biomart cache
  export BIOMART_CACHE="./cache/"

  ensembl.R \
    --species ${species} \
    --ensembl_version ${ensembl_version}
  """
}